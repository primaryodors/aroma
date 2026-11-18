#include "dock_worker.h"
#include "constants.h"

#if AROMA_OPT_GUARDED_THREADS && defined(__linux__)
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <system_error>
#include <thread>
#include <fcntl.h>
#include <glob.h>
#include <sched.h>
#include <sys/file.h>
#include <unistd.h>

namespace
{
using Clock = std::chrono::steady_clock;

struct CpuSample
{
    unsigned long long total = 0, idle = 0;
    unsigned cpus = 0;
};

CpuSample sample_cpu()
{
    std::ifstream stat("/proc/stat");
    std::string line, name;
    std::getline(stat, line);
    std::istringstream aggregate(line);
    CpuSample sample;
    if (!(aggregate >> name) || name != "cpu") return {};
    for (int i = 0; i < 8; ++i)
    {
        unsigned long long ticks;
        if (!(aggregate >> ticks)) return {};
        sample.total += ticks;
        if (i == 3 || i == 4) sample.idle += ticks;
    }
    while (std::getline(stat, line) && line.size() > 3 && line.compare(0, 3, "cpu") == 0
           && line[3] >= '0' && line[3] <= '9') ++sample.cpus;
    return sample;
}

bool has_cpu_capacity(int threads)
{
    cpu_set_t cpus;
    if (sched_getaffinity(0, sizeof(cpus), &cpus) != 0 || CPU_COUNT(&cpus) < threads) return false;
    std::ifstream groups("/proc/self/cgroup");
    std::string group, line;
    while (std::getline(groups, line))
        if (line.compare(0, 4, "0::/") == 0)
        {
            group = line.substr(3);
            break;
        }
    // A parent's CPU quota also limits every child group.
    for (;;)
    {
        std::ifstream quota("/sys/fs/cgroup" + group + "/cpu.max");
        if (quota.is_open())
        {
            std::string limit;
            double period;
            if (!(quota >> limit >> period) || !std::isfinite(period) || period <= 0) return false;
            if (limit != "max")
            {
                char* end = nullptr;
                const double value = std::strtod(limit.c_str(), &end);
                if (end == limit.c_str() || *end || !std::isfinite(value) || value < threads * period)
                    return false;
            }
        }
        if (group.empty()) return true;
        group.resize(group.rfind('/'));
    }
}

bool cpu_is_cool(double limit)
{
    glob_t sensors{};
    const int status = glob("/sys/class/hwmon/hwmon*/temp*_input", 0, nullptr, &sensors);
    bool found = false, cool = status == 0;
    for (size_t i = 0; cool && i < sensors.gl_pathc; ++i)
    {
        const std::string input = sensors.gl_pathv[i];
        const std::string directory = input.substr(0, input.rfind('/') + 1);
        const std::string prefix = input.substr(0, input.size() - 6);
        std::string driver, label;
        std::ifstream(directory + "name") >> driver;
        std::getline(std::ifstream(prefix + "_label"), label);
        const bool intel = driver == "coretemp" && label.compare(0, 11, "Package id ") == 0;
        const bool amd = (driver == "k10temp" || driver == "zenpower")
            && (label == "Tctl" || label == "Tdie" || (label.empty() && input == directory + "temp1_input"));
        if (!intel && !amd) continue;
        found = true;
        double reading;
        if (!(std::ifstream(input) >> reading) || !std::isfinite(reading) || reading <= 0)
        {
            cool = false;
            break;
        }
        double sensor_limit = limit, high;
        if ((std::ifstream(prefix + "_max") >> high) && high > 0)
            sensor_limit = std::min(sensor_limit, high / 1000 - 10);
        cool = reading / 1000 < sensor_limit;
    }
    globfree(&sensors);
    return found && cool;
}

class DockWorker
{
    static constexpr int max_workers = 3;
    bool requested = false, stopping = false;
    int admitted = 0, created = 0;
    int cool_samples = 0, lock_fd = -1;
    double temp_limit = 65;
    Clock::time_point sampled_at{};
    CpuSample previous{};
    std::mutex mutex;
    std::condition_variable changed;
    const std::function<void(int, int)>* task = nullptr;
    unsigned long long generation = 0;
    int batch_workers = 0, pending = 0;
    std::array<std::thread, max_workers> workers;

    void work(int index)
    {
        unsigned long long seen = 0;
        std::unique_lock<std::mutex> lock(mutex);
        for (;;)
        {
            changed.wait(lock, [&] {
                return stopping || (task && generation != seen && index <= batch_workers);
            });
            if (stopping) return;
            seen = generation;
            const auto* job = task;
            const int threads = batch_workers + 1;
            lock.unlock();
            (*job)(index, threads);
            lock.lock();
            if (--pending == 0)
            {
                task = nullptr;
                changed.notify_all();
            }
        }
    }

    bool start_worker()
    {
        if (lock_fd < 0)
        {
            lock_fd = open("/tmp/aromadock-worker.lock", O_RDONLY | O_CREAT | O_CLOEXEC | O_NOFOLLOW, 0644);
            if (lock_fd < 0) return false;
            if (flock(lock_fd, LOCK_EX | LOCK_NB) != 0)
            {
                close(lock_fd);
                lock_fd = -1;
                return false;
            }
        }
        try
        {
            workers[created] = std::thread(&DockWorker::work, this, created + 1);
            ++created;
            return true;
        }
        catch (const std::system_error&)
        {
            if (!created)
            {
                close(lock_fd);
                lock_fd = -1;
            }
            return false;
        }
    }

public:
    bool enable()
    {
        const char* setting = std::getenv("AROMA_CPU_MAX_C");
        if (setting)
        {
            char* end = nullptr;
            temp_limit = std::strtod(setting, &end);
            if (end == setting || *end || !std::isfinite(temp_limit) || temp_limit <= 0 || temp_limit > 90)
                return false;
        }
        requested = true;
        return true;
    }

    bool available()
    {
        if (!requested) return false;
        const auto now = Clock::now();
        if (now - sampled_at < std::chrono::seconds(1)) return admitted;
        const CpuSample current = sample_cpu();
        const bool recent = now - sampled_at <= std::chrono::seconds(3);
        sampled_at = now;
        const bool valid = recent && current.cpus && current.cpus == previous.cpus
            && current.total > previous.total && current.idle >= previous.idle
            && current.idle - previous.idle <= current.total - previous.total;
        double busy = 0;
        if (valid)
        {
            const double idle = double(current.idle - previous.idle) / (current.total - previous.total);
            busy = (1 - idle) * current.cpus;
        }
        previous = current;
        if (!valid || busy > admitted + 1.25 || !has_cpu_capacity(admitted + 1) || !cpu_is_cool(temp_limit))
        {
            if (admitted) std::cerr << "AromaDock: extra workers paused by CPU or temperature guard.\n";
            admitted = cool_samples = 0;
            return false;
        }

        // Observe a new window at each worker count before increasing concurrency again.
        if (admitted < max_workers && ++cool_samples >= 3)
        {
            cool_samples = 0;
            if (has_cpu_capacity(admitted + 2) && (admitted < created || start_worker()))
            {
                ++admitted;
                std::cerr << "AromaDock: " << admitted + 1 << " docking threads active.\n";
            }
        }
        return admitted > 0;
    }

    void run(const std::function<void(int, int)>& job)
    {
        if (!admitted)
        {
            job(0, 1);
            return;
        }
        {
            std::lock_guard<std::mutex> lock(mutex);
            task = &job;
            batch_workers = pending = admitted;
            ++generation;
        }
        changed.notify_all();
        job(0, admitted + 1);
        std::unique_lock<std::mutex> lock(mutex);
        changed.wait(lock, [&] { return !task; });
    }

    ~DockWorker()
    {
        {
            std::lock_guard<std::mutex> lock(mutex);
            stopping = true;
        }
        changed.notify_all();
        for (int i = 0; i < created; ++i) workers[i].join();
        if (lock_fd >= 0) close(lock_fd);
    }
};

DockWorker dock_worker;
}

bool enable_dock_worker() { return dock_worker.enable(); }
bool dock_worker_available() { return dock_worker.available(); }
void run_dock_worker(const std::function<void(int, int)>& job)
{
    dock_worker.run(job);
}

#else
bool enable_dock_worker() { return false; }
bool dock_worker_available() { return false; }
void run_dock_worker(const std::function<void(int, int)>& job)
{
    job(0, 1);
}
#endif
