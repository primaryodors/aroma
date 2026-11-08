/* Benchmark-only control of the wall-clock seeds used by native sampling. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdlib.h>
#include <time.h>

time_t time(time_t *destination)
{
    const char *setting = getenv("AROMA_BASELINE_TIME");
    if (setting && *setting)
    {
        char *end = NULL;
        long long parsed = strtoll(setting, &end, 10);
        if (!end || *end) abort();
        time_t value = (time_t)parsed;
        if (destination) *destination = value;
        return value;
    }
    time_t (*original)(time_t *) = dlsym(RTLD_NEXT, "time");
    if (!original) abort();
    return original(destination);
}
