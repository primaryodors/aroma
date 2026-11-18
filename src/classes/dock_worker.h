#ifndef AROMA_DOCK_WORKER_H
#define AROMA_DOCK_WORKER_H

#include <functional>

bool enable_dock_worker();
bool dock_worker_available();
void run_dock_worker(const std::function<void(int, int)>& job);

#endif
