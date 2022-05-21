#ifndef BENCHMARK_COUNTER_H
#define BENCHMARK_COUNTER_H
#include <future>
#include <thread>
#include <atomic>
using stopable_callable = void (*)(std::atomic<bool> token,
                                   std::promise<long> promise);

template <typename T> long start_benchmark(T callable, int milliseconds) {
  std::atomic<bool> should_stop;
  std::promise<long> promise;
  std::future<long> future{promise.get_future()};
  std::thread t(callable, std::ref(should_stop), std::move(promise));
  std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
  should_stop.store(true, std::memory_order_relaxed);
  t.join();
  return future.get();
}

#endif // BENCHMARK_COUNTER_H
