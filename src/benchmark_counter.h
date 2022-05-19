#ifndef BENCHMARK_COUNTER_H
#define BENCHMARK_COUNTER_H
#include <future>
#include <thread>
using stopable_callable = void (*)(std::stop_token token,
                                   std::promise<long> promise);

template <typename T> long start_benchmark(T callable, int sec) {
  std::stop_source should_stop;
  std::promise<long> promise;
  std::future<long> future{promise.get_future()};
  std::thread t(callable, should_stop.get_token(), std::move(promise));
  std::this_thread::sleep_for(std::chrono::seconds(sec));
  should_stop.request_stop();
  t.join();
  return future.get();
}

#endif // BENCHMARK_COUNTER_H
