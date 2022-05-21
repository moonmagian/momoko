#ifndef BACKGROUND_SAMPLER_H
#define BACKGROUND_SAMPLER_H
#include "gaussian_dist_sampler.hpp"
#include <semaphore>
#include <atomic>
#include <deque>
#include <thread>
#include <mutex>
namespace momoko::gaussian {

template <std::ptrdiff_t N>
class background_sampler : public gaussian_dist_sampler {
  private:
  std::counting_semaphore<N> current{0};
  std::counting_semaphore<N> empty{N};
  gaussian_dist_sampler &sampler;
  std::deque<long> cache;
  std::atomic<bool> stop{false};
  std::thread producer_thread;
  std::mutex mutex;
  void producer() {
    while (true) {
      // Block until the producer can read the value.
      empty.acquire();
      if (stop.load(std::memory_order_relaxed)) {
        break;
      }
      std::unique_lock lock(mutex);
      cache.push_back(sampler.sample_gaussian());
      // Notify consumers.
      current.release();
      lock.unlock();
    }
  }

  protected:
  long sample_halfside_gaussian() override {
    // Not used.
    return 0;
  }

  public:
  background_sampler(gaussian_dist_sampler &underlying_sampler,
                     base::ideal_lattice &latt)
      : sampler{underlying_sampler}, gaussian_dist_sampler{latt} {
    // Delay the start of the thread to ensure all things are ready.
    producer_thread = std::thread{&background_sampler::producer, this};
  }
  virtual ~background_sampler() {
    stop.store(true, std::memory_order_relaxed);
    if (current.try_acquire()) {
      empty.release();
    }
    producer_thread.join();
  }
  double get_std_var() override { return sampler.get_std_var(); }
  double get_tail() override { return sampler.get_tail(); }
  long sample_gaussian() override {
    // Block until there is a value to read.
    current.acquire();
    std::unique_lock lock(mutex);
    long result{cache.front()};
    cache.pop_front();
    lock.unlock();
    // Notify the producer.
    empty.release();
    return result;
  }

  base::ideal_lattice_element
  sample_lattice_element(bool pos_only = false) override {
    auto n = latt.getN();
    // Directly get the batch data to reduce lock waiting cost.
    std::vector<long> result(n);
    for (size_t i = 0; i < n; ++i) {
      current.acquire();
    }
    std::unique_lock lock(mutex);
    for (size_t i = 0; i < n; ++i) {
      result[i] = cache.front();
      if (pos_only && result[i] < 0) {
        result[i] = -i;
      }
      cache.pop_front();
    }
    empty.release(n);
    lock.unlock();
    return latt.make_element(result);
  }
};
} // namespace momoko::gaussian

#endif // BACKGROUND_SAMPLER_H
