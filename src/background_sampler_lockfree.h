#ifndef BACKGROUND_SAMPLER_LOCKFREE_H
#define BACKGROUND_SAMPLER_LOCKFREE_H
#ifdef BOOST_QUEUE_ENABLED
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#else
#include "lockfree_cache_queue.h"
#endif
#include <semaphore>
#include <atomic>
#include <thread>
#include "gaussian_dist_sampler.hpp"
#include "ideal_lattice_element.hpp"
namespace momoko::gaussian {
template <size_t N>
class background_sampler_lockfree : public gaussian_dist_sampler {
  private:
  std::counting_semaphore<N - 1> empty{N - 1};
#ifdef BOOST_QUEUE_ENABLED
  boost::lockfree::spsc_queue<long, boost::lockfree::capacity<N>> cache;
#else
  tools::lockfree_cache_queue<long, N> cache;
#endif

  gaussian_dist_sampler &sampler;
  std::atomic<bool> stop{false};
  std::thread producer_thread;
  void producer() {
    while (!stop.load(std::memory_order_relaxed)) {
      empty.acquire();
      cache.push(sampler.sample_gaussian());
    }
  }

  protected:
  long sample_halfside_gaussian() override {
    // Not used.
    return 0;
  }

  public:
  background_sampler_lockfree(gaussian_dist_sampler &sampler,
                              base::ideal_lattice &latt)
      : gaussian_dist_sampler{latt}, sampler{sampler} {
    // Delay the start of the thread to ensure all things are ready.
    producer_thread = std::thread{&background_sampler_lockfree::producer, this};
  }
  virtual ~background_sampler_lockfree() {
    stop.store(true, std::memory_order_relaxed);
    empty.release();
    producer_thread.join();
  }
  double get_std_var() override { return sampler.get_std_var(); }
  double get_tail() override { return sampler.get_tail(); }
  long sample_gaussian() override {
    long result;
    while (!cache.pop(result))
      ;
    empty.release();
    return result;
  }
  base::ideal_lattice_element
  sample_lattice_element(bool pos_only = false) override {
    size_t n = latt.getN();
    std::vector<long> result(n);
    for (size_t i = 0; i < n; ++i) {
      long pop_element;
      while (!cache.pop(pop_element))
        ;
      result[i] = pop_element;
      if (pos_only && result[i] < 0) {
        result[i] = -i;
      }
      empty.release();
    }
    return latt.make_element(result);
  }
};
} // namespace momoko::gaussian

#endif // BACKGROUND_SAMPLER_LOCKFREE_H
