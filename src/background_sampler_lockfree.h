#ifndef BACKGROUND_SAMPLER_LOCKFREE_H
#define BACKGROUND_SAMPLER_LOCKFREE_H
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <stop_token>
#include "gaussian_dist_sampler.hpp"
#include "ideal_lattice_element.hpp"
#include "lockfree_cache_queue.h"
namespace momoko::gaussian {
template <size_t N>
class background_sampler_lockfree : public gaussian_dist_sampler {
  private:
  std::counting_semaphore<N> empty{N};
  //  boost::lockfree::spsc_queue<long, boost::lockfree::capacity<N>> cache;
  tools::lockfree_cache_queue<long, N> cache;

  gaussian_dist_sampler &sampler;
  std::stop_source stop;
  std::thread producer_thread;
  void producer() {
    while (!stop.stop_requested()) {
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
    stop.request_stop();
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
