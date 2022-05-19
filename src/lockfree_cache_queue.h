#ifndef LOCKFREE_CACHE_QUEUE_H
#define LOCKFREE_CACHE_QUEUE_H
#include <array>
#include <atomic>
#include <concepts>
namespace momoko::tools {
template <typename T, size_t N> class lockfree_cache_queue {
  static_assert(N > 0, "size should be more than zero");

  private:
  std::atomic_size_t head_pos{0};
  std::atomic_size_t end_pos{0};
  std::array<T, N> ring;

  public:
  bool push(T v) {
    long end_pos_v = end_pos.load(std::memory_order_relaxed);
    // First push, then update the state.
    // Not doing full check, it's the semaphore's job.
    ring[end_pos_v] = v;
    if (end_pos_v == N - 1) {
      end_pos = 0;
    } else {
      ++end_pos;
    }
    return true;
  }
  bool pop(T &out) {
    long head_pos_v = head_pos.load(std::memory_order_relaxed);
    // First check, then pop, then update the state.

    // This line is not thread safe, but this is spsc, so consider:
    // 1. end_pos changed before the load of end_pos: there must have been a new
    // element to consume, so it's safe.
    // 2. end_pos changed after the load:
    //  2.1. returning false: in next loop, we are free to use the new value.
    //  2.2. not returning false: the element is available anyway.
    // All the assumptions work because of the semaphore in background sampler,
    // it blocks pushing if the queue is full.
    if (head_pos_v == end_pos.load(std::memory_order_relaxed)) {
      return false;
    }
    out = ring[head_pos_v];
    if (head_pos_v == N - 1) {
      head_pos = 0;
    } else {
      ++head_pos;
    }
    return true;
  }
};
} // namespace momoko::tools

#endif // LOCKFREE_CACHE_QUEUE_H
