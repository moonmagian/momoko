#include "bernoulli_sampler.hpp"
#include <cmath>
#include <stdexcept>
namespace momoko::gaussian {
long bernoulli_sampler::bin_std_var_halfside_sample() {
  for (;;) {
  algorithm_loop:
    auto zero_sign = rng();
    if ((zero_sign & 1) == 0) {
      return 0;
    }
    ulong i;
    for (i = 1;; ++i) {
      // Needed bit count.
      ulong k{2 * i - 1};
      while (k > 0) {
        // Use all bits of the next generated element.
        if (k > sizeof(decltype(rng)::result_type) * 8) {
          ulong bits{rng()};
          if (bits != 0) {
            // Restart the algorithm.
            goto algorithm_loop;
          }
          k -= sizeof(decltype(rng)::result_type) * 8;
        } else {
          // Only use part of the generated bits.
          ulong bits{rng()};
          // 1, 2, 3, ..., k - 1 th bits should be zero, otherwise restart the
          // algorithm.
          if ((bits & ((1ul << (k - 1)) - 1)) != 0) {
            goto algorithm_loop;
          }
          // And return i if the kth bit is zero.
          if ((bits & ((1ul) << (k - 1))) == 0) {
            return i;
          }
          // Else continue the loop.
          break;
        }
      }
    }
  }
}

bernoulli_sampler::bernoulli_sampler(ulong k, base::ideal_lattice &latt)
    : gaussian_dist_sampler{latt}, k{k}, udist_0_k{0, k - 1},
      sd{0.8493218002880191 * static_cast<double>(k)} {
  if (k == 0) {
    throw std::runtime_error("k should be larger than 0 in bernoulli sampler.");
  }
}

double bernoulli_sampler::get_std_var() { return sd; }

double bernoulli_sampler::get_tail() {
  return std::numeric_limits<double>::infinity();
}

long bernoulli_sampler::sample_halfside_gaussian() {
  switch (k) {
  case 1:
    return bin_std_var_halfside_sample();
  default:
    for (;;) {
      long x = bin_std_var_halfside_sample();
      long y = udist_0_k(rng);
      long z = k * x + y;
      auto cache_index = std::make_pair(x, y);
      if (exp_cache.find(cache_index) == exp_cache.end()) {
        exp_cache[cache_index] =
            std::exp(-(double)(y) * (y + 2 * k * x) / (double)(2 * sd * sd));
      }
      double refuse_rate = exp_cache[cache_index];
      std::bernoulli_distribution dist(refuse_rate);
      bool accept = dist(rng);
      if (!accept)
        continue;
      return z;
    }
  }
}
} // namespace momoko::gaussian
