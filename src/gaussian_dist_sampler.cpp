#include "gaussian_dist_sampler.hpp"
#include <vector>
namespace momoko::gaussian {

gaussian_dist_sampler::gaussian_dist_sampler(base::ideal_lattice &latt)
    : latt(latt), rng(std::random_device{}()) {}

long gaussian_dist_sampler::sample_gaussian() {
  bool sign;
  long half_result;
  do {
    sign = static_cast<bool>(rng() & 1);
    half_result = sample_halfside_gaussian();
    // Drop half of the zero samples.
  } while (half_result == 0 && sign);
  if (sign) {
    half_result = -half_result;
  }
  return half_result;
}

base::ideal_lattice_element
gaussian_dist_sampler::sample_lattice_element(bool pos_only) {
  auto n = latt.getN();
  std::vector<long> result(n);
  // Sum of a gaussian dist is also a gaussian dist.
  // In correct embedding, every axis is independent.
  for (ulong i = 0; i < n; ++i) {
    long r{sample_gaussian()};
    if (pos_only && r < 0) {
      r = -r;
    }
    result[i] = r;
  }
  return latt.make_element(result);
  // Make a vector from gaussian sample.
}

} // namespace momoko::gaussian
