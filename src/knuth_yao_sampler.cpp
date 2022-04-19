#include "knuth_yao_sampler.hpp"
#include <stdexcept>
#include <iostream>
long momoko::gaussian::knuth_yao_sampler::sample_halfside_gaussian() {
  // Knuth-Yao sampler doesn't need to sample half, directly sample the whole is
  // enough.
  return 0;
}

momoko::gaussian::knuth_yao_sampler::knuth_yao_sampler(
    double _std_var, double _tail, base::ideal_lattice &latt)
    : gaussian_dist_sampler(latt), std_var{_std_var}, tail{_tail} {
  size_t N{static_cast<size_t>(std_var * tail)};
  if (N == 0) {
    throw std::runtime_error("N is too small for Knuth-Yao sampler.");
  }
  mat.push(tools::discrete_gaussian(0, std_var));
  for (size_t i = 1; i < N - 2; ++i) {
    auto converted_result =
        tools::double_frac_to_integer(tools::discrete_gaussian(i, std_var) * 2);
    // Other parts are small enough, throw them directly.
    if ((converted_result & 0xFFFFFFFFFFFFF700UL) == 0) {
      break;
    }
    mat.push(tools::discrete_gaussian(i, std_var) * 2);
  }
  mat.complete();
}

double momoko::gaussian::knuth_yao_sampler::get_std_var() { return std_var; }

double momoko::gaussian::knuth_yao_sampler::get_tail() { return tail; }

long momoko::gaussian::knuth_yao_sampler::sample_gaussian() {
  long d{0};
  bool hit = false;
  size_t row{0};
  for (size_t col = 0; col < mat.precision(); ++col) {
    d = 2 * d + (bit_randomizer(rng) ? 1 : 0) -
        mat.column_hamming_distance(col);
    if (d < 0) {
      for (row = 0; row < mat.size(); ++row) {
        d = d + (mat.get(row, col) ? 1 : 0);
        if (d == 0) {
          hit = true;
          break;
        }
      }
      if (hit) {
        break;
      }
    }
  }
  if (bit_randomizer(rng)) {
    return -static_cast<long>(row);
  }
  return row;
}

long momoko::gaussian::knuth_yao_sampler::_check_carry() {
  size_t carry = 0;
  int i;
  std::cout << mat << std::endl;
  for (i = mat.precision() - 1; i >= 0; --i) {
    size_t now = mat.column_hamming_distance(i) + carry;
    if (now % 2 == 1) {
      std::cout << "die at " << i << std::endl;
    }
    carry = now / 2;
  }
  std::cout << carry << std::endl;
  return carry;
}
