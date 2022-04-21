#ifndef TOOLS_HPP
#define TOOLS_HPP
#include <cmath>
#include <iostream>
#include <vector>
namespace momoko::tools {
using element_type_T = long;
constexpr double SQRT_2PI{2.5066282746310002};
element_type_T mod_pow(element_type_T a, unsigned long b, unsigned long p);
element_type_T mod_inv(long a, long p);
std::vector<unsigned long> fatorize(unsigned long n);
///
/// \brief bit_reverse: Reverse the bits of the original number i.
/// \param i: The number to reverse the bits.
/// \param max_bit_n: Only reverse over the first n bits (LSBs).
/// \return: The reversed number.
///
unsigned long bit_reverse(unsigned long i, unsigned long max_bit_n);
///
/// \brief bitrevorder: Get the bit reversed order of the original vector in
/// place.
/// \param vec: Size must be pow of 2.
///
template <typename T> void bitrevorder(std::vector<T> &vec) {
  unsigned long max_bit = 0;
  size_t size = vec.size();
  if (size == 0) {
    return;
  }
  // Size is not power of 2.
  if ((size & (size - 1)) != 0) {
    throw std::runtime_error(
        "bitrevorder requires the vec to have size of power of 2.");
  }

  // Get max_bit_count.
  while (size != 1) {
    max_bit += 1;
    size >>= 1;
  }
  // Swap to get bit rev order.
  for (size_t i = 0; i < vec.size(); ++i) {
    size_t new_i = bit_reverse(i, max_bit);
    if (new_i > i) {
      std::swap(vec[i], vec[new_i]);
    }
  }
}
///
/// \brief mod_reduce: reduce the number i to the range [0, q) module to q.
///
inline element_type_T mod_reduce(long i, long q) {
  i %= q;
  if (i > 0) {
    return (i <= (q - 1) / 2) ? i : i - q;
  }
  if (i < 0) {
    return (i >= (1 - q) / 2) ? i : i + q;
  }
  return i;
};
element_type_T mod_reduce_lazy(long i, long q);
uint64_t double_frac_to_integer(double d);
// debug tools.
void _print_ulong_bit(ulong i);
template <typename T> void _print_vec(std::vector<T> vec) {
  for (const auto &x : vec) {
    std::cout << x << " ";
  }
}
constexpr double discrete_gaussian(double x, double SD) {
  // Use approx value: sqrt(2pi) * SD.
  // See also <Sampling from discrete Gaussians for lattice-based
  // cryptography on a constrained device>
  return std::exp(-(x * x) / (2.0 * SD * SD)) / (SQRT_2PI * SD);
};
} // namespace momoko::tools
#endif // TOOLS_HPP
