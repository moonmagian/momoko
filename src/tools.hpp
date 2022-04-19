#ifndef TOOLS_HPP
#define TOOLS_HPP
#include <cmath>
#include <iostream>
#include <vector>
namespace momoko::tools {
using ulong = unsigned long;
constexpr double SQRT_2PI{2.5066282746310002};
ulong mod_pow(ulong a, ulong b, ulong p);
ulong mod_inv(long a, long p);
std::vector<ulong> fatorize(ulong n);
///
/// \brief bit_reverse: Reverse the bits of the original number i.
/// \param i: The number to reverse the bits.
/// \param max_bit_n: Only reverse over the first n bits (LSBs).
/// \return: The reversed number.
///
ulong bit_reverse(ulong i, ulong max_bit_n);
///
/// \brief bitrevorder: Get the bit reversed order of the original vector in
/// place.
/// \param vec: Size must be pow of 2.
///
void bitrevorder(std::vector<ulong> &vec);
///
/// \brief mod_reduce: reduce the number i to the range [0, q) module to q.
///
ulong mod_reduce(long i, long q);
ulong mod_reduce_lazy(long i, long q);
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
