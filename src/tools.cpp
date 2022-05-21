#include "tools.hpp"
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <stdexcept>
momoko::tools::element_type_T
momoko::tools::mod_pow(momoko::tools::element_type_T a, unsigned long b,
                       unsigned long p) {
  element_type_T res = 1;
  while (b)
    if (b & 1) {
      res = mod_reduce(res * a, p);
      --b;
    } else {
      a = mod_reduce(a * a, p);
      b >>= 1;
    }
  return res;
}

std::vector<unsigned long> momoko::tools::fatorize(unsigned long n) {
  std::vector<unsigned long> result;
  for (unsigned long i = 2; i * i <= n; ++i) {
    if (n % i == 0) {
      result.push_back(i);
      while (n % i == 0) {
        n /= i;
      }
    }
  }
  if (n > 1) {
    result.push_back(n);
  }
  return result;
}

unsigned long momoko::tools::bit_reverse(unsigned long i,
                                         unsigned long max_bit_n) {
  unsigned long result{0ul};
  max_bit_n = std::min(sizeof(unsigned long) * 8, max_bit_n);
  for (unsigned int x = 0; i != 0 && x < max_bit_n; ++x) {
    if ((i & 1ul)) {
      result |= 1ul << (max_bit_n - x - 1);
    }
    i >>= 1;
  }
  return result;
}

void momoko::tools::_print_ulong_bit(unsigned long i) {
  for (unsigned int x = 0; x < sizeof(unsigned long) * 8; ++x) {
    if ((i & 1ul)) {
      std::cout << "1";
    } else {
      std::cout << "0";
    }
    i >>= 1;
  }
}
momoko::tools::element_type_T momoko::tools::mod_inv(long a, long p) {
  long b = p, u = 1, v = 0;
  while (b) {
    long t = a / b;
    a -= t * b;
    std::swap(a, b);
    u -= t * v;
    std::swap(u, v);
  }
  u = mod_reduce(u, p);
  return u;
}

uint64_t momoko::tools::double_frac_to_integer(double d) {
  uint64_t *ptr{reinterpret_cast<uint64_t *>(&d)};
  uint64_t result{0};
  // Copy 52 bits to result.
  result |= (*ptr) & (0x0FFFFFFFFFFFFFULL);
  result <<= 11;
  result |= (1ULL << 63);
  // Get exp part.
  uint16_t exp = ((*ptr) >> 52) & 0x07FF;
  // Bad condition, drop it.
  if (exp >= 1023) {
    throw std::runtime_error(
        "Convert error. Only numbers less or equal than 1 is supported.");
  }
  // Real offset.
  exp = 1022 - exp;
  // Right shift more than 64 bits will cause strange problem, return 0
  // directly.
  if (exp >= 64) {
    return 0;
  }
  result >>= exp;
  return result;
}

momoko::tools::element_type_T momoko::tools::mod_reduce_lazy(long i, long q) {
  if (i > (1L << (sizeof(long) * 8 / 4)) ||
      i < -(1L << (sizeof(long) * 8 / 4))) {
    return mod_reduce(i, q);
  } else {
    return i;
  }
}
