#include "pksystem.hpp"

momoko::pks::pksystem::pksystem(base::ideal_lattice &_latt)
    : latt(_latt), rng(std::random_device{}()) {}

momoko::base::ideal_lattice_element
momoko::pks::pksystem::make_uniform_element(long a, long b,
                                            bool disable_NTT_cache) {
  std::vector<long> factors(latt.getN());
  std::uniform_int_distribution<long> dist(a, b);
  for (auto &x : factors) {
    x = dist(dev);
  }
  return latt.make_element(factors, disable_NTT_cache);
}
