#ifndef PKSYSTEM_HPP
#define PKSYSTEM_HPP
#include "ideal_lattice.hpp"
#include "gaussian_dist_sampler.hpp"
#include <random>
#include <iostream>
namespace momoko::pks {
class pksystem {
  protected:
  base::ideal_lattice &latt;
  std::random_device dev;
  std::mt19937 rng;

  public:
  pksystem(base::ideal_lattice &_latt);
  base::ideal_lattice_element
  make_uniform_element(long a, long b, bool disable_NTT_cache = false);
  virtual void export_sk(std::ostream &os) = 0;
  virtual void import_sk(std::istream &is) = 0;
  virtual void export_pk(std::ostream &os) = 0;
  virtual void import_pk(std::istream &is) = 0;
  virtual bool pk_ready() = 0;
  virtual bool sk_ready() = 0;
};
} // namespace momoko::pks

#endif // PKSYSTEM_HPP
