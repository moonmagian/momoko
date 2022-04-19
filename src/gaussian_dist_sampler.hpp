#ifndef GAUSSIAN_DIST_SAMPLER_HPP
#define GAUSSIAN_DIST_SAMPLER_HPP

#include "ideal_lattice.hpp"
#include <random>
namespace momoko::gaussian {

class gaussian_dist_sampler {
  protected:
  momoko::base::ideal_lattice &latt;
  std::mt19937 rng;
  gaussian_dist_sampler(momoko::base::ideal_lattice &latt);
  virtual long sample_halfside_gaussian() = 0;

  public:
  virtual double get_std_var() = 0;
  virtual double get_tail() = 0;
  virtual long sample_gaussian();
  virtual momoko::base::ideal_lattice_element
  sample_lattice_element(bool pos_only = false);
};

} // namespace momoko::gaussian

#endif // GAUSSIAN_DIST_SAMPLER_HPP
