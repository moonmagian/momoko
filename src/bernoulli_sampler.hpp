#ifndef BERNOULLI_SAMPLER_HPP
#define BERNOULLI_SAMPLER_HPP
#include "gaussian_dist_sampler.hpp"
#include <map>
namespace momoko::gaussian {
class bernoulli_sampler : public gaussian_dist_sampler {
private:
  unsigned long k;
  long bin_std_var_halfside_sample();
  std::uniform_int_distribution<unsigned long> udist_0_k;
  std::map<std::pair<unsigned long, unsigned long>, double> exp_cache;

public:
  const double sd;
  bernoulli_sampler(unsigned long k, base::ideal_lattice &latt);
  double get_std_var() override;
  double get_tail() override;

protected:
  long sample_halfside_gaussian() override;
};
} // namespace momoko::gaussian

#endif // BERNOULLI_SAMPLER_HPP
