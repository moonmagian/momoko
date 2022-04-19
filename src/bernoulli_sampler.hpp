#ifndef BERNOULLI_SAMPLER_HPP
#define BERNOULLI_SAMPLER_HPP
#include "gaussian_dist_sampler.hpp"
#include <map>
namespace momoko::gaussian {
class bernoulli_sampler : public gaussian_dist_sampler {
private:
  ulong k;
  long bin_std_var_halfside_sample();
  std::uniform_int_distribution<ulong> udist_0_k;
  std::map<std::pair<ulong, ulong>, double> exp_cache;

public:
  const double sd;
  bernoulli_sampler(ulong k, base::ideal_lattice &latt);
  double get_std_var() override;
  double get_tail() override;

protected:
  long sample_halfside_gaussian() override;
};
} // namespace momoko::gaussian

#endif // BERNOULLI_SAMPLER_HPP
