#ifndef KNUTH_YAO_SAMPLER_HPP
#define KNUTH_YAO_SAMPLER_HPP
#include "gaussian_dist_sampler.hpp"
#include "bit_matrix.hpp"
#include <random>
namespace momoko::gaussian {
class knuth_yao_sampler : public gaussian_dist_sampler {
  private:
  double std_var;
  double tail;
  tools::bit_matrix mat;
  std::bernoulli_distribution bit_randomizer;

  protected:
  long sample_halfside_gaussian() override;

  public:
  knuth_yao_sampler(double _std_var, double _tail, base::ideal_lattice &latt);
  double get_std_var() override;
  double get_tail() override;
  long sample_gaussian() override;
  long _check_carry();
};
} // namespace momoko::gaussian

#endif // KNUTH_YAO_SAMPLER_HPP
