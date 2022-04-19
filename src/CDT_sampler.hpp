#ifndef CDT_SAMPLER_HPP
#define CDT_SAMPLER_HPP
#include "gaussian_dist_sampler.hpp"
#include <iostream>
#include <random>
namespace momoko::gaussian {
class CDT_sampler : public gaussian_dist_sampler {
  private:
  double SD;
  double tail;
  std::vector<double> CDT;
  std::uniform_real_distribution<double> real_dist;
  std::bernoulli_distribution sign_dist;

  protected:
  long sample_halfside_gaussian() final override;
  void read_param_from_stream(std::istream &s);

  public:
  // CDT!
  static constexpr char header[4]{'C', 'D', 'T', '!'};
  CDT_sampler(double SD, double tail, base::ideal_lattice &latt);
  CDT_sampler(std::istream &s, base::ideal_lattice &latt);
  CDT_sampler(const std::string &path, base::ideal_lattice &latt);
  void export_param(std::ostream &s);
  void export_param(const std::string &path);
  double get_std_var() final override;
  double get_tail() final override;
  long sample_gaussian() final override;
  bool operator==(const CDT_sampler &other) const;
};
} // namespace momoko::gaussian

#endif // CDT_SAMPLER_HPP
