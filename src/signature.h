#ifndef SIGNATURE_H
#define SIGNATURE_H
#include "pksystem.hpp"
#include <tuple>
namespace momoko::pks {
class signature : public pksystem {
  private:
  gaussian::gaussian_dist_sampler &sampler;
  std::optional<base::ideal_lattice_element> s1;
  std::optional<base::ideal_lattice_element> s2;
  std::optional<base::ideal_lattice_element> a;
  std::optional<base::ideal_lattice_element> t;
  std::uniform_int_distribution<long> distn;
  std::uniform_int_distribution<long> dist1;
  base::ideal_lattice_element sample_uniform();
  base::ideal_lattice_element sample_short_uniform();

  public:
  signature(base::ideal_lattice &latt,
            gaussian::gaussian_dist_sampler &sampler);
  void export_sk(std::ostream &os) override;
  void import_sk(std::istream &is) override;
  void export_pk(std::ostream &os) override;
  void import_pk(std::istream &is) override;
  bool pk_ready() override;
  bool sk_ready() override;
  std::tuple<base::ideal_lattice_element, base::ideal_lattice_element,
             base::ideal_lattice_element>
  sign_latt_element(const base::ideal_lattice_element &message);
  bool verify_latt_element(base::ideal_lattice_element &message,
                           std::pair<base::ideal_lattice_element,
                                     base::ideal_lattice_element> &sign);
};
} // namespace momoko::pks

#endif // SIGNATURE_H
