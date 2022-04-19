#ifndef SIGNATURE_H
#define SIGNATURE_H
#include "pksystem.hpp"
namespace momoko::pks {
class signature : public pksystem {
  private:
  std::optional<base::ideal_lattice_element> s;
  std::optional<base::ideal_lattice_element> e;
  std::optional<base::ideal_lattice_element> a;
  std::optional<base::ideal_lattice_element> neg_a;
  std::optional<base::ideal_lattice_element> neg_b;
  std::optional<base::ideal_lattice_element> b;
  std::uniform_int_distribution<ulong> dist;
  base::ideal_lattice_element sample_uniform();

  public:
  signature(base::ideal_lattice &latt,
            gaussian::gaussian_dist_sampler &sampler);
  void export_sk(std::ostream &os) override;
  void import_sk(std::istream &is) override;
  void export_pk(std::ostream &os) override;
  void import_pk(std::istream &is) override;
  bool pk_ready() override;
  bool sk_ready() override;
  std::pair<base::ideal_lattice_element, base::ideal_lattice_element>
  sign_latt_element(const base::ideal_lattice_element &message);
  bool verify_latt_element(base::ideal_lattice_element &message,
                           std::pair<base::ideal_lattice_element,
                                     base::ideal_lattice_element> &sign);
};
} // namespace momoko::pks

#endif // SIGNATURE_H
