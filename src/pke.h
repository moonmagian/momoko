#ifndef PKE_H
#define PKE_H
#include "pksystem.hpp"
#include "ideal_lattice_element.hpp"
#include <optional>
#include <memory>
namespace momoko::pks {
class pke : public pksystem {
  private:
  gaussian::gaussian_dist_sampler &sampler;
  std::optional<base::ideal_lattice_element> s;
  std::optional<base::ideal_lattice_element> a;
  std::optional<base::ideal_lattice_element> b;
  std::optional<base::ideal_lattice_element> e;
  static constexpr char header_pk[4]{'P', 'K', 'E', '-'};
  static constexpr char header_sk[4]{'P', 'K', 'E', '+'};

  public:
  pke(base::ideal_lattice &latt, gaussian::gaussian_dist_sampler &sampler);
  pke(std::istream &is);
  void generate_keys();
  void export_sk(std::ostream &os) override;
  void import_sk(std::istream &is) override;
  void export_pk(std::ostream &os) override;
  void import_pk(std::istream &is) override;
  bool pk_ready() override;
  bool sk_ready() override;
  std::pair<base::ideal_lattice_element, base::ideal_lattice_element>
  encrypt_latt_element(base::ideal_lattice_element message);
  base::ideal_lattice_element decrypt_latt_element(
      const std::pair<base::ideal_lattice_element, base::ideal_lattice_element>
          &message);
};
} // namespace momoko::pks

#endif // PKE_H
