#ifndef PKE_BLWE_H
#define PKE_BLWE_H
#include "pksystem.hpp"
#include "ideal_lattice_element.hpp"
#include <optional>
#include <memory>
namespace momoko::pks {
class pke_blwe : public pksystem {
  private:
  std::optional<base::ideal_lattice_element> s;
  std::optional<base::ideal_lattice_element> a;
  std::optional<base::ideal_lattice_element> b;
  std::optional<base::ideal_lattice_element> e;
  static constexpr char header_pk[4]{'P', 'K', 'E', '-'};
  static constexpr char header_sk[4]{'P', 'K', 'E', '+'};
  std::uniform_int_distribution<long> k_uniform_sampler;
  std::uniform_int_distribution<long> q_uniform_sampler;
  base::ideal_lattice_element
  sample_uniform(std::uniform_int_distribution<long> &dist);

  public:
  pke_blwe(base::ideal_lattice &latt);
  pke_blwe(std::istream &is);
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

#endif // PKE_BLWE_H
