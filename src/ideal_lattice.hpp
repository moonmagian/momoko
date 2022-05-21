#ifndef IDEAL_LATTICE_HPP
#define IDEAL_LATTICE_HPP
#include "ideal_lattice_element.hpp"
#include <optional>
#include <string>
#include <vector>
#include "tools.hpp"
using std::optional;
namespace momoko::base {
class ideal_lattice_element;
///
/// \brief An ideal lattice: Z[x] / (x^(2^n) + 1).
/// n should be 2^k and q should be 1 module 2n.
/// By default no check is made for performance, turn on the check with
/// LATT_PARAM_CHECK flag.
///
class ideal_lattice {
  friend class ideal_lattice_element;

  private:
  static constexpr char header_cache[4]{'L', 'A', 'T', '+'};
  static constexpr char header_no_cache[4]{'L', 'A', 'T', '-'};
  unsigned long n;
  unsigned long q;
  tools::element_type_T n_inv;
  std::optional<tools::element_type_T> psi;
  std::optional<std::vector<tools::element_type_T>> psi_list;
  std::optional<std::vector<tools::element_type_T>> psi_inv_list;
  bool use_NTT_cache = true;

  public:
  ideal_lattice(unsigned long n, unsigned long q, bool use_NTT_cache = true);
  ideal_lattice(std::istream &is, bool use_NTT_cache = true);
  tools::element_type_T get_2nth_root_of_unity();
  ideal_lattice_element make_element();
  ideal_lattice_element import_element(std::istream &is);
  ideal_lattice_element make_element(const std::vector<long> &factor,
                                     bool disable_NTT_cache = false);
  ideal_lattice_element make_element_from_NTT(const std::vector<long> &factor);
  std::vector<long> INTT_GS(std::vector<long> a);
  void export_to_stream(std::ostream &os, bool include_psi_cache = false);

  const std::vector<tools::element_type_T> &get_psi_list();
  const std::vector<tools::element_type_T> &get_psi_inv_list();
  unsigned long getN() const;
  unsigned long getQ() const;
  bool operator==(const ideal_lattice &other) const;
};
} // namespace momoko::base
#endif // IDEAL_LATTICE_HPP
