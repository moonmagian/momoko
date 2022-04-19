#ifndef IDEAL_LATTICE_ELEMENT_HPP
#define IDEAL_LATTICE_ELEMENT_HPP
#include "ideal_lattice.hpp"
#include <iostream>
#include <vector>
#include <optional>
namespace momoko::base {
class ideal_lattice;
using ulong = unsigned long;

class ideal_lattice_element {
  private:
  friend std::ostream &operator<<(std::ostream &os,
                                  const ideal_lattice_element &element);
  friend ideal_lattice_element operator-(ideal_lattice_element a);
  friend class ideal_lattice;
  friend ideal_lattice_element SPM_product(ideal_lattice_element a,
                                           const ideal_lattice_element &b);

  ideal_lattice &latt;
  std::vector<long> factors;
  bool use_NTT_cache = true;
  std::optional<std::vector<long>> NTT_cache;
  const ulong q;
  const ulong n;
  void normalize_factors();
  ideal_lattice_element(ideal_lattice &lattice,
                        const std::vector<long> &_factors = std::vector<long>{},
                        bool use_NTT_cache = true);
  ideal_lattice_element(ideal_lattice &lattice, std::vector<long> &&_factors,
                        bool use_NTT_cache = true);

  public:
  ideal_lattice_element(const ideal_lattice_element &other) = default;
  ideal_lattice_element &operator=(const ideal_lattice_element &other);
  ideal_lattice_element &operator*=(ideal_lattice_element &other);
  ideal_lattice_element &operator*=(const ulong &other);
  ideal_lattice_element &operator+=(const ideal_lattice_element &other);
  ideal_lattice_element &operator-=(const ideal_lattice_element &other);
  void export_to_stream(std::ostream &os) const;
  bool operator==(const ideal_lattice_element &other) const;
  long get_factor(size_t i) const;
  void set_factor(size_t i, long v);
  const std::vector<long> &NTT_CT();
};
ideal_lattice_element operator*(ideal_lattice_element a,
                                ideal_lattice_element &b);
ideal_lattice_element operator*(ideal_lattice_element a, ulong b);
ideal_lattice_element operator*(ulong b, ideal_lattice_element a);
ideal_lattice_element SPM_product(ideal_lattice_element a,
                                  const ideal_lattice_element &b);
ideal_lattice_element operator+(ideal_lattice_element a,
                                const ideal_lattice_element &b);
ideal_lattice_element operator-(ideal_lattice_element a);
ideal_lattice_element operator-(ideal_lattice_element a,
                                const ideal_lattice_element &b);
std::ostream &operator<<(std::ostream &os,
                         const ideal_lattice_element &element);
} // namespace momoko::base

#endif // IDEAL_LATTICE_ELEMENT_HPP
