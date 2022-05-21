#include "ideal_lattice_element.hpp"
#include "tools.hpp"
#include <utility>

void momoko::base::ideal_lattice_element::normalize_factors() {
  // Drop extra values.
  if (factors.size() > n) {
    factors.resize(n);
  }
  //  factors.resize(n);
  for (auto &x : factors) {
    x = tools::mod_reduce(x, q);
  }
}

momoko::base::ideal_lattice_element::ideal_lattice_element(
    ideal_lattice &lattice, const std::vector<long> &_factors,
    bool use_NTT_cache)
    : latt{lattice}, factors{_factors},
      use_NTT_cache{use_NTT_cache}, q{lattice.getQ()}, n{lattice.getN()} {
  normalize_factors();
}

momoko::base::ideal_lattice_element::ideal_lattice_element(
    ideal_lattice &lattice, std::vector<long> &&_factors, bool use_NTT_cache)
    : latt{lattice}, factors{std::forward<std::vector<long>>(_factors)},
      use_NTT_cache{use_NTT_cache}, q{lattice.getQ()}, n{lattice.getN()} {
  normalize_factors();
}

momoko::base::ideal_lattice_element &
momoko::base::ideal_lattice_element::operator=(
    const ideal_lattice_element &other) {
  if (latt != other.latt || n != other.n || q != other.q) {
    throw std::runtime_error(
        "Unable to assign an element from another lattice.");
  }
  factors = other.factors;
  NTT_cache = other.NTT_cache;
  return *this;
}

momoko::base::ideal_lattice_element &
momoko::base::ideal_lattice_element::operator*=(ideal_lattice_element &other) {
  auto ntt1 = NTT_CT();
  auto ntt2 = other.NTT_CT();
  for (size_t i = 0; i < ntt1.size(); ++i) {
    ntt1[i] = tools::mod_reduce(ntt1[i] * ntt2[i], q);
  }
  factors = latt.INTT_GS(ntt1);
  NTT_cache.value() = std::move(ntt1);
  return *this;
}

momoko::base::ideal_lattice_element &
momoko::base::ideal_lattice_element::operator*=(
    const tools::element_type_T &other) {
  for (size_t i = 0; i < factors.size(); ++i) {
    factors[i] = tools::mod_reduce(factors[i] * other, q);
  }
  if (NTT_cache.has_value()) {
    // Update the NTT cache at the same time.
    for (size_t i = 0; i < NTT_cache.value().size(); ++i) {
      NTT_cache.value()[i] = tools::mod_reduce(NTT_cache.value()[i] * other, q);
    }
  }
  return *this;
}

momoko::base::ideal_lattice_element &
momoko::base::ideal_lattice_element::operator+=(
    const ideal_lattice_element &other) {
  if (factors.size() < other.factors.size()) {
    factors.resize(other.factors.size(), 0);
  }
  for (size_t i = 0; i < other.factors.size(); ++i) {
    factors[i] = tools::mod_reduce(factors[i] + other.get_factor(i), q);
  }
  if (NTT_cache.has_value() && other.NTT_cache.has_value()) {
    for (size_t i = 0; i < NTT_cache->size(); ++i) {
      NTT_cache.value()[i] = tools::mod_reduce(
          NTT_cache.value()[i] + other.NTT_cache.value()[i], q);
    }
  } else if (NTT_cache.has_value()) {
    NTT_cache.reset();
  }
  return *this;
}

momoko::base::ideal_lattice_element &
momoko::base::ideal_lattice_element::operator-=(
    const ideal_lattice_element &other) {
  *this += -other;
  return *this;
}

void momoko::base::ideal_lattice_element::export_to_stream(
    std::ostream &os) const {
  os.write(reinterpret_cast<const char *>(&n), sizeof(n));
  os.write(reinterpret_cast<const char *>(&q), sizeof(q));
  size_t factor_count{factors.size()};
  os.write(reinterpret_cast<const char *>(&factor_count), sizeof(factor_count));
  for (size_t i = 0; i < factor_count; ++i) {
    os.write(reinterpret_cast<const char *>(&factors[i]), sizeof(factors[i]));
  }
}

bool momoko::base::ideal_lattice_element::operator==(
    const ideal_lattice_element &other) const {
  if (n != other.n || q != other.q) {
    return false;
  }
  for (size_t i = 0; i < std::max(factors.size(), other.factors.size()); ++i) {
    if (i >= factors.size() && other.factors[i] != 0) {
      return false;
    }
    if (i >= other.factors.size() && factors[i] != 0) {
      return false;
    }
    if (factors[i] != other.factors[i]) {
      return false;
    }
  }
  return true;
}

long momoko::base::ideal_lattice_element::get_factor(size_t i) const {
  if (i >= factors.size()) {
    return 0;
  }
  return factors[i];
}

void momoko::base::ideal_lattice_element::set_factor(size_t i, long v) {
  if (i >= n) {
    return;
  }
  if (i >= factors.size()) {
    factors.resize(i + 1);
  }
  factors[i] = v;
  NTT_cache.reset();
}

const std::vector<long> &momoko::base::ideal_lattice_element::NTT_CT() {
  if (NTT_cache.has_value() && use_NTT_cache) {
    return NTT_cache.value();
  }
  unsigned long t = n;
  auto phi_list = latt.get_psi_list();
  NTT_cache.emplace(factors);
  std::vector<long> &result = NTT_cache.value();
  result.resize(n, 0);
  for (unsigned long m = 1; m < n; m <<= 1) {
    t >>= 1;
    for (unsigned long i = 0; i < m; ++i) {
      unsigned long j1 = 2 * i * t;
      unsigned long j2 = j1 + t - 1;
      tools::element_type_T S = phi_list[m + i];
      for (unsigned long j = j1; j <= j2; ++j) {
        long U = result[j];
        long V = result[j + t] * S;
        result[j] = tools::mod_reduce(U + V, q);
        result[j + t] = tools::mod_reduce(U - V, q);
      }
    }
  }
  return result;
}

std::ostream &momoko::base::operator<<(std::ostream &os,
                                       const ideal_lattice_element &element) {
  for (auto i : element.factors) {
    os << i << " ";
  }
  return os;
}

momoko::base::ideal_lattice_element
momoko::base::operator*(ideal_lattice_element a, ideal_lattice_element &b) {
  a *= b;
  return a;
}
momoko::base::ideal_lattice_element
momoko::base::operator*(tools::element_type_T b, ideal_lattice_element a) {
  a *= b;
  return a;
}
momoko::base::ideal_lattice_element
momoko::base::operator*(ideal_lattice_element a, tools::element_type_T b) {
  a *= b;
  return a;
}
momoko::base::ideal_lattice_element
momoko::base::operator+(ideal_lattice_element a,
                        const ideal_lattice_element &b) {
  a += b;
  return a;
}

momoko::base::ideal_lattice_element
momoko::base::operator-(ideal_lattice_element a) {
  for (size_t i = 0; i < a.factors.size(); ++i) {
    a.factors[i] = tools::mod_reduce(-a.factors[i], a.q);
  }
  a.NTT_cache.reset();
  return a;
}
momoko::base::ideal_lattice_element
momoko::base::operator-(ideal_lattice_element a,
                        const ideal_lattice_element &b) {
  a -= b;
  return a;
}

momoko::base::ideal_lattice_element
momoko::base::SPM_product(ideal_lattice_element a,
                          const ideal_lattice_element &b) {
  std::vector<long> result(a.n, 0);
  for (size_t i = 0; i < a.factors.size(); ++i) {
    for (size_t j = 0; j < b.factors.size(); ++j) {
      if ((i + j) < a.n) {
        result[i + j] += (a.factors[i] * b.factors[j]);
        result[i + j] = tools::mod_reduce(result[i + j], a.q);
      } else {
        result[(i + j) % a.n] +=
            tools::mod_reduce(-a.factors[i] * b.factors[j], a.q);
        result[(i + j) % a.n] = tools::mod_reduce(result[(i + j) % a.n], a.q);
      }
    }
  }
  return a.latt.make_element(result);
}
