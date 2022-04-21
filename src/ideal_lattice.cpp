#include "ideal_lattice.hpp"
#include "tools.hpp"
#include <type_traits>
namespace momoko {
namespace base {
ideal_lattice::ideal_lattice(ulong n, ulong q, bool use_NTT_cache)
    : n{n}, q(q), n_inv(tools::mod_inv(n, q)), use_NTT_cache{use_NTT_cache} {}

ideal_lattice::ideal_lattice(std::istream &is, bool use_NTT_cache)
    : use_NTT_cache{use_NTT_cache} {
  char read_header[4];
  is.read(read_header, 4);
  bool has_cache{false};
  if (std::equal(std::begin(read_header), std::end(read_header),
                 std::begin(header_cache), std::end(header_cache))) {
    has_cache = true;
  } else if (std::equal(std::begin(read_header), std::end(read_header),
                        std::begin(header_no_cache),
                        std::end(header_no_cache))) {
    has_cache = false;
  } else {
    throw std::runtime_error("Invalid lattice stream.");
  }
  is.read(reinterpret_cast<char *>(&n), sizeof(n));
  is.read(reinterpret_cast<char *>(&q), sizeof(q));
  is.read(reinterpret_cast<char *>(&n_inv), sizeof(n_inv));
  if (has_cache) {
    psi.emplace();
    psi_list.emplace(n);
    psi_inv_list.emplace(n);
    is.read(reinterpret_cast<char *>(&psi.value()), sizeof(psi.value()));
    for (size_t i = 0; i < n; ++i) {
      is.read(reinterpret_cast<char *>(&(psi_list.value()[i])),
              sizeof(psi.value()));
      is.read(reinterpret_cast<char *>(&(psi_inv_list.value()[i])),
              sizeof(psi.value()));
    }
  }
}

tools::element_type_T ideal_lattice::get_2nth_root_of_unity() {
  // Use cached value.
  if (psi) {
    return psi.value();
  }
  auto x = tools::fatorize(2 * n);
  ulong res{};
  bool ok{};
  for (res = 2; res < q; ++res) {
    if (tools::mod_pow(res, 2 * n, q) != 1) {
      continue;
    }
    ok = true;
    for (auto i : x) {
      if (tools::mod_pow(res, 2 * n / i, q) == 1) {
        ok = false;
        break;
      }
    }
    if (ok) {
      break;
    }
  }
  // Can't find a root of unity.
  if (!ok) {
    psi = 0;
    return 0;
  }
  // Cache the result.
  res = tools::mod_reduce(res, q);
  psi = res;
  return res;
}

ideal_lattice_element ideal_lattice::make_element() {
  return ideal_lattice_element{*this};
}

ideal_lattice_element ideal_lattice::import_element(std::istream &is) {
  std::vector<long> result_factors;
  std::remove_const<decltype(ideal_lattice_element::n)>::type result_n;
  std::remove_const<decltype(ideal_lattice_element::q)>::type result_q;
  is.read(reinterpret_cast<char *>(&result_n),
          sizeof(ideal_lattice_element::n));
  is.read(reinterpret_cast<char *>(&result_q),
          sizeof(ideal_lattice_element::q));
  if (result_n != n || result_q != q) {
    throw std::runtime_error("Unmatched lattice and lattice element.");
  }
  size_t factor_count;
  is.read(reinterpret_cast<char *>(&factor_count), sizeof(factor_count));
  result_factors.resize(factor_count);
  for (size_t i = 0; i < factor_count; ++i) {
    is.read(reinterpret_cast<char *>(&result_factors[i]),
            sizeof(result_factors[i]));
  }
  return ideal_lattice_element{*this, result_factors};
}

ideal_lattice_element
ideal_lattice::make_element(const std::vector<long> &factor) {
  return ideal_lattice_element{*this, factor, use_NTT_cache};
}

ideal_lattice_element
ideal_lattice::make_element_from_NTT(const std::vector<long> &ntt) {
  return ideal_lattice_element{*this, INTT_GS(ntt), use_NTT_cache};
}

std::vector<long> ideal_lattice::INTT_GS(std::vector<long> a) {
  auto inv_list = get_psi_inv_list();
  a.resize(n, 0);
  ulong t = 1;
  for (ulong m = n; m > 1; m >>= 1) {
    ulong j1 = 0;
    ulong h = m >> 1;
    for (ulong i = 0; i < h; ++i) {
      ulong j2 = j1 + t - 1;
      long S = inv_list[h + i];
      for (ulong j = j1; j <= j2; ++j) {
        long U = a[j];
        long V = a[j + t];
        //                a[j] = tools::mod_reduce(U + V, q);
        //                a[j + t] = tools::mod_reduce((U - V) * S, q);
        a[j] = tools::mod_reduce_lazy(U + V, q);
        a[j + t] = tools::mod_reduce_lazy((U - V) * S, q);
      }
      j1 = j1 + (t << 1);
    }
    t <<= 1;
  }
  for (ulong j = 0; j < n; ++j) {
    a[j] = tools::mod_reduce(a[j] * n_inv, q);
  }
  return a;
}

void ideal_lattice::export_to_stream(std::ostream &os, bool include_psi_cache) {
  if (include_psi_cache) {
    os.write(header_cache, 4);
  } else {
    os.write(header_no_cache, 4);
  }
  os.write(reinterpret_cast<const char *>(&n), sizeof(n));
  os.write(reinterpret_cast<const char *>(&q), sizeof(q));
  os.write(reinterpret_cast<const char *>(&n_inv), sizeof(n_inv));
  if (include_psi_cache) {
    get_psi_list();
    get_psi_inv_list();
    os.write(reinterpret_cast<const char *>(&psi.value()), sizeof(psi.value()));
    for (size_t i = 0; i < n; ++i) {
      os.write(reinterpret_cast<const char *>(&(psi_list.value()[i])),
               sizeof(psi.value()));
      os.write(reinterpret_cast<const char *>(&(psi_inv_list.value()[i])),
               sizeof(psi.value()));
    }
  }
}

const std::vector<tools::element_type_T> &ideal_lattice::get_psi_list() {
  // Use cached result.
  if (psi_list) {
    return psi_list.value();
  }
  psi_list.emplace(n);
  tools::element_type_T psi_base{get_2nth_root_of_unity()};
  psi_list.value()[0] = 1;
  // Calculate the list.
  for (size_t i = 1; i < psi_list->size(); ++i) {
    psi_list.value()[i] =
        tools::mod_reduce(psi_list.value()[i - 1] * psi_base, q);
  }
  // Use bit reversed order (as required in CT and GS butterfly).
  tools::bitrevorder(psi_list.value());
  return psi_list.value();
}

const std::vector<tools::element_type_T> &ideal_lattice::get_psi_inv_list() {
  // Use cached result.
  if (psi_inv_list) {
    return psi_inv_list.value();
  }
  psi_inv_list.emplace(n);
  tools::element_type_T psi_inv_base{
      tools::mod_inv(get_2nth_root_of_unity(), q)};
  psi_inv_list.value()[0] = 1;
  // Calculate the list.
  for (size_t i = 1; i < psi_inv_list->size(); ++i) {
    psi_inv_list.value()[i] =
        tools::mod_reduce(psi_inv_list.value()[i - 1] * psi_inv_base, q);
  }
  // Use bit reversed order (as required in CT and GS butterfly).
  tools::bitrevorder(psi_inv_list.value());
  return psi_inv_list.value();
}

unsigned long ideal_lattice::getQ() const { return q; }

bool ideal_lattice::operator==(const ideal_lattice &other) const {
  return n == other.n && q == other.q;
}
unsigned long ideal_lattice::getN() const { return n; }

} // namespace base
} // namespace momoko
