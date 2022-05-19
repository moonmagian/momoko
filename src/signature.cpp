#include "signature.h"

momoko::base::ideal_lattice_element momoko::pks::signature::sample_uniform() {
  std::vector<long> factors(latt.getN());
  for (size_t i = 0; i < factors.size(); ++i) {
    factors[i] = distn(rng);
  }
  return latt.make_element(factors);
}

momoko::base::ideal_lattice_element
momoko::pks::signature::sample_short_uniform() {
  std::vector<long> factors(latt.getN());
  for (size_t i = 0; i < factors.size(); ++i) {
    factors[i] = dist1(rng);
  }
  return latt.make_element(factors);
}

momoko::pks::signature::signature(base::ideal_lattice &latt,
                                  gaussian::gaussian_dist_sampler &sampler)
    : pksystem{latt}, sampler{sampler},
      distn{-static_cast<long>(latt.getQ() - 1) / 2,
            static_cast<long>(latt.getQ() - 1) / 2},
      dist1{-1, 1} {
  s1.emplace(sample_short_uniform());
  s2.emplace(sample_short_uniform());
  a.emplace(sample_uniform());
  t.emplace(a.value() * s1.value() + s2.value());
}

void momoko::pks::signature::export_sk(std::ostream &os) {}

void momoko::pks::signature::import_sk(std::istream &is) {}

void momoko::pks::signature::export_pk(std::ostream &os) {}

void momoko::pks::signature::import_pk(std::istream &is) {}

bool momoko::pks::signature::pk_ready() { return true; }

bool momoko::pks::signature::sk_ready() { return true; }

std::tuple<momoko::base::ideal_lattice_element,
           momoko::base::ideal_lattice_element,
           momoko::base::ideal_lattice_element>
momoko::pks::signature::sign_latt_element(
    const base::ideal_lattice_element &message) {}

bool momoko::pks::signature::verify_latt_element(
    base::ideal_lattice_element &message,
    std::pair<base::ideal_lattice_element, base::ideal_lattice_element> &sign) {
  //  auto left(neg_a.value() * sign.second + b.value() * sign.first);
  //  auto right(neg_b.value() * message);
  //  for (size_t i = 0; i < latt.getN(); ++i) {
  //    if (tools::mod_reduce(left.get_factor(i), 3) !=
  //        tools::mod_reduce(right.get_factor(i), 3)) {
  //      std::cout << left.get_factor(i) << " " << right.get_factor(i)
  //                << std::endl;
  //      return false;
  //    }
  //  }
  return true;
}
