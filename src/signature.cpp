#include "signature.h"

momoko::base::ideal_lattice_element momoko::pks::signature::sample_uniform() {
  std::vector<long> factors(latt.getN());
  for (size_t i = 0; i < factors.size(); ++i) {
    factors[i] = dist(rng);
  }
  return latt.make_element(factors);
}

momoko::pks::signature::signature(base::ideal_lattice &latt,
                                  gaussian::gaussian_dist_sampler &sampler)
    : pksystem{latt}, dist{0, latt.getQ() - 1}, sampler{sampler} {
  s.emplace(sampler.sample_lattice_element());
  a.emplace(sampler.sample_lattice_element());
  e.emplace(sampler.sample_lattice_element());
  b.emplace(a.value() * s.value() + 3 * e.value());
  neg_a.emplace(-(a.value()));
  neg_b.emplace(-(b.value()));
}

void momoko::pks::signature::export_sk(std::ostream &os) {}

void momoko::pks::signature::import_sk(std::istream &is) {}

void momoko::pks::signature::export_pk(std::ostream &os) {}

void momoko::pks::signature::import_pk(std::istream &is) {}

bool momoko::pks::signature::pk_ready() {
  return a.has_value() && b.has_value() && neg_a.has_value();
}

bool momoko::pks::signature::sk_ready() {
  return s.has_value() && e.has_value();
}

std::pair<momoko::base::ideal_lattice_element,
          momoko::base::ideal_lattice_element>
momoko::pks::signature::sign_latt_element(
    const base::ideal_lattice_element &message) {
  auto v(sample_uniform());
  auto e1(sampler.sample_lattice_element());
  while (e1 == e.value()) {
    e1 = sampler.sample_lattice_element();
  }
  auto u((v + message) * s.value() + 3 * e1);
  return std::make_pair(v, u);
}

bool momoko::pks::signature::verify_latt_element(
    base::ideal_lattice_element &message,
    std::pair<base::ideal_lattice_element, base::ideal_lattice_element> &sign) {
  auto left(neg_a.value() * sign.second + b.value() * sign.first);
  auto right(neg_b.value() * message);
  for (size_t i = 0; i < latt.getN(); ++i) {
    if (left.get_factor(i) % 3 != right.get_factor(i) % 3) {
      return false;
    }
  }
  return true;
}
