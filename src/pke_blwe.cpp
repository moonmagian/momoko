#include "pke_blwe.h"
#include <iostream>

momoko::base::ideal_lattice_element momoko::pks::pke_blwe::sample_uniform(
    std::uniform_int_distribution<long> &dist) {
  std::vector<long> factors(latt.getN());
  for (size_t i = 0; i < factors.size(); ++i) {
    factors[i] = dist(rng);
  }
  return latt.make_element(factors);
}

momoko::pks::pke_blwe::pke_blwe(base::ideal_lattice &latt)
    : pksystem(latt), k_uniform_sampler{-1, 1},
      q_uniform_sampler{-static_cast<long>(latt.getQ() - 1) / 2,
                        static_cast<long>(latt.getQ() - 1) / 2} {}

void momoko::pks::pke_blwe::generate_keys() {

  s.emplace(sample_uniform(k_uniform_sampler));
  //  a.emplace(latt.make_element({1}));
  a.emplace(sample_uniform(q_uniform_sampler));
  e.emplace(sample_uniform(k_uniform_sampler));
  b.emplace(a.value() * s.value() + e.value());
}

void momoko::pks::pke_blwe::export_sk(std::ostream &os) {
  if (!sk_ready()) {
    throw std::runtime_error("Private key is not ready.");
  }
  os.write(header_sk, sizeof(header_sk));
  a->export_to_stream(os);
  b->export_to_stream(os);
  s->export_to_stream(os);
  e->export_to_stream(os);
}

void momoko::pks::pke_blwe::import_sk(std::istream &is) {
  char read_header[4];
  is.read(read_header, 4);
  if (!std::equal(std::begin(read_header), std::end(read_header),
                  std::begin(header_sk), std::end(header_sk))) {

    throw std::runtime_error("Invalid private key stream.");
  }
  a.emplace(latt.import_element(is));
  b.emplace(latt.import_element(is));
  s.emplace(latt.import_element(is));
  e.emplace(latt.import_element(is));
}

void momoko::pks::pke_blwe::export_pk(std::ostream &os) {
  if (!pk_ready()) {
    throw std::runtime_error("Public key is not ready.");
  }
  os.write(header_sk, sizeof(header_sk));
  a->export_to_stream(os);
  b->export_to_stream(os);
}
void momoko::pks::pke_blwe::import_pk(std::istream &is) {
  char read_header[4];
  is.read(read_header, 4);
  if (!std::equal(std::begin(read_header), std::end(read_header),
                  std::begin(header_pk), std::end(header_pk))) {

    throw std::runtime_error("Invalid public key stream.");
  }
  a.emplace(latt.import_element(is));
  b.emplace(latt.import_element(is));
}

bool momoko::pks::pke_blwe::sk_ready() {
  return s.has_value() && e.has_value() && a.has_value() && b.has_value();
}

bool momoko::pks::pke_blwe::pk_ready() {
  return a.has_value() && b.has_value();
}

std::pair<momoko::base::ideal_lattice_element,
          momoko::base::ideal_lattice_element>
momoko::pks::pke_blwe::encrypt_latt_element(
    base::ideal_lattice_element message) {
  if (!pk_ready()) {
    throw std::runtime_error("Public key is not ready.");
  }
  base::ideal_lattice_element r{sample_uniform(k_uniform_sampler)};
  base::ideal_lattice_element e1{sample_uniform(k_uniform_sampler)};
  base::ideal_lattice_element e2{sample_uniform(k_uniform_sampler)};
  message *= latt.getQ() / 2;
  return std::make_pair(a.value() * r + e1, b.value() * r + e2 + message);
}

momoko::base::ideal_lattice_element momoko::pks::pke_blwe::decrypt_latt_element(
    const std::pair<base::ideal_lattice_element, base::ideal_lattice_element>
        &message) {
  if (!sk_ready()) {
    throw std::runtime_error("Private key is not ready.");
  }
  auto result = message.second - (message.first * s.value());
  long n{static_cast<long>(latt.getN())};
  long q{static_cast<long>(latt.getQ())};
  // Begin from last to resize only once.
  for (int i = n - 1; i >= 0; --i) {
    long bit = result.get_factor(i);
    if (bit < 0) {
      bit = -bit;
    }
    if (bit > q / 4) {
      result.set_factor(i, 1);
    } else {
      result.set_factor(i, 0);
    }
  }
  return result;
}
