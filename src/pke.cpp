#include "pke.h"
#include <iostream>

momoko::pks::pke::pke(base::ideal_lattice &latt,
                      gaussian::gaussian_dist_sampler &sampler)
    : pksystem(latt), sampler{sampler} {}

void momoko::pks::pke::generate_keys() {

  s.emplace(sampler.sample_lattice_element());
  //  a.emplace(latt.make_element({1}));
  a.emplace(sampler.sample_lattice_element());
  e.emplace(sampler.sample_lattice_element());
  b.emplace(a.value() * s.value() + e.value());
}

void momoko::pks::pke::export_sk(std::ostream &os) {
  if (!sk_ready()) {
    throw std::runtime_error("Private key is not ready.");
  }
  os.write(header_sk, sizeof(header_sk));
  a->export_to_stream(os);
  b->export_to_stream(os);
  s->export_to_stream(os);
  e->export_to_stream(os);
}

void momoko::pks::pke::import_sk(std::istream &is) {
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

void momoko::pks::pke::export_pk(std::ostream &os) {
  if (!pk_ready()) {
    throw std::runtime_error("Public key is not ready.");
  }
  os.write(header_sk, sizeof(header_sk));
  a->export_to_stream(os);
  b->export_to_stream(os);
}
void momoko::pks::pke::import_pk(std::istream &is) {
  char read_header[4];
  is.read(read_header, 4);
  if (!std::equal(std::begin(read_header), std::end(read_header),
                  std::begin(header_pk), std::end(header_pk))) {

    throw std::runtime_error("Invalid public key stream.");
  }
  a.emplace(latt.import_element(is));
  b.emplace(latt.import_element(is));
}

bool momoko::pks::pke::sk_ready() {
  return s.has_value() && e.has_value() && a.has_value() && b.has_value();
}

bool momoko::pks::pke::pk_ready() { return a.has_value() && b.has_value(); }

std::pair<momoko::base::ideal_lattice_element,
          momoko::base::ideal_lattice_element>
momoko::pks::pke::encrypt_latt_element(base::ideal_lattice_element message) {
  if (!pk_ready()) {
    throw std::runtime_error("Public key is not ready.");
  }
  base::ideal_lattice_element r{sampler.sample_lattice_element()};
  base::ideal_lattice_element e1{sampler.sample_lattice_element()};
  base::ideal_lattice_element e2{sampler.sample_lattice_element()};
  message *= latt.getQ() / 2;
  return std::make_pair(a.value() * r + e1, b.value() * r + e2 + message);
}

momoko::base::ideal_lattice_element momoko::pks::pke::decrypt_latt_element(
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
