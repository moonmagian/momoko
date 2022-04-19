#include "CDT_sampler.hpp"
#include "tools.hpp"
#include <algorithm>
#include <fstream>
namespace momoko::gaussian {

long CDT_sampler::sample_halfside_gaussian() {
  // CDT doesn't use halfside rejection algorithm.
  // It can sample the whole gaussian directly.
  return 0;
}

void CDT_sampler::read_param_from_stream(std::istream &s) {
  char read_header[4];
  s.read(read_header, 4);
  if (!std::equal(std::begin(header), std::end(header), std::begin(read_header),
                  std::end(read_header))) {
    throw std::runtime_error("Bad CDT sampler file.");
  }
  s.read(reinterpret_cast<char *>(&SD), sizeof(SD));
  s.read(reinterpret_cast<char *>(&tail), sizeof(tail));
  size_t N{static_cast<size_t>(SD * tail) + 1};
  CDT.resize(N);
  for (size_t i = 0; i < N; ++i) {
    s.read(reinterpret_cast<char *>(&CDT[i]), sizeof(CDT[i]));
  }
}

CDT_sampler::CDT_sampler(double SD, double tail, base::ideal_lattice &latt)
    : gaussian_dist_sampler(latt), SD(SD), tail(tail),
      CDT(static_cast<size_t>(SD * tail) + 1) {
  size_t N{static_cast<size_t>(SD * tail) + 1};
  CDT[0] = 0;
  CDT[1] = tools::discrete_gaussian(0, SD);
  for (size_t i = 2; i < N - 1; ++i) {
    CDT[i] = CDT[i - 1] + tools::discrete_gaussian(i - 1, SD) * 2;
  }
  CDT[N - 1] = 1;
}

CDT_sampler::CDT_sampler(std::istream &s, base::ideal_lattice &latt)
    : gaussian_dist_sampler(latt) {
  read_param_from_stream(s);
}

CDT_sampler::CDT_sampler(const std::string &path, base::ideal_lattice &latt)
    : gaussian_dist_sampler(latt) {
  std::ifstream ifs(path);
  read_param_from_stream(ifs);
}

void CDT_sampler::export_param(std::ostream &s) {
  s.write(header, 4);
  s.write(reinterpret_cast<char *>(&SD), sizeof(SD));
  s.write(reinterpret_cast<char *>(&tail), sizeof(tail));
  for (auto x : CDT) {
    s.write(reinterpret_cast<char *>(&x), sizeof(x));
  }
}

void CDT_sampler::export_param(const std::string &path) {
  std::ofstream ofs(path);
  export_param(ofs);
}

double CDT_sampler::get_std_var() { return SD; }

double CDT_sampler::get_tail() { return tail; }

long CDT_sampler::sample_gaussian() {
  double p{real_dist(rng)};
  auto sample_bound = std::upper_bound(CDT.begin(), CDT.end(), p);
  long result = (sample_bound - CDT.begin()) - 1;
  if (result != 0 && sign_dist(rng)) {
    result = -result;
  }
  return result;
}

bool CDT_sampler::operator==(const CDT_sampler &other) const {
  return (SD == other.SD && tail == other.tail && CDT == other.CDT);
}

} // namespace momoko::gaussian
