#include "CDT_sampler.hpp"
#include "bernoulli_sampler.hpp"
#include "ideal_lattice.hpp"
#include "tools.hpp"
#include "bit_matrix.hpp"
#include "knuth_yao_sampler.hpp"
#include "background_sampler.h"
#include "background_sampler_lockfree.h"
#include "pke.h"
#include "pke_blwe.h"
#include "signature.h"
#include "benchmark_counter.h"
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <chrono>
#include <thread>
#include <fstream>
using namespace std;
void NTT_test() {
  momoko::base::ideal_lattice latt{512, 12289};
  std::vector<long> factors(512, 0);
  momoko::gaussian::CDT_sampler CDT{1.698644, 9.42, latt};
  for (int i = 0; i < 1000; ++i) {

    auto element1 = CDT.sample_lattice_element();
    auto element2 = CDT.sample_lattice_element();
    auto r1 = element1 * element2;
    auto r2 = momoko::base::SPM_product(element1, element2);
    if (!(r1 == r2)) {
      throw std::runtime_error("Error while testing NTT.");
    }
  }
  std::cout << "[PASS] NTT Correctness." << std::endl;
}
void sampler_test(momoko::gaussian::gaussian_dist_sampler &sampler) {
  int times = 1000000;
  //  int times = 1000;
  double mean = 0;
  double SD = 0;
  std::map<long, ulong> m;
  for (int i = 0; i < times; ++i) {
    long result = sampler.sample_gaussian();
    m[result] += 1;
  }

  for (const auto &pair : m) {
    std::cout << pair.first << "\t" << pair.second << std::endl;
    mean += static_cast<double>(pair.first) * pair.second / times;
  }
  for (const auto &pair : m) {
    double delta = static_cast<double>(pair.first) - mean;
    delta *= delta;
    delta *= pair.second;
    delta /= times;
    SD += delta;
  }
  SD = std::sqrt(SD);
  std::cout << "Mean: " << mean << std::endl;
  std::cout << "Desired mean: 0" << std::endl;
  std::cout << "SD: " << SD << std::endl;
  std::cout << "Desired SD: " << sampler.get_std_var() << std::endl;
}
void CDT_save_test(momoko::gaussian::CDT_sampler &sampler,
                   momoko::base::ideal_lattice &latt) {
  sampler.export_param("CDT.sampler");
  momoko::gaussian::CDT_sampler sampler2("CDT.sampler", latt);
  std::cout << (sampler2 == sampler) << std::endl;
}
void bit_matrix_test() {
  momoko::tools::bit_matrix mat;
  mat.push(0.2);
  mat.push(0.35);
  mat.push(0.33);
  mat.complete();
  std::cout << mat << std::endl;
}
void pke_benchmark() {
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::gaussian::CDT_sampler CDT{4.85, 9.42, latt};
  momoko::gaussian::background_sampler_lockfree<300000> bg{CDT, latt};
  momoko::pks::pke pke(latt, bg);
  pke.generate_keys();
  long result = start_benchmark(
      [&pke, &latt](std::stop_token token, std::promise<long> promise) {
        long result{0};
        while (!token.stop_requested()) {
          pke.encrypt_latt_element(latt.make_element({1, 1, 1, 0, 1, 1, 0, 1}));
          //          pke.decrypt_latt_element(
          //              std::make_pair(latt.make_element({1, 1, 1, 0, 1, 1, 0,
          //              1}),
          //                             latt.make_element({1, 1, 0})));
          result += 1;
        }
        promise.set_value(result);
      },
      10);
  std::cout << result << std::endl;
}
void pke_test_and_benchmark(bool warmup = true) {
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::gaussian::CDT_sampler CDT{4.85, 9.42, latt};
  //  momoko::gaussian::knuth_yao_sampler knuth_yao{1.698644, 9.42, latt};
  //  momoko::gaussian::background_sampler_lockfree<50000> bg{CDT, latt};
  momoko::gaussian::background_sampler_lockfree<300000> bg{CDT, latt};
  momoko::pks::pke pke(latt, bg);
  //  momoko::pks::pke pke2(latt, CDT);
  momoko::pks::pke_blwe pke2(latt);
  pke.generate_keys();
  pke2.generate_keys();
  if (warmup) {
    this_thread::sleep_for(std::chrono::seconds(10));
  }
  auto start = chrono::high_resolution_clock::now();
  for (long i = 0; i < 300; ++i) {
    auto r =
        pke.encrypt_latt_element(latt.make_element({1, 1, 1, 0, 1, 1, 0, 1}));
  }
  auto end = chrono::high_resolution_clock::now();
  std::cout << "[background sampling] Randomly encrypt 150 kbits(μs):"
            << chrono::duration_cast<chrono::microseconds>(end - start).count()
            << std::endl;
  start = chrono::high_resolution_clock::now();
  for (long i = 0; i < 300; ++i) {
    auto r =
        pke2.encrypt_latt_element(latt.make_element({1, 1, 1, 0, 1, 1, 0, 1}));
    //    auto result = pke2.decrypt_latt_element(r);
    //    std::cout << result << std::endl;
  }
  end = chrono::high_resolution_clock::now();
  std::cout << "[normal samping] Randomly encrypt 150 kbits(μs):"
            << chrono::duration_cast<chrono::microseconds>(end - start).count()
            << std::endl;
}
void NTT_benchmark() {
  momoko::base::ideal_lattice latt_NTT_cache{512, 12289, true};
  momoko::base::ideal_lattice latt{512, 12289, false};

  momoko::gaussian::CDT_sampler CDT_NTT_cache{4.85, 9.42, latt_NTT_cache};
  momoko::pks::pke pke_NTT_cache(latt_NTT_cache, CDT_NTT_cache);

  momoko::gaussian::CDT_sampler CDT{4.85, 9.42, latt};
  momoko::pks::pke pke(latt, CDT);

  pke.generate_keys();
  pke_NTT_cache.generate_keys();
  auto start = chrono::high_resolution_clock::now();
  for (long i = 0; i < 300; ++i) {
    auto r =
        pke.encrypt_latt_element(latt.make_element({1, 1, 1, 0, 1, 1, 0, 1}));
  }
  auto end = chrono::high_resolution_clock::now();
  std::cout << "[No NTT cache and lazy modular] Randomly encrypt 150 kbits(μs):"
            << chrono::duration_cast<chrono::microseconds>(end - start).count()
            << std::endl;
  start = chrono::high_resolution_clock::now();
  for (long i = 0; i < 300; ++i) {
    auto r = pke_NTT_cache.encrypt_latt_element(
        latt.make_element({1, 1, 1, 0, 1, 1, 0, 1}));
  }
  end = chrono::high_resolution_clock::now();
  std::cout << "[NTT cache + lazy modular] Randomly encrypt 150 kbits(μs):"
            << chrono::duration_cast<chrono::microseconds>(end - start).count()
            << std::endl;
}
void sign_test() {
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::gaussian::CDT_sampler CDT{4.85, 9.42, latt};
  momoko::pks::signature sign{latt, CDT};
  auto message(latt.make_element({1, 1, 1, 0, 1, 1, 0, 1}));
  auto sig = sign.sign_latt_element(message);
  //  std::cout << sign.verify_latt_element(message, sig) << std::endl;
}
void latt_export_test() {
  std::stringstream ss;
  momoko::base::ideal_lattice latt{512, 12289};
  latt.export_to_stream(ss, true);
  momoko::base::ideal_lattice latt2{ss};
  std::cout << latt2.getN() << std::endl;
  std::cout << latt2.getQ() << std::endl;
  momoko::tools::_print_vec(latt.get_psi_list());
  std::cout << std::endl << std::endl;
  momoko::tools::_print_vec(latt2.get_psi_list());
  std::cout << endl;
  auto element1 = latt.make_element({1, 1, 0, 0, 11, 4, 5, 14});
  std::cout << element1 << std::endl;
  element1.export_to_stream(ss);
  auto element2 = latt.import_element(ss);
  std::cout << element2 << std::endl;
}
int main() {
  //  latt_export_test();
  //  NTT_test();
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::gaussian::CDT_sampler CDT{1.698644, 9.42, latt};
  //  auto element = latt.make_element({1, 1, 4, 5, 1, 4});
  //  std::ofstream ofs{"out.element"};
  //  std::ofstream ofs_latt{"out.latt"};
  //  element.export_to_stream(ofs);
  //  latt.export_to_stream(ofs_latt, true);
  momoko::gaussian::background_sampler_lockfree<10000> bg{CDT, latt};
  sampler_test(bg);
  //  momoko::gaussian::bernoulli_sampler bernoulli{2, latt};
  //  momoko::gaussian::knuth_yao_sampler knuth_yao{1.698644, 9.42, latt};
  //  std::cout << "Sampler shape test:" << std::endl;
  //  sampler_test(CDT);
  //  std::cout << std::endl;
  //  sampler_test(bernoulli);
  //  sampler_test(knuth_yao);
  //  CDT_save_test(CDT, latt);
  //  bit_matrix_test();
  //  knuth_yao._check_carry();
  //  std::cout << "background sampling test (with warm up):" << std::endl;
  //  pke_test_and_benchmark(true);
  //  std::cout << std::endl;
  //  std::cout << "background sampling test (without warm up):" << std::endl;
  //  pke_test_and_benchmark(false);
  //  pke_benchmark();
  //  std::cout << std::endl;
  //  std::cout << "NTT cache and lazy modular test:" << std::endl;
  //  NTT_benchmark();
  //  sign_test();
  return 0;
}
