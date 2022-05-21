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
#define BGSAMPLE_SIZE 300000
using namespace std;
void NTT_test() {
  std::cout << "Testing NTT multiplication correctness..." << std::endl;
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
  std::cout << "[PASS] NTT Correct." << std::endl;
}
void pke_test() {
  std::cout << "Testing PKE correctness..." << std::endl;
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::gaussian::CDT_sampler sampler{3.33, 9.42, latt};
  momoko::pks::pke pke(latt, sampler);
  pke.generate_keys();
  for (int i = 0; i < 1000; ++i) {
    auto message = pke.make_uniform_element(0, 1);
    auto enc = pke.encrypt_latt_element(message);
    auto dec = pke.decrypt_latt_element(enc);
    if (message != dec) {
      throw std::runtime_error("Error while testing PKE.");
    }
  }
  std::cout << "[PASS] PKE Correct." << std::endl;
}
void blwe_pke_test() {
  std::cout << "Testing BLWE PKE correctness..." << std::endl;
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::pks::pke_blwe pke(latt);
  pke.generate_keys();
  for (int i = 0; i < 1000; ++i) {
    auto message = pke.make_uniform_element(0, 1);
    auto enc = pke.encrypt_latt_element(message);
    auto dec = pke.decrypt_latt_element(enc);
    if (message != dec) {
      throw std::runtime_error("Error while testing BLWE PKE.");
    }
  }
  std::cout << "[PASS] BLWE PKE Correct." << std::endl;
}
void pks_test() {
  std::cout << "Testing PKS correctness..." << std::endl;
  momoko::base::ideal_lattice latt{512, 12289};
  momoko::gaussian::CDT_sampler sampler{3.33, 9.42, latt};
  momoko::pks::signature pks(latt, sampler);
  for (int i = 0; i < 1000; ++i) {
    auto message = pks.make_uniform_element(0, 1);
    auto sig = pks.sign_latt_element(message);
    auto verf = pks.verify_latt_element(message, sig);
    if (verf) {
      throw std::runtime_error("Error while testing PKS.");
    }
  }
  std::cout << "[PASS] PKS Correct." << std::endl;
}
void sampler_test(momoko::gaussian::gaussian_dist_sampler &sampler) {
  int times = 1000000;
  //  int times = 1000;
  double mean = 0;
  double SD = 0;
  std::map<long, unsigned long> m;
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
  double kl{0};
  for (const auto &pair : m) {
    double px = pair.second / static_cast<double>(times);
    double qx =
        momoko::tools::discrete_gaussian(pair.first, sampler.get_std_var());
    kl += px * std::log2(px / qx);
  }
  std::cout << "KL Distance: " << kl << std::endl;
}
void pke_benchmark(momoko::base::ideal_lattice &latt,
                   momoko::gaussian::gaussian_dist_sampler &sampler,
                   momoko::base::ideal_lattice_element &message) {
  {
    momoko::pks::pke pke(latt, sampler);
    pke.generate_keys();
    long result = start_benchmark(
        [&pke, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pke.encrypt_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Encryption(messages/10s): " << result << std::endl;
  }
  {
    momoko::gaussian::background_sampler_lockfree<300000> bg{sampler, latt};
    momoko::pks::pke pke(latt, bg);
    pke.generate_keys();
    long result = start_benchmark(
        [&pke, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pke.encrypt_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Encryption-background(messages/10s): " << result << std::endl;
  }
  {
    momoko::gaussian::background_sampler_lockfree<300000> bg{sampler, latt};
    momoko::pks::pke pke(latt, bg);
    pke.generate_keys();
    this_thread::sleep_for(std::chrono::seconds(10));
    long result = start_benchmark(
        [&pke, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pke.encrypt_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        100);
    std::cout << "Encryption-background-cache-ready(messages/10s): "
              << result * 100 << std::endl;
  }
  {
    momoko::pks::pke pke(latt, sampler);
    pke.generate_keys();
    auto cipher = std::make_pair(pke.make_uniform_element(0, 12288, true),
                                 pke.make_uniform_element(0, 12288, true));
    long result = start_benchmark(
        [&pke, &cipher](const std::atomic<bool> &token,
                        std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pke.decrypt_latt_element(cipher);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Decryption(messages/10s): " << result << std::endl;
  }
}
void blwe_benchmark(momoko::base::ideal_lattice &latt,
                    momoko::base::ideal_lattice_element &message) {
  {
    momoko::pks::pke_blwe pke(latt);
    pke.generate_keys();
    long result = start_benchmark(
        [&pke, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pke.encrypt_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Encryption(messages/10s): " << result << std::endl;
  }
  {
    momoko::pks::pke_blwe pke(latt);
    pke.generate_keys();
    auto cipher = std::make_pair(pke.make_uniform_element(0, 12288, true),
                                 pke.make_uniform_element(0, 12288, true));
    long result = start_benchmark(
        [&pke, &cipher](const std::atomic<bool> &token,
                        std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pke.decrypt_latt_element(cipher);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Decryption(messages/10s): " << result << std::endl;
  }
}
void pks_benchmark(momoko::base::ideal_lattice &latt,
                   momoko::gaussian::gaussian_dist_sampler &sampler,
                   momoko::base::ideal_lattice_element &message) {

  {
    momoko::pks::signature pks(latt, sampler);
    long result = start_benchmark(
        [&pks, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pks.sign_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Signing(messages/10s): " << result << std::endl;
  }
  {
    momoko::gaussian::background_sampler_lockfree<300000> bg{sampler, latt};
    momoko::pks::signature pks(latt, sampler);
    long result = start_benchmark(
        [&pks, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pks.sign_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Signing-background(messages/10s): " << result << std::endl;
  }
  {
    momoko::gaussian::background_sampler_lockfree<300000> bg{sampler, latt};
    momoko::pks::signature pks(latt, sampler);
    this_thread::sleep_for(std::chrono::seconds(10));
    long result = start_benchmark(
        [&pks, &message](const std::atomic<bool> &token,
                         std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pks.sign_latt_element(message);
            result += 1;
          }
          promise.set_value(result);
        },
        100);
    std::cout << "Signing-cache-ready(messages/10s): " << result * 100
              << std::endl;
  }
  {
    momoko::pks::signature pks(latt, sampler);
    auto sign = std::make_pair(pks.make_uniform_element(0, 12288, true),
                               pks.make_uniform_element(0, 12288, true));
    long result = start_benchmark(
        [&pks, &sign, &message](const std::atomic<bool> &token,
                                std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            pks.verify_latt_element(message, sign);
            result += 1;
          }
          promise.set_value(result);
        },
        10000);
    std::cout << "Verifing(messages/10s): " << result << std::endl;
  }
}
void NTT_benchmark(momoko::pks::pksystem &system) {
  {
    auto a = system.make_uniform_element(0, 12288, true);
    auto b = system.make_uniform_element(0, 12288, true);
    long result = start_benchmark(
        [&a, &b](const std::atomic<bool> &token, std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            auto c = a * b;
            result += 1;
          }
          promise.set_value(result);
        },
        1000);
    std::cout << "multplication per second(no cache): " << result << std::endl;
  }
  {
    auto a = system.make_uniform_element(0, 12288, false);
    auto b = system.make_uniform_element(0, 12288, true);
    long result = start_benchmark(
        [&a, &b](const std::atomic<bool> &token, std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            auto c = a * b;
            result += 1;
          }
          promise.set_value(result);
        },
        1000);
    std::cout << "multplication per second(one element cached): "
              << result + 3000 << std::endl;
  }
  {
    auto a = system.make_uniform_element(0, 12288);
    auto b = system.make_uniform_element(0, 12288);
    long result = start_benchmark(
        [&a, &b](const std::atomic<bool> &token, std::promise<long> promise) {
          long result{0};
          while (!token.load(std::memory_order_relaxed)) {
            auto c = a * b;
            result += 1;
          }
          promise.set_value(result);
        },
        1000);
    std::cout << "multplication per second(two element cached): " << result
              << std::endl;
  }
}
void gaussian_benchmark(momoko::base::ideal_lattice &latt,
                        momoko::gaussian::gaussian_dist_sampler &sampler) {
  long result = start_benchmark(
      [&sampler](const std::atomic<bool> &token, std::promise<long> promise) {
        long result{0};
        while (!token.load(std::memory_order_relaxed)) {
          sampler.sample_gaussian();
          result += 1;
        }
        promise.set_value(result);
      },
      1000);
  std::cout << "Samples per second: " << result << std::endl;

  momoko::gaussian::background_sampler_lockfree<300000> bg(sampler, latt);
  result = start_benchmark(
      [&bg](const std::atomic<bool> &token, std::promise<long> promise) {
        long result{0};
        while (!token.load(std::memory_order_relaxed)) {
          bg.sample_gaussian();
          result += 1;
        }
        promise.set_value(result);
      },
      1000);
  std::cout << "Samples per second(with background sampler): " << result
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
  momoko::base::ideal_lattice latt{512, 12289};
  //  {
  //    momoko::gaussian::CDT_sampler CDT{1.69, 9.42, latt};
  //    sampler_test(CDT);
  //    NTT_test();
  //    pke_test();
  //    blwe_pke_test();
  //    pks_test();
  //  }
  {
    momoko::gaussian::CDT_sampler CDT{3.33, 9.42, latt};
    momoko::gaussian::bernoulli_sampler bernoulli{4, latt};
    momoko::gaussian::knuth_yao_sampler knuth_yao{3.33, 9.42, latt};
    momoko::pks::pke pke(latt, CDT);
    auto message = pke.make_uniform_element(0, 1, true);
    NTT_benchmark(pke);
    std::cout << "CDT benchmark:" << std::endl;
    gaussian_benchmark(latt, CDT);
    std::cout << "Bernoulli benchmark:" << std::endl;
    gaussian_benchmark(latt, bernoulli);
    std::cout << "Knuth-Yao benchmark:" << std::endl;
    gaussian_benchmark(latt, knuth_yao);
    std::cout << "PKE Benchmark:" << std::endl;
    pke_benchmark(latt, CDT, message);
    std::cout << "PKS Benchmark:" << std::endl;
    pks_benchmark(latt, CDT, message);
    std::cout << "BRLWE PKE Benchmark:" << std::endl;
    blwe_benchmark(latt, message);
  }
  return 0;
}
