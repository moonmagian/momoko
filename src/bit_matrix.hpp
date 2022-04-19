#ifndef BIT_MATRIX_HPP
#define BIT_MATRIX_HPP
#include <cinttypes>
#include <vector>
#include <array>
#include <iostream>
#include "tools.hpp"
namespace momoko::tools {
using std::size_t;
class bit_matrix {
  private:
  size_t prec;
  size_t row_vec_size;
  std::vector<std::vector<uint64_t>> elements;
  std::vector<uint64_t> col_hamming_sum;
  bool hamming_distance_cached = false;
  size_t sum_col(size_t col);
  long _check_carry();

  public:
  bit_matrix(size_t _precision = 53);
  bit_matrix(size_t _size, size_t _precision);
  bool get(size_t row, size_t col);
  void cache_hamming_distance();
  uint64_t column_hamming_distance(size_t col);
  void set(size_t row, size_t col, bool v);
  void set(size_t row, double d);
  void push(double d);
  void push();
  void push(const std::vector<uint64_t> &v);
  void complete();
  size_t size();
  size_t precision();
  void import_param(std::istream &is);
  void export_param(std::ostream &os);
};
std::ostream &operator<<(std::ostream &os, bit_matrix mat);
} // namespace momoko::tools

#endif // BIT_MATRIX_HPP
