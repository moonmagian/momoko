#include "bit_matrix.hpp"
#include <algorithm>
std::size_t momoko::tools::bit_matrix::sum_col(std::size_t col) {
  size_t result{0};
  for (size_t i = 0; i < elements.size(); ++i) {
    result += get(i, col);
  }
  return result;
}

long momoko::tools::bit_matrix::_check_carry() {
  size_t carry = 0;
  int i;
  for (i = prec - 1; i >= 0; --i) {
    size_t now = column_hamming_distance(i) + carry;
    carry = now / 2;
  }
  return carry;
}

momoko::tools::bit_matrix::bit_matrix(std::size_t _precision)
    : prec(_precision), row_vec_size{(_precision - 1) / (sizeof(uint64_t) * 8) +
                                     1},
      col_hamming_sum(_precision, 0) {}

momoko::tools::bit_matrix::bit_matrix(std::size_t _size, std::size_t _precision)
    : prec{_precision}, row_vec_size{(_precision - 1) / (sizeof(uint64_t) * 8) +
                                     1},
      elements(_size, std::vector<uint64_t>(row_vec_size, 0)),
      col_hamming_sum(_precision, 0) {}

bool momoko::tools::bit_matrix::get(size_t row, size_t col) {
  const std::vector<uint64_t> &row_elements = elements[row];
  size_t element_index{col / (sizeof(uint64_t) * 8)};
  size_t element_offset{sizeof(uint64_t) * 8 - col % (sizeof(uint64_t) * 8) -
                        1};
  auto element = row_elements[element_index];
  return element & (0x1ULL << element_offset);
}

void momoko::tools::bit_matrix::cache_hamming_distance() {
  if (hamming_distance_cached) {
    return;
  }
  for (size_t i = 0; i < prec; ++i) {
    col_hamming_sum[i] = sum_col(i);
  }
  hamming_distance_cached = true;
}

uint64_t momoko::tools::bit_matrix::column_hamming_distance(std::size_t col) {
  if (!hamming_distance_cached) {
    cache_hamming_distance();
  }
  return col_hamming_sum[col];
}

void momoko::tools::bit_matrix::set(std::size_t row, std::size_t col, bool v) {
  std::vector<uint64_t> &row_elements = elements[row];
  size_t element_index{col / (sizeof(uint64_t) * 8)};
  size_t element_offset{sizeof(uint64_t) * 8 - col % (sizeof(uint64_t) * 8) -
                        1};
  if (v) {
    row_elements[element_index] |= 0x1ULL << element_offset;
  } else {
    row_elements[element_index] &= ~(0x1ULL << element_offset);
  }
}

void momoko::tools::bit_matrix::set(std::size_t row, double d) {
  elements[row][0] = tools::double_frac_to_integer(d);
}

void momoko::tools::bit_matrix::push(double d) {
  hamming_distance_cached = false;
  elements.emplace_back(row_vec_size, 0);
  set(elements.size() - 1, d);
}

void momoko::tools::bit_matrix::push() {
  hamming_distance_cached = false;
  elements.emplace_back(row_vec_size, 0);
}

void momoko::tools::bit_matrix::push(const std::vector<uint64_t> &v) {
  hamming_distance_cached = false;
  elements.push_back(v);
  elements.back().resize(row_vec_size, 0);
}

void momoko::tools::bit_matrix::complete() {
  bool larger_than_1{_check_carry() >= 1};
  size_t carry = 0;
  std::vector<bool> result_row(prec, false);
  for (size_t i = 0; i < prec; ++i) {
    size_t col{prec - i - 1};
    auto col_sum_result{sum_col(col)};
    col_hamming_sum[col] = col_sum_result;
    auto result = col_sum_result + carry;
    if (result % 2 == 1) {
      if (larger_than_1) {
        // Reduce a bit.
        size_t bit_remaining{1};
        // If no bit to reduce in current bit, reduce two in next bit.
        for (size_t current_col{col}; bit_remaining && current_col < prec;
             ++current_col) {
          size_t reduced_bit{0};
          for (i = 0; i < elements.size(); ++i) {
            size_t row{elements.size() - i - 1};
            if (get(row, current_col)) {
              set(row, current_col, false);
              --bit_remaining;
              ++reduced_bit;
            }
            if (!bit_remaining) {
              break;
            }
          }
          col_hamming_sum[current_col] -= reduced_bit;
          bit_remaining *= 2;
        }
        result -= 1;
      } else {
        // Add a bit.
        result_row[col] = true;
        col_hamming_sum[col] += 1;
        result += 1;
      }
    }
    carry = result / 2;
  }
  if (!larger_than_1 && std::any_of(result_row.begin(), result_row.end(),
                                    [](bool x) { return x; })) {
    push();
    for (size_t i = 0; i < prec; ++i) {
      if (result_row[i]) {
        set(elements.size() - 1, i, true);
      }
    }
  }
  hamming_distance_cached = true;
}

std::size_t momoko::tools::bit_matrix::size() { return elements.size(); }

std::size_t momoko::tools::bit_matrix::precision() { return prec; }

std::ostream &momoko::tools::operator<<(std::ostream &os, bit_matrix mat) {
  for (size_t i = 0; i < mat.size(); ++i) {
    for (size_t j = 0; j < mat.precision(); ++j) {
      os << mat.get(i, j);
    }
    os << std::endl;
  }
  return os;
}
