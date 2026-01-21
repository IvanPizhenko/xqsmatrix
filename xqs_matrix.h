///////////////////////////////////////////////////////////////////////////////
// This is implementation of a 2D matrix class in C++,
// inspired by some ideas described in the following articles:
//
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Source-File
//
// Copyright (c) 2015-2018, 2020, 2024-2026 Ivan Pizhenko.
// All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// OR OTHER DEALINGS IN THE SOFTWARE.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENATION PLAN
///////////////////////////////////////////////////////////////////////////////
//
// 1. Create initial header file.
//    - [x] Create header file with include guards and license.
//
// 2. Class declaration:
//    - [x] Declare template class xqs_matrix<T, Alloc>.
//    - [x] Add public type aliases.
//    - [x] Add private member variables.
//
// 3. Constructors:
//    - [x] Default constructor.
//    - [x] Constructor with row and column counts.
//    - [x] Constructor with row and column counts and initial value.
//    - [x] Constructor with allocator.
//    - [x] Constructor with allocator, row and column counts.
//    - [x] Constructor with allocator, row and column counts
//          and initial value.
//    - [x] Constructor from range defined by iterators.
//    - [x] Constructor from range object.
//    - [x] Copy constructor.
//    - [x] Copy constructor with different allocator object.
//    - [x] Move constructor.
//    - [x] Move constructor with different allocator object.
//    - [x] Constructor from raw data pointer (view).
//    - [x] Constructor from initializer list.
//
// 4. Destructor.
//
// 5. Accessors:
//    - [x] row_count()
//    - [x] column_count()
//    - [x] size()
//    - [x] empty()
//    - [x] stride()
//    - [x] data()
//    - [x] capacity()
//    - [x] is_owner()
//    - [x] is_view()
//
// 7. Row and column access:
//    - [x] operator[](row)
//    - [x] row_at(row)
//    - [x] operator()(column)
//    - [x] column_at(column)
//
// 6. Element access:
//    - [x] operator()(row, column)
//    - [x] operator()(pair<row, column>)
//    - [x] operator[](pair<row, column>)
//    - [x] at(row, column)
//    - [x] at(pair<row, column>)
//
// 8. Assignment operators:
//    - [x] Copy assignment
//    - [x] Move assignment
//
// 9. Swap method.
//    - [x] swap() member function
//    - [x] swap() external function
//
// 10. Resize method.
//     - [x] resize() member function
//     - [x] resize_rows() external function
//     - [x] resize_columns() external function
//     - [x] reserve() member function
//
// 11. Member arithmetic operators:
//     - [x] operator*=(scalar)
//     - [x] operator/=(scalar)
//     - [x] operator+=(matrix)
//     - [x] operator-=(matrix)
//     - [x] operator*=(matrix)
//
// 12. External arithmetic operators:
//     - [x] friend operator*(matrix, scalar)
//     - [x] friend operator/(matrix, scalar)
//     - [x] friend operator+(matrix, matrix)
//     - [x] friend operator-(matrix, matrix)
//     - [x] friend operator*(matrix, matrix)
//
// 13. Special matrix factory functuons:
//     - [x] identity()
//     - [x] transposed_identity()
//
// 14. Matrix operations:
//     - [x] clear()
//     - [x] shrink_to_fit()
//     - [x] copy_full()
//     - [x] move_full()
//     - [x] window()
//     - [x] copy_window()
//     - [x] move_window()
//     - [x] transpose_copy()
//     - [x] transpose_move()
//     - [x] diag_to_hvec()
//     - [x] diag_to_vvec()
//     - [x] inverse_v1()
//     - [x] inverse_v2()
//     - [x] gaussian_reduction()
//
// 15. Stream operations:
//     - [x] operator<<()
//     - [x] operator>>()
//
// 16. Additional I/O operations:
//     - [ ] read_csv()
//
///////////////////////////////////////////////////////////////////////////////

#ifndef XQS_MATRIX_H__
#define XQS_MATRIX_H__

#ifndef XQS_MATRIX_NO_PRAGMA_ONCE
#pragma once
#endif

// CRT
#include <cmath>

// STL
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <type_traits>
#include <vector>

/// @brief Matrix class.
/// @tparam T Element type.
/// @tparam Alloc Allocator type.
template <typename T, class Alloc = std::allocator<T>>
class xqs_matrix {
public:

  // Types

  using value_type = T;
  using allocator_type = Alloc;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

private:

  // Tags for private constructors

  struct identity_matrix_tag {};
  struct transposed_identity_matrix_tag {};

  struct mult_by_scalar_tag {};
  struct div_by_scalar_tag {};

  struct addition_tag {};
  struct subtraction_tag {};
  struct multiplication_tag {};

  struct transpose_copy_tag {};
  struct transpose_move_tag {};

  struct diag_to_hvec_tag {};
  struct diag_to_vvec_tag {};

  struct window_copy_tag {};
  struct window_move_tag {};

public:

  // Public Constructors

  xqs_matrix() noexcept(std::is_nothrow_default_constructible_v<Alloc>) :
    m_capacity(0),
    m_data(nullptr),
    m_row_count(0),
    m_column_count(0),
    m_stride(0),
    m_is_owner(true)
  {
  }

  xqs_matrix(const Alloc& alloc) noexcept(
      std::is_nothrow_copy_constructible_v<Alloc>) :
    m_allocator(alloc),
    m_capacity(0),
    m_data(nullptr),
    m_row_count(0),
    m_column_count(0),
    m_stride(0),
    m_is_owner(true)
  {
  }

  explicit xqs_matrix(const size_type row_count, const size_type column_count) :
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;

    if constexpr (std::is_trivial_v<T>) {
      std::uninitialized_fill_n(m_data, m_capacity, T{});
      return;
    }

    auto p = m_data;
    try {
      const auto e = m_data + m_capacity;
      for (; p != e; ++p) {
        std::construct_at(p);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      const size_type row_count,
      const size_type column_count,
      const_reference v,
      const Alloc& allocator = Alloc()) :
    m_allocator(allocator),
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;

    if constexpr (std::is_trivial_v<T>) {
      std::uninitialized_fill_n(m_data, m_capacity, v);
      return;
    }

    auto p = m_data;
    try {
      const auto e = m_data + m_capacity;
      for (; p != e; ++p) {
        std::construct_at(p, v);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  template <typename InputIt>
  xqs_matrix(
      InputIt first,
      InputIt last,
      const size_type row_count,
      const size_type column_count,
      const Alloc& allocator = Alloc()) :
    m_allocator(allocator),
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;

    auto p = m_data;
    const auto e = m_data + m_capacity;

    if constexpr (std::is_trivial_v<T>) {
      for (; p != e && first != last; ++first, ++p) {
        *p = *first;
      }
      if (p != e) {
        const T a{};
        for (; p != e; ++p) {
          *p = a;
        }
      }
      return;
    }

    try {
      for (; p != e && first != last; ++first, ++p) {
        std::construct_at(p, *first);
      }
      for (; p != e; ++p) {
        std::construct_at(p, *first);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  template <typename Range>
  xqs_matrix(
      [[maybe_unused]] std::from_range_t tag,
      Range&& r,
      const size_type row_count,
      const size_type column_count,
      const Alloc& allocator = Alloc()) :
    m_allocator(allocator),
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;

    auto p = m_data;
    const auto e = m_data + m_capacity;

    if constexpr (std::is_trivial_v<T>) {
      auto first = std::begin(r);
      const auto last = std::end(r);
      for (; p != e && first != last; ++first, ++p) {
        *p = *first;
      }
      if (p != e) {
        const T a{};
        for (; p != e; ++p) {
          *p = a;
        }
      }
      return;
    }

    try {
      auto first = std::begin(r);
      const auto last = std::end(r);
      for (; p != e && first != last; ++first, ++p) {
        std::construct_at(p, *first);
      }
      for (; p != e; ++p) {
        std::construct_at(p, *first);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(const xqs_matrix& src) :
    m_capacity(src.size()),
    m_data(src.m_is_owner
      ? (m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr)
      : src.m_data),
    m_row_count(src.m_row_count),
    m_column_count(src.m_column_count),
    m_stride(m_column_count),
    m_is_owner(src.m_is_owner)
  {
    if (m_data == nullptr) [[unlikely]] return;

    if constexpr (std::is_trivial_v<T>) {
      std::uninitialized_copy_n(m_data, m_capacity, src.m_data);
      return;
    }

    auto p = m_data;
    try {
      const auto e = m_data + m_capacity;
      for (auto q = src.m_data; p != e; ++q, ++p) {
        std::construct_at(p, *q);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(const xqs_matrix& src, const Alloc& allocator) :
    m_allocator(allocator),
    m_capacity(src.size()),
    m_data(src.m_is_owner
      ? (m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr)
      : src.m_data),
    m_row_count(src.m_row_count),
    m_column_count(src.m_column_count),
    m_stride(m_column_count),
    m_is_owner(src.m_is_owner)
  {
    if (m_data == nullptr) [[unlikely]] return;

    if constexpr (std::is_trivial_v<T>) {
      std::uninitialized_copy_n(m_data, m_capacity, src.m_data);
      return;
    }

    auto p = m_data;
    try {
      const auto e = m_data + m_capacity;
      for (auto q = src.m_data; p != e; ++p, ++q) {
        std::construct_at(p, *q);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(xqs_matrix&& src)
      noexcept(std::is_nothrow_move_constructible_v<Alloc>) :
    m_allocator(std::move(src.m_allocator)),
    m_capacity(src.m_capacity),
    m_data(src.m_data),
    m_row_count(src.m_row_count),
    m_column_count(src.m_column_count),
    m_stride(src.m_stride),
    m_is_owner(src.m_is_owner)
  {
    src.m_data = nullptr;
    src.m_capacity = src.m_row_count = src.m_column_count = src.m_stride = 0;
    src.m_is_owner = true;
  }

  xqs_matrix(xqs_matrix&& src, const Alloc& allocator) noexcept(
      std::is_nothrow_copy_constructible_v<Alloc>) :
    m_allocator(allocator),
    m_capacity(src.m_capacity),
    m_data(src.m_data),
    m_row_count(src.m_row_count),
    m_column_count(src.m_column_count),
    m_stride(src.m_stride),
    m_is_owner(src.m_is_owner)
  {
    src.m_data = nullptr;
    src.m_capacity = src.m_row_count = src.m_column_count = src.m_stride = 0;
    src.m_is_owner = true;
  }

  xqs_matrix(
      T* const data,
      const size_type row_count,
      const size_type column_count,
      const size_type stride) noexcept :
    m_capacity(0),
    m_data(data),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(stride),
    m_is_owner(false)
  {
  }

  xqs_matrix(
      std::initializer_list<T> init,
      const size_type row_count,
      const size_type column_count,
      const Alloc& allocator = Alloc()) :
    m_allocator(allocator),
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;

    if constexpr (std::is_trivial_v<T>) {
      const auto n = std::min(init.size(), m_capacity);
      if (n > 0) {
        std::uninitialized_copy_n(init.begin(), n, m_data);
      }
      if (n < m_capacity) {
        const T a{};
        std::uninitialized_fill_n(m_data + n, m_capacity - n, a);
      }
      return;
    }

    auto p = m_data;
    try {
      const auto e = m_data + m_capacity;
      auto first = init.begin();
      const auto last = init.end();
      for (; p != e && first != last; ++first, ++p) {
        std::construct_at(p, *first);
      }
      for (; p != e; ++p) {
        std::construct_at(p);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

private:

// Private constructors

  xqs_matrix(
      [[maybe_unused]] identity_matrix_tag tag,
      const size_type dimension,
      const_reference value) :
    m_capacity(validate_dimensions(dimension, dimension)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(dimension),
    m_column_count(dimension),
    m_stride(dimension),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;
    auto p = m_data;
    if constexpr (std::is_trivial_v<T>) {
      if (m_column_count > 1) [[likely]] {
        const T t{};
        for (size_type i = 0, n = m_row_count - 1; i < n; ++i) {
          *p = value;
          ++p;
          std::uninitialized_fill_n(p, m_column_count, t);
          p += m_column_count;
        }
      }
      *p = value;
      return;
    }

    if (m_column_count > 1) [[likely]] {
      try {
        const auto step = m_column_count - 1;
        for (size_type i = 0, n = m_row_count - 1; i < n; ++i) {
          std::construct_at(p, value);
          ++p;
          for (const auto e = p + step; p != e; ++p) {
            std::construct_at(p);
          }
        }
        std::construct_at(p, value);
        return;
      } catch (...) {
        for (; p != m_data; --p) {
          std::destroy_at(p);
        }
        m_allocator.deallocate(m_data, m_capacity);
        throw;
      }
    }

    try {
      std::construct_at(m_data, value);
      return;
    } catch (...) {
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] transposed_identity_matrix_tag tag,
      const size_type dimension,
      const_reference value) :
    m_capacity(validate_dimensions(dimension, dimension)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(dimension),
    m_column_count(dimension),
    m_stride(dimension),
    m_is_owner(true)
  {
    if (m_data == nullptr) [[unlikely]] return;
    auto p = m_data;
    if constexpr (std::is_trivial_v<T>) {
      if (m_column_count > 1) [[likely]] {
        auto n = m_column_count - 1;
        const T t{};
        std::uninitialized_fill_n(p, n, t);
        p += n;
        --n;
        for (size_type i = 0; i < m_row_count; ++i) {
          *p = value;
          ++p;
          std::uninitialized_fill_n(p, n, t);
          p += n;
        }
        *p = t;
        return;
      }

      *m_data = value;
      return;
    }

    if (m_column_count > 1) [[likely]] {
      try {
        auto n = m_column_count - 1;
        for (const auto e = p + n; p != e; ++p) {
          std::construct_at(p);
        }
        --n;
        for (size_type i = 0; i < m_row_count; ++i) {
          std::construct_at(p);
          ++p;
          for (const auto e = p + n; p != e; ++p) {
            std::construct_at(p);
          }
        }
        std::construct_at(p);
        return;
      } catch (...) {
        for (; p != m_data; --p) {
          std::destroy_at(p);
        }
        m_allocator.deallocate(m_data, m_capacity);
        throw;
      }
    }

    try {
      std::construct_at(m_data, value);
    } catch (...) {
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] mult_by_scalar_tag tag,
      const xqs_matrix& lhs,
      const_reference rhs) :
    m_capacity(validate_dimensions(lhs.m_row_count, lhs.m_column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(m_row_count),
    m_column_count(m_column_count),
    m_stride(m_column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    try {
      auto l = lhs.m_data;
      const auto e = m_data + m_capacity;
      for (; p != e; ++l, ++p) {
        std::construct_at(p, (*l) * rhs);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] div_by_scalar_tag tag,
      const xqs_matrix& lhs,
      const_reference rhs) :
    m_capacity(validate_dimensions(lhs.m_row_count, lhs.m_column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(m_row_count),
    m_column_count(m_column_count),
    m_stride(m_column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    try {
      auto l = lhs.m_data;
      const auto e = m_data + m_capacity;
      for (; p != e; ++l, ++p) {
        std::construct_at(p, (*l) / rhs);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] addition_tag tag,
      const xqs_matrix& lhs,
      const xqs_matrix& rhs) :
    m_capacity(validate_dimensions(lhs.m_row_count, lhs.m_column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(lhs.m_row_count),
    m_column_count(rhs.m_column_count),
    m_stride(lhs.m_column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;

    auto p = m_data;
    auto l = lhs.m_data;
    auto r = rhs.m_data;

    if constexpr (std::is_trivial_v<T>) {
      if (lhs.m_stride == lhs.m_column_count) {
        const auto le = l + lhs.m_row_count * lhs.m_stride;
        if (rhs.m_stride == rhs.m_column_count) {
          for (; l != le; ++p, ++r, ++l) {
            *p = *l + *r;
          }
        } else {
          for (; l != le; r += rhs.m_stride) {
            auto rr = r;
            const auto e = r + rhs.m_column_count;
            for (; rr != e; ++p, ++l, ++rr) {
              *p = *l + *rr;
            }
          }
        }
      } else if (rhs.m_stride == rhs.m_column_count) {
        const auto re = r + rhs.m_row_count * rhs.m_stride;
        for (; r != re; l += m_stride) {
          auto ll = l;
          const auto le = l + lhs.m_column_count;
          for (; ll != le; ++p, ++r, ++ll) {
            *p = *ll + *r;
          }
        }
      } else {
        for (std::size_t i = 0; i < m_row_count;
             l += m_stride, r += rhs.m_stride, ++i) {
          auto rr = r;
          auto ll = l;
          const auto le = l + lhs.m_column_count;
          for (; ll != le; ++p, ++rr, ++ll) {
            *p = *ll + *rr;
          }
        }
      }
      return;
    }

    try {
      if (lhs.m_stride == lhs.m_column_count) {
        const auto le = l + lhs.m_row_count * lhs.m_stride;
        if (rhs.m_stride == rhs.m_column_count) {
          for (; l != le; ++p, ++r, ++l) {
            std::construct_at(p, *l + *r);
          }
        } else {
          for (; l != le; r += rhs.m_stride) {
            auto rr = r;
            const auto e = r + rhs.m_column_count;
            for (; rr != e; ++p, ++l, ++rr) {
              std::construct_at(p, *l + *rr);
            }
          }
        }
      } else if (rhs.m_stride == rhs.m_column_count) {
        const auto re = r + rhs.m_row_count * rhs.m_stride;
        for (; r != re; l += m_stride) {
          auto ll = l;
          const auto lend = l + lhs.m_column_count;
          for (; ll != lend; ++p, ++r, ++ll) {
            std::construct_at(p, *ll + *r);
          }
        }
      } else {
        for (std::size_t i = 0; i < m_row_count;
              l += m_stride, r += rhs.m_stride, ++i) {
          auto rr = r;
          auto ll = l;
          const auto le = l + lhs.m_column_count;
          for (; ll != le; ++p, ++rr, ++ll) {
            std::construct_at(p, *ll + *rr);
          }
        }
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] subtraction_tag tag,
      const xqs_matrix& lhs,
      const xqs_matrix& rhs) :
    m_capacity(validate_dimensions(lhs.m_row_count, lhs.m_column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(lhs.m_row_count),
    m_column_count(lhs.m_column_count),
    m_stride(lhs.m_column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    auto l = lhs.m_data;
    auto r = rhs.m_data;
    if constexpr (std::is_trivial_v<T>) {
      if (lhs.m_stride == lhs.m_column_count) {
        const auto le = l + lhs.m_row_count * lhs.m_stride;
        if (rhs.m_stride == rhs.m_column_count) {
          for (; l != le; ++p, ++r, ++l) {
            *p = *l - *r;
          }
        } else {
          for (; l != le; r += rhs.m_stride) {
            auto rr = r;
            const auto e = r + rhs.m_column_count;
            for (; rr != e; ++p, ++l, ++rr) {
              *p = *l - *rr;
            }
          }
        }
      } else if (rhs.m_stride == rhs.m_column_count) {
        const auto re = r + rhs.m_row_count * rhs.m_stride;
        for (; r != re; l += m_stride) {
          auto ll = l;
          const auto lend = l + lhs.m_column_count;
          for (; ll != lend; ++p, ++r, ++ll) {
            *p = *ll - *r;
          }
        }
      } else {
        for (std::size_t i = 0; i < m_row_count;
             l += m_stride, r += rhs.m_stride, ++i) {
          auto rr = r;
          auto ll = l;
          const auto le = l + lhs.m_column_count;
          for (; ll != le; ++p, ++rr, ++ll) {
            *p = *ll - *rr;
          }
        }
      }
      return;
    }

    try {
      if (lhs.m_stride == lhs.m_column_count) {
        const auto le = l + lhs.m_row_count * lhs.m_stride;
        if (rhs.m_stride == rhs.m_column_count) {
          for (; l != le; ++p, ++r, ++l) {
            std::construct_at(p, *l - *r);
          }
        } else {
          for (; l != le; r += rhs.m_stride) {
            auto rr = r;
            const auto e = r + rhs.m_column_count;
            for (; rr != e; ++p, ++l, ++rr) {
              std::construct_at(p, *l - *rr);
            }
          }
        }
      } else if (rhs.m_stride == rhs.m_column_count) {
        const auto re = r + rhs.m_row_count * rhs.m_stride;
        for (; r != re; l += m_stride) {
          auto ll = l;
          const auto le = l + lhs.m_column_count;
          for (; ll != le; ++p, ++r, ++ll) {
            std::construct_at(p, *ll - *r);
          }
        }
      } else {
        for (std::size_t i = 0; i < m_row_count;
              l += m_stride, r += rhs.m_stride, ++i) {
          auto rr = r;
          auto ll = l;
          const auto le = l + lhs.m_column_count;
          for (; ll != le; ++p, ++rr, ++ll) {
            std::construct_at(p, *ll - *rr);
          }
        }
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] multiplication_tag tag,
      const xqs_matrix& lhs,
      const xqs_matrix& rhs) :
    m_capacity(validate_dimensions(lhs.m_row_count, rhs.m_column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(lhs.m_row_count),
    m_column_count(rhs.m_column_count),
    m_stride(rhs.m_column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;

    auto p = m_data;

    try {
      auto l = lhs.m_data;
      const auto e = l + lhs.size();
      for (; l != e; l += lhs.m_column_count) {
        auto rr = rhs.m_data;
        auto const re = rr + rhs.m_column_count;
        for (; rr != re; ++rr, ++p) {
          auto ll = l;
          auto const le = l + lhs.m_column_count;
          auto r = rr;
          T v{};
          for (; ll != le; ++ll, r += rhs.m_stride) {
            v += (*ll) * (*r);
          }
          std::construct_at(p, std::move(v));
        }
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] transpose_copy_tag tag,
      const xqs_matrix& src) :
    m_capacity(validate_dimensions(src.m_column_count, src.m_row_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(src.m_column_count),
    m_column_count(src.m_row_count),
    m_stride(src.m_row_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    try {
      auto const e = m_data + size();
      auto s = src.m_data;
      const auto se = src.m_data + src.m_column_count;
      for (; s != se; ++s) {
        auto ss = s;
        for (; p != e; ++p, ss += src.m_stride) {
          std::construct_at(p, *ss);
        }
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] transpose_move_tag tag,
      xqs_matrix&& src) :
    m_capacity(validate_dimensions(src.m_column_count, src.m_row_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(src.m_column_count),
    m_column_count(src.m_row_count),
    m_stride(src.m_row_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    try {
      auto const e = m_data + size();
      auto s = src.m_data;
      const auto se = src.m_data + src.m_column_count;
      for (; s != se; ++s) {
        auto ss = s;
        for (; p != e; ++p, ss += src.m_stride) {
          std::construct_at(p, std::move(*ss));
        }
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] diag_to_hvec_tag tag,
      const xqs_matrix& src) :
    m_capacity(validate_dimensions(src.m_column_count, 1)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(src.m_row_count),
    m_column_count(1),
    m_stride(1),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    try {
      auto s = src.m_data;
      const auto step = src.m_stride + 1;
      const auto e = p + row_count;
      for (; p != e; ++p, s += step) {
        std::construct_at(p, *s);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix(
      [[maybe_unused]] diag_to_vvec_tag tag,
      const xqs_matrix& src) :
    m_capacity(validate_dimensions(1, src.m_column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(1),
    m_column_count(src.m_column_count),
    m_stride(src.m_column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    try {
      const auto s = src.m_data;
      const auto step = src.m_stride + 1;
      const auto e = p + row_count;
      for (; p != e; ++p, s += step) {
        std::construct_at(p, *s);
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

  xqs_matrix([[maybe_unused]] window_copy_tag tag,
      const xqs_matrix& src,
      const size_type row_offset,
      const size_type column_offset,
      const size_type row_count,
      const size_type column_count) :
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;
    auto p = m_data;
    auto s = src.m_data + row_offset * src.m_stride + column_offset;
    if (std::is_trivial_v<T>) {
      for (size_type i = 0; i < row_count; p += column_count, s += src.m_stride, ++i) {
        std::uninitialized_copy_n(s, 1, p);
      }
    } else {
      try {
        const auto e = m_data + m_capacity;
        for (size_type i = 0; i < row_count; s += src.m_stride, ++i) {
          auto ss = s;
          const auto se = ss + column_count;
          for (; ss != se; ++p, ++ss) {
            std::construct_at(p, *ss);
          }
        }
      } catch (...) {
        for (; p != m_data; --p) {
          std::destroy_at(p);
        }
        m_allocator.deallocate(m_data, m_capacity);
        throw;
      }
    }
  }

  xqs_matrix([[maybe_unused]] window_move_tag tag,
      const xqs_matrix& src,
      const size_type row_offset,
      const size_type column_offset,
      const size_type row_count,
      const size_type column_count) :
    m_capacity(validate_dimensions(row_count, column_count)),
    m_data(m_capacity != 0 ? m_allocator.allocate(m_capacity) : nullptr),
    m_row_count(row_count),
    m_column_count(column_count),
    m_stride(column_count),
    m_is_owner(true)
  {
    if (empty()) [[unlikely]] return;

    auto p = m_data;
    auto s = src.m_data + row_offset * src.m_stride + column_offset;
    if (std::is_trivial_v<T>) {
      for (size_type i = 0; i < row_count; p += column_count, s += src.m_stride, ++i) {
        std::uninitialized_copy_n(s, 1, p);
      }
      return;
    }

    try {
      const auto e = m_data + m_capacity;
      for (size_type i = 0; i < row_count; s += src.m_stride, ++i) {
        auto ss = s;
        const auto se = ss + column_count;
        for (; ss != se; ++p, ++ss) {
          std::construct_at(p, std::move(*ss));
        }
      }
    } catch (...) {
      for (; p != m_data; --p) {
        std::destroy_at(p);
      }
      m_allocator.deallocate(m_data, m_capacity);
      throw;
    }
  }

public:

  // Destructor

  ~xqs_matrix()
  {
    if (m_is_owner && m_data != nullptr) [[likely]] {
      if constexpr (!std::is_trivial_v<T>) {
        for (auto p = m_data, e = m_data + size(); p != e; ++p) {
          std::destroy_at(p);
        }
      }
      m_allocator.deallocate(m_data, m_capacity);
    }
  }


public:

  // Identity matrices

  static xqs_matrix identity(
      const size_type dimension,
      const_reference value = T(1))
  {
    return xqs_matrix(identity_matrix_tag{}, dimension, value);
  }

  static xqs_matrix transposed_identity(
      const size_type dimension,
      const_reference value = T(1))
  {
    return xqs_matrix(transposed_identity_matrix_tag{}, dimension, value);
  }

  // Helper functions

  void swap(xqs_matrix& other) noexcept
  {
    // m_allocator.swap(other.m_allocator);
    std::swap(m_capacity, other.m_capacity);
    std::swap(m_data, other.m_data);
    std::swap(m_row_count, other.m_row_count);
    std::swap(m_column_count, other.m_column_count);
    std::swap(m_stride, other.m_stride);
    std::swap(m_is_owner, other.m_is_owner);
  }

  // Assignment operators

  xqs_matrix& operator=(const xqs_matrix& rhs)
  {
    if (&rhs != this) [[likely]] {
      xqs_matrix tmp(rhs);
      do_move_assign(std::move(tmp));
    }
    return *this;
  }

  xqs_matrix& operator=(xqs_matrix&& rhs) noexcept
  {
    if (&rhs != this) [[likely]] {
      do_move_assign(std::move(rhs));
    }
    return *this;
  }

  // Math operations

  xqs_matrix& operator+=(const xqs_matrix& rhs)
  {
    check_equal_dimensions(rhs);
    if (!empty()) [[likely]] {
      auto l = m_data;
      auto r = rhs.m_data;
      if (m_stride == m_column_count) {
        const auto le = l + m_row_count * m_stride;
        if (rhs.m_stride == rhs.m_column_count) {
          for (; l != le; ++r, ++l) {
            *l += *r;
          }
        } else {
          for (; l != le; r += rhs.m_stride) {
            auto rr = r;
            const auto re = r + rhs.m_column_count;
            for (; rr != re; ++l, ++rr) {
              *l += *rr;
            }
          }
        }
      } else if (rhs.m_stride == rhs.m_column_count) {
        const auto re = r + rhs.m_row_count * rhs.m_stride;
        for (; r != re; l += m_stride) {
          auto ll = l;
          const auto le = l + rhs.m_column_count;
          for (; ll != le; ++r, ++ll) {
            *ll += *r;
          }
        }
      } else {
        for (size_type i = 0; i < m_row_count;
             l += m_stride, r += rhs.m_stride, ++i) {
          auto rr = r;
          auto ll = l;
          const auto le = l + rhs.m_column_count;
          for (; ll != le; ++rr, ++ll) {
            *ll += *rr;
          }
        }
      }
    }
    return *this;
  }

  xqs_matrix& operator-=(const xqs_matrix& rhs)
  {
    check_equal_dimensions(rhs);
    if (!empty()) [[likely]] {
      auto l = m_data;
      auto r = rhs.m_data;
      if (m_stride == m_column_count) {
        const auto le = l + m_row_count * m_stride;
        if (rhs.m_stride == rhs.m_column_count) {
          for (; l != le; ++r, ++l) {
            *l -= *r;
          }
        } else {
          for (; l != le; r += rhs.m_stride) {
            auto rr = r;
            const auto re = r + rhs.m_column_count;
            for (; rr != re; ++l, ++rr) {
              *l -= *rr;
            }
          }
        }
      } else if (rhs.m_stride == rhs.m_column_count) {
        const auto re = r + rhs.m_row_count * rhs.m_stride;
        for (; r != re; l += m_stride) {
          auto ll = l;
          const auto le = l + rhs.m_column_count;
          for (; ll != le; ++r, ++ll) {
            *ll -= *r;
          }
        }
      } else {
        for (size_type i = 0; i < m_row_count;
             l += m_stride, r += rhs.m_stride, ++i) {
          auto rr = r;
          auto ll = l;
          const auto le = l + rhs.m_column_count;
          for (; ll != le; ++rr, ++ll) {
            *ll -= *rr;
          }
        }
      }
    }
    return *this;
  }

  xqs_matrix& operator*=(const xqs_matrix& rhs)
  {
    *this = (*this) * rhs;
    return *this;
  }

  template <typename U, typename A>
  friend xqs_matrix<U, A> operator+(
      const xqs_matrix<U, A>& lhs,
      const xqs_matrix<U, A>& rhs)
  {
    lhs.check_equal_dimensions(rhs);
    return xqs_matrix(xqs_matrix::addition_tag{}, lhs, rhs);
  }

  template <typename U, typename A>
  friend xqs_matrix<U, A> operator-(
      const xqs_matrix<U, A>& lhs,
      const xqs_matrix<U, A>& rhs)
  {
    lhs.check_equal_dimensions(rhs);
    return xqs_matrix(xqs_matrix::subtraction_tag{}, lhs, rhs);
  }

  template <typename U, typename A>
  friend xqs_matrix<U, A> operator*(
    const xqs_matrix<U, A>& lhs,
    const xqs_matrix<U, A>& rhs)
  {
    lhs.check_suitable_for_product(rhs);
    return {xqs_matrix::multiplication_tag{}, lhs, rhs};
  }

  template <class U, class A>
  friend xqs_matrix<U, A> diag_to_hvec(const xqs_matrix<U, A>& m)
  {
    m.check_is_square(m);
    return xqs_matrix<U, A>(xqs_matrix::diag_to_hvec_tag{}, m);
  }

  template <class U, class A>
  friend xqs_matrix<U, A> diag_to_vvec(const xqs_matrix<U, A>& m)
  {
    m.check_is_square(m);
    return xqs_matrix<U, A>(xqs_matrix::diag_to_vvec_tag{}, m);
  }

  template <class U, class A>
  friend xqs_matrix<U, A> transpose_copy(const xqs_matrix<U, A>& m)
  {
    return xqs_matrix<U, A>(xqs_matrix::transpose_copy_tag{}, m);
  }

  template <class U, class A>
  friend xqs_matrix<U, A> transpose_move(xqs_matrix<U, A>& m)
  {
    return xqs_matrix<U, A>(xqs_matrix::transpose_move_tag{}, std::move(m));
  }

  template <class U, class A>
  friend xqs_matrix<U, A> transpose_move(xqs_matrix<U, A>&& m)
  {
    return xqs_matrix<U, A>(xqs_matrix::transpose_move_tag{}, std::move(m));
  }

  template <class U, class A>
  friend xqs_matrix<U, A> inverse_v1(const xqs_matrix<U, A>& m)
  {
    // Based on the ideas from
    // http://www.sanfoundry.com/java-program-find-inverse-matrix/

    m.check_is_square();

    const U zero{};
    const auto N = m.m_row_count;
    const auto N1 = N - 1;

    xqs_matrix<U, A> a(m);
    const auto index = a.gaussian_reduction();

    // Update the matrix b[i][j] with the ratios stored
    auto b = xqs_matrix<U, A>::identity(N);
    for (std::size_t i = 0; i < N1; ++i) {
      auto pi0 = b.m_data + index[i] * b.m_stride;
      for (std::size_t j = i + 1; j < N; ++j) {
        const auto& av = a.m_data[index[j] * a.m_stride + i];
        auto pi = pi0;
        auto pj = b.m_data + index[j] * b.m_stride;
        const auto pje = pj + N;
        for (; pj != pje; ++pi, ++pj) {
            *pj -= av * (*pi);
        }
      }
    }

    // Perform backward substitutions
    xqs_matrix<U, A> x(N, N);
    auto xr = x.m_data + (N - 1) * x.m_stride;
    auto ar = a.m_data + index[N1] * a.m_stride;
    auto br = b.m_data + index[N1] * b.m_stride;
    const auto& aa = ar[N - 1];
    if (aa == zero) {
      throw std::runtime_error("xqs_matrix: matrix can't be inverted 3");
    }

    for (std::size_t i = 0; i < N; ++i) {
      xr[i] = br[i] / ar[N1];
      for (auto jj = N1, j = N1 - 1; jj > 0; --j, --jj) {
        auto& xji = x.m_data[j * x.m_stride + i];
        xji = b.m_data[index[j] * b.m_stride + i];

        auto px = x.m_data + jj * x.m_stride + i;
        const auto ajrow = a.m_data + index[j] * a.m_stride;
        auto paj = ajrow + jj;
        const auto paje = ajrow + N;
        for (; paj != paje; px += x.m_stride, ++paj) {
          xji -= (*paj) * (*px);
        }

        if (ajrow[j] == zero) {
          throw std::runtime_error("xqs_matrix: matrix can't be inverted 4");
        }

        xji /= ajrow[j];
      }
    }
    return x;
  }

  template <class U, class A>
  friend xqs_matrix<U, A> inverse_v2(const xqs_matrix<U, A>& m)
  {
    m.check_is_square();

    const U zero{};
    const auto N = m.m_row_count;
    const auto N1 = N - 1;

    xqs_matrix<U, A> rm(m);
    auto im = xqs_matrix<U, A>::identity(N);

    for (std::size_t i = 0; i < N - 1; ++i) {
      auto ri = rm.m_data + i * rm.m_stride;
      auto d = ri + i;
      if (*d == zero) {
        throw std::logic_error("xqs_matrix: matrix can't be inverted 1");
      }
      auto ii = im.m_data + i * im.m_stride;
      for (std::size_t col = 0; col < N; ++col) {
        ri[col] /= *d;
        ii[col] /= *d;
      }
      for (std::size_t row = i + 1; row < N; ++row) {
        auto rr = rm.m_data + row * rm.m_stride;
        auto ir = im.m_data + row * im.m_stride;
        d = rr + i;
        for (std::size_t col = 0; col < N; ++col) {
          rr[col] -= ri[col] * (*d);
          ir[col] -= ii[col] * (*d);
        }
      }
    }

    const auto rrs = rm.m_stride - N;
    const auto irs = im.m_stride - N;

    auto ri = rm.m_data + N1 * rm.m_stride;
    auto ii = im.m_data + N1 * im.m_stride;
    for (auto i = N1; i > 0; ri -= rm.m_stride, ii -= im.m_stride, --i) {
      auto d = ri + i;
      if (*d == zero) {
        throw std::logic_error("xqs_matrix: matrix can't be inverted 2");
      }
      for (std::size_t col = 0; col < N; ++col) {
        ri[col] /= *d;
        ii[col] /= *d;
      }
      auto rr = rm.m_data;
      auto ir = im.m_data;
      for (std::size_t row = 0; row < i; rr += rrs, ir += irs, ++row) {
        d = rr + i;
        for (std::size_t col = 0; col < N; ++rr, ++ir, ++col) {
          rr[col] -= ri[col] * (*d);
          ir[col] -= ii[col] * (*d);
        }
      }
    }

    return im;
  }

  // Matrix/scalar operations

  template <class U, class A>
  friend xqs_matrix<U, A> operator*(const xqs_matrix<U, A>& lhs, const U& rhs)
  {
    return xqs_matrix<U, A>(xqs_matrix::mult_by_scalar_tag{}, lhs, rhs);
  }

  template <class U, class A>
  friend xqs_matrix<U, A> operator*(const U& lhs, const xqs_matrix<U, A>& rhs)
  {
    return xqs_matrix<U, A>(xqs_matrix::mult_by_scalar_tag{}, rhs, lhs);
  }

  template <class U, class A>
  friend xqs_matrix<U, A> operator/(const xqs_matrix<U, A>& lhs, const U& rhs)
  {
    return xqs_matrix<U, A>(xqs_matrix::div_by_scalar_tag{}, lhs, rhs);
  }

  template <class U>
  xqs_matrix& operator*=(const U& rhs)
  {
    if (!empty()) [[likely]] {
      auto p = m_data;
      if (rhs.m_stride == rhs.m_column_count) {
        const auto e = p + size();
        for (; p != e; ++p) {
          (*p) *= rhs;
        }
      } else {
        for (std::size_t i = 0; i < m_row_count; p += rhs.m_stride, ++i) {
          auto pp = p;
          const auto e = p + rhs.m_column_count;
          for (; pp != e; ++pp) {
            (*pp) *= *rhs;
          }
        }
      }
    }
    return *this;
  }

  template <class U>
  xqs_matrix& operator/=(const U& rhs)
  {
    if (!empty()) [[likely]] {
      auto p = m_data;
      if (rhs.m_stride == rhs.m_column_count) {
        const auto e = p + size();
        for (; p != e; ++p) {
          (*p) /= rhs;
        }
      } else {
        for (std::size_t i = 0; i < m_row_count; p += rhs.m_stride, ++i) {
          auto pp = p;
          const auto e = p + rhs.m_column_count;
          for (; pp != e; ++pp) {
            (*pp) /= *rhs;
          }
        }
      }
    }
    return *this;
  }

  // Add "count" columns at postion "pos" with inital value "v"
  void insert_columns(
      size_type pos,
      size_type count = 1,
      const_reference v = T());

  // Remove "count" columns at postion "pos"
  void remove_columns(size_type pos, size_type count = 1);

  // Fix elements to zero
  void fix_to_zero(const_reference threshold)
  {
    const T zero{};
    const auto p0e = m_data + m_row_count * m_stride;
    for (auto p0 = m_data; p0 != p0e; p0 += m_stride) {
      const auto pe = p0 + m_column_count;
      for (pointer p = p0; p != pe; ++p) {
        if (std::fabs(*p) < threshold) {
          *p = zero;
        }
      }
    }
  }

  // Access individual rows

  xqs_matrix operator[](const size_type i) noexcept
  {
    return xqs_matrix(m_data + i * m_stride, 1, m_column_count, m_stride);
  }

  const xqs_matrix operator[](const size_type i) const noexcept
  {
    return xqs_matrix(m_data + i * m_stride, 1, m_column_count, m_stride);
  }

  xqs_matrix row_at(const size_type i)
  {
    validate_row_index(i);
    return xqs_matrix(m_data + i * m_stride, 1, m_column_count, m_stride);
  }

  const xqs_matrix row_at(const size_type i) const
  {
    validate_row_index(i);
    return xqs_matrix(m_data + i * m_stride, 1, m_column_count, m_stride);
  }

  // Access individual columns

  xqs_matrix operator()(const size_type i) noexcept
  {
    return xqs_matrix(m_data + i, m_row_count, 1, m_stride);
  }

  const xqs_matrix operator()(const size_type i) const noexcept
  {
    return xqs_matrix(m_data + i, m_row_count, 1, m_stride);
  }

  xqs_matrix column_at(const size_type i)
  {
    validate_column_index(i);
    return xqs_matrix(m_data + i, m_row_count, 1, m_stride);
  }

  const xqs_matrix column_at(const size_type i) const
  {
    validate_column_index(i);
    return xqs_matrix(m_data + i, m_row_count, 1, m_stride);
  }

  // Access individual elements

  reference operator[](
      const std::pair<size_type, size_type> row_and_col) noexcept
  {
    return m_data[row_and_col.first * m_stride + row_and_col.second];
  }

  const_reference operator[](
      const std::pair<size_type, size_type> row_and_col) const noexcept
  {
    return m_data[row_and_col.first * m_stride + row_and_col.second];
  }

  reference operator()(
      const std::pair<size_type, size_type> row_and_col) noexcept
  {
    return m_data[row_and_col.first * m_stride + row_and_col.second];
  }

  const_reference operator()(
      const std::pair<size_type, size_type> row_and_col) const noexcept
  {
    return m_data[row_and_col.first * m_stride + row_and_col.second];
  }

  reference operator()(const size_type row, const size_type col) noexcept
  {
    return m_data[row * m_stride + col];
  }

  const_reference operator()(
      const size_type row, const size_type col) const noexcept
  {
    return m_data[row * m_stride + col];
  }

  reference at(const std::pair<size_type, size_type> row_and_col)
  {
    validate_row_index(row_and_col.first);
    validate_column_index(row_and_col.second);
    return m_data[row_and_col.first * m_stride + row_and_col.second];
  }

  const_reference at(const std::pair<size_type, size_type> row_and_col) const
  {
    validate_row_index(row_and_col.first);
    validate_column_index(row_and_col.second);
    return m_data[row_and_col.first * m_stride + row_and_col.second];
  }

  reference at(const size_type row, const size_type col)
  {
    validate_row_index(row);
    validate_column_index(col);
    return m_data[row * m_stride + col];
  }

  const_reference at(const size_type row, const size_type col) const
  {
    validate_row_index(row);
    validate_column_index(col);
    return m_data[row * m_stride + col];
  }

  // Access dimensions

  size_type row_count() const noexcept
  {
    return m_row_count;
  }

  size_type column_count() const noexcept
  {
    return m_column_count;
  }

  size_type size() const noexcept
  {
    return m_row_count * m_column_count;
  }

  static constexpr size_type max_size() noexcept
  {
    return (std::numeric_limits<size_type>::max() / 2) / sizeof(T);
  }

  bool empty() const noexcept
  {
    return m_row_count == 0 || m_column_count == 0;
  }

  size_type stride() const noexcept
  {
    return m_stride;
  }

  // Access raw data

  pointer data() noexcept
  {
    return m_data;
  }

  const_pointer data() const noexcept
  {
    return m_data;
  }

  // Access raw data buffer properties

  size_type capacity() const noexcept
  {
    return m_capacity;
  }

  bool is_owner() const noexcept
  {
    return m_is_owner;
  }

  bool is_view() const noexcept
  {
    return !m_is_owner;
  }

  // Resize

  void resize_rows(const size_type new_rows)
  {
    resize(new_rows, m_column_count);
  }

  void resize_columns(const size_type new_cols)
  {
    resize(m_row_count, new_cols);
  }

  void resize(const size_type new_rows, const size_type new_cols)
  {
    if (!m_is_owner) {
      throw std::logic_error("xqs_matrix: can't resize a non-owning matrix");
    }

    const auto new_capacity = validate_dimensions(new_rows, new_cols);
    if (new_capacity > m_capacity) {
      auto new_data = new_capacity != 0
          ? m_allocator.allocate(new_capacity)
          : nullptr;
      if (new_data != nullptr) {
        auto p = new_data;
        try {
          // Initialize new memory
          if constexpr (std::is_trivial_v<T>) {
            std::uninitialized_fill_n(new_data, new_capacity, T{});
          } else {
            const auto e = new_data + new_capacity;
            for (; p != e; ++p) {
              std::construct_at(p);
            }
          }

          // Copy existing data
          if (m_data != nullptr) {
            const auto min_rows = std::min(m_row_count, new_rows);
            const auto min_cols = std::min(m_column_count, new_cols);
            for (size_type r = 0; r < min_rows; ++r) {
              std::copy_n(
                m_data + r * m_stride,
                min_cols,
                new_data + r * new_cols);
            }
          }
        } catch (...) {
          if constexpr (!std::is_trivial_v<T>) {
            for (; p != new_data; --p) {
              std::destroy_at(p - 1);
            }
          }
          m_allocator.deallocate(new_data, new_capacity);
          throw;
        }
      }

      // Destroy and deallocate old data
      if (m_data != nullptr) {
        if constexpr (!std::is_trivial_v<T>) {
          for (auto p = m_data, e = m_data + size(); p != e; ++p) {
            std::destroy_at(p);
          }
        }
        m_allocator.deallocate(m_data, m_capacity);
      }

      m_data = new_data;
      m_capacity = new_capacity;
    } else if (new_cols != m_column_count) {
      // Adjust existing data to the new column count
      const auto min_rows = std::min(m_row_count, new_rows);
      const auto min_cols = std::min(m_column_count, new_cols);
      if (new_cols < m_column_count) {
        auto ps = m_data;
        const auto pse = m_data + min_rows * m_stride;
        auto pd = m_data;
        for (; ps != pse; ps += m_stride, pd += new_cols) {
          std::copy_n(ps, min_cols, pd);
        }
      } else {
        auto ps = m_data + (min_rows - 1) * m_stride;
        auto pd = m_data;
        // TODO: recheck carefully for loop condition correctness
        for (size_type r = min_rows; ps > m_data;
             ps -= m_stride, pd += new_cols) {
          std::copy_n(ps, min_cols, pd);
        }
      }
    }

    m_row_count = new_rows;
    m_column_count = new_cols;
    m_stride = new_cols;
  }

  void reserve(const size_type new_capacity)
  {
    if (!m_is_owner) {
      throw std::logic_error(
        "xqs_matrix: can't reserve capacity in a non-owning matrix");
    }

    if (new_capacity > m_capacity) {
      auto new_data = new_capacity != 0
        ? m_allocator.allocate(new_capacity)
        : nullptr;
      if (new_data != nullptr) {
        auto p = new_data;
        try {
          // Copy existing data
          if (m_data != nullptr) {
            const auto n = std::min(m_capacity, new_capacity);
            if constexpr (std::is_trivially_copyable_v<T>) {
              std::copy_n(m_data, n, new_data);
              p += n;
            } else {
              const auto e = m_data + n;
              for (auto s = m_data; s != e; ++s, ++p) {
                std::construct_at(p, *s);
              }
            }
          }

          // Default-initialize remaining elements
          if constexpr (!std::is_trivial_v<T>) {
            const auto e = new_data + new_capacity;
            for (; p != e; ++p) {
              std::construct_at(p);
            }
          }
        } catch (...) {
          if constexpr (!std::is_trivial_v<T>) {
            for (; p != new_data; --p) {
              std::destroy_at(p - 1);
            }
          }
          m_allocator.deallocate(new_data, new_capacity);
          throw;
        }
      }

      // Destroy and deallocate old data
      if (m_data != nullptr) {
        if constexpr (!std::is_trivial_v<T>) {
          for (auto p = m_data, e = m_data + m_capacity; p != e; ++p) {
            std::destroy_at(p);
          }
        }
        m_allocator.deallocate(m_data, m_capacity);
      }

      m_data = new_data;
      m_capacity = new_capacity;
    }
  }

  void clear() noexcept
  {
    if constexpr (!std::is_trivial_v<T>) {
      if (m_is_owner && m_data != nullptr) [[likely]] {
        for (auto p = m_data, e = m_data + size(); p != e; ++p) {
          std::destroy_at(p);
        }
        m_stride = 0; // TODO: check if this is needed
      }
    }
    m_row_count = 0;
    m_column_count = 0;
  }

  void shrink_to_fit()
  {
    if (!m_is_owner) {
      throw std::logic_error(
        "xqs_matrix: can't shrink to fit in a non-owning matrix");
    }

    const size_type new_capacity = m_row_count * m_column_count;

    if (new_capacity < m_capacity) {
      auto new_data = new_capacity != 0
          ? m_allocator.allocate(new_capacity)
          : nullptr;

      if (new_data != nullptr) {
        auto p = new_data;
        try {
          // Copy existing data
          if (m_data != nullptr) {
            const auto n = new_capacity;
            if constexpr (std::is_trivially_copyable_v<T>) {
              std::copy_n(m_data, n, new_data);
              p += n;
            } else {
              const auto e = m_data + n;
              for (auto s = m_data; s != e; ++s, ++p) {
                std::construct_at(p, *s);
              }
            }
          }
        } catch (...) {
          if constexpr (!std::is_trivial_v<T>) {
            for (; p != new_data; --p) {
              std::destroy_at(p - 1);
            }
          }
          m_allocator.deallocate(new_data, new_capacity);
          throw;
        }
      }

      // Destroy and deallocate old data
      if (m_data != nullptr) {
        if constexpr (!std::is_trivial_v<T>) {
          for (auto p = m_data, e = m_data + size(); p != e; ++p) {
            std::destroy_at(p);
          }
        }
        m_allocator.deallocate(m_data, m_capacity);
      }

      m_data = new_data;
      m_capacity = new_capacity;
    }
  }

  // Copy matrix

  template <typename OutputIt>
  OutputIt copy_full(OutputIt dest) const
  {
    if (!empty()) {
      auto src = m_data;
      for (size_type r = 0; r < m_row_count; src += m_stride, ++r) {
        std::copy_n(src, m_column_count, dest);
      }
    }
    return dest;
  }

  template <typename OutputIt>
  OutputIt move_full(OutputIt dest)
  {
    if (!empty()) {
      auto s = m_data;
      for (size_type r = 0; r < m_row_count; s += m_stride, ++r) {
        // Unfortunately, there is no std::move_n()
        if constexpr (std::is_trivially_copyable_v<T>) {
          dest = std::copy_n(s, m_column_count, dest);
        } else {
          const auto se = s + m_column_count;
          for (; s != se; ++s, ++dest) {
            *dest = std::move(*s);
          }
        }
      }
    }
    return dest;
  }

  xqs_matrix window(
    const size_type row,
    const size_type col,
    const size_type row_count,
    const size_type column_count) const
  {
    validate_row_index(row);
    validate_column_index(col);
    validate_row_index(row + row_count - 1);
    validate_column_index(col + column_count - 1);
    return xqs_matrix(
      m_data + row * m_stride + col,
      row_count,
      column_count,
      m_stride);
  }

  xqs_matrix copy_window(
    const size_type row,
    const size_type col,
    const size_type row_count,
    const size_type column_count) const
  {
    validate_row_index(row);
    validate_column_index(col);
    validate_row_index(row + row_count - 1);
    validate_column_index(col + column_count - 1);
    xqs_matrix result(row_count, column_count);
    if (!result.empty()) {
      auto s = m_data + row * m_stride + col;
      auto d = result.m_data;
      for (size_type r = 0; r < row_count;
           s += m_stride, d += column_count, ++r) {
        std::copy_n(s, column_count, d);
      }
    }
    return result;
  }

  xqs_matrix move_window(
    const size_type row,
    const size_type col,
    const size_type row_count,
    const size_type column_count)
  {
    validate_row_index(row);
    validate_column_index(col);
    validate_row_index(row + row_count - 1);
    validate_column_index(col + column_count - 1);
    xqs_matrix result(row_count, column_count);
    if (!result.empty()) {
      auto s = m_data + row * m_stride + col;
      auto d = result.m_data;
      for (size_type r = 0; r < row_count;
           s += m_stride, d += column_count, ++r) {
        // Unfortunately, there is no std::move_n()
        if constexpr (std::is_trivially_copyable_v<T>) {
          std::copy_n(s, column_count, d);
        } else {
          auto de = d + column_count;
          for (; d != de; ++s, ++d) {
            *d = std::move(*s);
          }
        }
      }
    }
    return result;
  }

private:

  // Check that matrices have equal dimensions
  void check_equal_dimensions(const xqs_matrix<T>& other) const
  {
    if (m_row_count == other.m_row_count &&
        m_column_count == other.m_column_count) [[likely]] return;
    std::ostringstream err;
    err << "xqs_matrix: dimensions of the other matrix differ "
        << m_row_count << '*' << m_column_count << " vs "
        << other.m_row_count << '*' << other.m_column_count << ')';
    throw std::invalid_argument(err.str());
  }

  // Check that matrix has dimensions that are suitable for product oeration
  void check_suitable_for_product(const xqs_matrix<T>& other) const
  {
    if (m_column_count == other.m_row_count) [[likely]] return;
    std::ostringstream err;
    err << "xqs_matrix: dimensions of the other matrix are not suitable"
            " for the product [this * other] ("
        << m_row_count << '*' << m_column_count << " vs "
        << other.m_row_count << '*' << other.m_column_count << ')';
    throw std::invalid_argument(err.str());
  }

  // Check that matrix is square
  void check_is_square() const
  {
    if (m_row_count == m_column_count) [[unlikely]] return;
    std::ostringstream err;
    err << "xqs_matrix: matrix is not square ("
        << m_row_count << '*' << m_column_count << ')';
    throw std::invalid_argument(err.str());
  }

  // Validate row index
  void validate_row_index(const size_type row) const
  {
    if (row < m_row_count) [[likely]] return;
    std::ostringstream err;
    err << "xqs_matrix: row index " << row << " is out of range ("
        << m_row_count << ')';
    throw std::invalid_argument(err.str());
  }

  // Validate column index
  void validate_column_index(std::size_t col) const
  {
    if (col < m_column_count) [[likely]] return;
    std::ostringstream err;
    err << "xqs_matrix: column index " << col << " is out of range ("
        << m_column_count << ')';
    throw std::invalid_argument(err.str());
  }

  // Initial validation of dimensions
  static size_type validate_dimensions(const size_type row_count, const size_type column_count)
  {
    if (row_count == 0 || column_count == 0) [[unlikely]] return 0;
    if (column_count > max_size() / row_count) [[unlikely]] {
      throw std::length_error("xqs_matrix: matrix size is too large");
    }
    return row_count * column_count;
  }

  void do_move_assign(xqs_matrix&& rhs) noexcept
  {
    if (m_is_owner) delete[] m_data;

    m_capacity = rhs.m_capacity;
    m_data = rhs.m_data;
    m_row_count = rhs.m_row_count;
    m_column_count = rhs.m_column_count;
    m_stride = rhs.m_stride;
    m_is_owner = rhs.m_is_owner;

    rhs.m_data = nullptr;
    rhs.m_capacity = rhs.m_row_count = rhs.m_column_count = rhs.m_stride = 0;
    rhs.m_is_owner = true;
  }

  // Helper function for gaussian reduction.
  // Used to find inverse matrix.
  std::unique_ptr<std::size_t[]> gaussian_reduction()
  {
    const auto N = m_row_count;
    std::unique_ptr<size_type[]> index(new size_type[N]);
    std::iota(index.get(), index.get() + N, 0);

    // Find the rescaling factors, one from each row
    const T zero{};
    auto c = m_allocator.allocate(N);
    auto pc = c;
    const auto pce = c + N;
    try {
      auto r = m_data;
      for (; pc != pce; ++pc, r += m_stride) {
        auto p = r;
        const auto e = r + m_column_count;
        std::construct_at(pc);
        for (; p != e; ++p) {
          auto c0 = std::abs(*p);
          if (c0 > *pc) *pc = std::move(c0);
        }
      }

      // Search the pivoting element from each column
      size_type k = 0;
      auto ppj = m_data;
      for (size_type j = 0; j < N - 1 ; ++ppj, ++j) {
        T pi1{};
        for (size_type i = j; i < N; ++i) {
          const auto ii = index[i];
          auto pi0 = std::abs(ppj[ii * m_stride]);
          if (c[ii] == zero) [[unlikely]] {
            throw std::runtime_error("xqs_matrix: matrix can't be inverted 5");
          }
          pi0 /= c[ii];
          if (pi0 > pi1) {
            pi1 = std::move(pi0);
            k = i;
          }
        }

        // Interchange rows according to the pivoting order
        std::swap(index[k], index[j]);
        const auto row0 = m_data + index[j] * m_stride;
        const auto& v = row0[j];
        if (v == zero) [[unlikely]] {
          throw std::runtime_error("xqs_matrix: matrix can't be inverted 6");
        }

        for (auto i = j + 1; i < N; ++i) {
          auto row = m_data + index[i]  * m_stride;
          auto& v2 = row[j];
          auto pj = v2 / v;

          // Record pivoting ratios below the diagonal
          v2 = pj;

          // Modify other elements accordingly
          for (auto l = j + 1; l < N; ++l) {
            row[l] -= pj * row0[l];
          }
        }
      }
    } catch (...) {
      if (!std::is_trivial_v<T>) {
        for (; pc != c; --pc) {
          std::destroy_at(pc);
        }
      }
      m_allocator.deallocate(c, N);
      throw;
    }

    if (!std::is_trivial_v<T>) {
      for (; pc != c; --pc) {
        std::destroy_at(pc);
      }
    }
    m_allocator.deallocate(c, N);

    return index;
  }

  Alloc m_allocator;
  size_type m_capacity;
  pointer m_data;
  size_type m_row_count;
  size_type m_column_count;
  size_type m_stride;
  bool m_is_owner;
};

// ***** end of class *****

template <typename T, class Alloc>
inline void swap(xqs_matrix<T, Alloc>& a, xqs_matrix<T, Alloc>& b) noexcept
{
  a.swap(b);
}

template<class T, class Alloc, class Ch, class Traits>
std::basic_ostream<Ch, Traits>& operator<<(
  std::basic_ostream<Ch, Traits>& os,
  const xqs_matrix<T, Alloc>& m)
{
  typename std::basic_ostream<Ch, Traits>::sentry sentry(os);

  const auto row_count = m.row_count();
  const auto column_count = m.column_count();

  os << row_count;
  if (!os) [[unlikely]] return os;
  os << Ch('\t');
  if (!os) [[unlikely]] return os;
  os << column_count;
  if (!os) [[unlikely]] return os;
  os << Ch('\n');
  if (!os) [[unlikely]] return os;

  auto p0 = m.data();
  const auto stride = m.stride();
  const auto p0e = p0 + row_count * stride;
  for (; p0 != p0e; p0 += stride) {
    auto p = p0;
    const auto pe = p + column_count;
    os << *p;
    if (!os) [[unlikely]] return os;
    for (++p; p != pe; ++p) {
      os << Ch('\t');
      if (!os) [[unlikely]] return os;
      os << *p;
      if (!os) [[unlikely]] return os;
    }
    os << Ch('\n');
    if (!os) [[unlikely]] return os;
  }

  return os;
}

template<class T, class Alloc, class Ch, class Traits>
std::basic_istream<Ch, Traits>& operator>>(
  std::basic_istream<Ch, Traits>& is,
  xqs_matrix<T, Alloc>& m)
{
  typename std::basic_istream<Ch, Traits>::sentry sentry(is);

  std::size_t row_count, column_count;
  is >> row_count;
  if (!is) [[unlikely]] return is;
  is >> column_count;
  if (!is) [[unlikely]] return is;
  m.resize(row_count, column_count);

  auto p0 = m.data();
  const auto stride = m.stride();
  const auto p0e = p0 + row_count * stride;
  for (; p0 != p0e; p0 += stride) {
    auto p = p0;
    const auto pe = p + column_count;
    for (; p != pe; ++p) {
      is >> *p;
      if (!is) [[unlikely]] return is;
    }
  }

  return is;
}

namespace detail {

void split()
{

}

} // namespace detail

template <
    typename T,
    typename Alloc,
    typename Converter,
    typename Ch,
    typename Traits,
    typename StrAlloc
>
xqs_matrix<T, Alloc> read_csv(
  const std::filesystem::path& path,
  const Ch line_delim,
  const std::basic_string_view<Ch>& field_delims,
  const Converter& conv,
  const bool has_header_line = true)
{
  xqs_matrix<U, Alloc> result;

  // Open input file
  std::basic_ifstream<Ch, Traits> in(path.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("read_csv: can't open input file");
  }

  std::string line;

  // Skip header lines
  if (has_header_line) {
    std::getlinr(in, line, line_delim);
    if (!in) throw std::runtime_error("read_csv: can't read file header");
  }

  // Parse data lines
  std::size_t data_line_count = 0;
  while (std::getline(in, line, line_delim)) {
    ++data_line_count;
    auto row = parse_vector(line, conv, field_delims);
    if (row.empty()) {
      throw std::runtime_error("read_csv: there is empty data line");
    }
    if (result.m_column_count != row.size()) {
      if (result.m_column_count < row.size()) {
        result.column_count(row.size());
      } else {
        row.resize(result.m_column_count);
      }
    }
    result.m_data.push_back(std::move(row));
    ++result.m_row_count;
  }

  // Ensure that at least one row have been successfully read
  if (data_line_count == 0) {
    throw std::runtime_error("read_csv: there is no data");
  }

  return result;
}

#endif // XQS_MATRIX_H__
