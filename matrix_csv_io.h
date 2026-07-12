///////////////////////////////////////////////////////////////////////////////
// This is implementation of the 2D matrix in C++,
// inspired by ideas described in the following articles:
//
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Source-File
//
// Copyright (c) 2015-2018, 2020, 2024-2026 Ivan Pizhenko.
// All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
// THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include "matrix.h"

// STL
#include <filesystem>
#include <fstream>

namespace stdx {

template <
    typename T,
    typename Alloc,
    typename Ch,
    typename Traits,
    typename Converter
>
matrix<T, Alloc> read_csv(
  const std::filesystem::path& path,
  const Ch line_delim,
  const std::basic_string_view<Ch, Traits>& field_delims,
  const Converter& conv,
  const bool has_header_line = true)
{
  matrix<T, Alloc> result;

  // Open input file
  std::basic_ifstream<Ch, Traits> in(path.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("read_csv: can't open input file");
  }

  std::string line;

  // Skip header lines
  std::size_t col_count = 0;
  if (has_header_line) {
    std::getline(in, line, line_delim);
    if (!in) throw std::runtime_error("read_csv: can't read file header");
    col_count = std::count_if(
      line.cbegin(),
      line.cend(),
      [&field_delims](const Ch c) noexcept
      {
         return field_delims.find(c) !=
                std::basic_string_view<Ch, Traits>::npos;
      }
    ) + 1;
  }

  // Parse data lines
  std::size_t row_count = 0;
  while (std::getline(in, line, line_delim)) {
    ++row_count;
    auto row = parse_vector(line, conv, field_delims);
    if (row.empty()) throw std::runtime_error("read_csv: empty data line");
    if (result.m_col_count != row.size()) {
      if (result.m_col_count < row.size()) {
        result.col_count(row.size());
      } else {
        row.resize(result.m_col_count);
      }
    }
    result.m_data.push_back(std::move(row));
    ++result.m_row_count;
  }

  return result;
}

} // namespace stdx
