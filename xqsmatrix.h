// This is slightly modified version of public available QSMatix class,
// described in these articles:
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Source-File
//
// License terms:
//
// Copyright © 2012-2017 Michael Halls-Moore
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// Changes from the original code presented in the article:
// - Renamed QSMatrix data members
// - Removed unnecesary "this->" before member access
// - Optimized constructors and assignment operators
// - Added move constructor and assignment operator
// - Using "std::size_t" instead of "unsigned" for the matrix dimensions
// - Modified meber operator() signature to eliminate passing primitive type 
//   by reference
// - Added new methods at() to enable access with checking dimension boundaries
// - Fixed error in the method transpose(), which worked incorrectly with
//   non-square matrices
// - Modified member operator*(const XQSMatrix&) to validate matrix dimensions.
// - Added new static method identity() for generation of identity matrix
// - added new method inverse() to compute of inverse matrix
// - Added new static method readCsv() to construct matrix from data in the
//   CSV file
// - Added new methods row_count() and col_count() to resize matrix
// - added method window() to create new matrix from rectangular window
//   in the current matrix
// - Added new method swap() and used it in some other methods to optimize
//   computation speed
// - Added new method add_columns() to inject new columns size_to matrix
// - Added new method remove_columns() to remove columns from matrix
// - Added new methods mul_by_row() to multiply this matrix by vector,
//   which represents single row matrix 
// - Added new methods mul_by_column() to multiply this matrix by vector,
//   which represents single column matrix 
// - Added new method row_scalar_product() to find scalar product of matix row
//   with a given vector
// - Added new method column_scalar_product() to find scalar product of matix
//   column with a given vector
// - Added new method row() to extract row as matrix
// - Added new method row_as_vector() to extract row as vector
// - Added new method col() to extract column as matrix
// - Added new method col_as_vector() to extract column as vector
// - Added stream output operator
// - Class renamed to XQSMatrix
// - Added operator[] and new versions of method at() to access individual rows
// - Added fix_to_zero() to quickly eliminate too small values
// - Reworked operators
//
// Changes are maintained on the GitHub:
// https://github.com/IvanPizhenko/xqsmatrix
//
// License terms for modifications:
//
// Copyright © 2015-2017, 2018, 2020, 2024 Ivan Pizhenko
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
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#ifndef XQSMATRIX_H__
#define XQSMATRIX_H__

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>


// Tokenize string and parse tokens as matrix cell values.
// This code is based on the public domain code taken from here:
// https://stackoverflow.com/a/1493195/1540501
template <typename T, typename Converter>
std::vector<T> tokenizeAndParse(const std::string& str,
  const Converter& converter, const std::string& delimiters = " ",
  bool trimEmpty = false)
{
  std::vector<T> result;
  std::string::size_type pos, lastPos = 0, length = str.length();
  while (lastPos < length + 1) {
    pos = str.find_first_of(delimiters, lastPos);
    if (pos == std::string::npos) {
      pos = length;
    }
    if (pos != lastPos || !trimEmpty) {
      result.push_back(converter(str.substr(lastPos, pos - lastPos)));
    }
    lastPos = pos + 1;
  }
  return result;
}

template <typename T>
class XQSMatrix {
public:
  // Constructors
  explicit XQSMatrix(std::size_t nrows = 1, std::size_t ncols = 1);
  XQSMatrix(std::size_t nrows, std::size_t ncols, const T& v);
  XQSMatrix(const XQSMatrix& src);
  XQSMatrix(XQSMatrix&& src);

  // Create identity matrix with K on the diaginal
  static XQSMatrix<T> identity(std::size_t n, const T& k = T(1));
  
  // Swap matrices
  void swap(XQSMatrix& other) noexcept;

  // Operator overloading, for "standard" mathematical matrix operations
  XQSMatrix& operator=(const XQSMatrix<T>& rhs);
  XQSMatrix& operator=(XQSMatrix<T>&& rhs);
  
  // Matrix mathematical operations

  template <class T1>
  friend XQSMatrix<T1> operator+(const XQSMatrix<T1>& lhs, const XQSMatrix<T1>& rhs);

  template <class T1>
  friend XQSMatrix<T1> operator-(const XQSMatrix<T1>& lhs, const XQSMatrix<T1>& rhs);

  template <class T1>
  friend XQSMatrix<T1> operator*(const XQSMatrix<T1>& lhs, const XQSMatrix<T1>& rhs);

  XQSMatrix& operator+=(const XQSMatrix<T>& rhs);
  XQSMatrix& operator-=(const XQSMatrix<T>& rhs);
  XQSMatrix&operator*=(const XQSMatrix<T>& rhs);

  template <class T1>
  friend XQSMatrix<T1> transpose(const XQSMatrix<T1>& m);

  template <class T1>
  friend XQSMatrix<T1> inverse_v1(const XQSMatrix<T1>& m);

  template <class T1>
  friend XQSMatrix<T1> inverse_v2(const XQSMatrix<T1>& m);

  // Matrix/scalar operations

  template <class T1>
  friend XQSMatrix<T1> operator*(const XQSMatrix<T1>& lhs, const T1& rhs);

  template <class T1>
  friend XQSMatrix<T1> operator/(const XQSMatrix<T1>& lhs, const T1& rhs);

  XQSMatrix& operator*=(const T& rhs);
  XQSMatrix& operator/=(const T& rhs);

  // Multiple by vector as row
  template <class T1>
  friend XQSMatrix<T1> mul_by_row(const XQSMatrix<T1>& lhs, const std::vector<T1>& rhs);

  // Multiple by vector as column
  template <class T1>
  friend std::vector<T1> mul_by_column(const XQSMatrix<T1>& lhs, const std::vector<T1>& rhs);

  // Scalar product of the row with a given vector
  template <class T1>
  friend T1 row_scalar_product(const XQSMatrix<T1>& lhs, std::size_t row_index, const std::vector<T>& rhs);

  // Scalar product of the column with a given vector
  template <class T1>
  friend T1 column_scalar_product(const XQSMatrix<T1>& lhs, std::size_t col_index, const std::vector<T1>& rhs);
  
  // Add "count" columns at postion "pos" with inital value "v"
  void add_columns(std::size_t pos, std::size_t count, const T& v);

  // Remove "count" columns at postion "pos"
  void remove_columns(std::size_t pos, std::size_t count);

  // Fix elements to zero
  void fix_to_zero(const T& threshold);

  // Access the rows
  std::vector<T>& operator[](std::size_t i) noexcept
  {
    return m_data[i];
  }

  const std::vector<T>& operator[](std::size_t i) const noexcept
  {
    return m_data[i];
  }

  std::vector<T>& at(std::size_t i)
  {
    return m_data.at(i);
  }

  const std::vector<T>& at(std::size_t i) const
  {
    return m_data.at(i);
  }

  // Access the individual elements
  T& operator()(std::size_t row, std::size_t col) noexcept
  {
    return m_data[row][col];
  }

  const T& operator() (std::size_t row, std::size_t col) const noexcept
  {
    return m_data[row][col];
  }

  T& at(std::size_t row, std::size_t col)
  {
    return m_data.at(row).at(col);
  }

  const T& at(std::size_t row, std::size_t col) const
  {
    return m_data.at(row).at(col);
  }

  // Access the row and column sizes
  std::size_t row_count() const noexcept
  {
    return m_nrows;
  }

  std::size_t col_count() const noexcept
  {
    return m_ncols;
  }

  // Change row and colums sizes
  void row_count(std::size_t new_rows);
  void col_count(std::size_t new_cols);

  // Extact rectangular window as new matrix
  XQSMatrix window(std::size_t row, std::size_t col,
    std::size_t nrows, std::size_t ncols) const;

  // Extract matrix row as matrix
  XQSMatrix<T> row(std::size_t index) const;

  // Extract matrix row as vector
  std::vector<T> row_as_vector(std::size_t index) const;

  // Extract matrix column as matrix
  XQSMatrix<T> col(std::size_t index) const;

  // Extract matrix column as vector
  std::vector<T> col_as_vector(std::size_t index) const;

  std::vector<std::vector<T>>& data() noexcept
  {
    return m_data;
  }

  const std::vector<std::vector<T>>& data() const noexcept
  {
    return m_data;
  }

  // Read from CSV file
  template <typename T1, typename Converter>
  friend XQSMatrix<T1> readCsv(const std::string& path, char lineEnding,
    const std::string& fieldDelimiters, const Converter& converter,
    std::size_t numberOfHeaderLines);

private:
  std::size_t m_nrows;
  std::size_t m_ncols;
  std::vector<std::vector<T>> m_data;

  // Check that matrix has equal dimensions
  void check_equal_dimensions(const XQSMatrix<T>& other) const;

  // Check that matrix has dimensions that are suitable for product oeration
  void check_suitable_for_product(const XQSMatrix<T>& other) const;

  // Validate row index
  void validate_row_index(std::size_t index) const;

  // Validate column index
  void validate_column_index(std::size_t index) const;

  // Helper function for gaussian reduction.
  // Used to find inverse matrix.
  std::vector<size_t> gaussian_reduction();
};

template <typename T>
XQSMatrix<T>::XQSMatrix(std::size_t nrows, std::size_t ncols) :
  m_nrows(nrows),
  m_ncols(ncols),
  m_data(m_nrows)
{
  for (std::size_t i = 0; i < nrows; ++i) {
    m_data[i].resize(ncols);
  }
}

template <typename T>
XQSMatrix<T>::XQSMatrix(std::size_t nrows, std::size_t ncols, const T& v) :
  m_nrows(nrows),
  m_ncols(ncols),
  m_data(m_nrows)
{
  for (std::size_t i = 0; i < nrows; ++i) {
    m_data[i].resize(ncols, v);
  }
}

template <typename T>
XQSMatrix<T>::XQSMatrix(const XQSMatrix<T>& src) :
  m_nrows(src.m_nrows),
  m_ncols(src.m_ncols),
  m_data(src.m_data)
{
}

template <typename T>
XQSMatrix<T>::XQSMatrix(XQSMatrix<T>&& src) :
  m_nrows(src.m_nrows),
  m_ncols(src.m_ncols),
  m_data(std::move( src.m_data))
{
}

template <typename T>
XQSMatrix<T> XQSMatrix<T>::identity(std::size_t n, const T& k)
{
  XQSMatrix result(n, n, 0);
  for (std::size_t i = 0; i < n; ++i)
    result.m_data[i][i] = k;
  return result;
}

template <typename T>
void XQSMatrix<T>::swap(XQSMatrix<T>& other) noexcept
{
  std::swap(m_nrows, other.m_nrows);
  std::swap(m_ncols, other.m_ncols);
  m_data.swap(other.m_data);
}

template <typename T>
inline void swap(XQSMatrix<T>& a, XQSMatrix<T>& b) noexcept
{
  a.swap(b);
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator=(const XQSMatrix<T>& rhs)
{
  if (&rhs != this) {
    m_data = rhs.m_data;
    m_nrows = rhs.m_nrows;
    m_ncols = rhs.m_ncols;
  }
  return *this;
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator=(XQSMatrix<T>&& rhs)
{
  if (&rhs != this) {
    std::vector<std::vector<T>> tmp(1);
    tmp[0].resize(1);
    tmp.swap(rhs.m_data);
    m_data.swap(tmp);
    m_nrows = rhs.m_nrows;
    m_ncols = rhs.m_ncols;
    rhs.m_nrows = 1;
    rhs.m_ncols = 1;
  }
  return *this;
}

template <typename T1>
XQSMatrix<T1> operator+(const XQSMatrix<T1>& lhs, const XQSMatrix<T1>&rhs)
{
  XQSMatrix<T1> result(lhs);
  result += rhs;
  return result;
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator+=(const XQSMatrix<T>& rhs)
{
  check_equal_dimensions(rhs);
  for (std::size_t i = 0; i < m_nrows; ++i) {
    auto& row = m_data[i];
    const auto& orow = rhs.m_data[i];
    for (std::size_t j = 0; j < m_ncols; ++j) {
      row[j] += orow[j];
    }
  }
  return *this;
}

template <typename T1>
XQSMatrix<T1> operator-(const XQSMatrix<T1>& lhs, const XQSMatrix<T1>& rhs)
{
  XQSMatrix<T1> result(lhs);
  result -= rhs;
  return result;
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator-=(const XQSMatrix<T>& rhs)
{
  check_equal_dimensions(rhs);
  for (std::size_t i = 0; i < m_nrows; ++i) {
    auto& row = m_data[i];
    const auto& orow = rhs.m_data[i];
    for (std::size_t j = 0; j < m_ncols; ++j) {
      row[j] -= orow[j];
    }
  }
  return *this;
}

template <typename T1>
XQSMatrix<T1> operator*(const XQSMatrix<T1>& lhs, const XQSMatrix<T1>& rhs)
{
  lhs.check_suitable_for_product(rhs);
  const auto ncols = rhs.col_count();
  XQSMatrix<T1> result(lhs.m_nrows, ncols, 0.0);
  for (std::size_t i = 0; i < lhs.m_nrows; ++i) {
    const auto& row = lhs.m_data[i];
    for (std::size_t j = 0; j < ncols; ++j) {
      auto& res = result.m_data[i][j];
      for (std::size_t k = 0; k < lhs.m_ncols; ++k) {
        res += row[k] * rhs.m_data[k][j];
      }
    }
  }
  return result;
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator*=(const XQSMatrix<T>& rhs)
{
  XQSMatrix result = (*this) * rhs;
  swap(result);
  return *this;
}

template <typename T1>
XQSMatrix<T1> transpose(const XQSMatrix<T1>& m)
{
  XQSMatrix<T1> result(m.m_ncols, m.m_nrows);
  for (std::size_t i = 0; i < m.m_nrows; ++i) {
    const auto& row = m.m_data[i];
    for (std::size_t j = 0; j < m.m_ncols; ++j) {
      result.m_data[j][i] = row[j];
    }
  }
  return result;
}


// based on the ideas from 
// http://www.sanfoundry.com/java-program-find-inverse-matrix/
template <typename T1>
XQSMatrix<T1> inverse_v1(const XQSMatrix<T1>& m)
{
  if (m.m_nrows != m.m_ncols) {
    throw std::logic_error("Can't invert non-square matrix");
  }

  const T1 zero = 0;
  const std::size_t N = m.m_nrows;
  XQSMatrix<T1> a(m);
  auto index = a.gaussian_reduction();

  // Update the matrix b[i][j] with the ratios stored
  auto b = XQSMatrix<T1>::identity(N);
  for (std::size_t i = 0; i < N - 1; ++i) {
    for (std::size_t j = i + 1; j < N; ++j) {
      const auto& av = a.m_data[index[j]][i];
      for (std::size_t k = 0; k < N; ++k) {
          b.m_data[index[j]][k] -= av * b.m_data[index[i]][k];
      }
    }
  }

  // Perform backward substitutions
  XQSMatrix<T1> x(N, N);
  auto& xrow = x.m_data[N-1];
  auto& arow = a.m_data[index[N-1]];
  auto& brow = b.m_data[index[N-1]];
  const auto& aa = arow[N-1];
  if (aa == zero) {
    throw std::runtime_error("Matrix can't be inverted 3");
  }
  for (std::size_t i = 0; i < N; ++i) {
    xrow[i] = brow[i] / arow[N-1];
    for (std::size_t jj = N-1; jj > 0; --jj) 
    {
      const auto j = jj - 1;
      const auto& ajrow = a[index[j]];
      auto& xji = x[j][i]; 
      xji = b[index[j]][i];
      for (std::size_t k = j + 1; k < N; ++k) 
        xji -= ajrow[k] * x[k][i];
      if (ajrow[j] == zero) {
        throw std::runtime_error("Matrix can't be inverted 4");
      }
      xji /= ajrow[j];
    }
  }
  return x;
}

// Calculate an inverse of this matrix (version #2)
template <typename T1>
XQSMatrix<T1> inverse_v2(const XQSMatrix<T1>& m)
{
  if (m.m_nrows != m.m_ncols) {
    throw std::logic_error("Can't invert non-square matrix");
  }

  const std::size_t N = m.m_nrows;
  XQSMatrix<T1> rm(m);
  auto im = XQSMatrix<T1>::identity(N);
  const T1 zero = 0;
  T1 d;

  for (std::size_t i = 0; i < N - 1; ++i) {
    auto& ri = rm.m_data[i];
    d = ri[i];
    if (d == zero) {
      throw std::logic_error("Matrix can't be inverted 1");
    }
    auto& ii = im.m_data[i];
    for (std::size_t col = 0; col < N; ++col) {
      ri[col] /= d;
      ii[col] /= d;
    }
    for (std::size_t row = i + 1; row < N; ++row) {
      auto& rr = rm.m_data[row];
      auto& ir = im.m_data[row];
      d = rr[i];
      for (std::size_t col = 0; col < N; ++col) {
        rr[col] -= ri[col] * d;
        ir[col] -= ii[col] * d;
      }
    }
  }

  for (std::size_t i = N - 1; i > 0; --i) {
    auto& ri = rm.m_data[i];
    d = ri[i];
    if (d == zero) {
      throw std::logic_error("Matrix can't be inverted 2");
    }
    auto& ii = im.m_data[i];
    for (size_t col = 0; col < N; ++col) {
      ri[col] /= d;
      ii[col] /= d;
    }
    for (size_t row = 0; row < i; ++row) {
      auto& rr = rm.m_data[row];
      auto& ir = im.m_data[row];
      d = rr[i];
      for (size_t col = 0; col < N; ++col) {
        rr[col] -= ri[col] * d;
        ir[col] -= ii[col] * d;
      }
    }
  }

  return im;
}

template <typename T>
std::vector<size_t> XQSMatrix<T>::gaussian_reduction()
{
  const std::size_t N = m_nrows;
  std::vector<double> c(N);

  // Initialize index
  std::vector<size_t> index(N);
  for (std::size_t i = 0; i < N; ++i) 
    index[i] = i;

  // Find the rescaling factors, one from each row
  for (std::size_t i = 0; i < N; ++i) 
  {
    const auto& row = m_data[i];
    double c1 = 0;
    for (std::size_t j = 0; j < N; ++j) 
    {
      auto c0 = std::abs(row[j]);
      if (c0 > c1) c1 = c0;
    }
    c[i] = c1;
  }

  // Search the pivoting element from each column
  std::size_t k = 0;
  for (std::size_t j = 0; j < N - 1 ; ++j) 
  {
    T pi1 = 0;
    for (std::size_t i = j; i < N; ++i) 
    {
      T pi0 = std::abs(m_data[index[i]][j]);
      if (c[index[i]] == 0) {
        throw std::runtime_error("Matrix can't be inverted 5");
      }
      pi0 /= c[index[i]];
      if (pi0 > pi1) {
        pi1 = pi0;
        k = i;
      }
    }

    // Interchange rows according to the pivoting order
    std::swap(index[k], index[j]);
    auto& row0 = m_data[index[j]];
    const auto& v = row0[j];
    if (v == 0) {
      throw std::runtime_error("Matrix can't be inverted 6");
    }

    for (std::size_t i = j + 1; i < N; ++i) 	
    {
      auto& row = m_data[index[i]];
      auto& v2 = row[j];
      auto pj = v2 / v;

      // Record pivoting ratios below the diagonal
      v2 = pj;

      // Modify other elements accordingly
      for (size_t l = j + 1; l < N; ++l) {
        row[l] -= pj * row0[l];
      }
    }
  }
  return index;
}

template <typename T1>
XQSMatrix<T1> operator*(const XQSMatrix<T1>& lhs, const T1& rhs)
{
  XQSMatrix<T1> result(lhs);
  result *= rhs;
  return result;
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator*=(const T& rhs)
{
  for (std::size_t i = 0; i < m_nrows; ++i) {
    auto& row = m_data[i];
    for (std::size_t j = 0; j < m_ncols; ++j) {
      row[j] *= rhs;
    }
  }
  return *this;
}

template <typename T1>
XQSMatrix<T1> operator/(const XQSMatrix<T1>& lhs, const T1& rhs)
{
  XQSMatrix<T1> result(lhs);
  result /= rhs;
  return result;
}

template <typename T>
XQSMatrix<T>& XQSMatrix<T>::operator/=(const T& rhs)
{
  for (std::size_t i = 0; i < m_nrows; ++i) {
    auto& row = m_data[i];
    for (std::size_t j = 0; j < m_ncols; ++j) {
      row[j] /= rhs;
    }
  }
  return *this;
}

template <typename T1>
XQSMatrix<T1> mul_by_row(const XQSMatrix<T1>& lhs, const std::vector<T1>& rhs)
{
  // Validate parameters
  if (rhs.empty()) {
    throw std::logic_error("Empty column data");
  }
  if (lhs.m_ncols != 1) {
    throw std::logic_error(
      "Matrix dimensions mismatch for product with vector row");
  }

  // Compute product
  const auto ncols = rhs.size();
  XQSMatrix<T1> result(lhs.m_nrows, ncols, 0.0);
  for (std::size_t i = 0; i < lhs.m_nrows; ++i) {
    const auto& row = lhs.m_data[i];
    for (std::size_t j = 0; j < ncols; ++j) {
      result.m_data[i][j] = row[j] * rhs[j];
    }
  }
  return result;
}

template <typename T1>
std::vector<T1> mul_by_column(const XQSMatrix<T1>& lhs, const std::vector<T1>& rhs)
{
  // Validate parameters
  if (lhs.m_ncols != rhs.size()) {
    throw std::invalid_argument(
      "Input vector size mismatch for product with vector column");
  }

  // Compute product
  std::vector<T1> result(lhs.m_nrows, 0.0);
  for (std::size_t i = 0; i < lhs.m_nrows; ++i) {
    const auto& row = lhs.m_data[i];
    for (std::size_t j = 0; j < lhs.m_ncols; ++j) {
      result[i] += row[j] * rhs[j];
    }
  }
  return result;
}

template <typename T1>
T1 row_scalar_product(const XQSMatrix<T1>& lhs, std::size_t row_index, const std::vector<T1>& v)
{
  lhs.validate_row_index(row_index);
  if (v.size() != lhs.m_ncols) {
    throw std::invalid_argument(
      "Input vector size mismatch for row scalar product");
  }
  T1 result = 0;
  const auto& row = lhs.m_data[row_index];
  for (std::size_t i = 0; i < lhs.m_ncols; ++i) {
    result += row[i] * v[i];
  }
  return result;
}

template <typename T1>
T1 column_scalar_product(const XQSMatrix<T1>& lhs, std::size_t col_index, const std::vector<T1>& v)
{
  lhs.validate_column_index(col_index);
  if (v.size() != lhs.m_nrows) {
    throw std::invalid_argument(
      "Input vector size mismatch for column scalar product");
  }
  T1 result = 0;
  for (std::size_t i = 0; i < lhs.m_ncols; ++i) {
    result += lhs.m_data[i][col_index] * v[i];
  }
  return result;
}

template <typename T>
void XQSMatrix<T>::add_columns(std::size_t pos, std::size_t count, const T& v)
{
  validate_column_index(pos);
  for (auto& row: m_data) {
    row.insert(row.begin() + pos, count, v);
  }
  m_ncols += count;
}

template <typename T>
void XQSMatrix<T>::remove_columns(std::size_t pos, std::size_t count)
{
  validate_column_index(pos);
  if (count > m_ncols - pos) {
    throw std::out_of_range("Removal count is out of range");
  }
  if (m_ncols == 1) {
    throw std::logic_error(
      "Can't remove column from matrix with single column");
  }
  for (auto& row: m_data) {
    auto it = row.begin() + pos;
    row.erase(it, it + count);
  }
  m_ncols -= count;
}

template <typename T>
void XQSMatrix<T>::fix_to_zero(const T& threshold)
{
  const T zero = 0;
  for (auto& row: m_data) {
    for (T* p = row.data(), *e = p + row.size(); p != e; ++p) {
      if (std::fabs(*p) < threshold) {
        *p = zero;
      }
    }
  }
}

template <typename T>
void XQSMatrix<T>::row_count(std::size_t new_rows)
{
  m_data.resize(new_rows);
  if (new_rows > m_nrows) {
    for (std::size_t i = m_nrows; i < new_rows; ++i) {
      m_data[i].resize(m_ncols);
    }
  }
  m_nrows = new_rows;
}

template <typename T>
void XQSMatrix<T>::col_count(std::size_t new_cols)
{
  for (std::size_t i = 0; i < m_data.size(); ++i)
    m_data[i].resize(new_cols);
  m_ncols = new_cols;
}

template <typename T>
XQSMatrix<T> XQSMatrix<T>::row(std::size_t index) const
{
  validate_row_index(index);
  XQSMatrix<T> result(1, m_ncols);
  result.m_data[0] = m_data[index];
  return result;
}

template <typename T>
std::vector<T> XQSMatrix<T>::row_as_vector(std::size_t index) const
{
  validate_row_index(index);
  return m_data[index];
}

template <typename T>
XQSMatrix<T> XQSMatrix<T>::col(std::size_t index) const
{
  validate_column_index(index);
  XQSMatrix<T> result(m_nrows, 1);
  for (std::size_t i = 0; i < m_nrows; ++i) {
    result.m_data[i][0] = m_data[i][index];
  }
  return result;
}

template <typename T>
std::vector<T> XQSMatrix<T>::col_as_vector(std::size_t index) const
{
  validate_column_index(index);
  std::vector<T> result(m_nrows);
  for (std::size_t i = 0; i < m_nrows; ++i) {
    result[i] = m_data[i][index];
  }
  return result;
}

template <typename T>
XQSMatrix<T> XQSMatrix<T>::window(
  std::size_t row, std::size_t col,
  std::size_t nrows, std::size_t ncols) const
{
  // Validate input parameters
  if (row >= m_nrows) {
    throw std::out_of_range("Row number is out of range");
  }
  if (col >= m_ncols) {
    throw std::out_of_range("Column number is out of range");
  }
  if (nrows == 0 || nrows > m_nrows - row) {
    throw std::out_of_range("Number of window rows is out of range");
  }
  if (ncols == 0 || ncols > m_ncols - col) {
    throw std::out_of_range("Number of window columns is out of range");
  }

  // Build new matrix
  XQSMatrix result(nrows, ncols);
  for (std::size_t i = 0; i < nrows; ++i) {
    auto& rrow = result.m_data[i];
    auto& r = m_data[row + i];
    for (std::size_t j = 0; j < ncols; ++j) {
      rrow[j] = r[col + j];
    }
  }
  return result;
}

// Check that matrix has equal dimensions
template <typename T>
void XQSMatrix<T>::check_equal_dimensions(const XQSMatrix<T>& other) const
{
  if (m_nrows != other.m_nrows && m_ncols != other.m_ncols) {
    std::ostringstream err;
    err << "Dimensions of the other matrix differ "
        "(this vs other (rows*cols): "
      << m_nrows << "*" << m_ncols << " vs " << other.m_nrows << "*" 
      << other.m_ncols << ")"; 
    throw std::invalid_argument(err.str());
  }
}

template <typename T>
void XQSMatrix<T>::check_suitable_for_product(const XQSMatrix<T>& other) const
{
  if (m_ncols != other.m_nrows) {
    std::ostringstream err;
    err << "Dimensions of the other matrix are not suitable for the"
        " product this*other (this vs other (rows*cols): "
      << m_nrows << "*" << m_ncols << " vs " << other.m_nrows << "*"
      << other.m_ncols << ")"; 
    throw std::invalid_argument(err.str());
  }
}

template <typename T>
void XQSMatrix<T>::validate_row_index(std::size_t index) const
{
  if (index >= m_nrows) {
    throw std::out_of_range("Row index is out of range");
  }
}

template <typename T>
void XQSMatrix<T>::validate_column_index(std::size_t index) const
{
  if (index >= m_ncols) {
    throw std::out_of_range("Column index is out of range");
  }
}

template <typename T1, typename Converter>
XQSMatrix<T1> readCsv(const std::string& path, char lineEnding,
  const std::string& fieldDelimiters, const Converter& converter,
  std::size_t numberOfHeaderLines)
{
  XQSMatrix<T1> result(0, 0);

  // Open input file
  std::ifstream in(path.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("Can't open input file");
  }

  // Skip header lines
  std::string line;
  size_t i = numberOfHeaderLines;
  while (i > 0 && std::getline(in, line, lineEnding)) {
    --i;
  }
  if (i > 0) {
    throw std::runtime_error("Missing some header m_nrows");
  }

  // Parse data lines
  size_t numberOfDataLines = 0;
  while (std::getline(in, line, lineEnding)) {
    ++numberOfDataLines;
    auto row = tokenizeAndParse(line, converter, fieldDelimiters);
    if (row.empty()) {
      throw std::
        runtime_error("There is empty data line");
    }
    if (result.m_ncols != row.size()) {
      if (result.m_ncols < row.size()) {
        result.col_count(row.size());
      } else {
        row.resize(result.m_ncols);
      }
    }
    result.m_data.push_back(std::move(row));
    ++result.m_nrows;
  }

  // Ensure that at least one row have been successfully read
  if (numberOfDataLines == 0) {
    throw std::runtime_error("There is no data");
  }

  return result;
}

template<class T, class CharT = char, class Traits>
std::basic_ostream<CharT, Traits>& operator<<(
  std::basic_ostream<CharT, Traits>& os,
  const XQSMatrix<T>& m)
{
  typename std::basic_ostream<CharT, Traits>::sentry sentry(os);
  const auto nrows = m.row_count();
  const auto ncols = m.col_count(); 
  const auto& v = m.data();
  for (std::size_t i = 0; i < nrows; ++i) {
    const auto& row = v[i];
    os << row[0];
    for (std::size_t j = 1; j < ncols; ++j) {
      os << '\t' << row[j];
    }
    os << '\n';
  }
  return os;
}

#endif // XQSMATRIX_H__
