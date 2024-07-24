// Copyright © 2018, 2024 Ivan Pizhenko. All rights reserved.
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

#include <iostream>
#include <random>
#include "xqsmatrix.h"

int main()
{
  std::mt19937 rng(111);
  std::uniform_real_distribution<double> d(0, 1.0);
  
  constexpr std::size_t N = 10;
  XQSMatrix<double> m(N, N);
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      m(i, j) = d(rng);
    }
  } 

  const auto mm1 = inverse_v1(m);
  const auto pr1 = m * mm1;

  const auto mm2 = inverse_v2(m);
  const auto pr2 = m * mm2;

  std::cout << pr1 << std::endl << std::endl;
  std::cout << pr2 << std::endl << std::endl;

  return 0;
}
