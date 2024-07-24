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

#ifndef TRANSFORM_IF_H__
#define TRANSFORM_IF_H__

namespace stdx {

template< class InputIt, class OutputIt, class UnaryOperation, class UnaryPredicate>
OutputIt transform_if(
  InputIt first, InputIt last, OutputIt d_first,
  UnaryOperation unary_op, UnaryPredicate pred)
{
  for (; first != last; ++first) {
    if (pred(*first)) {
      *d_first = unary_op(*first);
      ++d_first; 
    }
  }
  return d_first;
}

} // namespace stdx

#endif // TRANSFORM_IF_H__
