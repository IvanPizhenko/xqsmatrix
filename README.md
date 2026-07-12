# xqs_matrix

## Overview

`xqs_matrix` ("eXtended QSMatrix") is 2D matrix implementation,
inspired by ideas described in the following articles:

- [Matrix Classes in C++ - The Header File](https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File)
- [Matrix Classes in C++ - The Source File](https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Source-File)

Current implementation is completely reworked and share no common code with the
original `QSMatrix` class. The original class now serves only as "inspiration"
and source of some ideas.

## Implementation Progress

1. Create initial header file.
   - [x] Create header file with include guards and license.

2. Class declaration:
   - [x] Declare template class xqs_matrix<T, Alloc>.
   - [x] Add public type aliases.
   - [x] Add private member variables.

3. Constructors:
   - [x] Default constructor.
   - [x] Constructor with row and column counts.
   - [x] Constructor with row and column counts and initial value.
   - [x] Constructor with allocator.
   - [x] Constructor with allocator, row and column counts.
   - [x] Constructor with allocator, row and column counts
         and initial value.
   - [x] Constructor from range defined by iterators.
   - [x] Constructor from range object.
   - [x] Copy constructor.
   - [x] Copy constructor with different allocator object.
   - [x] Move constructor.
   - [x] Move constructor with different allocator object.
   - [x] Constructor from raw data pointer (view).
   - [x] Constructor from initializer list.

4. Destructor:
   - [x] ~xqs_matrix()

5. Accessors:
   - [x] row_count()
   - [x] column_count()
   - [x] size()
   - [x] empty()
   - [x] stride()
   - [x] data()
   - [x] capacity()
   - [x] is_owner()
   - [x] is_view()

6. Row and column access:
   - [x] operator[](row)
   - [x] row_at(row)
   - [x] operator()(column)
   - [x] column_at(column)

7. Element access:
   - [x] operator()(row, column)
   - [x] operator()(pair<row, column>)
   - [x] operator[](pair<row, column>)
   - [x] at(row, column)
   - [x] at(pair<row, column>)

8. Assignment operators:
   - [x] Copy assignment
   - [x] Move assignment

9. Swap method.
   - [x] swap() member function
   - [x] swap() external function

10. Resize method.
    - [x] resize() member function
    - [x] resize_rows() external function
    - [x] resize_columns() external function
    - [x] reserve() member function

11. Member arithmetic operators:
    - [x] operator*=(scalar)
    - [x] operator/=(scalar)
    - [x] operator+=(matrix)
    - [x] operator-=(matrix)
    - [x] operator*=(matrix)

12. External arithmetic operators:
    - [x] friend operator*(matrix, scalar)
    - [x] friend operator/(matrix, scalar)
    - [x] friend operator+(matrix, matrix)
    - [x] friend operator-(matrix, matrix)
    - [x] friend operator*(matrix, matrix)

13. Special matrix factory functuons:
    - [x] identity()
    - [x] transposed_identity()

14. Matrix operations:
    - [x] clear()
    - [x] shrink_to_fit()
    - [x] copy_full()
    - [x] move_full()
    - [x] window()
    - [x] copy_window()
    - [x] move_window()
    - [x] transpose_copy()
    - [x] transpose_move()
    - [x] diag_to_hvec()
    - [x] diag_to_vvec()
    - [x] inverse_v1()
    - [x] inverse_v2()
    - [x] gaussian_reduction()

15. Stream operations:
    - [x] operator<<()
    - [x] operator>>()

16. Additional I/O operations:
    - [x] read_csv()
