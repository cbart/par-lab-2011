// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#define DEBUGLEVEL 3

#include <iostream>
#include "matrix.h"
#include "random.h"
#include "fill.h"
#include "multiply.h"

const size_t SIZE = 1000L;

typedef double real_type;

typedef ::canon::square_matrix<real_type, ::canon::row_major, SIZE>::type r_matrix;
typedef ::canon::square_matrix<real_type, ::canon::col_major, SIZE>::type c_matrix;

inline void test()
{
    r_matrix a(SIZE, SIZE);
    c_matrix b(SIZE, SIZE);
    r_matrix c(SIZE, SIZE);
    ::canon::random_generator<real_type> generator;
    ::canon::fill<real_type, ::canon::row_major, SIZE>(a, generator);
    ::canon::fill<real_type, ::canon::col_major, SIZE>(b, generator);
    ::canon::prod<real_type, SIZE>(c, a, b);
    ::std::cerr << a(1,1) << ::std::endl;
}

int main(const int argc, char const * const * const argv)
{
    test();
    return 0;
}
