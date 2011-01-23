// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#include <iostream>
#include <boost/mpi.hpp>
#include "matrix.h"
#include "random.h"
#include "fill.h"
#include "multiply.h"
#include "constant.h"
#include "mpi.h"

const size_t SIZE = 1000L;
const size_t DIM_SIZE = 4;

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
    ::canon::fill<real_type, ::canon::row_major, SIZE>(c, & ::canon::constant<real_type, 0>);
    ::canon::prod<real_type, SIZE>(c, a, b);
    ::std::cerr << a(1,1) << ::std::endl;
}

int main(int argc, char * * argv)
{
    ::boost::mpi::environment env(argc, argv);
    ::boost::mpi::communicator cart_2d = ::canon::mpi::cart_square_sphere_create<DIM_SIZE>();
    ::canon::mpi::assert_processors<DIM_SIZE * DIM_SIZE>(cart_2d, env);
    test();
    return 0;
}
