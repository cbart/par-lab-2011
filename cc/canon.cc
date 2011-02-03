// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#define DEBUGLEVEL 3

#include <iostream>
#include <boost/mpi.hpp>
#include <boost/function.hpp>
#include "matrix.h"
#include "random.h"
#include "fill.h"
#include "multiply.h"
#include "constant.h"
#include "mpi.h"
#include "exceptions.h"
#include "debug.h"

const size_t SIZE = 1000L;
const size_t DIM_SIZE = 4;

typedef double real_type;

typedef ::canon::square_matrix<real_type, ::canon::row_major, SIZE>::type row_matrix_type;
typedef ::canon::square_matrix<real_type, ::canon::col_major, SIZE>::type col_matrix_type;
typedef row_matrix_type & row_matrix_ref;
typedef col_matrix_type & col_matrix_ref;


typedef ::boost::function<void (row_matrix_ref product_result, row_matrix_ref product_first_argument, col_matrix_ref product_second_argument) throw()> local_product_function_type;

inline void canon_mpi_algorithm(
        const ::boost::mpi::communicator & cart_2d,
        const local_product_function_type local_product)
{
    using namespace ::canon;
    row_matrix_type a(SIZE, SIZE);
    col_matrix_type b(SIZE, SIZE);
    row_matrix_type c(SIZE, SIZE);
    random_generator<real_type> generator;
    fill<real_type, row_major, SIZE>(a, generator);
    fill<real_type, col_major, SIZE>(b, generator);
    fill<real_type, row_major, SIZE>(c, & constant<real_type, 0>);
    mpi::ranks_array_type vertical_ranks = mpi::shift<mpi::DIRECTION_VERTICAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d);
    mpi::ranks_array_type horizontal_ranks = mpi::shift<mpi::DIRECTION_HORIZONTAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d);
    local_product(c, a, b);
}

int main(int argc, char * * argv)
{
    using namespace ::canon;
    ::boost::mpi::environment env(argc, argv);
    ::boost::mpi::communicator cart_2d = ::canon::mpi::cart_square_sphere_create<DIM_SIZE>();
    mpi::assert_processors<DIM_SIZE * DIM_SIZE>(cart_2d, env);
    canon_mpi_algorithm(cart_2d, prod<real_type, SIZE>);
    return 0;
}
