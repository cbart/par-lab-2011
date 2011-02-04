// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef DEBUGLEVEL
#define DEBUGLEVEL 3
#endif

#include <boost/mpi.hpp>
#include "matrix.h"
#include "random.h"
#include "fill.h"
#include "multiply.h"
#include "constant.h"
#include "mpi.h"
#include "exceptions.h"
#include "debug.h"
#include "canon_algorithm.h"

// The size of a single parial matrix.
const size_t SIZE = 500L;

// The MPI cart will be a `DIM_SIZE` x `DIM_SIZE` square.
const size_t DIM_SIZE = 4;

// The type we work with.
typedef double real_type;

// The algorithm we work with.
typedef ::canon::algorithm::canon_prod<real_type, SIZE, DIM_SIZE> canon_prod_type;

// The maintenance function.
int run_product(
        const ::boost::mpi::communicator & cart_2d,
        const canon_prod_type::product_function_type local_product);


int main(int argc, char * * argv)
{
    using namespace ::canon;
    ::boost::mpi::environment env(argc, argv);
    ::boost::mpi::communicator cart_2d = ::canon::mpi::cart_square_sphere_create<DIM_SIZE>();
    mpi::assert_processors<DIM_SIZE * DIM_SIZE>(cart_2d, env);
    int error_code = run_product(cart_2d, prod<real_type, SIZE>);
    return error_code;
}


inline int run_product(
        const ::boost::mpi::communicator & cart_2d,
        const canon_prod_type::product_function_type local_product)
{
    using namespace ::canon;
    // Create matrices.
    canon_prod_type::row_matrix_type left(SIZE, SIZE);
    canon_prod_type::col_matrix_type right(SIZE, SIZE);
    canon_prod_type::row_matrix_type result(SIZE, SIZE);
    canon_prod_type::row_matrix_type row_temp(SIZE, SIZE);
    canon_prod_type::col_matrix_type col_temp(SIZE, SIZE);
    // Generate pseudo-random values.
    random_generator<real_type> generator;
    fill<real_type, row_major, SIZE>(left, generator);
    fill<real_type, col_major, SIZE>(right, generator);
    fill<real_type, row_major, SIZE>(result, & constant<real_type, 0>);
    // Initiate the algorithm.
    canon_prod_type canon_product(cart_2d, local_product, row_temp, col_temp);
    // Run the algorithm.
    canon_product(result, left, right);
    return 0;
}
