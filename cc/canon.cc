// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#define DEBUGLEVEL 3

#include <iostream>
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
const size_t SIZE = 1000L;

// The MPI cart will be a `DIM_SIZE` x `DIM_SIZE` square.
const size_t DIM_SIZE = 4;

// The type we work with.
typedef double real_type;

// The algorithm we work with.
typedef ::canon::algorithm::canon_prod<real_type, SIZE> canon_prod_type;

// The maintenance function.
void run_product(
        const ::boost::mpi::communicator & cart_2d,
        const canon_prod_type::product_function_type local_product);


int main(int argc, char * * argv)
{
    using namespace ::canon;
    ::boost::mpi::environment env(argc, argv);
    ::boost::mpi::communicator cart_2d = ::canon::mpi::cart_square_sphere_create<DIM_SIZE>();
    mpi::assert_processors<DIM_SIZE * DIM_SIZE>(cart_2d, env);
    run_product(cart_2d, prod<real_type, SIZE>);
    return 0;
}


inline void run_product(
        const ::boost::mpi::communicator & cart_2d,
        const canon_prod_type::product_function_type local_product)
{
    using namespace ::canon;
    // Create matrices.
    canon_prod_type::row_matrix_type a(SIZE, SIZE);
    canon_prod_type::col_matrix_type b(SIZE, SIZE);
    canon_prod_type::row_matrix_type c(SIZE, SIZE);
    // Generate pseudo-random values.
    random_generator<real_type> generator;
    fill<real_type, row_major, SIZE>(a, generator);
    fill<real_type, col_major, SIZE>(b, generator);
    fill<real_type, row_major, SIZE>(c, & constant<real_type, 0>);
    // Initiate the algorithm.
    algorithm::canon_prod<real_type, SIZE> canon_product(cart_2d, local_product);
    // Run the algorithm.
    canon_product(c, a, b);
}
