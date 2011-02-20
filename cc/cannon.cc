// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl


#ifndef DEBUGLEVEL
#define DEBUGLEVEL 0
#endif


#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include "matrix.h"
#include "random.h"
#include "fill.h"
#include "multiply.h"
#include "constant.h"
#include "mpi.h"
#include "exceptions.h"
#include "algorithm.h"
#include "debug.h"


#ifndef MATRIXSIZE
#define MATRIXSIZE 65536l
#endif

// The size of a single parial matrix.
const size_t SIZE = MATRIXSIZE/CARTSIZE;

#ifndef CARTSIZE
#define CARTSIZE 8l
#endif

// The MPI cart will be a `CART_SIZE` x `CART_SIZE` square.
const size_t CART_SIZE = CARTSIZE;


// The type we work with.
typedef double real_type;

// Storage type - we use unbounded_array because of performance reasons
// typedef ::boost::numeric::ublas::unbounded_array<real_type> storage_type;
typedef ::std::vector<real_type> storage_type;

// The algorithm we work with.
typedef ::cannon::algorithm::cannon_prod<real_type, storage_type, SIZE, CART_SIZE> cannon_prod_type;


// The maintenance function.
int run_product(
        const ::boost::mpi::communicator & cart_2d,
        const cannon_prod_type::product_function_type local_product);


int main(int argc, char * * argv)
{
    ::debug::info << "Setting MPI environment..." << ::std::endl;
    ::boost::mpi::environment env(argc, argv);
    ::debug::info << "Setting the cartesian communicator..." << ::std::endl;
    ::boost::mpi::communicator cart_2d = ::cannon::mpi::cart_square_sphere_create<CART_SIZE>();
    ::debug::info << "Checking amount of processors..." << ::std::endl;
    ::cannon::mpi::assert_processors<CART_SIZE * CART_SIZE>(cart_2d, env);
    int error_code = run_product(cart_2d, ::cannon::prod<real_type, storage_type, SIZE>);
    return error_code;
}


inline int run_product(
        const ::boost::mpi::communicator & cart_2d,
        const cannon_prod_type::product_function_type local_product)
{
    using namespace ::cannon;
    ::debug::info << "Creating matrices..." << ::std::endl;
    cannon_prod_type::row_matrix_type left(SIZE, SIZE);
    cannon_prod_type::col_matrix_type right(SIZE, SIZE);
    cannon_prod_type::row_matrix_type result(SIZE, SIZE);
    cannon_prod_type::row_matrix_type row_temp(SIZE, SIZE);
    cannon_prod_type::col_matrix_type col_temp(SIZE, SIZE);
    // Generate pseudo-random values.
    ::debug::info << "Filling matrices with random values..." << ::std::endl;
    random_generator<real_type> generator;
    fill(left, generator);
    fill(right, generator);
    fill(result, & constant<real_type, 0>);
    // Initiate the algorithm.
    ::debug::info << "Initiating the algorithm..." << ::std::endl;
    cannon_prod_type cannon_product(cart_2d, local_product, row_temp, col_temp);
    // Run the algorithm.
    ::debug::info << "Running the algorithm..." << ::std::endl;
    cannon_product(result, left, right);
    return 0;
}
