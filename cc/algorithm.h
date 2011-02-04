// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__ALGORITHM__H__
#define __CANNON__ALGORITHM__H__


#include <algorithm>
#include <boost/function.hpp>
#include <boost/array.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>
#include <boost/mpi/nonblocking.hpp>
#include "matrix.h"
#include "mpi.h"


namespace cannon
{
namespace algorithm
{


// The Cannon's multiply algorithm signature
//   `real_t` The actual matrix element type
//   `SIZE` The size of a single partial, e.g.
//       Cannon algorithm will shuffle SIZE*SIZE
//       square matrices.
//   `CART_SIZE` The algorithm will work in
//       cartesian topology (the grid will have
//       CART_SIZE*CART_SIZE nodes) with periods
//       in both dimensions.
template<typename real_t, size_t SIZE, size_t CART_SIZE>
class cannon_prod
{
public:
    typedef real_t real_type;
    // Row-major matrix type
    typedef square_matrix<real_type, row_major, SIZE> row_matrix_concept;
    typedef typename row_matrix_concept::type row_matrix_type;
    // Column-major matrix type
    typedef square_matrix<real_type, col_major, SIZE> col_matrix_concept;
    typedef typename col_matrix_concept::type col_matrix_type;
    // Local multiplication function
    typedef ::boost::function<
        void (
                row_matrix_type & product_result,
                row_matrix_type & product_first_argument,
                col_matrix_type & product_second_argument)
        throw()> product_function_type;
    // MPI Communicator - boost wrapped
    typedef ::boost::mpi::communicator communicator_type;
private:
    typedef mpi::ranks_array_type ranks_array_type;
    typedef ::boost::mpi::request mpi_request_type;
private:
    const communicator_type & cart_2d;
    const product_function_type local_product;
    const ranks_array_type vertical_ranks;
    const ranks_array_type horizontal_ranks;
    row_matrix_type * result;
    row_matrix_type * left_current;
    col_matrix_type * right_current;
    row_matrix_type * left_temp;
    col_matrix_type * right_temp;
    row_matrix_type * row_temp;
    col_matrix_type * col_temp;
    mutable ::boost::array<mpi_request_type, 4> mpi_requests;
public:
    // Creates the algorithm framework, can reuse
    // the `row_temp` and `col_temp`
    // through sequential runs.
    cannon_prod(
            const communicator_type & cart_2d,
            product_function_type local_product,
            row_matrix_type & row_temp,
            col_matrix_type & col_temp)
        throw();
    ~cannon_prod()
        throw();
    // Performs the multiplication.
    //   result += first * second
    void operator()(
            row_matrix_type & result,
            row_matrix_type & left,
            col_matrix_type & right)
        throw();
private:
    void init_partials(
            row_matrix_type & result,
            row_matrix_type & left,
            col_matrix_type & right)
        throw();
    void swap_partials()
        throw();
    void ishift_partials()
        throw();
    void wait()
        throw();
};


template<typename real_t, size_t SIZE, size_t CART_SIZE>
cannon_prod<real_t, SIZE, CART_SIZE>::cannon_prod(
        const communicator_type & cart_2d,
        product_function_type local_product,
        row_matrix_type & row_temp,
        col_matrix_type & col_temp)
  : cart_2d(cart_2d),
    local_product(local_product),
    vertical_ranks(mpi::shift<mpi::DIRECTION_VERTICAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d)),
    horizontal_ranks(mpi::shift<mpi::DIRECTION_HORIZONTAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d)),
    result(NULL),
    left_current(NULL),
    right_current(NULL),
    left_temp(NULL),
    right_temp(NULL),
    row_temp(& row_temp),
    col_temp(& col_temp)
{
}


template<typename real_t, size_t SIZE, size_t CART_SIZE>
cannon_prod<real_t, SIZE, CART_SIZE>::~cannon()
{
}


template<typename real_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, SIZE, CART_SIZE>::operator()(
        row_matrix_type & result,
        row_matrix_type & left,
        col_matrix_type & right)
    throw()
{
    init_partials(result, left, right);
    for(uint32_t step = 0; step + 1 < CART_SIZE; ++step)
    {
        ishift_partials();
        local_product(* this->result, * left_current, * right_current);
        wait();
        swap_partials();
    }
    local_product(* this->result, * left_current, * right_current);
}


namespace
{


static const int CANNON_ALGORITHM_MPI_TAG = 42;


}  // namespace (unnamed)


template<typename real_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, SIZE, CART_SIZE>::ishift_partials()
    throw()
{
    mpi_requests[0] = cart_2d.isend(
            vertical_ranks[mpi::DESTINATION_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            row_matrix_concept::begin(left_current),
            SIZE * SIZE);
    mpi_requests[1] = cart_2d.irecv(
            vertical_ranks[mpi::SOURCE_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            row_matrix_concept::begin(left_temp),
            SIZE * SIZE);
    mpi_requests[2] = cart_2d.isend(
            horizontal_ranks[mpi::DESTINATION_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            col_matrix_concept::begin(right_current),
            SIZE * SIZE);
    mpi_requests[3] = cart_2d.irecv(
            vertical_ranks[mpi::SOURCE_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            col_matrix_concept::begin(right_temp),
            SIZE * SIZE);
}


template<typename real_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, SIZE, CART_SIZE>::wait()
    throw()
{
    ::boost::mpi::wait_all(mpi_requests.begin(), mpi_requests.end());
}


template<typename real_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, SIZE, CART_SIZE>::init_partials(
        row_matrix_type & result,
        row_matrix_type & left,
        col_matrix_type & right)
    throw()
{
    this->result = & result;
    left_current = & left;
    right_current = & right;
    left_temp = row_temp;
    right_temp = col_temp;
}


template<typename real_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, SIZE, CART_SIZE>::swap_partials()
    throw()
{
    ::std::swap(left_current, left_temp);
    ::std::swap(right_current, right_temp);
}


}  // namespace algorithm
}  // namespace cannon


#endif
