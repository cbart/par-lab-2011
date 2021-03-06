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
#include <boost/numeric/ublas/storage.hpp>
#include "matrix.h"
#include "mpi.h"
#include "debug.h"


namespace cannon
{
namespace algorithm
{


namespace
{


static const int CANNON_ALGORITHM_MPI_TAG = 42;


}  // namespace (unnamed)


// The Cannon's multiply algorithm signature
//   `real_t` The actual matrix element type
//   `SIZE` The size of a single partial, e.g.
//       Cannon algorithm will shuffle SIZE*SIZE
//       square matrices.
//   `CART_SIZE` The algorithm will work in
//       cartesian topology (the grid will have
//       CART_SIZE*CART_SIZE nodes) with periods
//       in both dimensions.
template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
class cannon_prod
{
public:
    typedef real_t real_type;
    typedef storage_t storage_type;
    // Row-major matrix type
    typedef square_matrix_concept<real_type, storage_type, row_major, SIZE> row_matrix_concept;
    typedef typename row_matrix_concept::type row_matrix_type;
    // Column-major matrix type
    typedef square_matrix_concept<real_type, storage_type, col_major, SIZE> col_matrix_concept;
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
    typedef ::boost::array<mpi_request_type, 2 * mpi::DIMS> mpi_request_array_type;
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
    mutable mpi_request_array_type mpi_requests;
public:
    // Creates the algorithm framework,
    // reuses `row_temp` and `col_temp`
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
    // Assigns current and temp inner
    // algorithm's pointers.
    void init_partials(
            row_matrix_type & result,
            row_matrix_type & left,
            col_matrix_type & right)
        throw();
    // Performs initial realignment of the partials
    // according to the Cannon's algorithm.
    void align_partials()
        throw();
    // Realigns the partials after the actual algorithm
    // so that the data is on its correct place.
    void realign_partials()
        throw();
    // Swaps temp with current pointers.
    void swap_partials()
        throw();
    // Performs nonblocking send on current pointers
    // and nonblocing receive on temp pointers
    // according to Cannon's algorithm.
    void ishift_partials()
        throw();
    // Waits for all requests performed by `ishift_partials`.
    void wait()
        throw();
    // Partial send.
    template<typename matrix_t>
    mpi_request_type isend(mpi::rank_type destination, matrix_t * matrix)
        throw()
    {
        return cart_2d.isend(destination, CANNON_ALGORITHM_MPI_TAG, begin(matrix), SIZE * SIZE);
    }
    // Patrial receive.
    template<typename matrix_t>
    mpi_request_type irecv(mpi::rank_type source, matrix_t * matrix)
        throw()
    {
        return cart_2d.irecv(source, CANNON_ALGORITHM_MPI_TAG, begin(matrix), SIZE * SIZE);
    }
    // Returns pointer to the memory a matrix is stored.
    // Placeholder
    template<typename matrix_t>
    real_type * begin(matrix_t * matrix)
        throw()
    {
        ::debug::info << "DUPA" << ::std::endl;
        return NULL;
    }
    // "Specializations" (overloaded actually)
    real_type * begin(row_matrix_type * matrix)
        throw();
    real_type * begin(col_matrix_type * matrix)
        throw();
};




template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::cannon_prod(
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


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::~cannon()
{
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::operator()(
        row_matrix_type & result,
        row_matrix_type & left,
        col_matrix_type & right)
    throw()
{
    init_partials(result, left, right);
    //align_partials();
    for(uint32_t step = 0; step + 1 < CART_SIZE; ++step)
    {
        ::debug::info << "Begin iteration " << step + 1 << ".\n" << ::std::flush;
        ishift_partials();
        ::debug::info << "Begin product.\n" << ::std::flush;
        local_product(* this->result, * left_current, * right_current);
        ::debug::info << "Waiting for exchange.\n" << ::std::flush;
        wait();
        ::debug::info << "Swapping partials.\n" << ::std::flush;
        swap_partials();
    }
    local_product(* this->result, * left_current, * right_current);
    //realign_partials();
}




template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::align_partials()
{
    mpi::coords_type coords = mpi::coords(cart_2d);
    if(coords[mpi::DIRECTION_VERTICAL] != 0)
    {
        const int step = coords[mpi::DIRECTION_VERTICAL];
        ranks_array_type align_ranks =
            mpi::shift<mpi::DIRECTION_VERTICAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d, step);
        mpi_requests[mpi::DIMS * mpi::DIRECTION_VERTICAL + 0] =
            isend(align_ranks[mpi::DESTINATION_RANK_INDEX], left_current);
        mpi_requests[mpi::DIMS * mpi::DIRECTION_VERTICAL + 1] =
            irecv(align_ranks[mpi::SOURCE_RANK_INDEX], left_temp);
    }
    if(coords[mpi::DIRECTION_HORIZONTAL] != 0)
    {
        const int step = coords[mpi::DIRECTION_HORIZONTAL];
        ranks_array_type align_ranks =
            mpi::shift<mpi::DIRECTION_HORIZONTAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d, step);
        mpi_requests[mpi::DIMS * mpi::DIRECTION_HORIZONTAL + 0] =
            isend(align_ranks[mpi::DESTINATION_RANK_INDEX], right_current);
        mpi_requests[mpi::DIMS * mpi::DIRECTION_HORIZONTAL + 1] =
            irecv(align_ranks[mpi::SOURCE_RANK_INDEX], right_temp);
    }
    if(coords[mpi::DIRECTION_VERTICAL] != 0)
    {
        const size_t offset = mpi::DIMS * mpi::DIRECTION_VERTICAL;
        ::boost::mpi::wait_all(
                mpi_requests.begin() + offset,
                mpi_requests.begin() + offset + mpi::DIMS);
        ::std::swap(left_current, left_temp);
    }
    if(coords[mpi::DIRECTION_HORIZONTAL] != 0)
    {
        const size_t offset = mpi::DIMS * mpi::DIRECTION_HORIZONTAL;
        ::boost::mpi::wait_all(
                mpi_requests.begin() + offset,
                mpi_requests.begin() + offset + mpi::DIMS);
        ::std::swap(right_current, right_temp);
    }
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::realign_partials()
{
    // TODO
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::ishift_partials()
    throw()
{
    mpi_requests[0] = isend(vertical_ranks[mpi::DESTINATION_RANK_INDEX], left_current);
    mpi_requests[1] = irecv(vertical_ranks[mpi::SOURCE_RANK_INDEX], left_temp);
    mpi_requests[2] = isend(horizontal_ranks[mpi::DESTINATION_RANK_INDEX], right_current);
    mpi_requests[3] = irecv(horizontal_ranks[mpi::SOURCE_RANK_INDEX], right_temp);
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::wait()
    throw()
{
    ::boost::mpi::wait_all(mpi_requests.begin(), mpi_requests.end());
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::init_partials(
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


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
inline void cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::swap_partials()
    throw()
{
    ::std::swap(left_current, left_temp);
    ::std::swap(right_current, right_temp);
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
typename cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::real_type *
cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::begin(row_matrix_type * matrix)
    throw()
{
    return row_matrix_concept::begin(matrix);
}


template<typename real_t, typename storage_t, size_t SIZE, size_t CART_SIZE>
typename cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::real_type *
cannon_prod<real_t, storage_t, SIZE, CART_SIZE>::begin(col_matrix_type * matrix)
    throw()
{
    return col_matrix_concept::begin(matrix);
}


}  // namespace algorithm
}  // namespace cannon


#endif
