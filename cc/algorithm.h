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


const int CANNON_ALGORITHM_MPI_TAG = 42;


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
class cannon_prod
{
public:
    typedef real_t real_type;
    typedef typename square_matrix<real_type, row_major, SIZE>::type row_matrix_type;
    typedef typename square_matrix<real_type, col_major, SIZE>::type col_matrix_type;
    typedef ::boost::function<
        void (
                row_matrix_type & product_result,
                row_matrix_type & product_first_argument,
                col_matrix_type & product_second_argument)
        throw()> product_function_type;
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


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
cannon_prod<real_t, SIZE, DIM_SIZE>::cannon_prod(
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


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
cannon_prod<real_t, SIZE, DIM_SIZE>::~cannon()
{
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void cannon_prod<real_t, SIZE, DIM_SIZE>::operator()(
        row_matrix_type & result,
        row_matrix_type & left,
        col_matrix_type & right)
    throw()
{
    init_partials(result, left, right);
    for(uint32_t step = 0; step + 1 < DIM_SIZE; ++step)
    {
        ishift_partials();
        local_product(* this->result, * left_current, * right_current);
        wait();
        swap_partials();
    }
    local_product(* this->result, * left_current, * right_current);
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void cannon_prod<real_t, SIZE, DIM_SIZE>::ishift_partials()
    throw()
{
    mpi_requests[0] = cart_2d.isend(
            vertical_ranks[mpi::DESTINATION_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            left_current->data().begin(),
            SIZE * SIZE);
    mpi_requests[1] = cart_2d.irecv(
            vertical_ranks[mpi::SOURCE_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            left_temp->data().begin(),
            SIZE * SIZE);
    mpi_requests[2] = cart_2d.isend(
            horizontal_ranks[mpi::DESTINATION_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            right_current->data().begin(),
            SIZE * SIZE);
    mpi_requests[3] = cart_2d.irecv(
            vertical_ranks[mpi::SOURCE_RANK_INDEX],
            CANNON_ALGORITHM_MPI_TAG,
            right_temp->data().begin(),
            SIZE * SIZE);
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void cannon_prod<real_t, SIZE, DIM_SIZE>::wait()
    throw()
{
    ::boost::mpi::wait_all(mpi_requests.begin(), mpi_requests.end());
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void cannon_prod<real_t, SIZE, DIM_SIZE>::init_partials(
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


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void cannon_prod<real_t, SIZE, DIM_SIZE>::swap_partials()
    throw()
{
    ::std::swap(left_current, left_temp);
    ::std::swap(right_current, right_temp);
}


}  // namespace algorithm
}  // namespace cannon


#endif
