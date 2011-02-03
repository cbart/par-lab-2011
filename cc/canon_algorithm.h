// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__ALGORITHM__H__
#define __CANON__ALGORITHM__H__


#include <algorithm>
#include <boost/function.hpp>
#include "matrix.h"
#include "mpi.h"


namespace canon
{
namespace algorithm
{


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
class canon_prod
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
private:
    const communicator_type & cart_2d;
    product_function_type local_product;
    ranks_array_type vertical_ranks;
    ranks_array_type horizontal_ranks;
    row_matrix_type * result_partial;
    row_matrix_type * first_argument_partial_current;
    col_matrix_type * second_argument_partial_current;
    row_matrix_type * first_argument_partial_temp;
    col_matrix_type * second_argument_partial_temp;
public:
    canon_prod(const communicator_type & cart_2d, product_function_type local_product)
        throw();
    ~canon_prod()
        throw();
    void operator()(
            row_matrix_type & product_result_partial,
            row_matrix_type & product_first_argument_partial,
            col_matrix_type & product_second_argument_partial)
        throw();
private:
    void init_partials(
            row_matrix_type & product_result_partial,
            row_matrix_type & product_first_argument_partial,
            col_matrix_type & product_second_argument_partial)
        throw();
    void clear_partials()
        throw();
    void swap_partials()
        throw();
};


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
canon_prod<real_t, SIZE, DIM_SIZE>::canon_prod(
        const communicator_type & cart_2d,
        product_function_type local_product)
  : cart_2d(cart_2d),
    local_product(local_product),
    vertical_ranks(
            mpi::shift<mpi::DIRECTION_VERTICAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d)),
    horizontal_ranks(
            mpi::shift<mpi::DIRECTION_HORIZONTAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d)),
    result_partial(NULL),
    first_argument_partial_current(NULL),
    second_argument_partial_current(NULL),
    first_argument_partial_temp(NULL),
    second_argument_partial_temp(NULL)
{
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
canon_prod<real_t, SIZE, DIM_SIZE>::~canon()
{
    delete result_partial;
    delete first_argument_partial_current;
    delete second_argument_partial_current;
    delete first_argument_partial_temp;
    delete second_argument_partial_temp;
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void canon_prod<real_t, SIZE, DIM_SIZE>::operator()(
        row_matrix_type & product_result_partial,
        row_matrix_type & product_first_argument_partial,
        col_matrix_type & product_second_argument_partial)
    throw()
{
    init_partials(
            product_result_partial,
            product_first_argument_partial,
            product_second_argument_partial);
    local_product(
            product_result_partial,
            product_first_argument_partial,
            product_second_argument_partial);
    clear_partials();
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void canon_prod<real_t, SIZE, DIM_SIZE>::init_partials(
        row_matrix_type & product_result_partial,
        row_matrix_type & product_first_argument_partial,
        col_matrix_type & product_second_argument_partial)
    throw()
{
    result_partial = & product_result_partial;
    first_argument_partial_current = & product_first_argument_partial;
    second_argument_partial_current = & product_second_argument_partial;
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void canon_prod<real_t, SIZE, DIM_SIZE>::clear_partials()
    throw()
{
    result_partial = NULL;
    first_argument_partial_current = NULL;
    second_argument_partial_current = NULL;
}


template<typename real_t, size_t SIZE, size_t DIM_SIZE>
inline void canon_prod<real_t, SIZE, DIM_SIZE>::swap_partials()
    throw()
{
    ::std::swap(first_argument_partial_current, first_argument_partial_temp);
    ::std::swap(second_argument_partial_current, second_argument_partial_temp);
}


}  // namespace algorithm
}  // namespace canon


#endif
