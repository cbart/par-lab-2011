// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__ALGORITHM__H__
#define __CANON__ALGORITHM__H__


#include <boost/function.hpp>
#include "matrix.h"
#include "mpi.h"


namespace canon
{
namespace algorithm
{


template<typename real_t, size_t SIZE>
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
public:
    canon_prod(const communicator_type & cart_2d, product_function_type local_product);
    ~canon_prod();
    void operator()(
            row_matrix_type & product_result_partial,
            row_matrix_type & product_first_argument_partial,
            col_matrix_type & product_second_argument_partial)
        throw();
};


template<typename real_t, size_t SIZE>
canon_prod<real_t, SIZE>::canon_prod(
        const communicator_type & cart_2d,
        product_function_type local_product)
  : cart_2d(cart_2d),
    local_product(local_product),
    vertical_ranks(
            mpi::shift<mpi::DIRECTION_VERTICAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d)),
    horizontal_ranks(
            mpi::shift<mpi::DIRECTION_HORIZONTAL, mpi::DISPLACEMENT_DOWNWARD>(cart_2d))
{
}


template<typename real_t, size_t SIZE>
canon_prod<real_t, SIZE>::~canon()
{
}


template<typename real_t, size_t SIZE>
void canon_prod<real_t, SIZE>::operator()(
        row_matrix_type & product_result_partial,
        row_matrix_type & product_first_argument_partial,
        col_matrix_type & product_second_argument_partial)
    throw()
{
    local_product(
            product_result_partial,
            product_first_argument_partial,
            product_second_argument_partial);
}


}  // namespace algorithm
}  // namespace canon


#endif
