// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__MULTIPLY__H__
#define __CANNON__MULTIPLY__H__


#include <boost/numeric/ublas/matrix_expression.hpp>
#include "matrix.h"
#include "exceptions.h"


namespace cannon
{


// Runs ublas FORTRAN procedure to perform a matrix product
// it's actually equal to:
//   result += left * right
template<typename element_t, typename storage_t, size_t size>
void prod(
        typename square_matrix_concept<element_t, storage_t, row_major, size>::type & result,
        typename square_matrix_concept<element_t, storage_t, row_major, size>::type & left,
        typename square_matrix_concept<element_t, storage_t, col_major, size>::type & right)
    throw()
{
    // `noalias` uses proxy in order not to copy additional memory.
    ::boost::numeric::ublas::noalias(result) += ::boost::numeric::ublas::prod(left, right);
}


}  // namespace cannon


#endif
