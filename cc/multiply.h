// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__MULTIPLY__H__
#define __CANON__MULTIPLY__H__


#include <boost/numeric/ublas/matrix_expression.hpp>
#include "matrix.h"
#include "exceptions.h"


namespace canon
{


// performs result += left * right;
template<typename element_t, size_t size>
void prod(
        typename square_matrix<element_t, row_major, size>::type & result,
        typename square_matrix<element_t, row_major, size>::type & left,
        typename square_matrix<element_t, col_major, size>::type & right)
    throw()
{
    ::boost::numeric::ublas::noalias(result) += ::boost::numeric::ublas::prod(left, right);
}


}  // namespace canon


#endif
