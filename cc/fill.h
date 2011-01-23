// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__FILL__H__
#define __CANON__FILL__H__


#include <algorithm>
#include "matrix.h"


namespace canon
{


template<typename element_t, typename layout_t, size_t size, typename generator_t>
void fill(typename square_matrix<element_t, layout_t, size>::type & matrix, generator_t generator)
    throw()
{
    ::std::generate(matrix.data().begin(), matrix.data().end(), generator);
}


}  // namespace canon


#endif
