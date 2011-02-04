// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__FILL__H__
#define __CANNON__FILL__H__


#include <algorithm>


namespace cannon
{


// Fills given matrix with elements returned by `generator()`.
template<typename matrix_t, typename generator_t>
void fill(matrix_t & matrix, generator_t generator)
    throw()
{
    ::std::generate(matrix.data().begin(), matrix.data().end(), generator);
}


}  // namespace cannon


#endif
