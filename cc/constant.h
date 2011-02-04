// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__CONSTANT__H__
#define __CANNON__CONSTANT__H__


namespace cannon
{


// Constant generator.
template<typename result_t, long value>
result_t constant()
{
    return static_cast<result_t>(value);
}


}  // namespace cannon


#endif
