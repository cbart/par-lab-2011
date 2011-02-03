// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__CONSTANT__H__
#define __CANON__CONSTANT__H__


namespace canon
{


// Constant generator.
template<typename result_t, long value>
result_t constant()
{
    return static_cast<result_t>(value);
}


}  // namespace canon


#endif
