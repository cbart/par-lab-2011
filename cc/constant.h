// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__CONSTANT__H__
#define __CANON__CONSTANT__H__


namespace canon
{


template<typename result_t, result_t value>
result_t constant()
{
    return value;
}


}  // namespace canon


#endif
