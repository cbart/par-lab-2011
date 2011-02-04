// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__MATRIX__H__
#define __CANNON__MATRIX__H__


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <vector>


namespace cannon
{


// Possible matrix layouts
typedef ::boost::numeric::ublas::row_major row_major;
typedef ::boost::numeric::ublas::column_major col_major;


// Just typedef for boost's bounded_matrix
template<typename element_t, typename layout_t, size_t size>
class square_matrix
{
public:
    typedef element_t element_type;
    typedef layout_t layout_type;
    typedef ::std::vector<element_type> storage_type;
    //typedef ::boost::numeric::ublas::matrix<element_type, layout_type, storage_type> type;
    typedef ::boost::numeric::ublas::bounded_matrix<element_type, size, size, layout_type> type;
};


}  // namespace cannon


#endif
