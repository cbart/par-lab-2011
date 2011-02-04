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

    // ::std::vector based allocation - considered faster and more parallel-able

    typedef ::std::vector<element_type> storage_type;
    typedef ::boost::numeric::ublas::matrix<element_type, layout_type, storage_type> type;

    // ::boost::numeric::bounded_array based allocation - enables full-stack allocation

    // typedef ::boost::numeric::ublas::bounded_matrix<element_type, size, size, layout_type> type;

    // Returns raw pointer to container begin (needed because ::std::vector uses fancy-iterators).
    static element_type * begin(type * matrix);

private:
    square_matrix();
    square_matrix(const square_matrix & copy);
};


template<typename element_t, typename layout_t, size_t size>
inline typename square_matrix<element_t, layout_t, size>::element_type * square_matrix<element_t, layout_t, size>::begin(type * matrix)
{
    return & matrix->data()[0];
}


// Bounded array implementation

// template<typename element_t, typename layout_t, size_t size>
// inline typename square_matrix<element_t, layout_t, size>::element_type * square_matrix<element_t, layout_t, size>::begin(type * matrix)
// {
//    return & matrix->data().begin();
// }



}  // namespace cannon


#endif
