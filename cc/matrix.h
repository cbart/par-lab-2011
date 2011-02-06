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


// Square matrix concept
template<typename element_t, typename storage_t, typename layout_t, size_t SIZE>
class square_matrix_concept
{
public:
    typedef element_t element_type;
    typedef layout_t layout_type;
    typedef storage_t storage_type;
    typedef ::boost::numeric::ublas::matrix<element_type, layout_type, storage_type> type;
public:
    static element_type * begin(type * matrix);
};


template<typename element_t, typename storage_t, typename layout_t, size_t SIZE>
typename square_matrix_concept<element_t, storage_t, layout_t, SIZE>::element_type *
square_matrix_concept<element_t, storage_t, layout_t, SIZE>::begin(type * matrix)
{
    return matrix->data().begin();
}


// Swuare matrix concept - vector specialization
template<typename element_t, template<typename a_t> class allocator_template, typename layout_t, size_t SIZE>
class square_matrix_concept<element_t, ::std::vector<element_t, allocator_template<element_t> >, layout_t, SIZE>
{
public:
    typedef element_t element_type;
    typedef layout_t layout_type;
    typedef allocator_template<element_t> allocator_type;
    typedef ::std::vector<element_type, allocator_type> storage_type;
    typedef ::boost::numeric::ublas::matrix<element_type, layout_type, storage_type> type;
public:
    static element_type * begin(type * matrix);
};


template<typename element_t, template<typename a_t> class allocator_template, typename layout_t, size_t SIZE>
typename square_matrix_concept<element_t, ::std::vector<element_t, allocator_template<element_t> >, layout_t, SIZE>::element_type *
square_matrix_concept<element_t, ::std::vector<element_t, allocator_template<element_t> >, layout_t, SIZE>::begin(type * matrix)
{
    // Specialization because of vector's fancy iterators
    return & matrix->data()[0];
}


}  // namespace cannon


#endif
