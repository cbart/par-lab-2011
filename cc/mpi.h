// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__MPI__H__
#define __CANON__MPI__H__


#include <boost/mpi/communicator.hpp>


namespace canon
{
namespace mpi
{


template<size_t DIM_SIZE>
inline ::boost::mpi::communicator cart_square_sphere_create()
{
    MPI_Comm comm_cart;
    int dims = 2;  // two dimensional
    int dim_size[2] = {DIM_SIZE, DIM_SIZE};  // square
    int periods[2] = {1, 1};  // periods in both dimensions
    int reorder = true;
    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_size, periods, reorder, & comm_cart);
    return ::boost::mpi::communicator(comm_cart, ::boost::mpi::comm_take_ownership);
}


template<size_t PROCESSORS>
inline void assert_processors(const ::boost::mpi::communicator & comm, const ::boost::mpi::environment & env)
{
    if(comm.size() != PROCESSORS)
    {
        ::std::cerr << "Please run with " << PROCESSORS << " processors!\n";
        ::std::cerr << "Aborting..." << ::std::endl;
        env.abort(-1);
    }
}


}  // namespace mpi
}  // namespace canon



#endif
