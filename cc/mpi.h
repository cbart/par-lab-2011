// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__MPI__H__
#define __CANNON__MPI__H__


#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/shared_array.hpp>
#include "debug.h"


namespace cannon
{
namespace mpi
{


typedef int rank_type;

typedef ::boost::shared_array<rank_type> ranks_array_type;

const size_t RANKS_ARRAY_SIZE = 2;

// ranks_array[`SOURCE_RANK_INDEX`] is the source rank
const size_t SOURCE_RANK_INDEX = 0;

// ranks_array[`DESTINATION_RANK_INDEX`] is the destination rank
const size_t DESTINATION_RANK_INDEX = 1;

// Possible direction arguments for `shift`
const int DIRECTION_VERTICAL = 0;
const int DIRECTION_HORIZONTAL = 1;

// Possible displacement arguments for `shift`
const int DISPLACEMENT_UPWARD = 1;
const int DISPLACEMENT_DOWNWARD = -1;


// Creates cartesian square sphere
// The "square sphere" means that the topology is a cartesian square
// of `DIM_SIZE` * `DIM_SIZE` nodes and there is a period in both dimensions.
template<size_t DIM_SIZE>
inline ::boost::mpi::communicator cart_square_sphere_create()
{
    MPI_Comm comm_cart;
    int dims = 2;  // two dimensional
    int dim_size[2] = {DIM_SIZE, DIM_SIZE};  // square
    int periods[2] = {true, true};  // periods in both dimensions
    int reorder = true;
    MPI_Cart_create(MPI_COMM_WORLD, dims, dim_size, periods, reorder, & comm_cart);
    return ::boost::mpi::communicator(comm_cart, ::boost::mpi::comm_take_ownership);
}


// Fails if mpi wasn't run with `PROCESSORS` amount of processors
template<size_t PROCESSORS>
inline void assert_processors(
        const ::boost::mpi::communicator & comm,
        const ::boost::mpi::environment & env)
{
    if(comm.size() != PROCESSORS)
    {
        ::debug::err << "Please run with " << PROCESSORS << " processors!\n";
        ::debug::err << "Aborting..." << ::std::endl;
        env.abort(-1);
    }
}


// Returns a ranks array which `SOURCE_RANK_INDEX`st element
// is the source rank and the `DESTINATION_RANK_INDEX`st element
// is the destination rank.
template<int DIRECTION, int DISPL>
inline ranks_array_type shift(const ::boost::mpi::communicator & comm)
{
    ranks_array_type ranks(new rank_type[RANKS_ARRAY_SIZE]);
    MPI_Cart_shift(comm, DIRECTION, DISPL, & ranks[SOURCE_RANK_INDEX], & ranks[DESTINATION_RANK_INDEX]);
    return ranks;
}


}  // namespace mpi
}  // namespace cannon



#endif
