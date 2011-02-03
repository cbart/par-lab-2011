// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__BOOST_EXCEPTIONS__H__
#define __CANON__BOOST_EXCEPTIONS__H__


#include <stdexcept>


namespace canon
{
namespace exception
{


const int EXCEPTION_ERROR = -1;


}  // namespace exception
}  // namespace canon


namespace boost
{


#ifdef BOOST_NO_EXCEPTIONS

void throw_exception(const ::std::exception & exception)
{
    ::boost::mpi::communicator comm(MPI_COMM_WORLD, ::boost::mpi::comm_attach);
    comm.abort(::canon::exception::EXCEPTION_ERROR);
}

#endif


}  // namespace boost


#endif
