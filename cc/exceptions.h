// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANNON__BOOST_EXCEPTIONS__H__
#define __CANNON__BOOST_EXCEPTIONS__H__


#include <stdexcept>


namespace cannon
{
namespace exception
{


const int EXCEPTION_ERROR = -1;


}  // namespace exception
}  // namespace cannon


#ifdef BOOST_NO_EXCEPTIONS


namespace boost
{


void throw_exception(const ::std::exception & exception)
{
    ::boost::mpi::communicator comm(MPI_COMM_WORLD, ::boost::mpi::comm_attach);
    comm.abort(::cannon::exception::EXCEPTION_ERROR);
}


}  // namespace boost


#endif


#endif
