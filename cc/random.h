// Author: Cezary Bartoszuk
// Email: cbart@students.mimuw.edu.pl

#ifndef __CANON__RANDOM__H__
#define __CANON__RANDOM__H__


#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>


namespace canon
{


// Random numbers generator. To run just create an instance yourself
// and run the operator() to get next random number.
template<typename result_t>
class random_generator
{
public:
    typedef result_t result_type;
    typedef ::boost::rand48 engine_type;
    typedef engine_type & engine_ref;
    typedef ::boost::uniform_01<result_t> distribution_type;
    typedef ::boost::variate_generator<engine_ref, distribution_type> generator_type;
private:
    engine_type engine;
    distribution_type distribution;
    generator_type generator;
public:
    random_generator()
        throw();
    ~random_generator()
        throw();
    result_t operator()();
};


template<typename result_t>
random_generator<result_t>::random_generator()
    throw()
  : engine(),
    distribution(),
    generator(engine, distribution)
{
}


template<typename result_t>
random_generator<result_t>::~random_generator()
    throw()
{
}


template<typename result_t>
typename random_generator<result_t>::result_type random_generator<result_t>::operator()()
{
    return generator();
}


}  // namespace canon


#endif
