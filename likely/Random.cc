// Created 08-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Random.h"

#include "boost/random/uniform_01.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"

namespace local = likely;

local::Random::Random() :
_uniform(boost::variate_generator<boost::mt19937&, boost::uniform_01<> >
    (_generator, boost::uniform_01<>())),
_gauss(boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    (_generator, boost::normal_distribution<>(0,1)))
{
}

local::Random &local::Random::instance() {
    static Random *_instance = new Random();
    return *_instance;
}

float local::Random::getFastUniform() {
    return 0;
}