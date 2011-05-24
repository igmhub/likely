// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TYPES
#define LIKELY_TYPES

#include "boost/function.hpp"

#include <vector>

namespace likely {

    // Represents a vector of parameter values.
    typedef std::vector<double> Parameters;
    
    // Represents a likelihood function that, by convention, returns -logL(p).
    // The likelihood L(p) is not required to be normalized with respect to its
    // parameters. In case a function cannot be evaluated for its input parameters,
    // it should return +-inf, nan, or throw an exception.
    typedef boost::function<double (Parameters const &)> Function;

} // likely

#endif // LIKELY_TYPES
