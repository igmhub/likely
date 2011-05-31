// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TYPES
#define LIKELY_TYPES

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/numeric/ublas/symmetric.hpp"

#include <vector>

namespace likely {

    // Represents a vector of parameter values.
    typedef std::vector<double> Parameters;
    
    // Represents a likelihood function that, by convention, returns -logL(p).
    // The likelihood L(p) is not required to be normalized with respect to its
    // parameters. In case a function cannot be evaluated for its input parameters,
    // it should return +-inf, nan, or throw an exception.
    typedef boost::function<double (Parameters const &pValues)> Function;
    
    // Represents a vector of function gradients with respect to each parameter.
    typedef std::vector<double> Gradients;

    typedef boost::function<void (Parameters const &pValues, Gradients &gValues)>
        FunctionGradientCalculator;

    // Represents a Function covariance matrix near a local minimum.
    typedef boost::numeric::ublas::symmetric_matrix<double> Covariance;

    typedef boost::function<double (Parameters const &pInitial, Parameters const &pErrors,
    Parameters &pFinal, Covariance &covariance)>
        MinimumAndCovarianceFinder;

} // likely

#endif // LIKELY_TYPES
