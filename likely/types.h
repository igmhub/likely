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
    
    // Represents a vector of function gradients with respect to each parameter.
    typedef std::vector<double> Gradients;

    // Represents a Function covariance matrix near a local minimum.
    typedef boost::numeric::ublas::symmetric_matrix<double> Covariance;

    // Encapsulates a minimization objective function.
    typedef boost::function<double (Parameters const &pValues)> Function;

    // Declares a smart pointer to an objective function.
    typedef boost::shared_ptr<Function> FunctionPtr;

    // Represents a smart pointer to a function minimum object.
    class FunctionMinimum;
    typedef boost::shared_ptr<FunctionMinimum> FunctionMinimumPtr;

} // likely

#endif // LIKELY_TYPES
