// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TYPES
#define LIKELY_TYPES

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

#include <vector>

namespace likely {

    // Represents a vector of parameter values.
    typedef std::vector<double> Parameters;
    
    // Represents a gradient vector of function partial derivatives.
    typedef std::vector<double> Gradient;

    // Represents a column-wise packed covariance matrix near a local minimum.
    typedef std::vector<double> PackedCovariance;
    
    // Declares a smart pointer to a packed covariance matrix.
    typedef boost::shared_ptr<PackedCovariance> PackedCovariancePtr;
    
    // Declares a smart pointer to a const binning object.
    class AbsBinning;
    typedef boost::shared_ptr<const AbsBinning> AbsBinningCPtr;

    // Declares a smart pointer to a const covariance matrix.
    class CovarianceMatrix;
    typedef boost::shared_ptr<const CovarianceMatrix> CovarianceMatrixCPtr;

    // Represents a smart pointer to a minimization engine.
    class AbsEngine;
    typedef boost::shared_ptr<AbsEngine> AbsEnginePtr;

    // Encapsulates a minimization objective function.
    typedef boost::function<double (Parameters const &pValues)> Function;

    // Declares a smart pointer to an objective function.
    typedef boost::shared_ptr<Function> FunctionPtr;
    
    // Encapsulates a gradient calculator for a minimization objective function.
    typedef boost::function<void (Parameters const &pValues, Gradient &gValues)>
        GradientCalculator;
        
    // Declares a smart pointer to a gradient calculator.
    typedef boost::shared_ptr<GradientCalculator> GradientCalculatorPtr;

    // Represents a smart pointer to a function minimum object.
    class FunctionMinimum;
    typedef boost::shared_ptr<FunctionMinimum> FunctionMinimumPtr;
    
    // Represents a smart pointer to an interpolator object.
    class Interpolator;
    typedef boost::shared_ptr<Interpolator> InterpolatorPtr;
    
    // Represents a smart pointer to a weighted accumulator object.
    class AbsAccumulator;
    typedef boost::shared_ptr<AbsAccumulator> AbsAccumulatorPtr;
    
} // likely

#endif // LIKELY_TYPES
