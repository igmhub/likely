// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FUNCTION_MINIMUM
#define LIKELY_FUNCTION_MINIMUM

#include "likely/types.h"

#include "boost/utility.hpp"

#include <string>
#include <iosfwd>

namespace likely {
    class Random;
	class FunctionMinimum : public boost::noncopyable {
	public:
	    // Represents the information known about an approximate function minimum.
		FunctionMinimum(double minValue, Parameters const &where);
		// The next form of the constructor adds an estimate of the convariance matrix
		// at the function minimum, which must be provided as a column-wise packed vector
		// of length npar*(npar+1)/2:
		//
		// m00 m01 m02 ... 
		//     m11 m12 ...  ==> { m00, m01, m11, m02, m12, m22, ... }
		//         m22 ...
		//             ...
		//
		// The corresponding index calculation is m(i,j) = array[i+j*(j+1)/2] for i<=j.
		// If j>i, then use m(i,j) = m(j,i). Set errorsOnly = true if the input covariance
		// vector should be interpreted as a list of diagonal errors of length npar.
		// Throws a RuntimeError if the input covariance is not (numerically) positive
		// definite. Use the updateCovariance method for more flexibility in handling
		// this error condition.
        FunctionMinimum(double minValue, Parameters const &where,
            PackedCovariance const &covar, bool errorsOnly = false);
		virtual ~FunctionMinimum();
		// Returns the function value at the minimum.
        double getMinValue() const;
		// Returns a copy of the parameter values at this minimum.
        Parameters getParameters() const;
        // Updates the location of the minimum and the function value at that point.
        void updateParameters(Parameters const &updatedParams, double updatedMinValue);
        // Returns true if a covariance matrix is available.
        bool haveCovariance() const;
        // Returns a vector of parameter error estimates or throws a RuntimeError if
        // no covariance matrix is available.
        Parameters getErrors() const;
        // Returns a smart pointer to the packed covariance matrix at this minimum.
        PackedCovariancePtr getCovariance() const;
        // Updates the covariance matrix associated with this minimum, if possible.
        // Refer to the constructor for details on the parameters. This method can be
        // used to add a covariance matrix to a minimum that did not originally have one.
        // Returns true if the covariance provided is (numerically) positive definite,
        // otherwise the covariance associated with this minimum is not changed.
        bool updateCovariance(PackedCovariance const &covar, bool errorsOnly = false);
        // Sets parameter values that are randomly sampled from this minimum and
        // returns the -log(weight) associated with the chosen parameters.
        double setRandomParameters(Parameters &params) const;
        // Ouptuts a multiline description of this minimum to the specified stream using
        // the specified printf format for floating point values.
        void printToStream(std::ostream &os, std::string formatSpec = "%.6f") const;
	private:
        double _minValue;
        Parameters _where;
        PackedCovariancePtr _covar;
        mutable PackedCovariancePtr _cholesky;
        mutable Random &_random;
	}; // FunctionMinimum
	
    inline double FunctionMinimum::getMinValue() const { return _minValue; }
    inline Parameters FunctionMinimum::getParameters() const { return Parameters(_where); }
    inline bool FunctionMinimum::haveCovariance() const { return bool(_covar); }
    inline PackedCovariancePtr FunctionMinimum::getCovariance() const { return _covar; }    
	
	// Returns a smart pointer to the Cholesky decomposition of a covariance matrix.
	// The return value will be null if the input covariance is not (numerically)
	// positive definite.
    PackedCovariancePtr choleskyDecomposition(PackedCovariance const &covar);
	
} // likely

#endif // LIKELY_FUNCTION_MINIMUM
