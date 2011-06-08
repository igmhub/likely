// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FUNCTION_MINIMUM
#define LIKELY_FUNCTION_MINIMUM

#include "likely/types.h"

#include <string>
#include <iosfwd>

namespace likely {
	class FunctionMinimum {
	public:
	    // Represents the information known about an approximate function minimum.
		FunctionMinimum(double minValue, Parameters const &where);
		// The second form of the constructor adds an estimate of the convariance matrix
		// at the function minimum, which must be provided as a column-wise packed vector:
		//
		// m00 m01 m02 ... 
		//     m11 m12 ...  ==> { m00, m01, m11, m02, m12, m22, ... }
		//         m22 ...
		//             ...
		//
		// The corresponding index calculation is m(i,j) = array[i+j*(j+1)/2] for i<=j.
		// If j>i, then use m(i,j) = m(j,i).
        FunctionMinimum(double minValue, Parameters const &where,
            PackedCovariance const &covar);
		virtual ~FunctionMinimum();
		// Returns the function value at the minimum.
        double getMinValue() const;
		// Returns a copy of the parameter values at this minimum.
        Parameters getParameters() const;
        // Returns true if a covariance matrix is available.
        bool haveCovariance() const;
        // Returns a vector of parameter error estimates or throws a RuntimeError if
        // no covariance matrix is available.
        Parameters getErrors() const;
        // Returns parameter values that are randomly sampled from this minimum.
        Parameters getRandomParameters() const;
        // Ouptuts a multiline description of this minimum to the specified stream using
        // the specified printf format for floating point values.
        void printToStream(std::ostream &os, std::string formatSpec = "%.6f") const;
	private:
        double _minValue;
        Parameters _where;
        bool _haveCovariance;
        PackedCovariance _covar;
        mutable bool _haveCholesky;
        mutable PackedCovariance _cholesky;
	}; // FunctionMinimum
	
    inline double FunctionMinimum::getMinValue() const { return _minValue; }
    inline Parameters FunctionMinimum::getParameters() const { return Parameters(_where); }
    inline bool FunctionMinimum::haveCovariance() const { return _haveCovariance; }
	
} // likely

#endif // LIKELY_FUNCTION_MINIMUM
