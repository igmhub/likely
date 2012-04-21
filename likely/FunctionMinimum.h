// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FUNCTION_MINIMUM
#define LIKELY_FUNCTION_MINIMUM

#include "likely/types.h"

#include "boost/utility.hpp"

#include <string>
#include <iosfwd>

namespace likely {
    // Represents the information known about an approximate function minimum.
	class FunctionMinimum : public boost::noncopyable {
	public:
	    // Creates a function minimum object for the specified function minimum value and
	    // estimated location of the minimum in parameter space.
		FunctionMinimum(double minValue, Parameters const &where);
		// Creates a function minimum object that also specifies an estimated covariance matrix.
        FunctionMinimum(double minValue, Parameters const &where, CovarianceMatrixCPtr covariance);
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
        CovarianceMatrixCPtr getCovariance() const;
        // Updates the covariance matrix associated with this minimum, if possible.
        // This method can be used to add a covariance matrix to a minimum that did not
        // originally have one.
        void updateCovariance(CovarianceMatrixCPtr covariance);
        // Sets parameter values that are randomly sampled from this minimum and
        // returns the -log(weight) associated with the chosen parameters.
        double setRandomParameters(Parameters &params) const;
        // Ouptuts a multiline description of this minimum to the specified stream using
        // the specified printf format for floating point values.
        void printToStream(std::ostream &os, std::string formatSpec = "%+.6f") const;
	private:
        double _minValue;
        Parameters _where;
        CovarianceMatrixCPtr _covar;
	}; // FunctionMinimum
	
    inline double FunctionMinimum::getMinValue() const { return _minValue; }
    inline Parameters FunctionMinimum::getParameters() const { return Parameters(_where); }
    inline bool FunctionMinimum::haveCovariance() const { return bool(_covar); }
    inline CovarianceMatrixCPtr FunctionMinimum::getCovariance() const { return _covar; }    
	
} // likely

#endif // LIKELY_FUNCTION_MINIMUM
