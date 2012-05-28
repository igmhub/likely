// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FUNCTION_MINIMUM
#define LIKELY_FUNCTION_MINIMUM

#include "likely/types.h"
#include "likely/FitParameter.h"

#include "boost/utility.hpp"

#include <string>
#include <iosfwd>

namespace likely {
    // Represents the information known about an approximate function minimum.
	class FunctionMinimum : public boost::noncopyable {
	public:
	    // Creates a function minimum object for the specified function minimum value and
	    // estimated location of the minimum in parameter space.
		FunctionMinimum(double minValue, FitParameters const &parameters);
		// Creates a function minimum object that also specifies an estimated covariance matrix
		// for the subset of floating parameters. Agreement between parameter errors and covariance
		// matrix diagonal elements is not checked or required.
        FunctionMinimum(double minValue, FitParameters const &parameters, CovarianceMatrixCPtr covariance);
		virtual ~FunctionMinimum();
		// Returns the function value at the minimum.
        double getMinValue() const;
		// Returns a vector of parameter values at this minimum. If onlyFloating is true, only
		// the values of floating parameters are included in the returned vector.
        Parameters getParameters(bool onlyFloating = false) const;
		// Filters the input parameters to only include floating parameters.
        void filterParameterValues(Parameters const &allValues, Parameters &floatingValues) const;
		// Returns a vector of parameter errors at this minimum. If onlyFloating is true, only
		// the errors of floating parameters are included in the returned vector.
        Parameters getErrors(bool onlyFloating = false) const;
        // Updates the fit parameters and function value at the minimum.
        void updateParameters(double minValue, FitParameters const &parameters);
        // Updates the location of the minimum and the function value at that point. If a covariance
        // matrix is available, its diagonal elements will be used to update fit parameter errors.
        void updateParameterValues(double minValue, Parameters const &values);
        // Returns true if a covariance matrix is available.
        bool hasCovariance() const;
        // Returns a pointer to the estimated covariance matrix of floating parameters at this minimum.
        CovarianceMatrixCPtr getCovariance() const;
        // Updates the covariance matrix associated with this minimum, if possible.
        // This method can be used to add a covariance matrix to a minimum that did not
        // originally have one. This method does not update the errors associated with our
        // parameters, but calling it before updateParameterValues() will have this effect.
        void updateCovariance(CovarianceMatrixCPtr covariance);
        // Sets parameter values that are randomly sampled from this minimum and
        // returns the -log(weight) associated with the chosen parameters. Parameters that
        // are not floating will be included in the result, but not varied.
        double setRandomParameters(Parameters &params) const;
        // Ouptuts a multiline description of this minimum to the specified stream using
        // the specified printf format for floating point values.
        void printToStream(std::ostream &os, std::string formatSpec = "%12.6f") const;
	private:
        double _minValue;
        int _nFloating;
        FitParameters _parameters;
        CovarianceMatrixCPtr _covar;
	}; // FunctionMinimum
	
    inline double FunctionMinimum::getMinValue() const { return _minValue; }
    inline bool FunctionMinimum::hasCovariance() const { return bool(_covar); }
    inline CovarianceMatrixCPtr FunctionMinimum::getCovariance() const { return _covar; }    
	
} // likely

#endif // LIKELY_FUNCTION_MINIMUM
