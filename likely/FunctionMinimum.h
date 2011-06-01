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
		virtual ~FunctionMinimum();
		// Returns the function value at the minimum.
        double getMinValue() const;
		// Returns a copy of the parameter values at this minimum.
        Parameters getParameters() const;
        // Ouptuts a multiline description of this minimum to the specified stream using
        // the specified printf format for floating point values.
        void printToStream(std::ostream &os, std::string formatSpec = "%.6f") const;
	private:
        double _minValue;
        Parameters _where;
	}; // FunctionMinimum
	
    inline double FunctionMinimum::getMinValue() const { return _minValue; }
    inline Parameters FunctionMinimum::getParameters() const { return Parameters(_where); }
	
} // likely

#endif // LIKELY_FUNCTION_MINIMUM
