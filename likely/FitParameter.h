// Created 27-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FIT_PARAMETER
#define LIKELY_FIT_PARAMETER

#include "likely/types.h"

#include <string>
#include <vector>

namespace likely {
	class FitParameter {
	// Represents a lightweight snapshot of a fit parameter, specified by its name,
	// current value and error.
	public:
		FitParameter(std::string const &name, double value, double error = 0);
		// No virtual destructor since this is plain-old-data
        std::string getName() const;
        double getValue() const;
        double getError() const;
        bool isFloating() const;
	private:
        std::string _name;
        double _value, _error;
	}; // FitParameter

    inline FitParameter::FitParameter(std::string const &name, double value, double error)
    : _name(name), _value(value), _error(error) { }

    inline std::string FitParameter::getName() const { return _name; }
    inline double FitParameter::getValue() const { return _value; }
    inline double FitParameter::getError() const { return _error; }
    inline bool FitParameter::isFloating() const { return (0 != _error); }
    
    // Defines a vector of fit parameters.
    typedef std::vector<FitParameter> FitParameters;
    
    // Returns a vector of parameter values, errors or names. Use the optional onlyFloating
    // parameter to only include floating parameters in the result.
    void getFitParameterValues(FitParameters const &parameters, Parameters &values,
        bool onlyFloating = false);
    void getFitParameterErrors(FitParameters const &parameters, Parameters &values,
        bool onlyFloating = false);
    void getFitParameterNames(FitParameters const &parameters, std::vector<std::string> &names,
        bool onlyFloating = false);
    int countFloatingFitParameters(FitParameters const &parameters);

} // likely

#endif // LIKELY_FIT_PARAMETER
