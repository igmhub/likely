// Created 27-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FIT_PARAMETER
#define LIKELY_FIT_PARAMETER

#include <string>

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
	private:
        std::string _name;
        double _value, _error;
	}; // FitParameter

    inline FitParameter::FitParameter(std::string const &name, double value, double error)
    : _name(name), _value(value), _error(error) { }

    inline std::string FitParameter::getName() const { return _name; }
    inline double FitParameter::getValue() const { return _value; }
    inline double FitParameter::getError() const { return _error; }
    
    typedef std::vector<FitParameter> FitParameters;

} // likely

#endif // LIKELY_FIT_PARAMETER
