// Created 27-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FIT_PARAMETER
#define LIKELY_FIT_PARAMETER

#include "likely/types.h"

#include <string>
#include <vector>
#include <iosfwd>

namespace likely {
	class FitParameter {
	// Represents a fit parameter, specified by its name, value and estimated error.
	// Class has no virtual destructor since this is plain-old-data. Do not inherit
	// from this class.
	public:
	    // Creates a new FitParameter or throws a RuntimeError if the optional error
	    // is negative or the name uses invalid characters. The valid characters are
	    // " a-zA-Z0-9()*/+-", and a name cannot end with an asterisk. A zero error
	    // implies that the parameter is fixed.
		FitParameter(std::string const &name, double value, double error = 0);
        // Returns a copy of this parameter's name.
        std::string getName() const;
        // Returns this parameter's value.
        double getValue() const;
        // Sets a new value for this parameter.
        void setValue(double value);
        // Returns this parameter's estimated error, or zero if the parameter is fixed.
        double getError() const;
        // Sets a new error for this parameter, or throws a RuntimeError if the specified
        // error is < 0. Set an error of zero to permanently fix a parameter. Use the fix()
        // method to temporarily fix a parameter, so that release() will restore its original
        // estimated error.
        void setError(double error);
        void fix();
        void release();
        // Returns true if this parameter is floating, in which case it will have an error > 0.
        bool isFloating() const;
        // Returns the type of prior imposed on this parameter.
        enum PriorType { None, Box, Gaussian };
        PriorType getPriorType() const;
        // Returns limits of an imposed prior or zero if no prior is imposed.
        double getPriorMin() const;
        double getPriorMax() const;
        // Sets (or resets) the prior on this parameter.
        void setPrior(double priorMin, double priorMax, PriorType type);
        // Removes any prior on this parameter.
        void removePrior();
        // Returns the set of valid characters in a FitParameter name.
        static std::string const &getValidNameCharacters();
	private:
        std::string _name;
        double _value, _error, _priorMin, _priorMax;
        PriorType _priorType;
        
	}; // FitParameter

    inline std::string FitParameter::getName() const { return _name; }
    inline double FitParameter::getValue() const { return _value; }
    inline void FitParameter::setValue(double value) { _value = value; }
    inline double FitParameter::getError() const { return _error < 0 ? 0 : _error; }
    inline void FitParameter::fix() { if(_error > 0) _error = -_error; }
    inline void FitParameter::release() { if(_error < 0) _error = -_error; }
    inline bool FitParameter::isFloating() const { return (_error > 0); }
    inline FitParameter::PriorType FitParameter::getPriorType() const { return _priorType; }
    inline double FitParameter::getPriorMin() const { return _priorType == None ? 0 : _priorMin; }
    inline double FitParameter::getPriorMax() const { return _priorType == None ? 0 : _priorMax; }
    inline void FitParameter::removePrior() { _priorType = None; }
    
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
        
    // Sets parameter values. Use the optional onlyFloating parameter to only update
    // floating parameter values.
    void setFitParameterValues(FitParameters &parameters, Parameters const &values,
        bool onlyFloating = false);
    void setFitParameterValues(FitParameters &parameters, Parameters::const_iterator first,
        Parameters::const_iterator last, bool onlyFloating = false);
        
    // Returns the number of total or floating parameters.
    int countFitParameters(FitParameters const &parameters, bool onlyFloating = false);
    
    // Prints a multi-line description of FitParameters to the specified output stream.
    void printFitParametersToStream(FitParameters const &parameters, std::ostream &out,
        std::string const &formatSpec = "%12.6f");

    // Returns the index of the parameter with the specified name or throws a RuntimeError.
    int findFitParameterByName(FitParameters const &parameters, std::string const &name);
        
    // Modifies FitParameters using instructions in the specified script or throws a
    // RuntimeError in case an error in the script is detected (in which case the input
    // parameters will not be modified). The following script commands are supported:
    //
    //  value [<name>] = <newvalue>
    //  error [<name>] = <newerror>
    //  fix [<name>] = <newvalue>
    //  fix [<name>]
    //  release [<hame>]
    //
    // Multiple commands separated by semicolons are executed in the order they appear.
    // Command verbs (value,error,fix,release) are case sensitive. Arbitrary whitespace
    // is allowed between tokens, except within names delimited by [ ], where whitespace
    // is considered part of the name. If <name> ends with an asterisk, it is interpreted
    // as a wildcard that matches zero or more trailing characters, and must match at
    // least one parameter name.
    void modifyFitParameters(FitParameters &parameters, std::string const &script);

    // A formatted string of a value and its error(s) is returned. The value and its 
    // error(s) are rounded to a precision that matches that of the smallest error. 
    // The precision of the smallest error is determined as follows: if the three highest 
    // order digits of the error lie between 100 and 354, round to two signiﬁcant digits. 
    // If they lie between 355 and 949, round to one signiﬁcant digit. Finally, if they 
    // lie between 950 and 999, round up to 1000 and keep two signiﬁcant digits. If the 
    // errors vector is empty, the value is returned as a string.
    std::string roundValueWithError(double value, std::vector<double> const &errors, 
        std::string const &seperator = "+/-");

} // likely

#endif // LIKELY_FIT_PARAMETER
