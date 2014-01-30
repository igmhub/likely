// Created 08-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FIT_MODEL
#define LIKELY_FIT_MODEL

#include "likely/FitParameter.h"

#include <string>
#include <iosfwd>
#include <map>

namespace likely {
	class FitModel {
	// Represents a model that is parameterized by FitParameters, with support for "lazy"
	// implementations that only do expensive operations when specific parameter values change,
	// and convenient parameter configuration via a simple scripting language.
	public:
	    // Creates a new FitModel with the specified name. The definition of model parameters is
	    // delegated to subclass constructors.
		FitModel(std::string const &name);
		virtual ~FitModel();
		// Returns a copy of our model name.
        std::string getName() const;
        // Returns the number of model parameters.
        int getNParameters(bool onlyFloating = false) const;
        // Prints a multi-line description of this object to the specified output stream,
        // including information about our parameters and their initial values.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
        // Prints a multi-line summary of our current parameter values to the specified
        // output stream. Values that are flagged as changed will be marked with an asterisk.
        void printCurrentValues(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
        // Configures our fit parameters using the specified script. See the documentation for
        // modifyFitParameters in FitParameters.h for details. Note that this method modifies our
        // parameter configuration, i.e., their default values and errors, but not their current values
        // as reported by setParameterValue() or returned by getParameterValue().
        virtual void configureFitParameters(std::string const &script);
        // Finds the minimum of the specified function with respect to our parameters,
        // using the specified method. Use the configureFitParameters method to change the initial
        // parameter values and errors for this and subsequent minimizations, or provide an optional
        // config script to make one-time changes that only apply to this minimization. Note that,
        // depending on how fptr is implemented, our current parameter values might change as a result
        // of this operation, but our parameter configuration is guaranteed not to change. If priors
        // should be imposed, the fptr provided should arrange to call our evaluatePrior() method.
        // Note that any modifications to priors in oneTimeConfig will have no effect on evaluatePrior().
        FunctionMinimumPtr findMinimum(FunctionPtr fptr, std::string const &method,
            std::string const &oneTimeConfig = "");
        // Returns the current value of the named parameter, or throws a RuntimeError for an
        // invalid parameter name. If no value has ever been set, the default value specified
        // when this parameter was first defined is returned. Calls to configureFitParmeters()
        // will not change the current value.
        double getParameterValue(std::string const &name) const;
        // Sets the current value of the named parameter, or throws a RuntimeError for an
        // invalid parameter name.
        void setParameterValue(std::string const &name, double value);
        // Returns the -log(priors-likelihood) for the current parameter values and fit parameter
        // configuration. Priors are not imposed on fixed parameters, even if they are defined.
        // Use configureFitParameters() to add priors.
        double evaluatePriors() const;
        // Returns a guess at the specified function's minimum based on our fit parameter
        // configuration and assuming a diagonal covariance matrix.
        FunctionMinimumPtr guessMinimum(FunctionPtr fptr) const;
    protected:
        // Subclasses use this method to define their parameters. Parameters should generally
        // be specified with a reasonable error > 0 since the configureFitParameters() method
        // provides a convenient way to fix a parameter before a fit. Throws a RuntimeError for an
        // invalid parameter name or error < 0. Immediately after a parameter has been defined,
        // getParameterValue() returns the value specified here, and isParameterValueChanged()
        // returns true. This method returns the internal index associated with the newly
        // defined parameter, which can be used with the get/setParameterValue methods below
        // (but the named forms should be generally preferred for readability, except in
        // performance critical code).
        int defineParameter(std::string const &name, double value, double error = 0);
        // Gets/sets the current value of the named parameter, or throws a RuntimeError for an
        // invalid parameter name. Subclasses should generally prefer the public named versions
        // of these methods for more readable code, unless there is a specific
        // reason not to (e.g., a demonstrated performance issue).
        double getParameterValue(int index) const;
        void setParameterValue(int index, double value);
        // Updates our parameter values to those in the specified vector. Throws a RuntimeError
        // if the input vector has the wrong size. Registers and records any changes to our
        // current parameter values. Returns true if any parameters are changed after this operation,
        // including changes that were present before this operation and since the last call to
        // resetParameterValuesChanged().
        bool updateParameterValues(Parameters const &values);
        // Returns true if the specified parameter has changed since the last call to
        // resetParameterValuesChanged(). Parameter value changes are registered by calls to
        // updateParameterValues() and setParameterValue(). Note that values of our FitParameters
        // represent defaults (e.g., for setting initial values in a fit) that are independent
        // of the current values tracked here. In particular, configureFitParameters() modifies
        // the defaults but not the current values.
        bool isParameterValueChanged(std::string const &name) const;
        bool isParameterValueChanged(int index) const;
        // Resets our tracking of parameter value changes so that all parameters are
        // considered unchanged, until future changes.
        void resetParameterValuesChanged();
	private:
        std::string _name;
        FitParameters _parameters;
        Parameters _parameterValue;
        std::vector<bool> _parameterValueChanged;
        typedef std::map<std::string,int> NameIndexMap;
        NameIndexMap _nameIndexMap;
        // Sets the value of the parameter with the specified index, which is assumed to be valid.
        void _setParameterValue(int validIndex, double value);
        // Throws a RuntimeError if the specified index is out of range, or else returns the index.
        int _checkIndex(int index) const;
        // Returns the index of the named parameter in our vectors, or throws a RuntimeError
        // for an invalid name.
        int _getIndex(std::string const &name) const;
	}; // FitModel
	
    inline std::string FitModel::getName() const { return _name; }

    inline double FitModel::getParameterValue(std::string const &name) const {
        return _parameterValue[_getIndex(name)];
    }
    inline double FitModel::getParameterValue(int index) const {
        return _parameterValue[_checkIndex(index)];
    }
    
    inline void FitModel::setParameterValue(std::string const &name, double value) {
        _setParameterValue(_getIndex(name), value);
    }
    inline void FitModel::setParameterValue(int index, double value) {
        _setParameterValue(_checkIndex(index), value);
    }
    inline void FitModel::_setParameterValue(int validIndex, double value) {
        if(value != _parameterValue[validIndex]) {
            _parameterValue[validIndex] = value;
            _parameterValueChanged[validIndex] = true;
        }
    }
    
    inline bool FitModel::isParameterValueChanged(std::string const &name) const {
        return _parameterValueChanged[_getIndex(name)];
    }
    inline bool FitModel::isParameterValueChanged(int index) const {
        return _parameterValueChanged[_checkIndex(index)];
    }
    inline void FitModel::resetParameterValuesChanged() {
        _parameterValueChanged.assign(_parameterValueChanged.size(),false);
    }
	
} // likely

#endif // LIKELY_FIT_MODEL
