// Created 08-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FIT_MODEL
#define LIKELY_FIT_MODEL

#include "likely/FitParameter.h"

#include <string>
#include <iosfwd>

namespace likely {
	class FitModel {
	public:
		FitModel(std::string const &name);
		virtual ~FitModel();
		// Returns a copy of our model name.
        std::string getName() const;
        // Returns the number of model parameters.
        int getNParameters(bool onlyFloating = false) const;
        // Prints a multi-line description of this object to the specified output stream.
        virtual void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
        // Configures our fit parameters using the specified script. See the documentation for
        // modifyFitParameters in FitParameters.h for details.
        virtual void configure(std::string const &script);
        // Finds the minimum of the specified function with respect to our parameters,
        // using the specified method. Use the configure method to change the initial
        // parameter values and errors.
        FunctionMinimumPtr findMinimum(FunctionPtr fptr, std::string const &method) const;
    protected:
        // Subclasses use this method to define their parameters. Parameters should generally
        // be specified with a reasonable error > 0 since the configure() method provides a
        // convenient way to fix a parameter before a fit. Throws a RuntimeError for an
        // invalid parameter name or error < 0.
        void defineParameter(std::string const &name, double value, double error);
	private:
        std::string _name;
        FitParameters _parameters;
	}; // FitModel
	
    inline std::string FitModel::getName() const { return _name; }
	
} // likely

#endif // LIKELY_FIT_MODEL
