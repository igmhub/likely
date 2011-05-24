// Created 22-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MINUIT_ENGINE
#define LIKELY_MINUIT_ENGINE

#include "likely/types.h"

#include "Minuit2/FCNBase.h"

#include "boost/smart_ptr.hpp"

#include <vector>

namespace ROOT {
namespace Minuit2 {
    class SimplexMinimizer;
}} // ROOT::Minuit2

namespace likely {
    // Implements minimization and error analysis using the Minuit2 library.
	class MinuitEngine : public ROOT::Minuit2::FCNBase {
	public:
	    // Creates a new engine for the specified function.
		MinuitEngine(Function f);
		virtual ~MinuitEngine();
		// Returns the value of -2log(Likelihood) for the specified input parameter values.
        virtual double operator()(std::vector<double> const& pValues) const;
        // Returns the change in function value corresponding to one unit of error.
        // Can be changed to calculate different confidence intervals. For 1-sigma errors,
        // this value should be 1 for both chi-square and -2log(L) functions.
        virtual double Up() const;
        // Runs a simplex minimization using the specified initial parameter values
        // and error estimates, and returns the parameters at the minimum.
        Parameters simplex(Parameters const &initial, Parameters const &errors);
	private:
        Function _f;
        boost::scoped_ptr<ROOT::Minuit2::SimplexMinimizer> _simplex;
	}; // MinuitEngine
} // likely

#endif // LIKELY_MINUIT_ENGINE
