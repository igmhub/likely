// Created 22-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MINUIT_ENGINE
#define LIKELY_MINUIT_ENGINE

#include "likely/types.h"

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"

#include "boost/smart_ptr.hpp"

#include <vector>

namespace ROOT {
namespace Minuit2 {
    class MnUserParameterState;
    class SimplexMinimizer;
    class VariableMetricMinimizer;
}} // ROOT::Minuit2

namespace likely {
    // Implements minimization and error analysis using the Minuit2 library.
	class MinuitEngine : public ROOT::Minuit2::FCNBase {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters. If names are not specified, they will be P0,P1,P2,...
		MinuitEngine(Function f, int nPar);
		MinuitEngine(Function f, std::vector<std::string> const &parNames);
		virtual ~MinuitEngine();
		// Evaluates the engine's function for the specified input parameter values.
        virtual double operator()(Parameters const& pValues) const;
        // Returns the change in function value corresponding to one unit of error.
        // Can be changed to calculate different confidence intervals. For 1-sigma errors,
        // this value should be 1 for both chi-square and -2log(L) functions.
        virtual double Up() const;
        // Returns a smart pointer to the initial parameter state.
        typedef boost::shared_ptr<ROOT::Minuit2::MnUserParameterState> StatePtr;
        StatePtr getInitialState();
        // Runs a simplex minimization using the specified initial parameter values
        // and error estimates.
        ROOT::Minuit2::FunctionMinimum
            simplex(Parameters const &initial, Parameters const &errors);
        // Runs a variable-metric minimization using the specified initial parameter values
        // and error estimates.
        ROOT::Minuit2::FunctionMinimum
            variableMetric(Parameters const &initial, Parameters const &errors);
	private:
        int _nPar;
        Function _f;
        StatePtr _initialState;
        boost::scoped_ptr<ROOT::Minuit2::SimplexMinimizer> _simplex;
        boost::scoped_ptr<ROOT::Minuit2::VariableMetricMinimizer> _variableMetric;
        void _setInitialState(Parameters const &initial, Parameters const &errors);
	}; // MinuitEngine
	
	inline MinuitEngine::StatePtr MinuitEngine::getInitialState() {
        return _initialState;
	}
	
} // likely

#endif // LIKELY_MINUIT_ENGINE
