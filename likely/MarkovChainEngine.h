// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MARKOV_CHAIN_ENGINE
#define LIKELY_MARKOV_CHAIN_ENGINE

#include "likely/types.h"
#include "likely/AbsEngine.h"

#include "boost/function.hpp"

namespace likely {
	class MarkovChainEngine : public AbsEngine {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters.
		MarkovChainEngine(FunctionPtr f, GradientCalculatorPtr gc, FitParameters const &parameters,
            std::string const &algorithm, RandomPtr random = RandomPtr());
		virtual ~MarkovChainEngine();
        // Searches for a minimum by taking a sequence of random steps.
        void minimize(FunctionMinimumPtr fmin,
            double prec, int maxSteps, int acceptsPerParam, int maxTrialsPerParam);
	    // Generates samples using a FunctionMinimum's covariance to specify the proposal
	    // function until the specified number of trials have been accepted or the specified
	    // maximum number of trials has been generated. Returns the total number of trials
	    // generated (including duplicates after a trial is rejected). Use maxTrials=0 for
	    // unlimited trials. Updates the function minimum with an improved estimate.
	    // Provide an optional callback to see all trial steps. The callback parameters
	    // are the current and trial parameters, the function values at these parameters,
	    // and a boolean to flag if the trial is accepted. The callback is invoked every
	    // callbackInterval steps, with the first call after callbackInterval steps.
        typedef boost::function<void (Parameters const&, Parameters const&, double, double, bool)> Callback;
        int generate(FunctionMinimumPtr fmin, int nAccepts, int maxTrials,
            Callback callback = Callback(), int callbackInterval = 1) const;
	private:
        int _nParam,_nFloating;
        FunctionPtr _f;
        mutable RandomPtr _random;
	}; // MarkovChainEngine

    // Registers our named methods.
    void registerMarkovChainEngineMethods();

} // likely

#endif // LIKELY_MARKOV_CHAIN_ENGINE
