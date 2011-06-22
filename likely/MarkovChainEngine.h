// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MARKOV_CHAIN_ENGINE
#define LIKELY_MARKOV_CHAIN_ENGINE

#include "likely/types.h"
#include "likely/AbsEngine.h"

#include "boost/function.hpp"

namespace likely {
    class Random;
	class MarkovChainEngine : public AbsEngine {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters.
		MarkovChainEngine(FunctionPtr f, GradientCalculatorPtr gc, int nPar,
            std::string const &algorithm);
		virtual ~MarkovChainEngine();
        // Searches for a minimum by taking a sequence of random steps.
        FunctionMinimumPtr minimize(Parameters const &initial, Parameters const &errors,
            double prec, int maxSteps, int acceptsPerParam, int maxTrialsPerParam);
	    // Generates samples using a FunctionMinimum's covariance to specify the proposal
	    // function until the specified number of trials have been accepted or the specified
	    // maximum number of trials has been generated. Use maxTrials=0 for unlimited trials.
	    // Updates the function minimum with an improved estimate and returns the total
	    // number of trials generated (including duplicates after a trial is rejected).
	    // Provide an optional callback to see all trial steps. The callback parameters
	    // are the trial parameters, the function value at these parameters, and a boolean
	    // to flag if the trial is accepted.
        typedef boost::function<void (Parameters const&, double, bool)> Callback;
        int generate(FunctionMinimumPtr fmin, int nAccepts, int maxTrials,
            Callback callback = Callback());
        // Registers our named methods.
        static void registerMarkovChainEngineMethods();
	private:
        int _nPar;
        FunctionPtr _f;
        bool _haveMinimum;
        double _minNLL;
        Parameters _current, _trial, _minParams;
        PackedCovariance _covariance;
        Random &_random;
	}; // MarkovChainEngine
} // likely

#endif // LIKELY_MARKOV_CHAIN_ENGINE
