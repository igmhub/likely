// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MARKOV_CHAIN_ENGINE
#define LIKELY_MARKOV_CHAIN_ENGINE

#include "likely/types.h"
#include "likely/AbsEngine.h"

namespace likely {
    class Random;
	class MarkovChainEngine : public AbsEngine {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters.
		MarkovChainEngine(FunctionPtr f, int nPar);
		virtual ~MarkovChainEngine();
        // Searches for a minimum by taking a sequence of random steps.
        FunctionMinimumPtr minimize(Parameters const &initial, Parameters const &errors,
            double prec, int maxSteps);
	private:
	    // Generates the specified number of samples using a FunctionMinimum's
	    // covariance to specify the trial function. Stores the final sample in the
	    // parameter vector provided and returns the function value at this point.
        double _advance(FunctionMinimumPtr fmin, Parameters &params, double fVal,
            int nSamples = 1);
        int _nPar;
        FunctionPtr _f;
        Random &_random;
	}; // MarkovChainEngine
} // likely

#endif // LIKELY_MARKOV_CHAIN_ENGINE
