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
        Parameters advance(FunctionMinimumPtr fmin, int nSamples = 1);
        // Searches for a minimum by taking a sequence of random steps.
        FunctionMinimumPtr minimize(Parameters const &initial, Parameters const &errors,
            double prec, int maxSteps);
	private:
        int _nPar;
        FunctionPtr _f;
        Random &_random;
	}; // MarkovChainEngine
} // likely

#endif // LIKELY_MARKOV_CHAIN_ENGINE
