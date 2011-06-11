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
		MarkovChainEngine(FunctionPtr f, int nPar);
		virtual ~MarkovChainEngine();
        // Searches for a minimum by taking a sequence of random steps.
        FunctionMinimumPtr minimize(Parameters const &initial, Parameters const &errors,
            double prec, int maxSteps);
	    // Generates the specified number of samples using a FunctionMinimum's
	    // covariance to specify the trial function. Updates the function minimum with
	    // an improved estimate, if possible, and returns the number of samples accepted.
        typedef boost::function<void (Parameters const&, double, bool)> Callback;
        int generate(FunctionMinimumPtr fmin, Callback callback, int nSamples = 1);
	private:
        int _nPar;
        FunctionPtr _f;
        bool _haveMinimum;
        double _minNLL;
        Parameters _current, _trial, _genSum, _minParams;
        PackedCovariance _genPairSum;
        Random &_random;
	}; // MarkovChainEngine
} // likely

#endif // LIKELY_MARKOV_CHAIN_ENGINE
