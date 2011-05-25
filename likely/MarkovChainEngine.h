// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MARKOV_CHAIN_ENGINE
#define LIKELY_MARKOV_CHAIN_ENGINE

#include "likely/types.h"

#include "boost/random/mersenne_twister.hpp"

#include "boost/function.hpp"

namespace likely {
	class MarkovChainEngine {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters.
		MarkovChainEngine(Function f, int nPar);
		virtual ~MarkovChainEngine();
        Parameters advance(Parameters const &initial, Parameters const &errors,
            int nSamples = 1);
        // Sets the random seed for generating subsequent samples.
        void setSeed(int seedValue);
	private:
        int _nPar;
        Function _f;
        boost::mt19937 _uniform;
        boost::function<double ()> _gauss;
	}; // MarkovChainEngine
} // likely

#endif // LIKELY_MARKOV_CHAIN_ENGINE
