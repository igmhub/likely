// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MarkovChainEngine.h"
#include "likely/RuntimeError.h"

#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include <cmath>

namespace local = likely;

local::MarkovChainEngine::MarkovChainEngine(Function f, int nPar)
: _f(f), _nPar(nPar)
{
    if(_nPar <= 0) {
        throw RuntimeError("MarkovChainEngine: number of parameters must be > 0.");
    }
    _gauss = boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
        (_uniform, boost::normal_distribution<>(0,1));
}

local::MarkovChainEngine::~MarkovChainEngine() { }

void local::MarkovChainEngine::setSeed(int seedValue) {
    _uniform.seed(seedValue);
}

local::Parameters local::MarkovChainEngine::advance(
Parameters const &initial, Parameters const &errors,int nSamples) {
    Parameters current(initial),trial(_nPar);
    double currentNLL(_f(current));
    while(nSamples--) {
        bool accepted(false);
        while(!accepted) {
            // Take a trial step by adding Gaussian offsets (scaled by the input errors)
            // to the current parameter values.
            for(int i = 0; i < _nPar; ++i) {
                trial[i] = current[i] + _gauss()*errors[i];
            }
            double trialNLL(_f(trial));
            // calculate log( L(trial)/L(current) )
            double logProbRatio(currentNLL-trialNLL);
            if(logProbRatio >= 0 || _uniform() < std::exp(logProbRatio)) {
                // Accept the trial step.
                current = trial;
                currentNLL = trialNLL;
                accepted = true;
            }
        }
    }
    return current;
}

local::FunctionMinimumPtr local::MarkovChainEngine::minimize(
Parameters const &initial, Parameters const &errors, double prec, int maxSteps) {
    /*
    - initialize covar as diag of input errors
    - advance should use covar matrix as input
    - advance takes currentNLL value as input
    - advance keeps track of min NLL value seen so far
    - advance has options to update covar
    */
}