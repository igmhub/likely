// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MarkovChainEngine.h"
#include "likely/FunctionMinimum.h"
#include "likely/Random.h"
#include "likely/RuntimeError.h"

#include <cmath>

namespace local = likely;

local::MarkovChainEngine::MarkovChainEngine(FunctionPtr f, int nPar)
: _f(f), _nPar(nPar), _random(Random::instance())
{
    if(_nPar <= 0) {
        throw RuntimeError("MarkovChainEngine: number of parameters must be > 0.");
    }
}

local::MarkovChainEngine::~MarkovChainEngine() { }

local::Parameters local::MarkovChainEngine::advance(FunctionMinimumPtr fmin, int nSamples) {
    Parameters current(fmin->getParameters()),trial(_nPar);
    double currentNLL((*_f)(current));
    while(nSamples--) {
        bool accepted(false);
        while(!accepted) {
            // Take a trial step sampled from the input minimum's covariance.
            fmin->setRandomParameters(trial);
            double trialNLL((*_f)(trial));
            // calculate log( L(trial)/L(current) )
            double logProbRatio(currentNLL-trialNLL);
            if(logProbRatio >= 0 || _random.getUniform() < std::exp(logProbRatio)) {
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