// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MarkovChainEngine.h"
#include "likely/FunctionMinimum.h"
#include "likely/Random.h"
#include "likely/RuntimeError.h"

#include <algorithm>
#include <cmath>

namespace local = likely;

local::MarkovChainEngine::MarkovChainEngine(FunctionPtr f, int nPar)
: _f(f), _nPar(nPar), _current(nPar), _trial(nPar), _minParams(nPar), _haveMinimum(false),
_genSum(nPar), _genPairSum(nPar*(nPar+1)/2), _random(Random::instance())
{
    if(_nPar <= 0) {
        throw RuntimeError("MarkovChainEngine: number of parameters must be > 0.");
    }
}

local::MarkovChainEngine::~MarkovChainEngine() { }

int local::MarkovChainEngine::generate(FunctionMinimumPtr fmin, int nAccepts,
Callback callback) {
    // Set our initial parameters to the estimated function minimum, where the
    // NLW = -log(weight) is zero, by definition.
    _current = fmin->getParameters();
    double currentNLL((*_f)(_current)), currentNLW(0);
    // Initialize our minimum tracker, if necessary.
    if(!_haveMinimum) {
        _minParams = _current;
        _minNLL = currentNLL;
        _haveMinimum = true;
    }
    // Zero our statistics.
    std::fill(_genSum.begin(),_genSum.end(),0);
    std::fill(_genPairSum.begin(),_genPairSum.end(),0);
    // Loop over the requested samples.
    int nSamples(0),remaining(nAccepts);
    while(remaining > 0) {
        // Take a trial step sampled from the estimated function minimum's covariance.
        double trialNLW(fmin->setRandomParameters(_trial));
        // Evaluate the true NLL at this trial point.
        double trialNLL((*_f)(_trial));
        // Is this a new minimum?
        if(trialNLL < _minNLL) {
            _minParams = _trial;
            _minNLL = trialNLL;
        }
        // Calculate log( L(trial)/L(current) W(current)/W(trial) )
        double logProbRatio(currentNLL-trialNLL-currentNLW+trialNLW);
        // Do we accept this trial step?
        if(logProbRatio >= 0 || _random.getUniform() < std::exp(logProbRatio)) {
            std::swap(_current,_trial);
            currentNLL = trialNLL;
            currentNLW = trialNLW;
            if(callback) callback(_trial, trialNLL, true);
            remaining--;
        }
        else {
            if(callback) callback(_trial, trialNLL, false);
        }
        // Accumulate covariance statistics.
        nSamples++;
        PackedCovariance::iterator pairSum(_genPairSum.begin());
        for(int j = 0; j < _nPar; ++j) {
            double jCurrent(_current[j]);
            _genSum[j] += jCurrent;
            for(int i = 0; i <= j; ++i) {
                *pairSum++ += _current[i]*jCurrent;
            }
        }
    }
    // Record the best minimum found so far (rather than the sample mean).
    fmin->updateParameters(_minParams, _minNLL);
    // Calculate the covariance of the samples we have generated.
    //int index(0);
    double nSamplesSq(nSamples*nSamples);
    PackedCovariance::iterator pairSum(_genPairSum.begin());
    for(int j = 0; j < _nPar; ++j) {
        for(int i = 0; i <= j; ++i) {
            *pairSum/= nSamples;
            *pairSum++ -= _genSum[i]*_genSum[j]/nSamplesSq;
        }
    }
    // Update our guess at the function minimum. 
    fmin->updateCovariance(_genPairSum);
    // Return the number of samples accepted.
    return nSamples;
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