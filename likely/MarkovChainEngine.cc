// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MarkovChainEngine.h"
#include "likely/FunctionMinimum.h"
#include "likely/Random.h"
#include "likely/RuntimeError.h"

//!!#include "boost/lambda/lambda.hpp"
#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/covariance.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/variates/covariate.hpp"

#include <vector>
#include <algorithm>
#include <cmath>

using namespace boost::accumulators;

typedef accumulator_set<double, stats<
    tag::covariance<double, tag::covariate1> > > Accumulator;
typedef std::vector<Accumulator> Accumulators;

namespace local = likely;

local::MarkovChainEngine::MarkovChainEngine(FunctionPtr f, int nPar)
: _f(f), _nPar(nPar), _current(nPar), _trial(nPar), _minParams(nPar), _haveMinimum(false),
_covariance(nPar*(nPar+1)/2), _random(Random::instance())
{
    if(_nPar <= 0) {
        throw RuntimeError("MarkovChainEngine: number of parameters must be > 0.");
    }
}

local::MarkovChainEngine::~MarkovChainEngine() { }

int local::MarkovChainEngine::generate(FunctionMinimumPtr fmin, int nAccepts,
int maxTrials, Callback callback) {
    //!!using namespace boost::lambda;
    // Set our initial parameters to the estimated function minimum, where the
    // NLW = -log(weight) is zero, by definition.
    Parameters const &initial(fmin->getParameters());
    _current = initial;
    double currentNLL((*_f)(_current)), currentNLW(0);
    // Initialize our minimum tracker, if necessary.
    if(!_haveMinimum) {
        _minParams = _current;
        _minNLL = currentNLL;
        _haveMinimum = true;
    }
    // Initialize our statistics accumulators.
    int nCov(_nPar*(_nPar+1)/2);
    Accumulators accumulators(nCov);
    Parameters residual(_nPar);
    // Loop over the requested samples.
    int nTrials(0),remaining(nAccepts);
    while(remaining > 0 && (maxTrials == 0 || nTrials < maxTrials)) {
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
        // Accumulate covariance statistics...
        nTrials++;
        // Use the initial guess at the minimum for caculating residuals that are
        // hopefully small, to minimize round-off error.
        Accumulators::iterator accumulate(accumulators.begin());
        for(int j = 0; j < _nPar; ++j) {
            double jResidual(_current[j]-initial[j]);
            residual[j] = jResidual; // save for the inner loop
            for(int i = 0; i <= j; ++i) {
                (*accumulate++)(jResidual, covariate1 = residual[i]);
            }
        }
    }
    // Record the best minimum found so far (rather than the sample mean).
    fmin->updateParameters(_minParams, _minNLL);
    // Calculate and record the covariance of the samples we have generated.
    for(int k = 0; k < nCov; ++k) _covariance[k] = covariance(accumulators[k]);
    fmin->updateCovariance(_covariance);
    // Return the number of samples generated.
    return nTrials;
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