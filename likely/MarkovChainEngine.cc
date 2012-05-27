// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MarkovChainEngine.h"
#include "likely/FunctionMinimum.h"
#include "likely/Random.h"
#include "likely/EngineRegistry.h"
#include "likely/CovarianceMatrix.h"
#include "likely/CovarianceAccumulator.h"
#include "likely/RuntimeError.h"

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/covariance.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/variates/covariate.hpp"
#include "boost/functional/factory.hpp"
#include "boost/bind.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace boost::accumulators;

typedef accumulator_set<double, stats<
    tag::covariance<double, tag::covariate1> > > Accumulator;
typedef std::vector<Accumulator> Accumulators;

namespace local = likely;

local::MarkovChainEngine::MarkovChainEngine(FunctionPtr f, GradientCalculatorPtr gc,
FitParameters const &parameters, std::string const &algorithm)
: _f(f), _nPar(parameters.size()), _current(_nPar), _trial(_nPar), _minParams(_nPar), _haveMinimum(false),
_random(Random::instance())
{
    if(_nPar <= 0) {
        throw RuntimeError("MarkovChainEngine: number of parameters must be > 0.");
    }
    if(algorithm == "saunter") {
        minimumFinder = boost::bind(&MarkovChainEngine::minimize,this,
            _1,_2,_3,_4,100,1000);
    }
    else if(algorithm == "stroll") {
        minimumFinder = boost::bind(&MarkovChainEngine::minimize,this,
            _1,_2,_3,_4,50,5000);
    }
    else {
        throw RuntimeError("MarkovChainEngine: unknown algorithm '" + algorithm + "'");
    }
}

local::MarkovChainEngine::~MarkovChainEngine() { }

int local::MarkovChainEngine::generate(FunctionMinimumPtr fmin, int nAccepts,
int maxTrials, Callback callback) {
    // We are using the standard Metropolis-Hastings algorithm here. The only subtlety is
    // that we don't generate trials by taking a random step from our current location,
    // so our proposal pdf Q(p',p) for moving from p (current) to p' (trial) is
    // W(p') = exp(-delta.Cinv.delta/2) and the ratio Q(p,p')/Q(p',p) is not 1, as
    // usually assumed, but W(current)/W(trial).

    // Set our initial parameters to the estimated function minimum, where the
    // NLW = -log(W(current)) is zero, by definition.
    Parameters const &initial(fmin->getParameters());
    _current = initial;
    double currentNLL(fmin->getMinValue()), currentNLW(0);
    // Initialize our minimum tracker, if necessary.
    if(!_haveMinimum) {
        _minParams = _current;
        _minNLL = currentNLL;
        _haveMinimum = true;
    }
    // Initialize our statistics accumulators.
    CovarianceAccumulator accumulator(_nPar);
    Parameters residual(_nPar);
    // Loop over the requested samples.
    int nTrials(0),remaining(nAccepts);
    while(remaining > 0 && (maxTrials == 0 || nTrials < maxTrials)) {
        // Take a trial step sampled from the estimated function minimum's covariance.
        // The setRandomParameters method returns the value -log(W(trial)).
        double trialNLW(fmin->setRandomParameters(_trial));
        // Evaluate the true NLL at this trial point.
        double trialNLL((*_f)(_trial));
        incrementEvalCount();
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
        for(int j = 0; j < _nPar; ++j) {
            residual[j] = _current[j]-initial[j];
        }
        accumulator.accumulate(residual);
    }
    // Record the best minimum found so far (rather than the sample mean).
    fmin->updateParameters(_minParams, _minNLL);
    // Record the covariance of the samples we have generated.
    fmin->updateCovariance(accumulator.getCovariance());
    // Return the number of samples generated.
    return nTrials;
}

local::FunctionMinimumPtr local::MarkovChainEngine::minimize(
Parameters const &initial, Parameters const &errors, double prec, int maxSteps,
int acceptsPerParam, int maxTrialsPerParam) {
    // Build an initial diagonal convariance using the errors provided.
    double fval((*_f)(initial));
    incrementEvalCount();
    boost::shared_ptr<CovarianceMatrix> covariance(new CovarianceMatrix(_nPar));
    for(int k = 0; k < _nPar; ++k) {
        covariance->setCovariance(k,k,errors[k]);
    }
    FunctionMinimumPtr fmin(new FunctionMinimum(fval,initial,covariance));
    // Configure our cycles.
    int nAccepts(acceptsPerParam*_nPar), maxTrials(maxTrialsPerParam*_nPar), trials(0);
    while(maxSteps == 0 || trials < maxSteps) {
        double initialFval(fmin->getMinValue());
        trials += generate(fmin, nAccepts, maxTrials);
        // Check if we have reached the requested "precision"
        fval = fmin->getMinValue();
        if(fval < initialFval && initialFval - fval < prec) break;
    }
    return fmin;
}

void local::registerMarkovChainEngineMethods() {
    static bool registered = false;
    if(registered) return;
    // Create a function object that constructs a MarkovChainEngine with parameters
    // (FunctionPtr f, GradientCalculatorPtr gc, int npar, std::string const &methodName).
    EngineFactory factory = boost::bind(boost::factory<MarkovChainEngine*>(),_1,_2,_3,_4);
    // Register our minimization methods.
    getEngineRegistry()["mc"] = factory;
    registered = true;
}
