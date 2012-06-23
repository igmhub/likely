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
FitParameters const &parameters, std::string const &algorithm, RandomPtr random)
: _f(f), _random(random)
{
    _nParam = countFitParameters(parameters,false);
    _nFloating = countFitParameters(parameters,true);
    if(0 == _nFloating) {
        throw RuntimeError("MarkovChainEngine: number of floating parameters must be > 0.");
    }
    if(algorithm == "saunter") {
        minimumFinder = boost::bind(&MarkovChainEngine::minimize,this,
            _1,_2,_3,100,1000);
    }
    else if(algorithm == "stroll") {
        minimumFinder = boost::bind(&MarkovChainEngine::minimize,this,
            _1,_2,_3,50,5000);
    }
    else {
        throw RuntimeError("MarkovChainEngine: unknown algorithm '" + algorithm + "'");
    }
    if(!_random) _random = Random::instance();
}

local::MarkovChainEngine::~MarkovChainEngine() { }

int local::MarkovChainEngine::generate(FunctionMinimumPtr fmin, int nAccepts,
int maxTrials, Callback callback, int callbackInterval) const {
    // We are using the standard Metropolis-Hastings algorithm here. The only subtlety is
    // that we don't generate trials by taking a random step from our current location,
    // so our proposal pdf Q(p',p) for moving from p (current) to p' (trial) is
    // W(p') = exp(-delta.Cinv.delta/2) and the ratio Q(p,p')/Q(p',p) is not 1, as
    // usually assumed, but W(current)/W(trial).

    // Set our initial parameters to the estimated function minimum, where the
    // NLW = -log(W(current)) is zero, by definition.
    Parameters current(fmin->getParameters());
    double currentNLL(fmin->getMinValue());
    double currentNLW(0);
    
    // Our starting point is our current best guess at the minimum.
    Parameters minParams(current);
    double minNLL(currentNLL);

    // Remember the initial floating parameters for calculating residuals later.
    Parameters initialFloating(fmin->getParameters(true)), residual;

    // Initialize our covariance accumulator.
    CovarianceAccumulator accumulator(_nFloating);

    // Loop over the requested samples.
    int nTrials(0),remaining(nAccepts);
    Parameters trial;
    while(remaining > 0 && (maxTrials == 0 || nTrials < maxTrials)) {
        nTrials++;
        // Take a trial step sampled from the estimated function minimum's covariance.
        // The setRandomParameters method returns the value of -log(W(trial)) and 
        // includes any fixed parameters in trial.
        double trialNLW(fmin->setRandomParameters(trial));
        // Evaluate the true NLL at this trial point.
        double trialNLL((*_f)(trial));
        incrementEvalCount();
        // Is this a new minimum?
        if(trialNLL < minNLL) {
            minParams = trial;
            minNLL = trialNLL;
        }
        // Calculate log( L(trial)/L(current) W(current)/W(trial) )
        double logProbRatio(currentNLL-trialNLL-currentNLW+trialNLW);
        // Do we accept this trial step?
        bool accepted(false);
        if(logProbRatio >= 0 || _random->getUniform() < std::exp(logProbRatio)) {
            current = trial;
            currentNLL = trialNLL;
            currentNLW = trialNLW;
            accepted = true;
            remaining--;
        }
        // Invoke the callback now, if any.
        if(callback && (0 == nTrials%callbackInterval)) {
            callback(current, trial, currentNLL, trialNLL, accepted);
        }
        // Accumulate covariance statistics...
        // Use the initial guess at the minimum for calculating residuals of our
        // floating parameters that are hopefully small, to minimize round-off error.
        fmin->filterParameterValues(current,residual);
        for(int j = 0; j < _nFloating; ++j) residual[j] -= initialFloating[j];
        accumulator.accumulate(residual);
    }
    // Record the covariance of the samples we have generated. Do this before updating
    // the parameter values, so that the updated errors are available.
    try {
        // Make sure we have a valid positive-definite matrix before we use it.
        CovarianceMatrixCPtr C = accumulator.getCovariance();
        C->getDeterminant();
        fmin->updateCovariance(C);
    }
    catch(RuntimeError const &e) {
        // Stick with our original covariance estimate for now.
    }
    // Record the best minimum found so far (rather than the sample mean).
    fmin->updateParameterValues(minNLL, minParams);
    // Return the number of samples generated.
    return nTrials;
}

void local::MarkovChainEngine::minimize(FunctionMinimumPtr fmin, double prec, int maxSteps,
int acceptsPerParam, int maxTrialsPerParam) {
    
    // Build an initial diagonal convariance using the fit parameter errors provided.
    Parameters initialErrors = fmin->getErrors(true);
    int nFloating(initialErrors.size());
    CovarianceMatrixPtr covariance(new CovarianceMatrix(nFloating));
    for(int k = 0; k < nFloating; ++k) {
        covariance->setCovariance(k,k,initialErrors[k]*initialErrors[k]);
    }
    fmin->updateCovariance(covariance);
    
    // Configure our cycles.
    int nAccepts(acceptsPerParam*_nFloating), maxTrials(maxTrialsPerParam*_nFloating), trials(0);
    while(maxSteps == 0 || trials < maxSteps) {
        double initialFval(fmin->getMinValue());
        trials += generate(fmin, nAccepts, maxTrials);
        // Check if we have reached the requested "precision"
        double fval = fmin->getMinValue();
        if(fval < initialFval && initialFval - fval < prec) break;
    }
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
