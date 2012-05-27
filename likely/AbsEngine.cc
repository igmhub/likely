// Created 28-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/AbsEngine.h"
#include "likely/RuntimeError.h"
#include "likely/FunctionMinimum.h"
#include "likely/EngineRegistry.h"

namespace local = likely;

local::AbsEngine::AbsEngine()
: _evalCount(0), _gradCount(0)
{ }

local::AbsEngine::~AbsEngine() { }

local::FunctionMinimumPtr local::findMinimum(FunctionPtr f, GradientCalculatorPtr gc,
FitParameters const &parameters, std::string const &methodName,
double precision, long maxIterations) {
    // Create a new engine for this function.
    AbsEnginePtr engine = getEngine(methodName,f,gc,parameters.size());
    // Run the algorithm.
    std::vector<double> initial,errors;
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        initial.push_back(iter->getValue());
        errors.push_back(iter->getError());
    }
    FunctionMinimumPtr fmin(engine->minimumFinder(initial,errors,precision,maxIterations));
    // Save the evaluation counts.
    lastMinEvalCount = engine->getEvalCount();
    lastMinGradCount = engine->getGradCount();
    return fmin;
}

local::FunctionMinimumPtr local::findMinimum(FunctionPtr f,
FitParameters const &parameters, std::string const &methodName,
double precision, long maxIterations) {
    // Use a null gradient calculator.
    GradientCalculatorPtr gc;
    return findMinimum(f,gc,parameters,methodName,precision,maxIterations);
}

// Initialize the global evaluation counters.
long local::lastMinEvalCount = 0, local::lastMinGradCount = 0;
