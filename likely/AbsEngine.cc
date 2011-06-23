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
Parameters const &initial, Parameters const &errors, std::string const &methodName,
double precision, long maxIterations) {
    // Check that the input vectors have the same length.
    int nPar(initial.size());
    if(errors.size() != nPar) {
        throw RuntimeError(
            "findMinimum: initial parameter and error vectors have different sizes.");
    }
    // Create a new engine for this function.
    AbsEnginePtr engine = getEngine(methodName,f,gc,nPar);
    // Run the algorithm.
    FunctionMinimumPtr fmin(engine->minimumFinder(initial,errors,precision,maxIterations));
    // Save the evaluation counts.
    lastMinEvalCount = engine->getEvalCount();
    lastMinGradCount = engine->getGradCount();
    return fmin;
}

local::FunctionMinimumPtr local::findMinimum(FunctionPtr f,
Parameters const &initial, Parameters const &errors, std::string const &methodName,
double precision, long maxIterations) {
    // Use a null gradient calculator.
    GradientCalculatorPtr gc;
    return findMinimum(f,gc,initial,errors,methodName,precision,maxIterations);
}

// Initialize the global evaluation counters.
long local::lastMinEvalCount = 0, local::lastMinGradCount = 0;
