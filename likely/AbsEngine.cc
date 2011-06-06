// Created 28-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/AbsEngine.h"
#include "likely/RuntimeError.h"
#include "likely/FunctionMinimum.h"

#include "boost/regex.hpp"

namespace local = likely;

local::AbsEngine::AbsEngine()
: _evalCount(0), _gradCount(0)
{ }

local::AbsEngine::~AbsEngine() { }

local::AbsEngine::Registry &local::AbsEngine::getRegistry() {
    static Registry *registry = new Registry();
    return *registry;
}

local::AbsEngine::RegistryWithGC &local::AbsEngine::getRegistryWithGC() {
    static RegistryWithGC *registry = new RegistryWithGC();
    return *registry;
}

local::ParsedMethodName local::parseMethodName(std::string const &methodName) {
    static boost::regex pattern("([a-z0-9]+)::([0-9a-z_]+)");
    boost::smatch found;
    if(!boost::regex_match(methodName,found,pattern)) {
        throw RuntimeError("parseMethodName: invalid method name '" + methodName + "'");
    }
    std::string engine(found[1].first,found[1].second);
    std::string algorithm(found[2].first,found[2].second);
    return ParsedMethodName(engine,algorithm);
}

local::FunctionMinimumPtr local::findMinimum(FunctionPtr f,
Parameters const &initial, Parameters const &errors, std::string const &methodName,
double precision, long maxIterations) {
    // Check that the input vectors have the same length.
    int nPar(initial.size());
    if(errors.size() != nPar) {
        throw RuntimeError(
            "findMinimum: initial parameter and error vectors have different sizes.");
    }
    // Parse the method name, which should have the form <engine>::<algorithm>
    ParsedMethodName parsed(parseMethodName(methodName));
    // Lookup the factory that creates this type of engine.
    AbsEngine::Registry::iterator found = AbsEngine::getRegistry().find(parsed.first);
    if(found == AbsEngine::getRegistry().end()) {
        throw RuntimeError("findMinimum: no such engine '" + methodName + "'");
    }
    // Create a new engine for this function.
    AbsEngine::Factory factory = found->second;
    boost::scoped_ptr<AbsEngine> engine(factory(f,nPar,parsed.second));
    // Run the algorithm.
    FunctionMinimumPtr fmin(engine->minimumFinder(initial,errors,precision,maxIterations));
    // Save the evaluation counts.
    lastMinEvalCount = engine->getEvalCount();
    lastMinGradCount = 0;
    return fmin;
}

local::FunctionMinimumPtr local::findMinimum(FunctionPtr f, GradientCalculatorPtr gc,
Parameters const &initial, Parameters const &errors, std::string const &methodName,
double precision, long maxIterations) {
    // Check that the input vectors have the same length.
    int nPar(initial.size());
    if(errors.size() != nPar) {
        throw RuntimeError(
            "findMinimum: initial parameter and error vectors have different sizes.");
    }
    // Parse the method name, which should have the form <engine>::<algorithm>
    ParsedMethodName parsed(parseMethodName(methodName));
    // Lookup the factory that creates this type of engine.
    AbsEngine::RegistryWithGC::iterator found =
        AbsEngine::getRegistryWithGC().find(parsed.first);
    if(found == AbsEngine::getRegistryWithGC().end()) {
        throw RuntimeError("findMinimum: no such engine '" + methodName + "'");
    }
    // Create a new engine for this function.
    AbsEngine::FactoryWithGC factory = found->second;
    boost::scoped_ptr<AbsEngine> engine(factory(f,gc,nPar,parsed.second));
    // Run the algorithm.
    FunctionMinimumPtr fmin(engine->minimumFinder(initial,errors,precision,maxIterations));
    // Save the evaluation counts.
    lastMinEvalCount = engine->getEvalCount();
    lastMinGradCount = engine->getGradCount();
    return fmin;
}

// Initialize the global evaluation counters.
long local::lastMinEvalCount = 0, local::lastMinGradCount = 0;
