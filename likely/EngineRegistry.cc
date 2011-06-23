// Created 23-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/EngineRegistry.h"
#include "likely/RuntimeError.h"

#include "config.h"
#ifdef HAVE_LIBGSL
#include "likely/GslEngine.h"
#endif
#ifdef HAVE_LIBMINUIT2
#include "likely/MinuitEngine.h"
#endif
#include "likely/MarkovChainEngine.h"

#include "boost/regex.hpp"

namespace local = likely;

local::EngineRegistry &local::getEngineRegistry() {
    static EngineRegistry *registry = new EngineRegistry();
    return *registry;
}

local::AbsEnginePtr local::getEngine(std::string const methodName,
FunctionPtr f, GradientCalculatorPtr gc, int nPar) {
    // Trigger first-time registration of all engine implementations.
#ifdef HAVE_LIBGSL
    registerGslEngineMethods();
#endif
#ifdef HAVE_LIBMINUIT2
    registerMinuitEngineMethods();
#endif
    registerMarkovChainEngineMethods();
    // Parse the method name to split out the fields of <engine>::<algorithm>
    static boost::regex pattern("([a-z0-9]+)::([0-9a-z_]+)");
    boost::smatch parsed;
    if(!boost::regex_match(methodName,parsed,pattern)) {
        throw RuntimeError("getEngine: invalid method name '" + methodName + "'");
    }
    std::string engineName(parsed[1].first,parsed[1].second);
    std::string algorithmName(parsed[2].first,parsed[2].second);
    // Lookup the factory that creates this type of engine.
    EngineRegistry::iterator found = getEngineRegistry().find(engineName);
    if(found == getEngineRegistry().end()) {
        throw RuntimeError("getEngine: no such engine '" + engineName + "'");
    }
    // Create and return a new engine for this function.
    EngineFactory factory = found->second;
    AbsEnginePtr engine(factory(f,gc,nPar,algorithmName));
    return engine;
}