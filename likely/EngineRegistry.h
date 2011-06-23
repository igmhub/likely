// Created 23-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ENGINE_REGISTRY
#define LIKELY_ENGINE_REGISTRY

#include "likely/types.h"

#include "boost/function.hpp"

#include <string>
#include <map>

namespace likely {

    // Declares a global registry for creating engines by name.
    class AbsEngine;
    typedef boost::function<AbsEngine* (FunctionPtr, GradientCalculatorPtr,
        int, std::string const&)> EngineFactory;
    typedef std::map<std::string, EngineFactory> EngineRegistry;
    
    // Returns the unique engine registry.
    EngineRegistry &getEngineRegistry();
    
    // Returns a smart pointer to a newly-created engine inferred from a method name
    // of the form <engine>::<algorithm> or throws a RuntimeError. The engine is
    // created using the specified function, gradient calculator and number of parameters.
    AbsEnginePtr getEngine(std::string const methodName,
        FunctionPtr f, GradientCalculatorPtr gc, int nPar);

} // likely

#endif // LIKELY_ENGINE_REGISTRY
