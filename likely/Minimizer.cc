// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Minimizer.h"
#include "likely/AbsEngine.h"
#include "likely/RuntimeError.h"

namespace local = likely;

local::Minimizer::Minimizer(Function f, int nPar, std::string const &methodName)
{
    // Look up the requested method by name.
    Registry::iterator found = getRegistry().find(methodName);
    if(found == getRegistry().end()) {
        throw RuntimeError("Minimizer: no such method '" + methodName + "'");
    }
    // Create our engine.
    MethodFactory methodFactory = found->second;
    _engine.reset(methodFactory(f,nPar));
}

local::Minimizer::~Minimizer() { }

local::Parameters local::Minimizer::minimize(
Parameters const& initial, Parameters const &errors) {
    Parameters final(initial);
    return final;
}

local::Minimizer::Registry &local::Minimizer::getRegistry() {
    static Registry *registry = new Registry();
    return *registry;
}
