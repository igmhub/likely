// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Minimizer.h"
#include "likely/AbsEngine.h"

namespace local = likely;

local::Minimizer::Minimizer(Function f, int nPar, std::string const &method)
{ }

local::Minimizer::~Minimizer() { }

local::Parameters local::Minimizer::minimize(
Parameters const& initial, Parameters const &errors) {
    Parameters final(initial);
    return final;
}

local::Minimizer::Registry &local::Minimizer::getRegistry() {
    static Registry *_registry = new Registry();
    return *_registry;
}
