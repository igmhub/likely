// Created 22-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MinuitEngine.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/SimplexMinimizer.h"

namespace local = likely;

local::MinuitEngine::MinuitEngine(Function f)
: _f(f)
{ }

local::MinuitEngine::~MinuitEngine() { }

double local::MinuitEngine::operator()(Parameters const &pValues) const {
    return _f(pValues);
}

double local::MinuitEngine::Up() const {
    return 1;
}

local::Parameters local::MinuitEngine::simplex(
Parameters const &initial, Parameters const &errors) {
    // Allocate a SimplexMinimizer if this is the first time we are called.
    if(!_simplex) _simplex.reset(new ROOT::Minuit2::SimplexMinimizer());
    // Run the mimizer.
    ROOT::Minuit2::FunctionMinimum min = _simplex->Minimize(*this, initial, errors, 1, 1000, 0.1);
    // Extract the parameters at the function minimum.
    Parameters result(min.UserParameters().Params());
    return result;
}