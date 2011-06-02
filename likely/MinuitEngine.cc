// Created 22-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MinuitEngine.h"
#include "likely/RuntimeError.h"
#include "likely/FunctionMinimum.h"

#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/FumiliMinimizer.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"

#include "boost/format.hpp"

#include <iostream>

namespace local = likely;
namespace mn = ROOT::Minuit2;

local::MinuitEngine::MinuitEngine(FunctionPtr f, int nPar)
: _nPar(nPar), _f(f), _initialState(new mn::MnUserParameterState())
{
    if(_nPar <= 0) {
        throw RuntimeError("MinuitEngine: number of parameters must be > 0.");
    }
    boost::format fmt("p%d");
    // Minuit2 crashes during the fit if a parameter is initially defined fixed and
    // then later released, so we always create the parameter as floating with
    // zero error below.
    for(int i = 0; i < _nPar; ++i) {
        _initialState->Add(boost::str(fmt % i),0,0);
    }
}

local::MinuitEngine::MinuitEngine(FunctionPtr f, std::vector<std::string> const &parNames)
: _nPar(parNames.size()), _f(f), _initialState(new mn::MnUserParameterState())
{
    if(_nPar == 0) {
        throw RuntimeError("MinuitEngine: no parameter names specified.");
    }
    // Minuit2 crashes during the fit if a parameter is initially defined fixed and
    // then later released, so we always create the parameter as floating with
    // zero error below.
    for(int i = 0; i < _nPar; ++i) {
        _initialState->Add(parNames[i],0,0);
    }   
}

local::MinuitEngine::~MinuitEngine() { }

double local::MinuitEngine::operator()(Parameters const &pValues) const {
    if(pValues.size() != _nPar) {
        throw RuntimeError("MinuitEngine: function evaluated with wrong number of parameters.");
    }
    return (*_f)(pValues);
}

double local::MinuitEngine::Up() const {
    // Assumes that Fcn returns -logL(p) and that we are interested in 1-sigma errors.
    return 0.5;
}

void local::MinuitEngine::_setInitialState(
Parameters const &initial, Parameters const &errors) {
    // Check for the expected input vector sizes.
    if(initial.size() != _nPar) {
        throw RuntimeError("MinuitEngine: got unexpected number of initial parameter values.");
    }
    if(errors.size() != _nPar) {
        throw RuntimeError("MinuitEngine: got unexpected number of initial parameter errors.");
    }
    // Set the parameter values and error estimates from the input vectors.
    for(int i = 0; i < _nPar; ++i) {
        _initialState->SetValue(i,initial[i]);
        if(errors[i] > 0) {
            _initialState->SetError(i,errors[i]);
            _initialState->Release(i);
        }
        else {
            _initialState->SetError(i,0);
            _initialState->Fix(i);
        }
    }
}

template <class T>
local::FunctionMinimumPtr local::MinuitEngine::minimize(
Parameters const &initial, Parameters const &errors, double toler, int maxfcn) {
    _setInitialState(initial,errors);
    T algorithm;
    mn::MnStrategy strategy(1);
    mn::FunctionMinimum mnmin = algorithm.Minimize(*this, *_initialState, strategy);
    FunctionMinimumPtr fmin(new FunctionMinimum(_nPar,initial));
    return fmin;
}

// Explicit template instantiations.
template local::FunctionMinimumPtr
    local::MinuitEngine::minimize<mn::SimplexMinimizer>
    (Parameters const&,Parameters const&,double,int);
template local::FunctionMinimumPtr
    local::MinuitEngine::minimize<mn::VariableMetricMinimizer>
    (Parameters const&,Parameters const&,double,int);
template local::FunctionMinimumPtr
    local::MinuitEngine::minimize<mn::FumiliMinimizer>
    (Parameters const&,Parameters const&,double,int);

mn::FunctionMinimum
local::MinuitEngine::simplex(Parameters const &initial, Parameters const &errors) {
    _setInitialState(initial,errors);
    // Allocate a SimplexMinimizer if this is the first time we are called.
    if(!_simplex) _simplex.reset(new mn::SimplexMinimizer());    
    // Run the mimizer.
    mn::MnStrategy strategy(1);
    return _simplex->Minimize(*this, *_initialState, strategy);
}

mn::FunctionMinimum
local::MinuitEngine::variableMetric(Parameters const &initial, Parameters const &errors) {
    _setInitialState(initial,errors);
    // Allocate a VariableMetricMinimizer if this is the first time we are called.
    if(!_variableMetric) _variableMetric.reset(new mn::VariableMetricMinimizer());
    // Run the mimizer.
    mn::MnStrategy strategy(1);
    return _variableMetric->Minimize(*this, *_initialState, strategy);
}