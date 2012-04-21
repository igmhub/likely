// Created 22-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/MinuitEngine.h"
#include "likely/FunctionMinimum.h"
#include "likely/EngineRegistry.h"
#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"

#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"

#include "boost/format.hpp"
#include "boost/functional/factory.hpp"
#include "boost/bind.hpp"

#include <cmath>
#include <iostream>

namespace local = likely;
namespace mn = ROOT::Minuit2;

local::MinuitEngine::MinuitEngine(FunctionPtr f, GradientCalculatorPtr gc,
int nPar, std::string const &algorithm)
: _nPar(nPar), _f(f), _gc(gc), _initialState(new mn::MnUserParameterState())
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
    // Select the requested algorithm (last parameter is the MnStrategy value)
    bool useGradient(false);
    int fast(0), normal(1);
    if(algorithm == "simplex") {
        minimumFinder = boost::bind(
            &MinuitEngine::minimize<mn::SimplexMinimizer>,this,_1,_2,_3,_4,normal);
        useGradient = false;
    }
    else if(algorithm == "vmetric") {
        minimumFinder = boost::bind(
            &MinuitEngine::minimize<mn::VariableMetricMinimizer>,this,_1,_2,_3,_4,normal);
        useGradient = false;
    }
    else if(algorithm == "vmetric_fast") {
        minimumFinder = boost::bind(
            &MinuitEngine::minimize<mn::VariableMetricMinimizer>,this,_1,_2,_3,_4,fast);
        useGradient = false;
    }
    else if(algorithm == "vmetric_grad") {
        minimumFinder = boost::bind(
            &MinuitEngine::minimize<mn::VariableMetricMinimizer>,this,_1,_2,_3,_4,normal);
        useGradient = true;
    }
    else if(algorithm == "vmetric_grad_fast") {
        minimumFinder = boost::bind(
            &MinuitEngine::minimize<mn::VariableMetricMinimizer>,this,_1,_2,_3,_4,fast);
        useGradient = true;
    }
    else {
        throw RuntimeError("MinuitEngine: unknown algorithm '" + algorithm + "'");
    }
    if(useGradient) {
        // Check that we have a gradient calculator to use.
        if(!_gc) {
            throw RuntimeError(
                "MinuitEngine: selected algorithm needs a gradient calculator.");
        }
    }
}

local::MinuitEngine::~MinuitEngine() { }

double local::MinuitEngine::operator()(Parameters const &pValues) const {
    if(pValues.size() != _nPar) {
        throw RuntimeError(
            "MinuitEngine: function evaluated with wrong number of parameters.");
    }
    incrementEvalCount();
    return (*_f)(pValues);
}

local::Gradient local::MinuitEngine::Gradient(Parameters const& pValues) const {
    local::Gradient grad(_nPar);
    incrementGradCount();
    (*_gc)(pValues,grad);
    return grad;
}

double local::MinuitEngine::Up() const {
    // Assumes that Fcn returns -logL(p) and that we are interested in 1-sigma errors.
    return 0.5;
}

void local::MinuitEngine::_setInitialState(
Parameters const &initial, Parameters const &errors) {
    // Check for the expected input vector sizes.
    if(initial.size() != _nPar) {
        throw RuntimeError(
            "MinuitEngine: got unexpected number of initial parameter values.");
    }
    if(errors.size() != _nPar) {
        throw RuntimeError(
            "MinuitEngine: got unexpected number of initial parameter errors.");
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
local::FunctionMinimumPtr local::MinuitEngine::minimize(Parameters const &initial,
Parameters const &errors, double prec, int maxfcn, int strategy) {
    _setInitialState(initial,errors);
    // Minuit converts maxfcn = 0 to 200 + 100*npar + 5*npar*npar, but we want to
    // interpret zero as effectively unlimited calls.
    if(maxfcn == 0) maxfcn = 100*_nPar*_nPar;
    // Downscale the requested precision to convert it to an EDM tolerance value.
    double edmTolerance(1*std::sqrt(prec));
    // Do the minimization using the templated class, which is assumed to have
    // a default constructor and provide a Minimize method, e.g., a subclass
    // of ROOT::Minuit2::FunctionMinimizer.
    T algorithm;
    mn::FunctionMinimum mnmin = _gc ?
        algorithm.Minimize(*this,
            *_initialState, mn::MnStrategy(strategy), maxfcn, edmTolerance) :
        algorithm.Minimize((ROOT::Minuit2::FCNBase const&)(*this),
            *_initialState, mn::MnStrategy(strategy), maxfcn, edmTolerance);
    // Transfer the minimization results from the Minuit-specific return object
    // to our engine-neutral return object.
    FunctionMinimumPtr fmin;
    if(mnmin.HasValidParameters()) {
        if(mnmin.HasValidCovariance()) {
            // The Minuit packing of covariance matrix elements is directly
            // compatible with what the FunctionMinimum ctor expects.
            CovarianceMatrixCPtr covariance(new CovarianceMatrix(mnmin.UserCovariance().Data()));
            fmin.reset(new FunctionMinimum(mnmin.Fval(),
                mnmin.UserParameters().Params(),covariance));
        }
        else {
            fmin.reset(new FunctionMinimum(mnmin.Fval(),
                mnmin.UserParameters().Params()));
        }
    }
    return fmin;
}

// Explicit template instantiations.
template local::FunctionMinimumPtr
    local::MinuitEngine::minimize<mn::SimplexMinimizer>
    (Parameters const&,Parameters const&,double,int,int);
template local::FunctionMinimumPtr
    local::MinuitEngine::minimize<mn::VariableMetricMinimizer>
    (Parameters const&,Parameters const&,double,int,int);

void local::registerMinuitEngineMethods() {
    static bool registered = false;
    if(registered) return;
    // Create a function object that constructs a MinuitEngine with parameters
    // (FunctionPtr f, GradientCalculatorPtr gc, int npar, std::string const &methodName).
    EngineFactory factory = boost::bind(boost::factory<MinuitEngine*>(),_1,_2,_3,_4);
    // Register our minimization methods.
    getEngineRegistry()["mn2"] = factory;
    // Return a dummy value so that we can be called at program startup.
    registered = true;
}
