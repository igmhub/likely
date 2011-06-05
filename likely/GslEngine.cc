// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslEngine.h"
#include "likely/GslErrorHandler.h"
#include "likely/FunctionMinimum.h"
#include "likely/RuntimeError.h"

#include "boost/functional/factory.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <cmath>

namespace local = likely;

local::GslEngine::GslEngine(FunctionPtr f, int nPar, std::string const &algorithm)
: _nPar(nPar), _f(f)
{
    if(_nPar <= 0) {
        throw RuntimeError("GslEngine: number of parameters must be > 0.");
    }
    // Bind this function to our GSL global function
    _func.n = nPar;
    _func.f = _evaluate;
    _func.params = 0;
    _params = Parameters(nPar);
    _getFunctionStack().push(Binding(_f,boost::ref(_params),_gc,boost::ref(_grad)));
    // Select the requested algorithm.
    if(algorithm == "simplex2") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex2,_1,_2,_3,_4);
    }
    else if(algorithm == "simplex2rand") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex2rand,_1,_2,_3,_4);
    }
    else {
        throw RuntimeError("GslEngine: unknown algorithm '" + algorithm + "'");
    }
}

local::GslEngine::GslEngine(FunctionPtr f, GradientCalculatorPtr gc, int nPar,
std::string const &algorithm)
: _nPar(nPar), _f(f), _gc(gc)
{
    if(_nPar <= 0) {
        throw RuntimeError("GslEngine: number of parameters must be > 0.");
    }
    // Bind this function and its gradient calculator to our GSL global function
    _funcWithGradient.n = nPar;
    _funcWithGradient.f = _evaluate;
    _funcWithGradient.df = _evaluateGradient;
    _funcWithGradient.fdf = _evaluateBoth;
    _funcWithGradient.params = 0;
    _params = Parameters(nPar);
    _grad = Gradient(nPar);
    _getFunctionStack().push(Binding(_f,boost::ref(_params),_gc,boost::ref(_grad)));
    // Select the requested algorithm.
    if(algorithm == "conjugate_fr") {
        // ...
    }
    else {
        throw RuntimeError("GslEngine: unknown algorithm '" + algorithm + "'");
    }
}

local::GslEngine::~GslEngine() {
    _getFunctionStack().pop();
}

local::FunctionMinimumPtr local::GslEngine::minimize(Method method,
Parameters const &initial, Parameters const &errors,
double prec, long maxIterations) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::minimize");
    // Copy the input initial values and errors to GSL vectors.
    gsl_vector *gsl_initial(gsl_vector_alloc(_nPar)), *gsl_errors(gsl_vector_alloc(_nPar));
    for(int i = 0; i < _nPar; ++i) {
        gsl_vector_set(gsl_initial,i,initial[i]);
        gsl_vector_set(gsl_errors,i,errors[i]);
    }
    // Initialize the minimizer
    gsl_multimin_fminimizer *state(gsl_multimin_fminimizer_alloc(method,_nPar));
    gsl_multimin_fminimizer_set(state, &_func, gsl_initial, gsl_errors);
    // Do the minimization...
    long nIterations(0);
    double minSize = (prec > 0) ? std::sqrt(prec) : 1e-3;
    while(0 == maxIterations || nIterations < maxIterations) {
        nIterations++;
        if(gsl_multimin_fminimizer_iterate(state)) break;
        double size(gsl_multimin_fminimizer_size(state));
        if(gsl_multimin_test_size(size,minSize) != GSL_CONTINUE) break;
    }
    // Copy the results into our result object.
    Parameters final(_nPar);
    for(int i = 0; i < _nPar; ++i) {
        final[i] = gsl_vector_get(state->x,i);
    }
    // Initialize our result object.
    FunctionMinimumPtr fmin(new FunctionMinimum(state->fval,final));
    // Clean up.
    gsl_vector_free(gsl_errors);
    gsl_vector_free(gsl_initial);
    gsl_multimin_fminimizer_free(state);

    return fmin;
}

double local::GslEngine::_evaluate(const gsl_vector *v, void *p) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::_evaluate");
    // Get the top function on the stack.
    Binding const& bound(_useTopBinding(v));
    FunctionPtr f(bound.get<0>());
    Parameters &values(bound.get<1>());
    // Call the function and return its value.
    return (*f)(values);
}

void local::GslEngine::_evaluateGradient(const gsl_vector *v, void *p, gsl_vector *g) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::_evaluateGradient");
    // Get the top function on the stack.
    Binding const& bound(_useTopBinding(v));
    Parameters &values(bound.get<1>());
    GradientCalculatorPtr gc(bound.get<2>());
    Gradient &grad(bound.get<3>());
    // Fill our cached gradient object.
    (*gc)(values,grad);
    // Copy the gradient components to the GSL vector provided.
    int nPar(values.size());
    for(int i = 0; i < nPar; ++i) gsl_vector_set(g,i,grad[i]);
}

void local::GslEngine::_evaluateBoth(const gsl_vector *v, void *p,
double *fval, gsl_vector *g) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::_evaluateBoth");
    // Use the top function on the stack.
    Binding const& bound(_useTopBinding(v));
    FunctionPtr f(bound.get<0>());
    Parameters &values(bound.get<1>());
    GradientCalculatorPtr gc(bound.get<2>());
    Gradient &grad(bound.get<3>());
    // Call the function and save the result.
    *fval = (*f)(values);
    // Fill our cached gradient object.
    (*gc)(values,grad);
    // Copy the gradient components to the GSL vector provided.
    int nPar(values.size());
    for(int i = 0; i < nPar; ++i) gsl_vector_set(g,i,grad[i]);    
}

local::GslEngine::Binding const& local::GslEngine::_useTopBinding(const gsl_vector *v) {
    // Get the top function on the stack.
    Binding &bound(_getFunctionStack().top());
    // Lookup the cached parameter vector for this function.
    Parameters &values(bound.get<1>());
    // Copy the input GSL vector to our cached Parameters object.
    int nPar(values.size());
    for(int i = 0; i < nPar; ++i) values[i] = gsl_vector_get(v,i);
    return bound;
}

std::stack<local::GslEngine::Binding> &local::GslEngine::_getFunctionStack() {
    static std::stack<Binding> *stack = new std::stack<Binding>();
    return *stack;
}

bool local::GslEngine::registerGslEngineMethods() {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::registerEngineMethods");
    // Create a function object that constructs a GslEngine with parameters
    // (Function f, int npar, std::string const &methodName).
    AbsEngine::Factory factory = boost::bind(boost::factory<GslEngine*>(),_1,_2,_3);
    // Register our minimization methods.
    AbsEngine::getRegistry()["gsl"] = factory;
    // Return a dummy value so that we can be called at program startup.
    return true;
}

bool local::GslEngine::_registered = local::GslEngine::registerGslEngineMethods();
