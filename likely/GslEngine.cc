// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslEngine.h"
#include "likely/GslErrorHandler.h"
#include "likely/FunctionMinimum.h"
#include "likely/RuntimeError.h"

#include "boost/functional/factory.hpp"
#include "boost/bind.hpp"

namespace local = likely;

local::GslEngine::GslEngine(Function f, int nPar, std::string const &algorithm)
: _nPar(nPar), _f(f)
{
    if(_nPar <= 0) {
        throw RuntimeError("GslEngine: number of parameters must be > 0.");
    }
    // Bind this function to our GSL global function
    _func.n = nPar;
    _func.f = _evaluate;
    _func.params = 0;
    getFunctionStack().push(Binding(f,Parameters(nPar)));
    // Select the requested algorithm.
    if(algorithm == "simplex2") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex2,_1,_2,1e-3,1000);
    }
    else if(algorithm == "simplex2rand") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex2rand,_1,_2,1e-3,1000);
    }
    else {
        throw RuntimeError("GslEngine: unknown algorithm '" + algorithm + "'");
    }
}

local::GslEngine::~GslEngine() {
    getFunctionStack().pop();
}

local::FunctionMinimumPtr local::GslEngine::minimize(Method method,
Parameters const &initial, Parameters const &errors,
double minSize, int maxIterations) {
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
    int nIterations(0);
    while(nIterations++ < maxIterations) {
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

/*
double local::GslEngine::operator()(Parameters const& pValues) const {
    // This method is not intended to be a streamlined way to call our function,
    // but rather a way to exercise and test the function stack machinery.
    
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::operator()");
    if(pValues.size() != _nPar) {
        throw RuntimeError(
            "GslEngine: function evaluated with wrong number of parameters.");
    }
    gsl_vector *v(gsl_vector_alloc(_nPar));
    for(int i = 0; i < _nPar; ++i) {
        gsl_vector_set(v,i,pValues[i]);
    }
    double result(_evaluate(v,0));
    gsl_vector_free(v);
    return result;
}
*/

double local::GslEngine::_evaluate(const gsl_vector *v, void *params) {
    Binding bound(getFunctionStack().top());
    Function &f(bound.first);
    Parameters &values(bound.second);
    int nPar(values.size());
    for(int i = 0; i < nPar; ++i) values[i] = gsl_vector_get(v,i);
    return f(values);
}

std::stack<local::GslEngine::Binding> &local::GslEngine::getFunctionStack() {
    static std::stack<Binding> *stack = new std::stack<Binding>();
    return *stack;
}

bool local::GslEngine::registerGslEngineMethods() {
    // Create a function object that constructs a GslEngine with parameters
    // (Function f, int npar, std::string const &methodName).
    AbsEngine::Factory factory = boost::bind(boost::factory<GslEngine*>(),_1,_2,_3);
    // Register our minimization methods.
    AbsEngine::getRegistry()["gsl"] = factory;
    // Return a dummy value so that we can be called at program startup.
    return true;
}

bool local::GslEngine::_registered = local::GslEngine::registerGslEngineMethods();
