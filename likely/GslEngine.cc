// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslEngine.h"
#include "likely/GslErrorHandler.h"
#include "likely/RuntimeError.h"

#include <cstdio>

namespace local = likely;

local::GslEngine::GslEngine(Function f, int nPar)
: _nPar(nPar), _f(f)
{
    _func.n = nPar;
    _func.f = _evaluate;
    _func.params = 0;
    _functionStack.push(Binding(f,Parameters(nPar)));
}

local::GslEngine::~GslEngine() {
    _functionStack.pop();
}

void local::GslEngine::minimize(Method method,
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
    //const gsl_multimin_fminimizer_type *T(gsl_multimin_fminimizer_nmsimplex2);
    gsl_multimin_fminimizer *state(gsl_multimin_fminimizer_alloc(method,_nPar));
    gsl_multimin_fminimizer_set(state, &_func, gsl_initial, gsl_errors);
    // Do the minimization...
    int nIterations(0);
    while(nIterations++ < maxIterations) {
        if(gsl_multimin_fminimizer_iterate(state)) break;
        double size(gsl_multimin_fminimizer_size(state));
        if(gsl_multimin_test_size(size,minSize) != GSL_CONTINUE) break;
    }
    // Clean up.
    gsl_vector_free(gsl_errors);
    gsl_vector_free(gsl_initial);
}

double local::GslEngine::operator()(Parameters const& pValues) const {
    // This method is not intended to be a streamlined way to call our function,
    // but rather a way to exercise and test the function stack machinery.
    
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::operator()");
    if(pValues.size() != _nPar) {
        throw RuntimeError("GslEngine: function evaluated with wrong number of parameters.");
    }
    gsl_vector *v(gsl_vector_alloc(_nPar));
    for(int i = 0; i < _nPar; ++i) {
        gsl_vector_set(v,i,pValues[i]);
    }
    double result(_evaluate(v,0));
    gsl_vector_free(v);
    return result;
}

double local::GslEngine::_evaluate(const gsl_vector *v, void *params) {
    Binding bound(_functionStack.top());
    Function &f(bound.first);
    Parameters &values(bound.second);
    int nPar(values.size());
    for(int i = 0; i < nPar; ++i) values[i] = gsl_vector_get(v,i);
    return f(values);
}

std::stack<local::GslEngine::Binding>
    local::GslEngine::_functionStack = std::stack<local::GslEngine::Binding>();