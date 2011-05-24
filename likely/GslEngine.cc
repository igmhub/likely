// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslEngine.h"
#include "likely/GslErrorHandler.h"
#include "likely/RuntimeError.h"

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

double local::GslEngine::operator()(Parameters const& pValues) const {
    // This method is not intended to be a streamlined way to call our function,
    // but rather a way to exercise and test the function stack machinery.
    if(pValues.size() != _nPar) {
        throw RuntimeError("GslEngine: function evaluated with wrong number of parameters.");
    }
    gsl_vector *v = gsl_vector_alloc(_nPar);
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