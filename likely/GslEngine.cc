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

local::GslEngine::GslEngine(FunctionPtr f, GradientCalculatorPtr gc,
int nPar, std::string const &algorithm)
: _nPar(nPar), _f(f), _gc(gc)
{
    if(_nPar <= 0) {
        throw RuntimeError("GslEngine: number of parameters must be > 0.");
    }
    // Use a fixed line minimization tolerance of 0.1 for methods that use derivatives:
    // http://www.gnu.org/software/gsl/manual/html_node/
    //   Initializing-the-Multidimensional-Minimizer.html
    double lineMinTol(0.1);
    // Select the requested algorithm.
    bool useGradient(true);
    if(algorithm == "nmsimplex") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex,_1,_2,_3,_4);
        useGradient = false;
    }
    else if(algorithm == "nmsimplex2") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex2,_1,_2,_3,_4);
        useGradient = false;
    }
    else if(algorithm == "nmsimplex2rand") {
        minimumFinder = boost::bind(&GslEngine::minimize,this,
            gsl_multimin_fminimizer_nmsimplex2rand,_1,_2,_3,_4);
        useGradient = false;
    }
    else if(algorithm == "conjugate_fr") {
        minimumFinder = boost::bind(&GslEngine::minimizeWithGradient,this,
            gsl_multimin_fdfminimizer_conjugate_fr,_1,_2,_3,_4,lineMinTol);
        useGradient = true;
    }
    else if(algorithm == "conjugate_pr") {
        minimumFinder = boost::bind(&GslEngine::minimizeWithGradient,this,
            gsl_multimin_fdfminimizer_conjugate_pr,_1,_2,_3,_4,lineMinTol);
        useGradient = true;
    }
    else if(algorithm == "vector_bfgs") {
        minimumFinder = boost::bind(&GslEngine::minimizeWithGradient,this,
            gsl_multimin_fdfminimizer_vector_bfgs,_1,_2,_3,_4,lineMinTol);
        useGradient = true;
    }
    else if(algorithm == "vector_bfgs2") {
        minimumFinder = boost::bind(&GslEngine::minimizeWithGradient,this,
            gsl_multimin_fdfminimizer_vector_bfgs2,_1,_2,_3,_4,lineMinTol);
        useGradient = true;
    }
    else if(algorithm == "steepest_descent") {
        minimumFinder = boost::bind(&GslEngine::minimizeWithGradient,this,
            gsl_multimin_fdfminimizer_steepest_descent,_1,_2,_3,_4,lineMinTol);
        useGradient = true;
    }
    else {
        throw RuntimeError("GslEngine: unknown algorithm '" + algorithm + "'");
    }
    if(useGradient) {
        // Check that we have a gradient calculator to use.
        if(!_gc) {
            throw RuntimeError(
                "GslEngine: selected algorithm needs a gradient calculator.");
        }
        // Bind this function and its gradient calculator to our GSL global function
        _funcWithGradient.n = nPar;
        _funcWithGradient.f = _evaluate;
        _funcWithGradient.df = _evaluateGradient;
        _funcWithGradient.fdf = _evaluateBoth;
        _funcWithGradient.params = 0;
        _grad = Gradient(nPar);        
    }
    else {
        // Bind this function to our GSL global function
        _func.n = nPar;
        _func.f = _evaluate;
        _func.params = 0;
    }
    _params = Parameters(nPar);
    _getEngineStack().push(this);        
}

local::GslEngine::~GslEngine() {
    _getEngineStack().pop();
}

local::FunctionMinimumPtr local::GslEngine::minimizeWithGradient(fdfMethod method,
Parameters const &initial, Parameters const &errors,
double prec, long maxIterations, double lineMinTol) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::minimizeWithGradient");
    // Copy the input initial values to a GSL vector.
    gsl_vector *gsl_initial(gsl_vector_alloc(_nPar));
    for(int i = 0; i < _nPar; ++i) {
        gsl_vector_set(gsl_initial,i,initial[i]);
    }
    // Initialize the minimizer.
    gsl_multimin_fdfminimizer *state(gsl_multimin_fdfminimizer_alloc(method,_nPar));
    // Calculate the RMS of the error vector components.
    double sumsq(0);
    for(int i = 0; i < _nPar; ++i) {
        double value(errors[i]);
        sumsq += value*value;
    }
    double rmsError(std::sqrt(sumsq));
    // Use an initial step size proportional to the RMS error.
    double stepSize = rmsError > 0 ? 0.1*rmsError : 1e-2;
    // Re-scale the requested precision by the initial step size so that it corresponds
    // to a target maximum |gradient| value.
    double maxGradient = prec > 0 ? prec/stepSize : 1e-3;
    gsl_multimin_fdfminimizer_set(state, &_funcWithGradient, gsl_initial,
        stepSize, lineMinTol);
    // Do the minimization...
    long nIterations(0);
    while(0 == maxIterations || nIterations < maxIterations) {
        nIterations++;
        if(gsl_multimin_fdfminimizer_iterate(state)) break;
        // The input precision specifies the 
        if(gsl_multimin_test_gradient(state->gradient, maxGradient) != GSL_CONTINUE) break;
    }
    // Copy the results into our result object.
    Parameters final(_nPar);
    for(int i = 0; i < _nPar; ++i) {
        final[i] = gsl_vector_get(state->x,i);
    }
    FunctionMinimumPtr fmin(new FunctionMinimum(state->f,final));
    // Clean up.
    gsl_vector_free(gsl_initial);
    gsl_multimin_fdfminimizer_free(state);

    return fmin;    
}

local::FunctionMinimumPtr local::GslEngine::minimize(fMethod method,
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
    // Get the top engine on the stack.
    GslEngine *top(_useTopEngine(v));
    // Call the function and return its value.
    top->incrementEvalCount();
    return (*(top->_f))(top->_params);
}

void local::GslEngine::_evaluateGradient(const gsl_vector *v, void *p, gsl_vector *g) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::_evaluateGradient");
    // Get the top engine on the stack.
    GslEngine *top(_useTopEngine(v));
    // Fill the engine's gradient vector.
    top->incrementGradCount();
    (*(top->_gc))(top->_params,top->_grad);
    // Copy the gradient components to the GSL vector provided.
    for(int i = 0; i < top->_nPar; ++i) gsl_vector_set(g,i,top->_grad[i]);
}

void local::GslEngine::_evaluateBoth(const gsl_vector *v, void *p,
double *fval, gsl_vector *g) {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::_evaluateBoth");
    // Get the top engine on the stack.
    GslEngine *top(_useTopEngine(v));
    // Call the function and save its value.
    top->incrementEvalCount();
    *fval = (*(top->_f))(top->_params);
    // Fill the engine's gradient vector.
    top->incrementGradCount();
    (*(top->_gc))(top->_params,top->_grad);
    // Copy the gradient components to the GSL vector provided.
    for(int i = 0; i < top->_nPar; ++i) gsl_vector_set(g,i,top->_grad[i]);
}

local::GslEngine* local::GslEngine::_useTopEngine(const gsl_vector *v) {
    // Get the top function on the stack.
    GslEngine *top(_getEngineStack().top());
    // Copy the input GSL vector to the top engine's _params.
    for(int i = 0; i < top->_nPar; ++i) top->_params[i] = gsl_vector_get(v,i);
    return top;
}

std::stack<local::GslEngine*> &local::GslEngine::_getEngineStack() {
    static std::stack<GslEngine*> *stack = new std::stack<GslEngine*>();
    return *stack;
}

bool local::GslEngine::registerGslEngineMethods() {
    // Declare our error-handling context.
    GslErrorHandler eh("GslEngine::registerEngineMethods");
    // Create a function object that constructs a GslEngine with parameters
    // (FunctionPtr f, GradientCalculatorPtr gc, int npar, std::string const &methodName).
    AbsEngine::EngineFactory factory =
        boost::bind(boost::factory<GslEngine*>(),_1,_2,_3,_4);
    // Register our factory method.
    AbsEngine::getEngineRegistry()["gsl"] = factory;
    // Return a dummy value so that we can be called at program startup.
    return true;
}

bool local::GslEngine::_registered = local::GslEngine::registerGslEngineMethods();
