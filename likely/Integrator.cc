// Created 07-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Integrator.h"
#include "likely/RuntimeError.h"

#include "config.h" // propagates HAVE_LIBGSL from configure
#ifdef HAVE_LIBGSL
#include "likely/GslErrorHandler.h"
#include "gsl/gsl_integration.h"
#endif

namespace local = likely;

// Declares our implementation data container.
namespace likely {
    struct Integrator::Implementation {
#ifdef HAVE_LIBGSL
        size_t workspaceSize;
        gsl_integration_workspace *workspace;
        gsl_function function;
#endif
    }; // Integrator::Implementation
} // likely::

local::Integrator::Integrator(IntegrandPtr integrand, double epsAbs, double epsRel)
: _integrand(integrand), _epsAbs(epsAbs), _epsRel(epsRel), _absError(0),
_pimpl(new Implementation())
{
    if(epsRel < 0) {
        throw RuntimeError("Integrator: bad epsRel < 0.");
    }
    if(epsAbs < 0) {
        throw RuntimeError("Integrator: bad epsAbs < 0.");
    }
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Integrator::Integrator");
    // Allocate an integration workspace.
    _pimpl->workspace = gsl_integration_workspace_alloc(_pimpl->workspaceSize = 1024);
    // Link the function wrapper to our static evaluator.
    _pimpl->function.function = &_evaluate;
    _pimpl->function.params = 0;
#else
    throw RuntimeError("Integrator: GSL required for all integration methods.");
#endif
}

local::Integrator::~Integrator() {
#ifdef HAVE_LIBGSL
    gsl_integration_workspace_free(_pimpl->workspace);
#endif
}

double local::Integrator::integrateSmooth(double a, double b) {
    double result(0);
    _getStack().push(this);
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Integrator::integrateSmooth");
    int status = gsl_integration_qag(&_pimpl->function,a,b,_epsAbs,_epsRel,
        _pimpl->workspaceSize,GSL_INTEG_GAUSS61,_pimpl->workspace,&result,&_absError);
#endif
    _getStack().pop();
    return result;
}

double local::Integrator::integrateSingular(double a, double b) {
    double result(0);
    _getStack().push(this);
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Integrator::integrateSingular");
    int status = gsl_integration_qags(&_pimpl->function,a,b,_epsAbs,_epsRel,
        _pimpl->workspaceSize,_pimpl->workspace,&result,&_absError);
#endif
    _getStack().pop();
    return result;
}

double local::Integrator::integrateUp(double a) {
        double result(0);
        _getStack().push(this);
    #ifdef HAVE_LIBGSL
        // Declare our error-handling context.
        GslErrorHandler eh("Integrator::integrateUp");
        int status = gsl_integration_qagiu(&_pimpl->function,a,_epsAbs,_epsRel,
            _pimpl->workspaceSize,_pimpl->workspace,&result,&_absError);
    #endif
        _getStack().pop();
        return result;    
}

double local::Integrator::integrateDown(double b) {
        double result(0);
        _getStack().push(this);
    #ifdef HAVE_LIBGSL
        // Declare our error-handling context.
        GslErrorHandler eh("Integrator::integrateDown");
        int status = gsl_integration_qagil(&_pimpl->function,b,_epsAbs,_epsRel,
            _pimpl->workspaceSize,_pimpl->workspace,&result,&_absError);
    #endif
        _getStack().pop();
        return result;    
}

double local::Integrator::integrateAll() {
        double result(0);
        _getStack().push(this);
    #ifdef HAVE_LIBGSL
        // Declare our error-handling context.
        GslErrorHandler eh("Integrator::integrateDown");
        int status = gsl_integration_qagi(&_pimpl->function,_epsAbs,_epsRel,
            _pimpl->workspaceSize,_pimpl->workspace,&result,&_absError);
    #endif
        _getStack().pop();
        return result;    
}

double local::Integrator::_evaluate(double x, void *params) {
    const Integrator *top(_getStack().top());
    return (*(top->_integrand))(x);
}

std::stack<const local::Integrator*> &local::Integrator::_getStack() {
    static std::stack<const Integrator*> *stack = new std::stack<const Integrator*>();
    return *stack;
}
