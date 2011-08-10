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
        gsl_integration_cquad_workspace *cquad_workspace;
        gsl_integration_qawo_table *qawo_table;
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
    // Integration workspaces will be allocated on demand.
    _pimpl->workspaceSize = 1024;
    _pimpl->workspace = 0;
    _pimpl->cquad_workspace = 0;
    _pimpl->qawo_table = 0;
    // Link the function wrapper to our static evaluator.
    _pimpl->function.function = &_evaluate;
    _pimpl->function.params = 0;
#else
    throw RuntimeError("Integrator: GSL required for all integration methods.");
#endif
}

local::Integrator::~Integrator() {
#ifdef HAVE_LIBGSL
    if(0 != _pimpl->workspace) gsl_integration_workspace_free(_pimpl->workspace);
    if(0 != _pimpl->cquad_workspace)
        gsl_integration_cquad_workspace_free(_pimpl->cquad_workspace);
    if(0 != _pimpl->qawo_table) gsl_integration_qawo_table_free(_pimpl->qawo_table);
#endif
}

double local::Integrator::integrateSmooth(double a, double b) {
    double result(0);
    _getStack().push(this);
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Integrator::integrateSmooth");
    if(0 == _pimpl->workspace) _pimpl->workspace =
        gsl_integration_workspace_alloc(_pimpl->workspaceSize);
    int status = gsl_integration_qag(&_pimpl->function,a,b,_epsAbs,_epsRel,
        _pimpl->workspaceSize,GSL_INTEG_GAUSS61,_pimpl->workspace,&result,&_absError);
#endif
    _getStack().pop();
    return result;
}

double local::Integrator::integrateRobust(double a, double b) {
    double result(0);
    _getStack().push(this);
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Integrator::integrateRobust");
    if(0 == _pimpl->cquad_workspace) _pimpl->cquad_workspace =
        gsl_integration_cquad_workspace_alloc(_pimpl->workspaceSize);
    size_t nEvals;
    int status = gsl_integration_cquad(&_pimpl->function,a,b,_epsAbs,_epsRel,
        _pimpl->cquad_workspace,&result,&_absError,&nEvals);
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
    if(0 == _pimpl->workspace) _pimpl->workspace =
        gsl_integration_workspace_alloc(_pimpl->workspaceSize);
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
        if(0 == _pimpl->workspace) _pimpl->workspace =
            gsl_integration_workspace_alloc(_pimpl->workspaceSize);
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
        if(0 == _pimpl->workspace) _pimpl->workspace =
            gsl_integration_workspace_alloc(_pimpl->workspaceSize);
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
        if(0 == _pimpl->workspace) _pimpl->workspace =
            gsl_integration_workspace_alloc(_pimpl->workspaceSize);
        int status = gsl_integration_qagi(&_pimpl->function,_epsAbs,_epsRel,
            _pimpl->workspaceSize,_pimpl->workspace,&result,&_absError);
    #endif
        _getStack().pop();
        return result;    
}

double local::Integrator::integrateOsc(double a, double b, double omega, bool useSin) {
        double result(0);
        _getStack().push(this);
    #ifdef HAVE_LIBGSL
        // Declare our error-handling context.
        GslErrorHandler eh("Integrator::integrateOsc");
        if(0 == _pimpl->workspace) _pimpl->workspace =
            gsl_integration_workspace_alloc(_pimpl->workspaceSize);
        gsl_integration_qawo_enum which(useSin ? GSL_INTEG_SINE : GSL_INTEG_COSINE);
        if(0 == _pimpl->qawo_table) {
            _pimpl->qawo_table =
                gsl_integration_qawo_table_alloc(omega,b-a,which,_pimpl->workspaceSize);
        }
        else {
            gsl_integration_qawo_table_set(_pimpl->qawo_table,omega,b-a,which);
        }
        int status = gsl_integration_qawo(&_pimpl->function,a,_epsAbs,_epsRel,
            _pimpl->workspaceSize,_pimpl->workspace,_pimpl->qawo_table,&result,&_absError);
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
