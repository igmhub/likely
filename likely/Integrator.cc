// Created 07-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Integrator.h"

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
#endif
    }; // Integrator::Implementation
} // likely::

local::Integrator::Integrator(IntegrandPtr integrand)
: _integrand(integrand), _pimpl(new Implementation())
{
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Integrator::Integrator");
    _pimpl->workspace = gsl_integration_workspace_alloc(_pimpl->workspaceSize = 1024);
#else
    throw RuntimeError("Integrator: GSL required for all integration methods.");
#endif
}

local::Integrator::~Integrator() {
#ifdef HAVE_LIBGSL
    gsl_integration_workspace_free(_pimpl->workspace);
#endif
}
