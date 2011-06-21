// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Interpolator.h"
#include "likely/RuntimeError.h"

#ifdef HAVE_LIBGSL
#include "likely/GslErrorHandler.h"
#endif

namespace local = likely;

local::Interpolator::Interpolator(CoordinateValues const &x, CoordinateValues const &y,
std::string const &algorithm)
: _x(x), _y(y), _nValues(x.size())
{
    // Check that the input vectors have the same length.
    if(x.size() != y.size()) {
        throw RuntimeError("Interpolator: input vectors must have the same length.");
    }
    // Check that x is sorted here, or does gsl_interp_init handle that?
#ifdef HAVE_LIBGSL
    // Declare our error-handling context.
    GslErrorHandler eh("Interpolator::Interpolator");
    // Create our accelerator.
    _accelerator = gsl_interp_accel_alloc();
    // Lookup the GSL engine for the requested algorithm.
    if(algorithm == "linear") _engine = gsl_interp_linear;
    else if(algorithm == "polynomial") _engine = gsl_interp_polynomial;
    else if(algorithm == "cspline") _engine = gsl_interp_cspline;
    else if(algorithm == "cspline_periodic") _engine = gsl_interp_cspline_periodic;
    else if(algorithm == "cspline_akima") _engine = gsl_interp_akima;
    else if(algorithm == "cspline_akima_periodic") _engine = gsl_interp_akima_periodic;
    else {
        throw RuntimeError("Interpolator: unknown algorithm '" + algorithm + "'.");
    }
    // Check that we have enough coordinate values for the requested scheme.
    if(_nValues < _engine->min_size) {
        throw RuntimeError("Interpolator: need more values for the requested algorithm.");
    }
    // Create the engine's data structure.
    _interpolator = gsl_interp_alloc(_engine, _nValues);
    // Initialize the interpolation using the coordinate values provided.
    gsl_interp_init(_interpolator, &_x[0], &_y[0], _nValues);
#else
    throw RuntimeError("Interpolator: GSL required for all interpolation methods.");
#endif
}

local::Interpolator::~Interpolator() {
#ifdef HAVE_LIBGSL
    gsl_interp_free(_interpolator);
    gsl_interp_accel_free(_accelerator);
#endif
}

double local::Interpolator::operator()(double x) const {
    // Declare our error-handling context.
    GslErrorHandler eh("Interpolator::operator()");
    // Check for an out-of-range x value.
    if(x <= _x.front()) return _y.front();
    if(x >= _x.back()) return _y.back();
    return gsl_interp_eval(_interpolator, &_x[0], &_y[0], x, _accelerator);
}