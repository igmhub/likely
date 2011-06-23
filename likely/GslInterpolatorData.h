// Created 23-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_GSL_INTERPOLATOR_DATA
#define LIKELY_GSL_INTERPOLATOR_DATA

#include "config.h" // propagates HAVE_LIBGSL from configure

#ifdef HAVE_LIBGSL
#include "gsl/gsl_interp.h"
#endif

namespace likely {
	struct GslInterpolatorData {
#ifdef HAVE_LIBGSL
        const gsl_interp_type *engine;
        gsl_interp_accel *accelerator;
        gsl_interp *interpolator;
#endif
	}; // GslInterpolatorData
} // likely

#endif // LIKELY_GSL_INTERPOLATOR_DATA
