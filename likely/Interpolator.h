// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_INTERPOLATOR
#define LIKELY_INTERPOLATOR

#include "config.h" // propagates HAVE_LIBGSL from configure

#ifdef HAVE_LIBGSL
#include "gsl/gsl_interp.h"
#endif

#include <vector>
#include <string>

namespace likely {
    // Implements interpolation algorithms.
	class Interpolator {
	public:
        typedef std::vector<double> CoordinateValues;
        Interpolator(CoordinateValues const &x, CoordinateValues const &y,
            std::string const &algorithm);
        virtual ~Interpolator();
        // Returns the interpolated y value for the specified x value. Returns the
        // appropriate endpoint y value if x is outside the interpolation domain.
        double operator()(double x) const;
	private:
        int _nValues;
        CoordinateValues _x, _y;
#ifdef HAVE_LIBGSL
        const gsl_interp_type *_engine;
        gsl_interp_accel *_accelerator;
        gsl_interp *_interpolator;
#endif
	}; // Interpolator
} // likely

#endif // LIKELY_INTERPOLATOR
