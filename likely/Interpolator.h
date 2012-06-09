// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_INTERPOLATOR
#define LIKELY_INTERPOLATOR

#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <vector>
#include <string>
#include <iosfwd>

namespace likely {
    // Implements interpolation algorithms.
	class Interpolator {
	public:
        typedef std::vector<double> CoordinateValues;
        // Creates a new interpolator from the specified x,y vectors. Supported algorithms are
        // described at http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html
        // (but with the gsl_interp_ prefix ommitted from the GSL function name).
        Interpolator(CoordinateValues const &x, CoordinateValues const &y,
            std::string const &algorithm);
        virtual ~Interpolator();
        // Returns the interpolated y value for the specified x value. Returns the
        // appropriate endpoint y value if x is outside the interpolation domain.
        double operator()(double x) const;
	private:
        int _nValues;
        CoordinateValues _x, _y;
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
	}; // Interpolator
	
    // Returns a smart pointer to an interpolator based on control points read
    // from the specified file name.
	InterpolatorPtr createInterpolator(std::string const &filename,
        std::string const &algorithm);

    // Fills the vectors provided from the columns of the specified input stream.
    // Returns the number of rows successfully read or throws a RuntimeError.
    // Any input beyond the required column values is silently ignore if ignoreExtra
    // is set or, otherwise, generates a RuntimeError.
    int readVectors(std::istream &input, std::vector<std::vector<double> > &vectors,
        bool ignoreExtra = true);
	
} // likely

#endif // LIKELY_INTERPOLATOR
