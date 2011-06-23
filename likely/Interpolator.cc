// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Interpolator.h"
#include "likely/RuntimeError.h"
#include "likely/GslInterpolatorData.h"

#ifdef HAVE_LIBGSL
#include "likely/GslErrorHandler.h"
#endif

#include "boost/lexical_cast.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

namespace local = likely;

local::Interpolator::Interpolator(CoordinateValues const &x, CoordinateValues const &y,
std::string const &algorithm)
: _x(x), _y(y), _nValues(x.size()), _gslData(new GslInterpolatorData())
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
    _gslData->accelerator = gsl_interp_accel_alloc();
    // Lookup the GSL engine for the requested algorithm.
    if(algorithm == "linear") _gslData->engine = gsl_interp_linear;
    else if(algorithm == "polynomial") _gslData->engine = gsl_interp_polynomial;
    else if(algorithm == "cspline") _gslData->engine = gsl_interp_cspline;
    else if(algorithm == "cspline_periodic") _gslData->engine = gsl_interp_cspline_periodic;
    else if(algorithm == "cspline_akima") _gslData->engine = gsl_interp_akima;
    else if(algorithm == "cspline_akima_periodic") _gslData->engine = gsl_interp_akima_periodic;
    else {
        throw RuntimeError("Interpolator: unknown algorithm '" + algorithm + "'.");
    }
    // Check that we have enough coordinate values for the requested scheme.
    if(_nValues < _gslData->engine->min_size) {
        throw RuntimeError("Interpolator: need more values for the requested algorithm.");
    }
    // Create the engine's data structure.
    _gslData->interpolator = gsl_interp_alloc(_gslData->engine, _nValues);
    // Initialize the interpolation using the coordinate values provided.
    gsl_interp_init(_gslData->interpolator, &_x[0], &_y[0], _nValues);
#else
    throw RuntimeError("Interpolator: GSL required for all interpolation methods.");
#endif
}

local::Interpolator::~Interpolator() {
#ifdef HAVE_LIBGSL
    gsl_interp_free(_gslData->interpolator);
    gsl_interp_accel_free(_gslData->accelerator);
#endif
}

double local::Interpolator::operator()(double x) const {
    // Declare our error-handling context.
    GslErrorHandler eh("Interpolator::operator()");
    // Check for an out-of-range x value.
    if(x <= _x.front()) return _y.front();
    if(x >= _x.back()) return _y.back();
    return gsl_interp_eval(_gslData->interpolator,
        &_x[0], &_y[0], x, _gslData->accelerator);
}

int local::readVectors(std::istream &input, std::vector<std::vector<double> > &vectors,
bool ignoreExtra) {
    // Loop over input lines.
    int lineNumber(0);
    std::string line;
    int nColumns(vectors.size());
    while(1) {
        // Try to read the next line.
        std::getline(input, line);
        if(!input.good() || input.eof()) break;
        lineNumber++;
        std::vector<double> values;
        std::istringstream iss(line);
        std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(),
            std::back_inserter(values));
        if(values.size() < nColumns || (!ignoreExtra && values.size() > nColumns)) {
            throw RuntimeError("Badly formed input on line " +
                boost::lexical_cast<std::string>(lineNumber));
        }
        for(int index = 0; index < nColumns; ++index) {
            vectors[index].push_back(values[index]);
        }
    }
    if(!input.eof()) {
        throw RuntimeError("Error reading input on line " +
            boost::lexical_cast<std::string>(lineNumber));
    }
    return lineNumber;
}
