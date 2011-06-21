// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_INTERPOLATOR
#define LIKELY_INTERPOLATOR

#include "likely/types.h"

#include <string>

namespace likely {
    // Implements interpolation algorithms.
	class Interpolator {
	public:
        Interpolator(Parameters const &x, Parameters const &y,
            std::string const &methodName);
		virtual ~Interpolator();
	private:
	}; // Interpolator
} // likely

#endif // LIKELY_INTERPOLATOR
