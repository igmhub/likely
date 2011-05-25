// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_MINIMIZER
#define LIKELY_ABS_MINIMIZER

#include "likely/types.h"

namespace likely {
	class AbsMinimizer {
	public:
		AbsMinimizer();
		virtual ~AbsMinimizer();
        virtual Parameters
            minimize(Parameters const& initial, Parameters const &errors) = 0;
	private:
	}; // AbsMinimizer
} // likely

#endif // LIKELY_ABS_MINIMIZER
