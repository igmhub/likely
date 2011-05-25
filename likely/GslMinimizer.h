// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_GSL_MINIMIZER
#define LIKELY_GSL_MINIMIZER

#include "likely/AbsMinimizer.h"

#include "boost/smart_ptr.hpp"

namespace likely {
    class GslEngine;
	class GslMinimizer : public AbsMinimizer {
	public:
		GslMinimizer(Function f, int nPar);
		virtual ~GslMinimizer();
        virtual Parameters minimize(Parameters const& initial, Parameters const &errors);
	private:
        boost::scoped_ptr<GslEngine> _engine;
	}; // GslMinimizer
} // likely

#endif // LIKELY_GSL_MINIMIZER
