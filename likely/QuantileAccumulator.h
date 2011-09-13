// Created 13-Sep-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_QUANTILE_ACCUMULATOR
#define LIKELY_QUANTILE_ACCUMULATOR

#include "boost/smart_ptr.hpp"

namespace likely {
	class QuantileAccumulator {
	public:
	    // Creates a new quantile accumulator for the specified probability level which
	    // defaults to 50% (median).
		QuantileAccumulator(double quantileProbability = 0.5);
		virtual ~QuantileAccumulator();
        // Accumulates one (possibly weighted) sample value.
        void accumulate(double value, double weight = 1);
        // Returns the number of weighted samples accumulated.
        int count() const;
        // Returns the quantile value based on the samples accumulated so far.
        double getQuantile() const;
	private:
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
	}; // QuantileAccumulator
} // likely

#endif // LIKELY_QUANTILE_ACCUMULATOR
