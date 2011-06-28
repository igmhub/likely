// Created 25-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_WEIGHTED_ACCUMULATOR
#define LIKELY_WEIGHTED_ACCUMULATOR

#include "likely/AbsAccumulator.h"

#include "boost/smart_ptr.hpp"

namespace likely {
	class WeightedAccumulator : public AbsAccumulator {
	public:
		WeightedAccumulator();
		virtual ~WeightedAccumulator();
        // Accumulates one weighted sample or throws a RuntimeError if weight <= 0.
        void accumulate(double value, double weight);
        // Returns the number of weighted samples accumulated.
        virtual int count() const;
        // Returns the weighted mean of the samples accumulated so far or zero if
        // no samples have been accumulated yet.        
        virtual double mean() const;
        // Returns the the weighted variance of the samples accumulated so far or zero.
        virtual double variance() const;
        // Returns the sum of weights accumulated so far.
        virtual double sumOfWeights() const;
	private:
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
	}; // WeightedAccumulator
} // likely

#endif // LIKELY_WEIGHTED_ACCUMULATOR
