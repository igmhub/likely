// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_WEIGHTED_COMBINER
#define LIKELY_WEIGHTED_COMBINER

#include "likely/AbsAccumulator.h"
#include "likely/WeightedAccumulator.h"

namespace likely {
	class WeightedCombiner : public AbsAccumulator {
	public:
		WeightedCombiner();
		virtual ~WeightedCombiner();
		// Combines the results of another accumulator with our results.
        void combine(int count, double sumOfWeights, double mean, double variance);
        void combine(AbsAccumulator const &other);
        // Returns the number of weighted samples accumulated.
        virtual int count() const;
        // Returns the weighted mean of the samples accumulated so far or zero if
        // no samples have been accumulated yet.        
        virtual double mean() const;
        // Returns the the weighted variance of the samples accumulated so far or zero.
        virtual double variance() const;
        // Returns the sum of weights accumulated so far.
        virtual double sumOfWeights() const;
        // Returns the estimated variance of the means that have been combined
        // so far.
        double meanVariance() const;
	private:
        int _count;
        WeightedAccumulator _combinedMean, _combinedSecondMoment;
	}; // WeightedCombiner
} // likely

#endif // LIKELY_WEIGHTED_COMBINER
