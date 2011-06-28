// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_ACCUMULATOR
#define LIKELY_ABS_ACCUMULATOR

namespace likely {
	class AbsAccumulator {
	// Delcares a read-only interface for an abstract sample statistics accumulator.
	public:
		AbsAccumulator();
		virtual ~AbsAccumulator();
        // Returns the number of weighted samples accumulated.
        virtual int count() const = 0;
        // Returns the weighted mean of the samples accumulated so far or zero if
        // no samples have been accumulated yet.        
        virtual double mean() const = 0;
        // Returns the the weighted variance of the samples accumulated so far or zero.
        virtual double variance() const = 0;
        // Returns sqrt(weightedVariance).
        virtual double error() const = 0;
        // Returns the sum of weights accumulated so far.
        virtual double sumOfWeights() const = 0;
        // Returns sqrt(1/weightSum) or zero. Note that this is only a valid estimate of
        // the error on the mean if each weight is calculated as 1/sigma^2.
        virtual double errorOnMean() const = 0;
	private:
	}; // AbsAccumulator
} // likely

#endif // LIKELY_ABS_ACCUMULATOR
