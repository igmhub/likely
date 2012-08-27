// Created 27-Aug-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef LIKELY_EXACT_QUANTILE_ACCUMULATOR
#define LIKELY_EXACT_QUANTILE_ACCUMULATOR

#include <set>
#include <utility>

namespace likely {
	class ExactQuantileAccumulator {

	public:
	    // Creates a new exact quantile accumulator.
		ExactQuantileAccumulator();
		~ExactQuantileAccumulator();
		// Accumulates one (possibly weighted) sample value;
		void accumulate(double value, double weight = 1);
		// Returns the number of weighted samples accumulated.
		int count() const;
		// Returns the quantile value to the specified probability level based on 
		// the samples accumulated so far.
		double getQuantile(double quantileProbability) const;
	private:
		typedef std::pair<double,double> ValueWeightPair;
		typedef std::multiset<ValueWeightPair> ValueWeightPairMultiSet;
		ValueWeightPairMultiSet _valueWeightPairs;
		double _weightedCount;

	}; // ExactQuantileAccumulator
} // likely

#endif // LIKELY_EXACT_QUANTILE_ACCUMULATOR
