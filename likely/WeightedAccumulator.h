// Created 25-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_WEIGHTED_ACCUMULATOR
#define LIKELY_WEIGHTED_ACCUMULATOR

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/count.hpp"
#include "boost/accumulators/statistics/weighted_mean.hpp"
#include "boost/accumulators/statistics/weighted_variance.hpp"
#include "boost/accumulators/statistics/sum.hpp"

#include <cmath>

namespace likely {
    
    // Uses the boost statistical accumulators library to calculate the weighted
    // mean and variance of a sample:
    //
    // likely::WeightedAccumulator bin;
    // accumulateWeighted(bin,1.2,0.5);
    // accumulateWeighted(bin,1.5,1.2);
    // std::cout << weightedMean(bin) << "+/-" << weightedError(bin) << std::endl;

    typedef boost::accumulators::accumulator_set<double,
        boost::accumulators::stats<
            boost::accumulators::features<boost::accumulators::tag::count>,
            boost::accumulators::tag::weighted_mean,
            boost::accumulators::tag::weighted_variance,
            boost::accumulators::tag::sum_of_weights
        >, double
    > WeightedAccumulator;
    
    // Accumulates one weighted sample. Samples with weight <= 0 are silently ignored.
    inline void accumulateWeighted(WeightedAccumulator &acc, double value, double weight) {
        if(weight > 0) acc(value, boost::accumulators::weight = weight);
    }

    // Returns the number of weighted samples accumulated.
    inline int weightedCount(WeightedAccumulator const &acc) {
        return boost::accumulators::count(acc);
    }
    
    // Returns the weighted mean of the samples accumulated so far or zero if
    // no samples have been accumulated yet.
    inline double weightedMean(WeightedAccumulator const &acc) {
        return boost::accumulators::count(acc) > 0 ?
            boost::accumulators::weighted_mean(acc) : 0;
    }

    // Returns the the weighted variance of the samples accumulated so far or zero.
    inline double weightedVariance(WeightedAccumulator const &acc) {
        return boost::accumulators::count(acc) > 0 ?
            boost::accumulators::weighted_variance(acc) : 0;
    }
    
    // Returns sqrt(weightedVariance).
    inline double weightedError(WeightedAccumulator const &acc) {
        return std::sqrt(weightedVariance(acc));
    }

    // Returns the sum of weights accumulated so far.
    inline double weightSum(WeightedAccumulator const &acc) {
        return boost::accumulators::sum_of_weights(acc);
    }
    
    // Returns sqrt(1/weightSum) or zero. Note that this is only a valid estimate of
    // the error on the mean if each weight is calculated as 1/sigma^2.
    inline double weightedErrorOnMean(WeightedAccumulator const &acc) {
        double wsum(boost::accumulators::sum_of_weights(acc));
        return wsum > 0 ? 1/std::sqrt(wsum) : 0;
    }

} // likely

#endif // LIKELY_WEIGHTED_ACCUMULATOR
