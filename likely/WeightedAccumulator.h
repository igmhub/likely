// Created 25-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_WEIGHTED_ACCUMULATOR
#define LIKELY_WEIGHTED_ACCUMULATOR

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/weighted_mean.hpp"
#include "boost/accumulators/statistics/weighted_variance.hpp"

#include <cmath>

namespace likely {
    
    // Uses the boost statistical accumulators library to calculate the weighted
    // mean and variance of a sample:
    //
    // likely::WeightedAccumulator bin;
    // accumulate(bin,1.2,0.5);
    // accumulate(bin,1.5,1.2);
    // std::cout << weightedMean(bin) << "+/-" << weightedError(bin) << std::endl;

    typedef boost::accumulators::accumulator_set<double,
        boost::accumulators::stats<
            boost::accumulators::tag::weighted_mean,
            boost::accumulators::tag::weighted_variance
        >, double
    > WeightedAccumulator;
    
    inline void accumulate(WeightedAccumulator &acc, double value, double weight) {
        acc(value, boost::accumulators::weight = weight);
    }
    
    inline double weightedMean(WeightedAccumulator const &acc) {
        return boost::accumulators::weighted_mean(acc);
    }

    inline double weightedVariance(WeightedAccumulator const &acc) {
        return boost::accumulators::weighted_variance(acc);
    }
    
    inline double weightedError(WeightedAccumulator const &acc) {
        return std::sqrt(weightedVariance(acc));
    }

} // likely

#endif // LIKELY_WEIGHTED_ACCUMULATOR
