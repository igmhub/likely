// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/WeightedAccumulator.h"
#include "likely/RuntimeError.h"

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/count.hpp"
#include "boost/accumulators/statistics/weighted_mean.hpp"
#include "boost/accumulators/statistics/weighted_variance.hpp"
#include "boost/accumulators/statistics/sum.hpp"

#include <cmath>

namespace likely {
    struct WeightedAccumulator::Implementation {
        // Uses the boost statistical accumulators library.
        boost::accumulators::accumulator_set<double,
            boost::accumulators::stats<
                boost::accumulators::features<boost::accumulators::tag::count>,
                boost::accumulators::tag::weighted_mean,
                boost::accumulators::tag::weighted_variance,
                boost::accumulators::tag::sum_of_weights
            >, double
        > data;
    };
} // likely::

namespace local = likely;

local::WeightedAccumulator::WeightedAccumulator()
: _pimpl(new Implementation())
{
}

local::WeightedAccumulator::~WeightedAccumulator() { }

void local::WeightedAccumulator::accumulate(double value, double weight) {
    if(weight <= 0) {
        throw RuntimeError("WeightedAccumulator::accumulate found weight <= 0.");
    }
    _pimpl->data(value, boost::accumulators::weight = weight);
}

int local::WeightedAccumulator::count() const {
    return boost::accumulators::count(_pimpl->data);
}

double local::WeightedAccumulator::mean() const {
    return boost::accumulators::count(_pimpl->data) > 0 ?
        boost::accumulators::weighted_mean(_pimpl->data) : 0;
}

double local::WeightedAccumulator::variance() const {
    return boost::accumulators::count(_pimpl->data) > 0 ?
        boost::accumulators::weighted_variance(_pimpl->data) : 0;
}

double local::WeightedAccumulator::error() const {
    return std::sqrt(variance());
}

double local::WeightedAccumulator::sumOfWeights() const {
    return boost::accumulators::sum_of_weights(_pimpl->data);
}

double local::WeightedAccumulator::errorOnMean() const {
    double wsum(boost::accumulators::sum_of_weights(_pimpl->data));
    return wsum > 0 ? 1/std::sqrt(wsum) : 0;
}
