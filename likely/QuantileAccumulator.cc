// Created 13-Sep-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/QuantileAccumulator.h"
#include "likely/RuntimeError.h"

//==================================================================================
// Need to add #include <cmath> to the top of
// /usr/local/include/boost/accumulators/statistics/weighted_p_square_quantile.hpp
// for this to compile on OS-X with boost 1.45.0
// See https://svn.boost.org/trac/boost/ticket/5894
// 
// The include below is a workaround kludge.
//==================================================================================
#include <cmath>

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/count.hpp"
#include "boost/accumulators/statistics/weighted_p_square_quantile.hpp"

namespace likely {
    struct QuantileAccumulator::Implementation {
        Implementation(double quantileProbability)
        : data(boost::accumulators::quantile_probability = quantileProbability) { }
        // Uses the boost statistical accumulators library.
        boost::accumulators::accumulator_set<double,
            boost::accumulators::stats<
                boost::accumulators::features<boost::accumulators::tag::count>,
                boost::accumulators::tag::weighted_p_square_quantile
            >, double
        > data;
    };
} // likely::

namespace local = likely;

local::QuantileAccumulator::QuantileAccumulator(double quantileProbability)
: _pimpl(new Implementation(quantileProbability))
{
}

local::QuantileAccumulator::~QuantileAccumulator() { }

void local::QuantileAccumulator::accumulate(double value, double weight) {
    if(weight <= 0) {
        throw RuntimeError("WeightedAccumulator::accumulate found weight <= 0.");
    }
    _pimpl->data(value, boost::accumulators::weight = weight);    
}

int local::QuantileAccumulator::count() const {
    return boost::accumulators::count(_pimpl->data);
}

double local::QuantileAccumulator::getQuantile() const {
    return boost::accumulators::weighted_p_square_quantile(_pimpl->data);
}