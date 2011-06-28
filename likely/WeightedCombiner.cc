// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/WeightedCombiner.h"

namespace local = likely;

local::WeightedCombiner::WeightedCombiner() { }

local::WeightedCombiner::~WeightedCombiner() { }

void local::WeightedCombiner::combine(AbsAccumulator const &other) {
    _count += other.count();
    double mu(other.mean()),wsum(other.sumOfWeights());
    _combinedMean.accumulate(mu, wsum);
    _combinedSecondMoment.accumulate(other.variance() + mu*mu, wsum);
}
