// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/WeightedCombiner.h"

namespace local = likely;

local::WeightedCombiner::WeightedCombiner() : _count(0) { }

local::WeightedCombiner::~WeightedCombiner() { }

void local::WeightedCombiner::combine(int count, double sumOfWeights,
double mean, double secondMoment) {
    if(count <= 0 || sumOfWeights <= 0) return;
    _count += count;
    _combinedMean.accumulate(mean,sumOfWeights);
    _combinedSecondMoment.accumulate(secondMoment,sumOfWeights);
}

void local::WeightedCombiner::combine(AbsAccumulator const &other) {
    double mu(other.mean()),wsum(other.sumOfWeights());
    combine(other.count(),wsum,mu,other.variance()+mu*mu);
}
