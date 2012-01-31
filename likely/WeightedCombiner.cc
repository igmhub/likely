// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/WeightedCombiner.h"

namespace local = likely;

local::WeightedCombiner::WeightedCombiner() : _count(0) { }

local::WeightedCombiner::~WeightedCombiner() { }

void local::WeightedCombiner::combine(int count, double sumOfWeights,
double mean, double variance) {
    if(count <= 0 || sumOfWeights <= 0) return;
    _count += count;
    _combinedMean.accumulate(mean,sumOfWeights);
    _combinedSecondMoment.accumulate(variance+mean*mean,sumOfWeights);
}

void local::WeightedCombiner::combine(AbsAccumulator const &other) {
    combine(other.count(),other.sumOfWeights(),other.mean(),other.variance());
}

int local::WeightedCombiner::count() const {
    return _count;
}

double local::WeightedCombiner::sum() const {
    return _combinedMean.sum();
}

double local::WeightedCombiner::mean() const {
    return _combinedMean.mean();
}

double local::WeightedCombiner::variance() const {
    double mu(mean());
    return _combinedSecondMoment.mean() - mu*mu;
}

double local::WeightedCombiner::sumOfWeights() const {
    return _combinedMean.sumOfWeights();
}

double local::WeightedCombiner::max() const {
	return _combinedMean.max();
}

double local::WeightedCombiner::min() const {
	return _combinedMean.min();
}

double local::WeightedCombiner::meanVariance() const {
    return _combinedMean.variance();
}
