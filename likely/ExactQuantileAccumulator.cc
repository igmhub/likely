// Created 27-Aug-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "likely/ExactQuantileAccumulator.h"
#include "likely/RuntimeError.h"

namespace local = likely;

local::ExactQuantileAccumulator::ExactQuantileAccumulator()
{
    _weightedCount = 0;
}

local::ExactQuantileAccumulator::~ExactQuantileAccumulator() { }

void local::ExactQuantileAccumulator::accumulate(double value, double weight) {
    if(weight <= 0) {
        throw RuntimeError("ExactQuantileAccumulator::accumulate found weight <= 0.");
    }
    _valueWeightPairs.insert(ValueWeightPair(value, weight));
    _weightedCount += weight;
}

int local::ExactQuantileAccumulator::count() const {
    return _valueWeightPairs.size();
}

double local::ExactQuantileAccumulator::getQuantile(double quantileProbability) const {
    if(_weightedCount == 0 ) {
        throw RuntimeError("ExactQuantileAccumulator::getQuantile : _weightedCount must be > 0, no values have been accumulated so far.");
    }
    if(quantileProbability < 0 || quantileProbability > 1) {
        throw RuntimeError("ExactQuantileAccumulator::getQuantile : quantileProbability should be in range [0,1].");
    }
    double _weightedSoFar = 0;
    ValueWeightPairMultiSet::iterator it;
    for(it = _valueWeightPairs.begin(); it != _valueWeightPairs.end(); ++it) {
        _weightedSoFar += it->second;
        if( _weightedSoFar/_weightedCount >= quantileProbability) break;
    }
    return it->first;
}