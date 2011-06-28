// Created 27-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/AbsAccumulator.h"

#include <cmath>

namespace local = likely;

local::AbsAccumulator::AbsAccumulator() { }

local::AbsAccumulator::~AbsAccumulator() { }

double local::AbsAccumulator::error() const {
    return std::sqrt(variance());
}

double local::AbsAccumulator::errorOnMean() const {
    double wsum(sumOfWeights());
    return wsum > 0 ? 1/std::sqrt(wsum) : 0;
}
