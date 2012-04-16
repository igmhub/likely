// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/UniformSampling.h"
#include "likely/BinningError.h"

namespace local = likely;

local::UniformSampling::UniformSampling(double minValue, double maxValue, int nSamples)
: _minValue(minValue), _maxValue(maxValue), _nSamples(nSamples)
{
    if(maxValue < minValue) {
        throw BinningError("UniformSampling: expected minValue <= maxValue.");
    }
    if(nSamples < 1) {
        throw BinningError("UniformSampling: expected nSamples > 1.");
    }
    if(nSamples == 1 && maxValue != minValue) {
        throw BinningError("UniformSampling: must have minValue==maxValue when nSamples==1.");
    }
    _sampleSpacing = (1==nSamples) ? 0 : (maxValue - minValue)/(nSamples-1.0);
}

local::UniformSampling::~UniformSampling() { }

int local::UniformSampling::getNBins() const {
    return _nSamples;
}

double local::UniformSampling::getBinLowEdge(int index) const {
    return getBinCenter(index);
}

double local::UniformSampling::getBinHighEdge(int index) const {
    return getBinCenter(index);
}

double local::UniformSampling::getBinWidth(int index) const {
    isValidBinIndex(index,"getBinWidth: invalid bin index %d.");
    return 0;
}

double local::UniformSampling::getBinCenter(int index) const {
    isValidBinIndex(index,"getBinCenter: invalid bin index %d.");
    return _minValue + index*_sampleSpacing;
}
