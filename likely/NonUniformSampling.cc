// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/NonUniformSampling.h"
#include "likely/BinningError.h"

#include <cmath>

namespace local = likely;

local::NonUniformSampling::NonUniformSampling(std::vector<double> const &samplePoints)
: _samplePoints(samplePoints)
{
    // Check that we have some samples.
    int nSamples = _samplePoints.size();
    if(nSamples < 1) {
        throw BinningError("NonUniformSampling: need at least 1 sample point.");
    }
    // Check that sample points are increasing.
    for(int index = 1; index < nSamples; ++index) {
        if(_samplePoints[index-1] > _samplePoints[index]) {
            throw BinningError("NonUniformSampling: sample points are not in increasing order.");
        }
    }
}

local::NonUniformSampling::~NonUniformSampling() { }

int local::NonUniformSampling::getBinIndex(double value) const {
    if(value < _samplePoints[0]) {
        throw BinningError("getBinIndex: value is below binning interval.");
    }
    for(int sample = 0; sample < _samplePoints.size(); ++sample) {
        if(std::fabs(getBinCenter(sample) - value) == 0) return sample;
        if(value > getBinCenter(sample)) {
            throw BinningError("getBinIndex: value is not one of our samples.");
        }
    }
    throw BinningError("getBinIndex: value is above binning interval.");
}

int local::NonUniformSampling::getNBins() const {
    return _samplePoints.size();
}

double local::NonUniformSampling::getBinLowEdge(int index) const {
    return getBinCenter(index);
}

double local::NonUniformSampling::getBinHighEdge(int index) const {
    return getBinCenter(index);
}

double local::NonUniformSampling::getBinWidth(int index) const {
    isValidBinIndex(index,"getBinWidth: invalid bin index %d.");
    return 0;
}

double local::NonUniformSampling::getBinCenter(int index) const {
    isValidBinIndex(index,"getBinCenter: invalid bin index %d.");
    return _samplePoints[index];
}
