// Created 14-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/UniformBinning.h"
#include "likely/BinningError.h"

namespace local = likely;

local::UniformBinning::UniformBinning(double minValue, double maxValue, int nBins)
: _minValue(minValue), _maxValue(maxValue), _nBins(nBins)
{
    if(maxValue <= minValue) {
        throw BinningError("UniformBinning: expected minValue < maxValue.");
    }
    if(nBins <= 0) {
        throw BinningError("UniformBinning: expected nBins > 0.");
    }
    _binWidth = (maxValue - minValue)/nBins;
}

local::UniformBinning::~UniformBinning() { }

int local::UniformBinning::getNBins() const {
    return _nBins;
}

double local::UniformBinning::getBinLowEdge(int index) const {
    isValidBinIndex(index,"getBinLowEdge: invalid bin index %d.");
    return _minValue + index*_binWidth;
}

double local::UniformBinning::getBinHighEdge(int index) const {
    isValidBinIndex(index,"getBinHighEdge: invalid bin index %d.");
    return _minValue + (index+1)*_binWidth;
}

double local::UniformBinning::getBinWidth(int index) const {
    isValidBinIndex(index,"getBinWidth: invalid bin index %d.");
    return _binWidth;
}

double local::UniformBinning::getBinCenter(int index) const {
    isValidBinIndex(index,"getBinCenter: invalid bin index %d.");
    return _minValue + (index+0.5)*_binWidth;
}
