// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/NonUniformBinning.h"
#include "likely/BinningError.h"

#include <ostream>

namespace local = likely;

local::NonUniformBinning::NonUniformBinning(std::vector<double> const &binEdges)
: _binEdges(binEdges)
{
    // Check that we have some bin edges.
    int nBins = _binEdges.size();
    if(nBins < 2) {
        throw BinningError("NonUniformBinning: need at least 2 bin edges.");
    }
    // Check that bin edges are increasing.
    for(int index = 1; index < nBins; ++index) {
        if(_binEdges[index-1] > _binEdges[index]) {
            throw BinningError("NonUniformBinning: bin edges are not in increasing order.");
        }
    }
}

local::NonUniformBinning::~NonUniformBinning() { }

int local::NonUniformBinning::getBinIndex(double value) const {
    // should use bisection for this, and cache the last bin found...
    if(value < _binEdges[0]) {
        throw BinningError("getBinIndex: value is below binning interval.");
    }
    for(int bin = 1; bin < _binEdges.size(); ++bin) {
        if(value < _binEdges[bin]) return bin-1;
    }
    throw BinningError("getBinIndex: value is above binning interval.");
}

int local::NonUniformBinning::getNBins() const {
    return _binEdges.size()-1;
}

double local::NonUniformBinning::getBinLowEdge(int index) const {
    isValidBinIndex(index,"getBinLowEdge: invalid bin index %d.");
    return _binEdges[index];
}

double local::NonUniformBinning::getBinHighEdge(int index) const {
    isValidBinIndex(index,"getBinHighEdge: invalid bin index %d.");
    return _binEdges[index+1];
}

double local::NonUniformBinning::getBinWidth(int index) const {
    isValidBinIndex(index,"getBinWidth: invalid bin index %d.");
    return _binEdges[index+1] - _binEdges[index];
}

double local::NonUniformBinning::getBinCenter(int index) const {
    isValidBinIndex(index,"getBinCenter: invalid bin index %d.");
    return 0.5*(_binEdges[index] + _binEdges[index+1]);
}

void local::NonUniformBinning::printToStream(std::ostream &os) const {
    os << '[';
    for(int bin = 0; bin < _binEdges.size(); ++bin) {
        if(bin) os << ',';
        os << _binEdges[bin];
    }
    os << ']';
}
