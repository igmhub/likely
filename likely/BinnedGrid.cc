// Created 10-May-2013 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/BinnedGrid.h"
#include "likely/AbsBinning.h"
#include "likely/RuntimeError.h"

#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"

namespace local = likely;

local::BinnedGrid::BinnedGrid(std::vector<AbsBinningCPtr> axes)
: _axisBinning(axes)
{
    if(0 == axes.size()) {
        throw RuntimeError("BinnedGrid: no axes provided.");
    }
    _nbins = 1;
    BOOST_FOREACH(AbsBinningCPtr binning, _axisBinning) {
        _nbins *= binning->getNBins();
    }
}

local::BinnedGrid::~BinnedGrid() { }

int local::BinnedGrid::getIndex(std::vector<int> const &binIndices) const {
    int nAxes(getNAxes());
    if(binIndices.size() != nAxes) {
        throw RuntimeError("BinnedGrid::getIndex: invalid input vector size.");
    }
    int index(0);
    for(int axis = 0; axis < nAxes; ++axis) {
        int binIndex(binIndices[axis]), nBins(_axisBinning[axis]->getNBins());
        if(binIndex < 0 || binIndex >= nBins) {
            throw RuntimeError("BinnedGrid::getIndex: invalid bin index.");
        }
        index = binIndex + index*nBins;
    }
    return index;
}

int local::BinnedGrid::getIndex(std::vector<double> const &values) const {
    int index(0), nAxes(getNAxes());
    if(values.size() != nAxes) {
        throw RuntimeError("BinnedGrid::getIndex: invalid input vector size.");
    }
    std::vector<int> binIndices;
    for(int axis = 0; axis < nAxes; ++axis) {
        binIndices.push_back(_axisBinning[axis]->getBinIndex(values[axis]));
    }
    return getIndex(binIndices);
}

void local::BinnedGrid::getBinIndices(int index, std::vector<int> &binIndices) const {
    _checkIndex(index);
    int nAxes(getNAxes());
    binIndices.resize(nAxes,0);
    int partial(index);
    for(int axis = nAxes-1; axis >= 0; --axis) {
        AbsBinningCPtr binning = _axisBinning[axis];
        int nBins(binning->getNBins()), binIndex(partial % nBins);
        binIndices[axis] = binIndex;
        partial = (partial - binIndex)/nBins;
    }
}

void local::BinnedGrid::getBinCenters(int index, std::vector<double> &binCenters) const {
    _checkIndex(index);
    binCenters.resize(0);
    binCenters.reserve(getNAxes());
    std::vector<int> binIndices;
    getBinIndices(index,binIndices);
    int nAxes(getNAxes());
    for(int axis = 0; axis < nAxes; ++axis) {
        AbsBinningCPtr binning = _axisBinning[axis];
        binCenters.push_back(binning->getBinCenter(binIndices[axis]));
    }
}

void local::BinnedGrid::getBinWidths(int index, std::vector<double> &binWidths) const {
    _checkIndex(index);
    binWidths.resize(0);
    binWidths.reserve(getNAxes());
    std::vector<int> binIndices;
    getBinIndices(index,binIndices);
    int nAxes(getNAxes());
    for(int axis = 0; axis < nAxes; ++axis) {
        AbsBinningCPtr binning = _axisBinning[axis];
        binWidths.push_back(binning->getBinWidth(binIndices[axis]));
    }
}

void local::BinnedGrid::_checkIndex(int index) const {
    if(index < 0 || index >= _nbins) {
        throw RuntimeError("BinnedGrid: invalid index " +
            boost::lexical_cast<std::string>(index));
    }
}
