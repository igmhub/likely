// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/BinnedData.h"
#include "likely/RuntimeError.h"
#include "likely/AbsBinning.h"

#include "boost/foreach.hpp"

namespace local = likely;

local::BinnedData::BinnedData(std::vector<AbsBinningCPtr> axes)
: _axisBinning(axes)
{
    if(0 == axes.size()) {
        throw RuntimeError("BinnedData: no axes provided.");
    }
    _initialize();
}

local::BinnedData::BinnedData(AbsBinningCPtr axis1)
{
    if(!axis1) {
        throw RuntimeError("BinnedData: missing axis1.");
    }
    _axisBinning.push_back(axis1);
    _initialize();
}

local::BinnedData::BinnedData(AbsBinningCPtr axis1, AbsBinningCPtr axis2)
{
    if(!axis1 || !axis2) {
        throw RuntimeError("BinnedData: missing axis data.");
    }
    _axisBinning.push_back(axis1);    
    _axisBinning.push_back(axis2);
    _initialize();
}

local::BinnedData::BinnedData(AbsBinningCPtr axis1, AbsBinningCPtr axis2, AbsBinningCPtr axis3)
{
    if(!axis1 || !axis2 || !axis3) {
        throw RuntimeError("BinnedData: missing axis data.");
    }
    _axisBinning.push_back(axis1);
    _axisBinning.push_back(axis2);
    _axisBinning.push_back(axis3);
    _initialize();
}

void local::BinnedData::_initialize() {
    _ndata = 0;
    _nbins = 1;
    BOOST_FOREACH(AbsBinningCPtr binning, _axisBinning) {
        _nbins *= binning->getNBins();
    }
}

local::BinnedData::~BinnedData() { }
