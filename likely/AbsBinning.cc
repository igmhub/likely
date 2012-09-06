// Created 14-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/AbsBinning.h"
#include "likely/BinningError.h"

#include "boost/format.hpp"

#include <iostream>

namespace local = likely;

local::AbsBinning::AbsBinning() { }

local::AbsBinning::~AbsBinning() { }

double local::AbsBinning::getBinCenter(int index) const {
    isValidBinIndex(index,"getBinCenter: invalid bin index %d.");
    return 0.5*(getBinLowEdge(index) + getBinHighEdge(index));
}

double local::AbsBinning::getBinWidth(int index) const {
    isValidBinIndex(index,"getBinWidth: invalid bin index %d.");
    return getBinHighEdge(index) - getBinLowEdge(index);
}

void local::AbsBinning::dump(std::ostream &os) const {
    int nbins(getNBins());
    os << nbins << ' ' << getBinLowEdge(0);
    for(int bin = 0; bin < nbins; ++bin) os << ' ' << getBinHighEdge(bin);
    os << std::endl;
}

bool local::AbsBinning::isValidBinIndex(int index, std::string const &errorFormat) const {
    if(index >= 0 && index < getNBins()) return true;
    if(0 < errorFormat.length()) {
        std::string errorMessage = boost::str(boost::format(errorFormat) % index);
        throw BinningError(errorMessage);
    }
    else {
        return false;
    }
}
