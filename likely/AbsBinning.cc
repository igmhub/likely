// Created 14-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/AbsBinning.h"

#include <iostream>

namespace local = likely;

local::AbsBinning::AbsBinning() { }

local::AbsBinning::~AbsBinning() { }

double local::AbsBinning::getBinCenter(int index) const {
    return 0.5*(getBinLowEdge(index) + getBinHighEdge(index));
}

double local::AbsBinning::getBinWidth(int index) const {
    return getBinHighEdge(index) - getBinLowEdge(index);
}

void local::AbsBinning::dump(std::ostream &os) const {
    int nbins(getNBins());
    os << nbins;
    for(int bin = 0; bin <= nbins; ++bin) os << ' ' << getBinLowEdge(bin);
    os << std::endl;
}
