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

#include "likely/UniformSampling.h"
#include "likely/NonUniformSampling.h"

#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

local::AbsBinningCPtr local::createBinning(std::string const &binningSpec) {
    // import boost spirit parser symbols
    using qi::double_;
    using qi::_1;
    using qi::lit;
    using phoenix::ref;
    using phoenix::push_back;
    
    std::vector<double> centers;

    // Parse the points string into a vector of doubles.
    std::string::const_iterator iter = binningSpec.begin();
    bool ok = qi::phrase_parse(iter,binningSpec.end(),
        (
            double_[push_back(ref(centers),_1)] % ','
        ),
        ascii::space);
    if(!ok || iter != binningSpec.end()) {
        throw RuntimeError("createBinning: badly formatted specification string.");
    }
    AbsBinningCPtr binning;
    if(centers.size() > 2) {
        binning.reset(new NonUniformSampling(centers));
    }
    else {
        binning.reset(new UniformSampling(centers.front(),centers.back(),centers.size()));
    }
    return binning;
}
