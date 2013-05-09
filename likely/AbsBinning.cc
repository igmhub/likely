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

void local::AbsBinning::printToStream(std::ostream &os) const {
    // The default format is a comma-separated list of bin centers, which is one of
    // the formats supported by createBinning. Subclasses can implement alternate formats.
    int nbins(getNBins());
    for(int bin = 0; bin < nbins; ++bin) {
        if(bin) os << ',';
        os << getBinCenter(bin);
    }
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
#include "likely/UniformBinning.h"

#include "boost/bind.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

namespace likely {
namespace binning {
    // Declare our parser grammar (no skip parser since whitespace is not allowed)
    struct Grammar : qi::grammar<std::string::const_iterator> {

        Grammar() : base_type(bspec) {

            using qi::double_;
            using qi::int_;
            using qi::_1;
            using qi::lit;
            using phoenix::ref;
            using phoenix::push_back;

            bspec = plist | brange;
            
            // Parse the format x1,x2,...,xn and push values into the centers vector
            plist = ( double_[push_back(ref(centers),_1)] % ',' )[boost::bind(&Grammar::createWithCenters,this)];
            
            // Parse the format [lo,hi]*n and fill the corresponding data members
            brange = ( '[' >> double_[ref(lo)=_1] >> ',' >> double_[ref(hi)=_1] >> "]*" >> int_[ref(nbins)=_1] )[
                boost::bind(&Grammar::createWithRange,this)];

        }
        qi::rule<std::string::const_iterator> bspec,plist,brange;
        likely::AbsBinningCPtr binning;

        std::vector<double> centers;
        double lo,hi;
        int nbins;
        
        void createWithCenters() {
            if(centers.size() > 2) {
                binning.reset(new likely::NonUniformSampling(centers));
            }
            else {
                binning.reset(new likely::UniformSampling(centers.front(),centers.back(),centers.size()));
            }
        }
        
        void createWithRange() {
            binning.reset(new likely::UniformBinning(lo,hi,nbins));
        }
    };
} // binning
} // baofit

local::AbsBinningCPtr local::createBinning(std::string const &binningSpec) {
    // Parse the binning specification string and create the corresponding binning object.
    binning::Grammar parser;
    std::string::const_iterator iter = binningSpec.begin();
    bool ok = qi::parse(iter,binningSpec.end(),parser);
    if(!ok || iter != binningSpec.end()) {
        throw RuntimeError("createBinning: badly formatted specification string.");
    }
    return parser.binning;
}
