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
#include "likely/NonUniformBinning.h"

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

            bspec = ( blist | slist | brange | srange );
            
            // Parse the formats [x1,x2,...,xn] and {x1,x2,...,xn}
            blist = ( '[' >> valuesList >> ']' )[boost::bind(&Grammar::createBinsWithEdges,this)];
            slist = ( '{' >> valuesList >> '}' )[boost::bind(&Grammar::createSamplesWithCenters,this)];
                
            valuesList = ( double_[push_back(ref(values),_1)] % ',' );
            
            // Parse the formats [lo:hi]*n and {lo:hi}*n
            brange = ( '[' >> double_[ref(lo)=_1] >> ':' >> double_[ref(hi)=_1] >> "]*" >> int_[ref(nbins)=_1] )[
                boost::bind(&Grammar::createBinsWithRange,this)];
            srange = ( '{' >> double_[ref(lo)=_1] >> ':' >> double_[ref(hi)=_1] >> "}*" >> int_[ref(nbins)=_1] )[
                boost::bind(&Grammar::createSamplesWithRange,this)];

        }
        qi::rule<std::string::const_iterator> bspec,blist,slist,valuesList,brange,srange;
        likely::AbsBinningCPtr binning;

        std::vector<double> values;
        double lo,hi;
        int nbins;
        
        void createBinsWithEdges() {
            binning.reset(new likely::NonUniformBinning(values));
        }
        void createSamplesWithCenters() {
            if(values.size() > 2) {
                binning.reset(new likely::NonUniformSampling(values));
            }
            else {
                binning.reset(new likely::UniformSampling(values.front(),values.back(),values.size()));
            }
        }        
        void createBinsWithRange() {
            binning.reset(new likely::UniformBinning(lo,hi,nbins));
        }
        void createSamplesWithRange() {
            binning.reset(new likely::UniformSampling(lo,hi,nbins));
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
