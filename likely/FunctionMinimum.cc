// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FunctionMinimum.h"

#include "boost/format.hpp"

#include <iostream>

namespace local = likely;

local::FunctionMinimum::FunctionMinimum(double minValue, Parameters const& where)
: _minValue(minValue), _where(where)
{
}

local::FunctionMinimum::~FunctionMinimum() { }

void local::FunctionMinimum::printToStream(std::ostream &os, std::string formatSpec) const {
    boost::format formatter(formatSpec);
    os << "F(" << formatter % _where[0];
    for(int i = 1; i < _where.size(); ++i) {
        os << ',' << formatter % _where[i];
    }
    os << ") = " << formatter % _minValue << std::endl;
}