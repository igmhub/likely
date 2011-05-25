// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslMinimizer.h"
#include "likely/GslEngine.h"

namespace local = likely;

local::GslMinimizer::GslMinimizer(Function f, int nPar)
: _engine(new GslEngine(f,nPar))
{
}

local::GslMinimizer::~GslMinimizer() { }

local::Parameters local::GslMinimizer::minimize(
Parameters const& initial, Parameters const &errors) {
    Parameters final(initial);
    return final;
}