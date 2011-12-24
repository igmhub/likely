// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/TriCubicInterpolator.h"
#include "likely/RuntimeError.h"

namespace local = likely;

local::TriCubicInterpolator::TriCubicInterpolator(DataCube data, int n1, int n2, int n3)
: _data(data), _n1(n1), _n2(n2), _n3(n3)
{
    if(_n2 == 0 && _n3 == 0) {
        _n3 = _n2 = _n1;
    }
    if(_n1 <= 0 || _n2 <= 0 || _n3 <= 0) throw RuntimeError("Bad datacube dimensions.");
}

local::TriCubicInterpolator::~TriCubicInterpolator() { }
