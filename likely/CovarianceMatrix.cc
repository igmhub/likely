// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"

namespace local = likely;

local::CovarianceMatrix::CovarianceMatrix(int size)
: _size(size), _compressed(false)
{
    if(size <= 0) {
        throw RuntimeError("CovarianceMatrix: expected size > 0.");
    }
}

local::CovarianceMatrix::~CovarianceMatrix() { }
