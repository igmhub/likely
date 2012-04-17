// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/AbsCovarianceMatrix.h"
#include "likely/RuntimeError.h"

namespace local = likely;

local::AbsCovarianceMatrix::AbsCovarianceMatrix(int size)
: _size(size), _compressed(false)
{
    if(size <= 0) {
        throw RuntimeError("AbsCovarianceMatrix: expected size > 0.");
    }
}

local::AbsCovarianceMatrix::~AbsCovarianceMatrix() { }
