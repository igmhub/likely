// Created 20-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/CovarianceAccumulator.h"
#include "likely/RuntimeError.h"
#include "likely/CovarianceMatrix.h"
#include "likely/BinnedData.h"

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/covariance.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/variates/covariate.hpp"

using namespace boost::accumulators;

typedef accumulator_set<double, stats<
    tag::covariance<double, tag::covariate1> > > Accumulator;
typedef std::vector<Accumulator> Accumulators;

namespace local = likely;

namespace likely {
    struct CovarianceAccumulator::Implementation {
        std::vector<Accumulator> accumulators;
    }; // CovarianceAccumulator::Implementation
} // likely::

local::CovarianceAccumulator::CovarianceAccumulator(int size)
: _size(size), _pimpl(new Implementation())
{
    if(size <= 0) {
        throw RuntimeError("CovarianceAccumulator: expected size > 0.");
    }
    _pimpl->accumulators.resize((size*(size+1))/2);
}

local::CovarianceAccumulator::~CovarianceAccumulator() { }

void local::CovarianceAccumulator::accumulate(std::vector<double> const &vector) {
    if(vector.size() != _size) {
        throw RuntimeError("CovarianceAccumulator::accumulate: invalid vector size.");
    }
    accumulate(&vector[0]);
}

void local::CovarianceAccumulator::accumulate(double const *vector) {
    int index(0);
    for(int i = 0; i < _size; ++i) {
        double xi(vector[i]);
        for(int j = 0; j <= i; ++j) {
            _pimpl->accumulators[index++](xi, covariate1 = vector[j]);
        }
    }
}

void local::CovarianceAccumulator::accumulate(BinnedDataCPtr data) {
    int index(0);
    for(BinnedData::IndexIterator row = data->begin(); row != data->end(); ++row) {
        double xi(data->getData(*row));
        for(BinnedData::IndexIterator col = data->begin(); col <= row; ++col) {
            _pimpl->accumulators[index++](xi, covariate1 = data->getData(*col));
        }
    }
}

local::CovarianceMatrixPtr local::CovarianceAccumulator::getCovariance() const {
    CovarianceMatrixPtr ptr(new CovarianceMatrix(_size));
    int index(0);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            ptr->setCovariance(row,col,covariance(_pimpl->accumulators[index++]));
        }
    }
    return ptr;
}