// Created 20-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/CovarianceAccumulator.h"
#include "likely/RuntimeError.h"
#include "likely/CovarianceMatrix.h"
#include "likely/BinnedData.h"

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/weighted_covariance.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/count.hpp"
#include "boost/accumulators/statistics/variates/covariate.hpp"
#include "boost/lexical_cast.hpp"

#include <iostream>

using namespace boost::accumulators;

typedef accumulator_set<double, stats<
    tag::weighted_covariance<double, tag::covariate1> >, double > Accumulator;
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

void local::CovarianceAccumulator::accumulate(std::vector<double> const &vector, double wgt) {
    if(vector.size() != _size) {
        throw RuntimeError("CovarianceAccumulator::accumulate: invalid vector size.");
    }
    accumulate(&vector[0],wgt);
}

void local::CovarianceAccumulator::accumulate(double const *vector, double wgt) {
    int index(0);
    for(int i = 0; i < _size; ++i) {
        double xi(vector[i]);
        for(int j = 0; j <= i; ++j) {
            _pimpl->accumulators[index++](xi, weight = wgt, covariate1 = vector[j]);
        }
    }
}

void local::CovarianceAccumulator::accumulate(BinnedDataCPtr data, double wgt) {
    if(data->getNBinsWithData() != _size) {
        throw RuntimeError("CovarianceAccumulator::accumulate: invalid data size.");
    }
    int index(0);
    bool weighted(false);
    for(BinnedData::IndexIterator row = data->begin(); row != data->end(); ++row) {
        double xi(data->getData(*row,weighted));
        for(BinnedData::IndexIterator col = data->begin(); col <= row; ++col) {
            _pimpl->accumulators[index++](xi, weight = 1, covariate1 = data->getData(*col,weighted));
        }
    }
}

int local::CovarianceAccumulator::count() const {
    return boost::accumulators::count(_pimpl->accumulators[0]);
}

local::CovarianceMatrixPtr local::CovarianceAccumulator::getCovariance() const {
    CovarianceMatrixPtr cov(new CovarianceMatrix(_size));
    int index(0);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            double value(weighted_covariance(_pimpl->accumulators[index++]));
            cov->setCovariance(row,col,value);
        }
    }
    return cov;
}

void local::CovarianceAccumulator::dump(std::ostream &out) const {
    // matrix dimension
    out << _size << std::endl;
    // number of samples accumulated
    out << count() << std::endl;
    // total weight of accumulated samples (use lexical_cast to get full precision)
    out << boost::lexical_cast<std::string>(
        sum_of_weights(_pimpl->accumulators[0])) << std::endl;
    // weighted means
    for(int col = 0; col < _size; ++col) {
        int index = symmetricMatrixIndex(col,col,_size);
        out << col << ' ' << boost::lexical_cast<std::string>(
            weighted_mean(_pimpl->accumulators[index])) << std::endl;
    }
    // weighted second moments
    int index(0);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            out << row << ' ' << col << ' ' << boost::lexical_cast<std::string>(
                weighted_covariance(_pimpl->accumulators[index++])) << std::endl;
        }
    }
}
