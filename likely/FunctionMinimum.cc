// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FunctionMinimum.h"
#include "likely/RuntimeError.h"
#include "likely/CovarianceMatrix.h"

#include "boost/format.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace local = likely;

local::FunctionMinimum::FunctionMinimum(double minValue, Parameters const& where)
: _minValue(minValue), _where(where)
{
}

local::FunctionMinimum::FunctionMinimum(double minValue, Parameters const& where,
CovarianceMatrixCPtr covariance)
: _minValue(minValue), _where(where)
{
    updateCovariance(covariance);
}

local::FunctionMinimum::~FunctionMinimum() { }

void local::FunctionMinimum::updateParameters(Parameters const &params, double fval) {
    _minValue = fval;
    _where = params;
}

void local::FunctionMinimum::updateCovariance(CovarianceMatrixCPtr covariance) {
    if(!covariance) return;
    if(_where.size() != covariance->getSize()) {
        throw RuntimeError("FunctionMinimum: incompatible sizes for parameters and covariance.");
    }
    _covar = covariance;
}

local::Parameters local::FunctionMinimum::getErrors() const {
    if(!haveCovariance()) {
        throw RuntimeError("FunctionMinimum::getErrors: no covariance matrix available.");
    }
    int nPar(_where.size());
    Parameters errors(nPar);
    for(int i = 0; i < nPar; ++i) {
        double sigsq(_covar->getCovariance(i,i));
        errors[i] = sigsq > 0 ? std::sqrt(sigsq) : 0;
    }
    return errors;
}

double local::FunctionMinimum::setRandomParameters(Parameters &params) const {
    if(!haveCovariance()) {
        throw RuntimeError(
            "FunctionMinimum::getRandomParameters: no covariance matrix available.");
    }
    int nPar(_where.size());
    Parameters gauss(nPar);
    double nlWeight = _covar->sample(params);
    for(int i = 0; i < nPar; ++i) {
        params[i] += _where[i];
    }
    return nlWeight;
}

void local::FunctionMinimum::printToStream(std::ostream &os,
std::string formatSpec) const {
    boost::format formatter(formatSpec);
    os << "F(" << formatter % _where[0];
    int nPar(_where.size());
    for(int i = 1; i < nPar; ++i) {
        os << ',' << formatter % _where[i];
    }
    os << ") = " << formatter % _minValue << std::endl;
    if(haveCovariance()) {
        Parameters errors(getErrors());
        os << "ERRORS:";
        for(int i = 0; i < nPar; ++i) {
            os << ' ' << formatter % errors[i];
        }
        os << std::endl << "COVARIANCE:" << std::endl;
        _covar->printToStream(os,formatSpec);
    }
}
