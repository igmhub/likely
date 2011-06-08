// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FunctionMinimum.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"

#include <iostream>
#include <cmath>

// Declare a binding to this LAPACK Cholesky decomposition routine:
// http://www.netlib.org/lapack/double/dpptrf.f
extern "C" {
    void dpptrf_(char* uplo, int* n, double* ap, int *info);
}

namespace local = likely;

local::FunctionMinimum::FunctionMinimum(double minValue, Parameters const& where)
: _minValue(minValue), _where(where), _haveCovariance(false)
{
}

local::FunctionMinimum::FunctionMinimum(double minValue, Parameters const& where,
PackedCovariance const &covar)
: _minValue(minValue), _where(where), _haveCovariance(true), _covar(covar)
{
    int nPar(_where.size());
    if(_covar.size() != nPar*(nPar+1)/2) {
        throw RuntimeError(
            "FunctionMinimum: parameter and covariance vectors have incompatible sizes.");
    }
}

local::FunctionMinimum::~FunctionMinimum() { }

local::Parameters local::FunctionMinimum::getErrors() const {
    if(!haveCovariance()) {
        throw RuntimeError("FunctionMinimum::getErrors: no covariance matrix available.");
    }
    int nPar(_where.size());
    Parameters errors(nPar);
    for(int i = 0; i < nPar; ++i) {
        double sigsq(_covar[i*(i+3)/2]);
        errors[i] = sigsq > 0 ? std::sqrt(sigsq) : 0;
    }
    return errors;
}

local::Parameters local::FunctionMinimum::getRandomParameters() const {
    Parameters params(_where.size());
    setRandomParameters(params);
    return params;
}

void local::FunctionMinimum::setRandomParameters(Parameters &params) const {
    if(!haveCovariance()) {
        throw RuntimeError(
            "FunctionMinimum::getRandomParameters: no covariance matrix available.");
    }
    int nPar(_where.size());
    // Compute the Cholesky decomposition of our covariance matrix if necessary.
    if(!_haveCholesky) {
        // Copy our covariance matrix.
        _cholesky = _covar;
        // Use LAPACK to perform the decomposition.
        char uplo('U');
        int info(0);
        dpptrf_(&uplo,&nPar,&_cholesky[0],&info);
        std::cout << "info = " << info << std::endl;
        for(int i = 0; i < _cholesky.size(); ++i) {
            std::cout << i << ' ' << _cholesky[i] << std::endl;
        }
        _haveCholesky = true;
    }
    // Fill a vector of random Gaussian variables.
    Parameters gauss(nPar);
    for(int i = 0; i < nPar; ++i) gauss[i] = 0;
}

void local::FunctionMinimum::printToStream(std::ostream &os, std::string formatSpec) const {
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
        for(int i = 0; i < nPar; ++i) {
            for(int j = 0; j < nPar; ++j) {
                int index = (i <= j) ? i + j*(j+1)/2 : j + i*(i+1)/2;
                os << ' ' << formatter % _covar[index];
            }
            os << std::endl;
        }
    }
}