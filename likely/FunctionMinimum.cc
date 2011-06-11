// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FunctionMinimum.h"
#include "likely/Random.h"
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
: _minValue(minValue), _where(where), _random(Random::instance())
{
}

local::FunctionMinimum::FunctionMinimum(double minValue, Parameters const& where,
PackedCovariance const &covar, bool errorsOnly)
: _minValue(minValue), _where(where), _random(Random::instance())
{
    updateCovariance(covar,errorsOnly);
}

local::FunctionMinimum::~FunctionMinimum() { }

void local::FunctionMinimum::updateParameters(Parameters const &params, double fval) {
    _minValue = fval;
    _where = params;
}

void local::FunctionMinimum::updateCovariance(PackedCovariance const &covar,
bool errorsOnly) {
    int nPar(_where.size()),nCovar(nPar*(nPar+1)/2);
    if(errorsOnly) {
        if(covar.size() != nPar) throw RuntimeError(
            "FunctionMinimum: parameter and error vectors have incompatible sizes.");
        // We have a vector of errors, instead of a full covariance matrix.
        _covar.reset(new PackedCovariance(nCovar,0));
        for(int i = 0; i < nPar; ++i) {
            double error(covar[i]);
            if(error <= 0) {
                throw RuntimeError("FunctionMinimum: errors must be > 0.");
            }
            (*_covar)[i*(i+3)/2] = error*error;
        }
    }
    else {
        if(covar.size() != nCovar) throw RuntimeError(
            "FunctionMinimum: parameter and covariance vectors have incompatible sizes.");
        _covar.reset(new PackedCovariance(covar));
    }
    // Forget any previously calculated Cholesky decomposition.
    _cholesky.reset();
    // Should check that covariance is positive definite here...    
}

local::Parameters local::FunctionMinimum::getErrors() const {
    if(!haveCovariance()) {
        throw RuntimeError("FunctionMinimum::getErrors: no covariance matrix available.");
    }
    int nPar(_where.size());
    Parameters errors(nPar);
    for(int i = 0; i < nPar; ++i) {
        double sigsq((*_covar)[i*(i+3)/2]);
        errors[i] = sigsq > 0 ? std::sqrt(sigsq) : 0;
    }
    return errors;
}

local::PackedCovariancePtr local::FunctionMinimum::getCholesky() const {
    if(!_cholesky) {
        // Copy our covariance matrix.
        _cholesky.reset(new PackedCovariance(*_covar));
        // Use LAPACK to perform the decomposition.
        char uplo('U');
        int info(0),nPar(_where.size());
        dpptrf_(&uplo,&nPar,&(*_cholesky)[0],&info);
        if(0 != info) {
            throw RuntimeError(
                "FunctionMinimum::setRandomParameters: Cholesky decomposition failed.");
        }        
    }
    return _cholesky;
}

double local::FunctionMinimum::setRandomParameters(Parameters &params) const {
    if(!haveCovariance()) {
        throw RuntimeError(
            "FunctionMinimum::getRandomParameters: no covariance matrix available.");
    }
    int nPar(_where.size());
    Parameters gauss(nPar);
    double nlWeight(0);
    for(int i = 0; i < nPar; ++i) {
        // Initialize the generated parameters to the function minimum.
        params[i] = _where[i];
        // Fill a vector of random Gaussian variables.
        double r(_random.getNormal());
        gauss[i] = r;
        nlWeight += r*r;
    }
    // Multiply by the Cholesky decomposition matrix.
    PackedCovariance::const_iterator next(getCholesky()->begin());
    for(int j = 0; j < nPar; ++j) {
        for(int i = 0; i <= j; ++i) {
            params[j] += (*next++)*gauss[i];
        }
    }
    return nlWeight/2;
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
        for(int i = 0; i < nPar; ++i) {
            for(int j = 0; j < nPar; ++j) {
                int index = (i <= j) ? i + j*(j+1)/2 : j + i*(i+1)/2;
                os << ' ' << formatter % (*_covar)[index];
            }
            os << std::endl;
        }
    }
}