// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/test/TestLikelihood.h"
#include "likely/RuntimeError.h"

#include "boost/lexical_cast.hpp"
#include "boost/foreach.hpp"

#include <cmath>

#include <iostream> // !!!

namespace local = likely::test;

local::TestLikelihood::TestLikelihood(int npar, double sigma, double rho, double alpha)
: _npar(npar), _sigma(sigma), _rho(rho), _alpha(alpha)
{
    // Check for valid inputs.
    if(npar < 0) {
        throw RuntimeError("TestLikelihood: invalid npar = " +
            boost::lexical_cast<std::string>(npar));
    }
    if(rho <= -1 || rho >= +1) {
        throw RuntimeError("TestLikelihood: invalid rho = " +
            boost::lexical_cast<std::string>(rho));
    }
    double tmp(1 + (npar-1)*rho),sigmasq(sigma*sigma);
    if(tmp == 0) {
        throw RuntimeError("TestLikelihood: determinant is zero.");
    }
    // Calculate the two distinct elements of the inverse covariance matrix.
    double denom(tmp*(1-rho)*sigmasq);
    _inverseDiagonal = (1+(npar-2)*rho)/denom;
    _inverseOffDiagonal = -rho/denom;
    // Calculate the covariance matrix determinant.
    double determinant(std::pow(1-rho,npar-1)*tmp*std::pow(sigmasq,npar));
    // Calculate the normalization factor.
    double twopi(8*std::atan(1));
    _norm = std::pow(twopi,0.5*npar)*std::sqrt(std::fabs(determinant));
    std::cout << _inverseDiagonal << ", " << _inverseOffDiagonal << " det="
        << determinant << std::endl;
}

local::TestLikelihood::~TestLikelihood() { }

double local::TestLikelihood::operator()(std::vector<double> const &params) const {
    if(params.size() != _npar) {
        throw RuntimeError("TestLikelihood() called with wrong number of parameters.");
    }
    double arg(0);
    for(int i = 0; i < _npar; ++i) {
        arg += (_inverseDiagonal - _inverseOffDiagonal)*params[i];
        for(int j = 0; j < _npar; ++j) {
            arg += _inverseOffDiagonal*params[i]*params[j];
        }
    }
    return std::exp(-arg/2)/_norm;
}
