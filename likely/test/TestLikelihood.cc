// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/test/TestLikelihood.h"
#include "likely/RuntimeError.h"

#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"

#include <cmath>
#include <iostream>

namespace local = likely::test;

local::TestLikelihood::TestLikelihood(int npar, double sigma, double rho, double alpha)
: _npar(npar), _sigma(sigma), _rho(rho), _alpha(alpha), _trace(false), _counts(0,0)
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
    double determinant(std::pow(1-rho,npar-1)*tmp*std::pow(sigmasq,npar));
    if(determinant == 0) {
        throw RuntimeError("TestLikelihood: determinant is zero.");
    }
    // Calculate the two distinct elements of the inverse covariance matrix.
    double denom((1-rho)*tmp*sigmasq);
    _inverseDiagonal = (1+(npar-2)*rho)/denom;
    _inverseOffDiagonal = -rho/denom;
}

local::TestLikelihood::~TestLikelihood() { }

double local::TestLikelihood::evaluate(Parameters const &params) const {
    if(params.size() != _npar) {
        throw RuntimeError("TestLikelihood() called with wrong number of parameters: "
            + boost::lexical_cast<std::string>(params.size()));
    }
    // Make a non-linear transformation of the input parameters to internal
    // Gaussian parameters.
    Parameters internalParams(params);
    if(_alpha != 0 && _npar > 1) {
        double normSq(0);
        for(int i = 1; i < _npar; ++i) {
            normSq += params[i]*params[i];
        }
        internalParams[0] -= _alpha*normSq/(_npar-1);
    }
    // Evaluate the correlated Gaussian probability density using internal parameters.
    double arg1(0),arg2(0);
    for(int i = 0; i < _npar; ++i) {
        arg1 += internalParams[i]*internalParams[i];
        for(int j = i+1; j < _npar; ++j) {
            arg2 += internalParams[i]*internalParams[j];
        }
    }
    arg1 *= _inverseDiagonal/2;
    arg2 *= _inverseOffDiagonal;
    double result(arg1+arg2);
    
    _counts.first++;    
    if(_trace) {
        boost::format pFormat("%.5f");
        std::cout << '[' << _counts.first << "] TestLikelihood(" << pFormat % params[0];
        for(int i = 1; i < _npar; ++i) {
            std::cout << ',' << pFormat % params[i];
        }
        std::cout << ") = " << result << std::endl;
    }
    return result;
}

void local::TestLikelihood::evaluateGradient(Parameters const &params,
Gradient &grad) const {
    _counts.second++;
}