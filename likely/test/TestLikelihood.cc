// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/test/TestLikelihood.h"
#include "likely/RuntimeError.h"

#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"

#include <cmath>
#include <iostream>

namespace local = likely::test;

local::TestLikelihood::TestLikelihood(int npar, double sigma, double rho, double alpha)
: _npar(npar), _sigma(sigma), _rho(rho), _alpha(alpha), _trace(false)
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
    _delta = (1-rho)*tmp*sigmasq;
    if(_delta <= 0) {
        throw RuntimeError("TestLikelihood: invalid determinant.");
    }
    // Calculate the two distinct elements of the inverse covariance matrix.
    _inverseDiagonal = (1+(npar-2)*rho)/_delta;
    _inverseOffDiagonal = -rho/_delta;
}

local::TestLikelihood::~TestLikelihood() { }

likely::Parameters local::TestLikelihood::_transform(Parameters const &x) const {
    Parameters y(x);
    double sumsq(0);
    for(int i = 0; i < _npar; ++i) {
        double value(x[i]),valuesq(_alpha*value*value);
        y[i] += valuesq;
        sumsq += valuesq;
    }
    for(int i = 0; i < _npar; ++i) y[i] -= sumsq;
    return y;
}

double local::TestLikelihood::evaluate(Parameters const &x) const {
    if(x.size() != _npar) {
        throw RuntimeError("evaluate() called with wrong number of parameters: "
            + boost::lexical_cast<std::string>(x.size()));
    }
    // Apply a coordinate transform to add non-parabolic behavior
    Parameters y = _transform(x);
    // Evaluate the correlated Gaussian NLL in y.
    double arg1(0),arg2(0);
    for(int i = 0; i < _npar; ++i) {
        arg1 += y[i]*y[i];
        for(int j = i+1; j < _npar; ++j) {
            arg2 += y[i]*y[j];
        }
    }
    arg1 *= _inverseDiagonal/2;
    arg2 *= _inverseOffDiagonal;
    double result(arg1+arg2);
    // Print an optional trace.
    if(_trace) {
        boost::format pFormat("%.5f");
        std::cout << "TestLikelihood(" << pFormat % x[0];
        for(int i = 1; i < _npar; ++i) {
            std::cout << ',' << pFormat % x[i];
        }
        std::cout << ") = " << result << std::endl;
    }
    return result;
}

void local::TestLikelihood::evaluateGradient(Parameters const &x, Gradient &grad) const {
    if(x.size() != _npar) {
        throw RuntimeError("evaluteGradient() called with wrong number of parameters: "
            + boost::lexical_cast<std::string>(x.size()));
    }
    // Apply a coordinate transform to add non-parabolic behavior.
    Parameters y = _transform(x);
    double ysum(0);
    for(int i = 0; i < _npar; ++i) ysum += y[i];
    // Calculate the gradient components.
    double tmp((1-_rho)*_sigma*_sigma),delta(tmp*(1+(_npar-1)*_rho));
    for(int i = 0; i < _npar; ++i) {
        grad[i] = y[i]*(1 + 2*_alpha*x[i])/tmp - (_rho + 2*_alpha*x[i])*ysum/delta;
    }
    // Print an optional trace.
    if(_trace) {
        boost::format pFormat("%.5f");
        std::cout << "TestLikelihood(" << pFormat % x[0];
        for(int i = 1; i < _npar; ++i) {
            std::cout << ',' << pFormat % x[i];
        }
        std::cout << ") GRADIENT = (" << pFormat % grad[0];
        for(int i = 1; i < _npar; ++i) {
            std::cout << ',' << pFormat % grad[i];
        }
        std::cout << ')' << std::endl;
    }
}