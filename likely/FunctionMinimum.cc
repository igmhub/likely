// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FunctionMinimum.h"
#include "likely/RuntimeError.h"
#include "likely/CovarianceMatrix.h"

#include "boost/format.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

#include <cassert>

namespace local = likely;

local::FunctionMinimum::FunctionMinimum(double minValue, FitParameters const &parameters)
: _status(OK)
{
    updateParameters(minValue, parameters);
}

local::FunctionMinimum::FunctionMinimum(double minValue, FitParameters const &parameters,
CovarianceMatrixCPtr covariance)
: _status(OK)
{
    updateParameters(minValue, parameters);
    updateCovariance(covariance);
}

local::FunctionMinimum::~FunctionMinimum() { }

void local::FunctionMinimum::updateParameters(double minValue, FitParameters const &parameters) {
    if(_parameters.size() > 0) {
        if(parameters.size() != _parameters.size()) {
            throw RuntimeError("FunctionMinimum::updateParameters: unexpected parameters size.");
        }
        if(countFloatingFitParameters(parameters) != _nFloating) {
            throw RuntimeError("FunctionMinimum::updateParameters: wrong number of floating parameters.");
        }
    }
    else {
        _nFloating = countFloatingFitParameters(parameters);
    }
    _parameters = parameters;
    _minValue = minValue;
}

void local::FunctionMinimum::updateParameterValues(double minValue, Parameters const &values) {
    if(_parameters.size() == 0) {
        throw RuntimeError("FunctionMinimum::updateParameters: no parameters set yet.");
    }
    if(values.size() != _parameters.size()) {
        throw RuntimeError("FunctionMinimum::updateParameters: unexpected number of values.");
    }
    FitParameters updated;
    int nextFloatIndex(0);
    updated.reserve(_parameters.size());
    for(int k = 0; k < _parameters.size(); ++k) {
        double error(_parameters[k].getError());
        if(0 != error && hasCovariance()) {
            error = std::sqrt(_covar->getCovariance(nextFloatIndex,nextFloatIndex));
            nextFloatIndex++;
        }
        updated.push_back(FitParameter(_parameters[k].getName(),values[k],error));
    }
    _parameters = updated;
    _minValue = minValue;
}

void local::FunctionMinimum::setParameterValue(std::string const &name, double value) {
    int index = findName(name);
    double error(_parameters[index].getError());
    if(0 != error && hasCovariance()) {
        error = std::sqrt(_covar->getCovariance(index,index));
    }
    _parameters[index] = FitParameter(name,value,error);
}

void local::FunctionMinimum::filterParameterValues(
Parameters const &allValues, Parameters &floatingValues) const {
    floatingValues.resize(0);
    floatingValues.reserve(_nFloating);
    for(int k = 0; k < _parameters.size(); ++k) {
        if(_parameters[k].isFloating()) floatingValues.push_back(allValues[k]);
    }
}

void local::FunctionMinimum::updateCovariance(CovarianceMatrixCPtr covariance) {
    if(covariance->getSize() != _nFloating) {
        assert(0);
        throw RuntimeError("FunctionMinimum: covariance size != number of floating parameters.");
    }
    _covar = covariance;
}

local::Parameters local::FunctionMinimum::getParameters(bool onlyFloating) const {
    Parameters params;
    getFitParameterValues(_parameters,params,onlyFloating);
    return params;
}

local::Parameters local::FunctionMinimum::getErrors(bool onlyFloating) const {
    Parameters errors;
    getFitParameterErrors(_parameters,errors,onlyFloating);
    return errors;
}

std::vector<std::string> local::FunctionMinimum::getNames(bool onlyFloating) const {
    std::vector<std::string> names;
    getFitParameterNames(_parameters,names,onlyFloating);
    return names;
}

int local::FunctionMinimum::findName(std::string const &name) const {
    return findFitParameterByName(_parameters,name);
}

double local::FunctionMinimum::setRandomParameters(Parameters &params) const {
    if(!hasCovariance()) {
        throw RuntimeError(
            "FunctionMinimum::getRandomParameters: no covariance matrix available.");
    }
    // Generate random offsets for our floating parameters.
    std::vector<double> floating;
    double nlWeight = _covar->sample(floating);
    std::vector<double>::const_iterator nextOffset(floating.begin());
    // Prepare to fill the parameter values vector we are provided.
    params.resize(0);
    params.reserve(_parameters.size());
    for(FitParameters::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter) {
        double value(iter->getValue());
        if(iter->isFloating()) value += *nextOffset++;
        params.push_back(value);
    }
    return nlWeight;
}

void local::FunctionMinimum::printToStream(std::ostream &os, std::string formatSpec) const {
    boost::format formatter(formatSpec), label("%20s = ");
    std::vector<std::string> labels;
    if(getStatus() == ERROR) {
        os << "FMIN ERROR: " << getStatusMessage() << std::endl;
        // Don't print out any more since it is presumably wrong.
        return;
    }
    else if(getStatus() == WARNING) {
        os << "FMIN WARNING: " << getStatusMessage() << std::endl;
        // Keep going, after this warning...
    }
    os << "FMIN Value = " << formatter % _minValue << " at:" << std::endl;
    for(FitParameters::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter) {
        os << (label % iter->getName()) << (formatter % iter->getValue());
        if(iter->isFloating()) {
            os << " +/- " << formatter % iter->getError();
            labels.push_back(iter->getName());
        }
        os << std::endl;
    }
    os << std::endl << "Number of function evaluations used: " << getNEvalCount() << std::endl;
    if(getNGradCount() > 0) {
        os << "Number of gradient evaluations used: " << getNGradCount() << std::endl;
    }
    if(hasCovariance()) {
        os << std::endl << "FMIN Errors & Correlations =" << std::endl;
        _covar->printToStream(os,true,formatSpec,labels);
    }
}
