// Created 30-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FunctionMinimum.h"
#include "likely/RuntimeError.h"
#include "likely/CovarianceMatrix.h"

#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace local = likely;

local::FunctionMinimum::FunctionMinimum(double minValue, FitParameters const &parameters)
: _status(OK)
{
    updateParameters(minValue, parameters);
    setCounts(0,0);
}

local::FunctionMinimum::FunctionMinimum(double minValue, FitParameters const &parameters,
CovarianceMatrixCPtr covariance)
: _status(OK)
{
    updateParameters(minValue, parameters);
    updateCovariance(covariance);
    setCounts(0,0);
}

local::FunctionMinimum::~FunctionMinimum() { }

void local::FunctionMinimum::updateParameters(double minValue, FitParameters const &parameters) {
    if(_parameters.size() > 0) {
        if(parameters.size() != _parameters.size()) {
            throw RuntimeError("FunctionMinimum::updateParameters: unexpected parameters size.");
        }
        if(countFitParameters(parameters,true) != _nFloating) {
            throw RuntimeError("FunctionMinimum::updateParameters: wrong number of floating parameters.");
        }
    }
    else {
        _nFloating = countFitParameters(parameters,true);
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
    int nextFloatIndex(0);
    for(int k = 0; k < _parameters.size(); ++k) {
        // Assign the new value to this parameter.
        _parameters[k].setValue(values[k]);
        // Assign an error from our covariance if we have one and if this is a floating parameter.
        double error(_parameters[k].getError());
        if(0 != error && hasCovariance()) {
            error = std::sqrt(_covar->getCovariance(nextFloatIndex,nextFloatIndex));
            nextFloatIndex++;
            _parameters[k].setError(error);
        }
    }
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

double local::FunctionMinimum::setRandomParameters(Parameters const &fromParams,
Parameters &toParams) const {
    if(!hasCovariance()) {
        throw RuntimeError(
            "FunctionMinimum::getRandomParameters: no covariance matrix available.");
    }
    // Generate random offsets for our floating parameters.
    std::vector<double> floating;
    double nlWeight = _covar->sample(floating);
    std::vector<double>::const_iterator nextOffset(floating.begin());
    // Prepare to fill the parameter values vector we are provided.
    toParams=fromParams;
    int index(0);
    for(FitParameters::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter) {
      if(iter->isFloating()) toParams[index] += *nextOffset++;
      index++;
    }
    return nlWeight;
}

void local::FunctionMinimum::printToStream(std::ostream &os, std::string const &formatSpec) const {
    boost::format formatter(formatSpec);
    std::vector<std::string> labels;
    getFitParameterNames(_parameters,labels,true);
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
    printFitParametersToStream(_parameters,os,formatSpec);
    os << std::endl << "Number of function evaluations used: " << getNEvalCount() << std::endl;
    if(getNGradCount() > 0) {
        os << "Number of gradient evaluations used: " << getNGradCount() << std::endl;
    }
    if(hasCovariance()) {
        os << std::endl << "FMIN Errors & Correlations =" << std::endl;
        _covar->printToStream(os,true,formatSpec,labels);
    }
}

void local::FunctionMinimum::saveParameters(std::ostream &os, bool onlyFloating) const {    
    int index(0);
    for(FitParameters::const_iterator iter = _parameters.begin(); iter != _parameters.end(); ++iter,++index) {
        if(onlyFloating && !iter->isFloating()) continue;
        // Use lexical_cast to ensure that the full double precision is saved.
        os << index << ' ' << boost::lexical_cast<std::string>(iter->getValue())
            << ' ' << boost::lexical_cast<std::string>(iter->getError()) << std::endl;
    }
}

void local::FunctionMinimum::saveFloatingParameterCovariance(std::ostream &os, double scale) const {
    if(!_covar->isPositiveDefinite()) {
        throw RuntimeError("FunctionMinimum::saveFloatingParameterCovariance: matrix is not positive definite.");
    }
    int index1(0),floatingIndex1(0);
    for(FitParameters::const_iterator iter1 = _parameters.begin(); iter1 != _parameters.end(); ++iter1,++index1) {
        if(!iter1->isFloating()) continue;
        // Save all diagonal elements.
        double value = scale*_covar->getCovariance(floatingIndex1,floatingIndex1);
        os << index1 << ' ' << index1 << ' '
            << boost::lexical_cast<std::string>(value) << std::endl;        
        // Loop over pairs with index2 > index1
        int index2(index1+1),floatingIndex2(floatingIndex1+1);
        for(FitParameters::const_iterator iter2 = iter1; ++iter2 != _parameters.end(); ++index2) {
            if(!iter2->isFloating()) continue;
            value = scale*_covar->getCovariance(floatingIndex1,floatingIndex2);
            // Use lexical_cast to ensure that the full double precision is saved.
            os << index1 << ' ' << index2 << ' ' << boost::lexical_cast<std::string>(value) << std::endl;
            floatingIndex2++;
        }
        floatingIndex1++;
    }
}
