// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FitParameter.h"
#include "likely/RuntimeError.h"

#include <iterator>

namespace local = likely;

local::FitParameter::FitParameter(std::string const &name, double value, double error)
: _name(name)
{
    setValue(value);
    setError(error);
}

void local::FitParameter::setError(double error) {
    if(error < 0) {
        throw RuntimeError("FitParameter::setError: error must be >= 0.");
    }
    _error = error;
}

void local::getFitParameterValues(FitParameters const &parameters, Parameters &values, bool onlyFloating) {
    values.resize(0);
    if(!onlyFloating) values.reserve(parameters.size());
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(!onlyFloating || iter->isFloating()) values.push_back(iter->getValue());
    }
}

void local::getFitParameterErrors(FitParameters const &parameters, Parameters &errors, bool onlyFloating) {
    errors.resize(0);
    if(!onlyFloating) errors.reserve(parameters.size());
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(!onlyFloating || iter->isFloating()) errors.push_back(iter->getError());
    }
}

void local::getFitParameterNames(FitParameters const &parameters, std::vector<std::string> &names,
bool onlyFloating) {
    names.resize(0);
    if(!onlyFloating) names.reserve(parameters.size());
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(!onlyFloating || iter->isFloating()) names.push_back(iter->getName());
    }
}

int local::countFloatingFitParameters(FitParameters const &parameters) {
    int count(0);
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(iter->isFloating()) count++;
    }
    return count;
}

int local::findFitParameterByName(FitParameters const &parameters, std::string const &name) {
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(name == iter->getName()) return std::distance(parameters.begin(),iter);
    }
    throw RuntimeError("findFitParameterByName: name not found: " + name);
}