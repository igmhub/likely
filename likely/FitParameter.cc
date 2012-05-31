// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FitParameter.h"

namespace local = likely;

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
