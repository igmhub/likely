// Created 08-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FitModel.h"
#include "likely/RuntimeError.h"
#include "likely/AbsEngine.h"

#include <iostream>

namespace local = likely;

local::FitModel::FitModel(std::string const &name)
: _name(name)
{ }

local::FitModel::~FitModel() { }

void local::FitModel::defineParameter(std::string const &name, double value, double error) {
    _parameters.push_back(FitParameter(name,value,error));
}

int local::FitModel::getNParameters(bool onlyFloating) const {
    return onlyFloating ? countFloatingFitParameters(_parameters) : _parameters.size();
}

void  local::FitModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    out << "Fit Model \"" << _name << "\" has initial parameters:" << std::endl;
    printFitParametersToStream(_parameters,out,formatSpec);
}

void local::FitModel::configure(std::string const &script) {
    modifyFitParameters(_parameters,script);
}

local::FunctionMinimumPtr local::FitModel::findMinimum(FunctionPtr fptr, std::string const &method) const {
    return local::findMinimum(fptr, _parameters, method);
}