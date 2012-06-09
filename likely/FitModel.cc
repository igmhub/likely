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
    _nameIndexMap.insert(NameIndexMap::value_type(name,_parameters.size()));
    _parameters.push_back(FitParameter(name,value,error));
    _parameterValue.push_back(value);
    _parameterValueChanged.push_back(true);
}

int local::FitModel::getNParameters(bool onlyFloating) const {
    return onlyFloating ? countFloatingFitParameters(_parameters) : _parameters.size();
}

bool local::FitModel::updateParameterValues(Parameters const &values) {
    int nValues(_parameterValue.size());
    if(values.size() != nValues) {
        throw RuntimeError("FitModel::updateParameterValues: invalid values size.");
    }
    bool anyChanged(false);
    for(int index = 0; index < nValues; ++index) {
        _setParameterValue(index,values[index]);
        if(_parameterValueChanged[index]) anyChanged = true;
    }
    return anyChanged;
}

void  local::FitModel::printToStream(std::ostream &out, std::string const &formatSpec) const {
    out << "Fit Model \"" << _name << "\" has initial parameters:" << std::endl;
    printFitParametersToStream(_parameters,out,formatSpec);
}

void local::FitModel::configureFitParameters(std::string const &script) {
    modifyFitParameters(_parameters,script);
}

local::FunctionMinimumPtr local::FitModel::findMinimum(FunctionPtr fptr, std::string const &method) {
    return local::findMinimum(fptr, _parameters, method);
}

int local::FitModel::_checkIndex(int index) const {
    if(index < 0 || index >= _parameters.size()) {
        throw RuntimeError("FitModel: invalid parameter index.");
    }
    return index;
}

int local::FitModel::_getIndex(std::string const &name) const {
    // Could remember the last find result to speed this up if necessary.
    NameIndexMap::const_iterator where = _nameIndexMap.find(name);
    if(where == _nameIndexMap.end()) {
        throw RuntimeError("FitModel: unknown parameter \"" + name + "\"");
    }
    return where->second;
}