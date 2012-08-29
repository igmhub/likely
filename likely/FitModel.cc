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
    return countFitParameters(_parameters, onlyFloating);
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

local::FunctionMinimumPtr local::FitModel::findMinimum(FunctionPtr fptr, std::string const &method,
std::string const &oneTimeConfig) {
    if(0 < oneTimeConfig.length()) {
        // Apply the config script to a copy of our parameters.
        FitParameters modified(_parameters);
        modifyFitParameters(modified,oneTimeConfig);
        // Minimize using the modified parameters.
        return local::findMinimum(fptr, modified, method);
    }
    else {
        // Minimize using un-modified parameters.
        return local::findMinimum(fptr, _parameters, method);
    }
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

double local::FitModel::evaluatePriors() const {
    double penalty(0);
    for(int index = 0; index < _parameters.size(); ++index) {
        FitParameter const &param(_parameters[index]);
        if(!param.isFloating()) continue;
        FitParameter::PriorType priorType = param.getPriorType();
        if(priorType == FitParameter::NoPrior) continue;
        double value = _parameterValue[index];
        double priorMin = param.getPriorMin();
        double priorMax = param.getPriorMax();
        double range = priorMax - priorMin;
        if(priorType == FitParameter::BoxPrior) {
            if(value < priorMin) {
                double diff((value - priorMin)/(1e-3*range));
                penalty += diff*diff/2;
            }
            else if(value > priorMax) {
                double diff((priorMax - value)/(1e-3*range));
                penalty += diff*diff/2;
            }
        }
        else if(priorType == FitParameter::GaussPrior) {
            double diff((value - 0.5*(priorMin+priorMax))/(0.5*range));
            penalty += diff*diff/2;
        }
    }
    return penalty;
}
