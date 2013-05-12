// Created 05-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the FitParameter class.

#include "likely/FitModel.h"
#include "likely/FitParameter.h"
#include "likely/RuntimeError.h"

#include <iostream>
#include <cassert>

namespace lk = likely;

class TestModel : public lk::FitModel {
public:
    TestModel() : lk::FitModel("Test Model") {
        defineParameter("p1",1);
        defineParameter("p2",2);
        defineParameter("p3",3);
    }
    virtual ~TestModel() { }
    void printChangedAndReset() {
        for(int index = 0; index < getNParameters(); ++index) {
            std::cout << "index " << index << " value " << getParameterValue(index)
                << " changed? " << isParameterValueChanged(index) << std::endl;
        }
        resetParameterValuesChanged();
    }
    bool update(lk::Parameters const &values) {
        return updateParameterValues(values);
    }
};

int main(int argc, char *argv[]) {
    lk::FitParameters params;
    params.push_back(lk::FitParameter("param1",1,0.1));
    params.push_back(lk::FitParameter("param2",2,0.2));
    params.push_back(lk::FitParameter("param3",3,0.3));
    params.push_back(lk::FitParameter("(1-beta)*bias",1.23));
    lk::printFitParametersToStream(params,std::cout);

    // This should be ok
    lk::modifyFitParameters(params,"");
    // This should fail since whitespace by itself is not a valid config script
    try {
        lk::modifyFitParameters(params," ");
    }
    catch(lk::RuntimeError const &e) {
        std::cout << "Got expected error with config = \" \"" << std::endl;
    }

    lk::modifyFitParameters(params," fix [param2] = -2; fix [param1]");
    lk::printFitParametersToStream(params,std::cout);    

    lk::modifyFitParameters(params,"release [par*] ;");
    lk::printFitParametersToStream(params,std::cout);    

    lk::modifyFitParameters(params,"value[param3*]=-123; error[param1]=1e-2");
    lk::printFitParametersToStream(params,std::cout);

    lk::modifyFitParameters(params,"error [(1-beta)*bias] = 0.5");
    lk::printFitParametersToStream(params,std::cout);
    
    lk::modifyFitParameters(params,"boxprior[param2] @ (-1,2.3)");
    lk::modifyFitParameters(params,"gaussprior[param1] @ (0.5,1.5;0.5)");
    
    lk::modifyFitParameters(params,"binning[param1] = [-1,+1]*5");
    lk::modifyFitParameters(params,"binning[param3] = 0.1,0.2,0.4");

    try {
        lk::modifyFitParameters(params,"value [param3]=0;error [param3] = -123");
    }
    catch(lk::RuntimeError const &e) {
        // We expect this since error < 0 is not allowed.
    }
    // Check that the parameters were not actually modified.
    lk::printFitParametersToStream(params,std::cout);
    std::cout
        << "-- script begin" << std::endl
        << lk::fitParametersToScript(params)
        << "-- script end" << std::endl;
    
    try {
        lk::FitParameter badName("a,b",123);
    }
    catch(lk::RuntimeError const &e) {
        // We expect this since commas are not allowed in names.
        std::cout << "Got expected RuntimeError" << std::endl;
    }    
    try {
        lk::modifyFitParameters(params,"release [aram*]");
    }
    catch(lk::RuntimeError const &e) {
        // We expect this since the pattern does not match any names.
        std::cout << "Got expected RuntimeError" << std::endl;
    }
    
    // Test grid iteration
    lk::BinnedGrid grid = getFitParametersGrid(params);
    for(lk::BinnedGrid::Iterator iter = grid.begin(); iter != grid.end(); ++iter) {
        std::string config = getFitParametersGridConfig(params,grid,iter);
        std::cout << "config[" << (*iter) << "]: " << config << std::endl;
    }
    
    TestModel model;
    model.printToStream(std::cout);
    model.printChangedAndReset();
    model.printChangedAndReset();
    model.setParameterValue("p2",22);
    model.printChangedAndReset();
    lk::Parameters pvalues(model.getNParameters());
    pvalues[0] = -1;
    pvalues[1] = model.getParameterValue("p2");
    pvalues[2] = model.getParameterValue("p3");
    std::cout << "update any changes? " << model.update(pvalues) << std::endl;
    model.printChangedAndReset();
    pvalues[1] = -2;
    pvalues[2] = -3;
    std::cout << "update any changes? " << model.update(pvalues) << std::endl;
    model.printChangedAndReset();
    std::cout << "update any changes? " << model.update(pvalues) << std::endl;
    model.printChangedAndReset();
}
