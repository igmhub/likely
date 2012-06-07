// Created 05-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the FitParameter class.

#include "likely/likely.h"

#include <iostream>
#include <cassert>

namespace lk = likely;

int main(int argc, char *argv[]) {
    lk::FitParameters params;
    params.push_back(lk::FitParameter("param1",1,0.1));
    params.push_back(lk::FitParameter("param2",2,0.2));
    params.push_back(lk::FitParameter("param3",3,0.3));
    params.push_back(lk::FitParameter("(1-beta)*bias",1.23));
    lk::printFitParametersToStream(params,std::cout);    

    lk::modifyFitParameters(params," fix [param2] = -2; fix [param1]");
    lk::printFitParametersToStream(params,std::cout);    

    lk::modifyFitParameters(params,"release [param1] ;");
    lk::printFitParametersToStream(params,std::cout);    

    lk::modifyFitParameters(params,"value[param3]=-123; error[param1]=1e-2");
    lk::printFitParametersToStream(params,std::cout);

    lk::modifyFitParameters(params,"error [(1-beta)*bias] = 0.5");
    lk::printFitParametersToStream(params,std::cout);

    try {
        lk::modifyFitParameters(params,"value [param3]=0;error [param3] = -123");
    }
    catch(lk::RuntimeError const &e) {
        // We expect this since error < 0 is not allowed.
    }
    // Check that the parameters were not actually modified.
    lk::printFitParametersToStream(params,std::cout);
    
    try {
        lk::FitParameter badName("a,b",123);
    }
    catch(lk::RuntimeError const &e) {
        // We expect this since commas are not allowed in names.
    }
}
