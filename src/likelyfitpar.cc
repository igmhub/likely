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
    lk::printFitParametersToStream(params,std::cout);
}
