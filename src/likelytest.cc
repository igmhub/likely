// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/likely.h"

#include <iostream>

namespace test = likely::test;

int main(int argc, char *argv) {
    test::TestLikelihood likelyfn(3,1,-0.75,0.25);
    std::vector<double> params(3,0);
    params[1] = 1;
    std::cout << likelyfn(params) << std::endl;
}