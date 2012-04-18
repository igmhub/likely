// Created 18-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the CovarianceMatrix class.

#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

int main(int argc, char **argv) {
    int size(3);
    lk::CovarianceMatrix cov(size);
    for(int k = 0; k < size; ++k) {
        cov.setCovariance(k,k,1);
    }
}
