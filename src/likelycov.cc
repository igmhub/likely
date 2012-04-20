// Created 18-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the CovarianceMatrix class.

#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

int main(int argc, char **argv) {
    int size(3);
    lk::CovarianceMatrix cov(size);
    std::cout << cov.getMemoryState() << std::endl;
    for(int k = 0; k < size; ++k) {
        cov.setCovariance(k,k,k+1);
    }
    cov.setCovariance(0,1,0.1);
    cov.setCovariance(1,2,-0.2);
    std::cout << cov.getMemoryState() << std::endl;
    cov.setInverseCovariance(2,2,0.3);
    std::cout << cov.getMemoryState() << std::endl;
    for(int row = 0; row < size; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << row << ',' << col << " = " << cov.getCovariance(row,col) << std::endl;
        }
    }
    std::cout << cov.getMemoryState() << std::endl;
    for(int row = 0; row < size; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << row << ',' << col << " = " << cov.getInverseCovariance(row,col) << std::endl;
        }
    }
    std::cout << cov.getMemoryState() << std::endl;
    cov.compress();
    std::cout << cov.getMemoryState() << std::endl;
    for(int row = 0; row < size; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << row << ',' << col << " = " << cov.getCovariance(row,col) << std::endl;
        }
    }
    std::cout << cov.getMemoryState() << std::endl;
    cov.setCovariance(0,0,1.5);
    std::cout << cov.getMemoryState() << std::endl;
    cov.getInverseCovariance(0,0);
    std::cout << cov.getMemoryState() << std::endl;
    cov.setCovariance(0,0,2.5);
    std::cout << cov.getMemoryState() << std::endl;
    cov.compress();
    std::cout << cov.getMemoryState() << std::endl;
    
    std::vector<double> vec(3);
    vec[0] = 1; vec[1] = 2; vec[2] = 3;
    cov.multiplyByInverseCovariance(vec);
    std::cout << cov.getMemoryState() << std::endl;

    cov.chiSquare(vec);
    std::cout << cov.getMemoryState() << std::endl;

    cov.compress();
    std::cout << cov.getMemoryState() << std::endl;

    // Test random sampling...
    cov.setCovariance(0,0,1);
    cov.setCovariance(1,1,1);
    cov.setCovariance(2,2,1);
    cov.setCovariance(0,1,0.5);
    cov.setCovariance(0,2,0.5);
    cov.setCovariance(1,2,0.5);
    std::cout << cov.getMemoryState() << std::endl;

    int nsample(10);
    boost::shared_array<double> residuals = cov.sample(nsample);
    for(int row = 0; row < nsample; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << ' ' << residuals[size*row + col];
        }
        std::cout << std::endl;
    }
    std::cout << cov.getMemoryState() << std::endl;
}
