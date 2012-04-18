// Created 18-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the CovarianceMatrix class.

#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

int main(int argc, char **argv) {
    int size(3);
    lk::CovarianceMatrix cov(size);
    for(int k = 0; k < size; ++k) {
        cov.setCovariance(k,k,k+1);
    }
    cov.setCovariance(0,1,0.1);
    cov.setCovariance(1,2,-0.2);
    cov.setInverseCovariance(2,2,0.3);
    for(int row = 0; row < size; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << row << ',' << col << " = " << cov.getCovariance(row,col) << std::endl;
        }
    }    
    for(int row = 0; row < size; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << row << ',' << col << " = " << cov.getInverseCovariance(row,col) << std::endl;
        }
    }    
    for(int row = 0; row < size; ++row) {
        for(int col = 0; col < size; ++col) {
            std::cout << row << ',' << col << " = " << cov.getCovariance(row,col) << std::endl;
        }
    }    
}
