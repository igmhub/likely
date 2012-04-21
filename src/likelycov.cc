// Created 18-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the CovarianceMatrix class.

#include "likely/likely.h"

#include <iostream>
#include <sys/resource.h>

namespace lk = likely;

// Returns the number of elapsed microseconds from before to after.
double elapsed(struct timeval const &before, struct timeval const &after) {
    return (after.tv_sec - before.tv_sec)*1e6 + (after.tv_usec - before.tv_usec);
}
double elapsed(struct rusage const &before, struct rusage const &after) {
    return elapsed(before.ru_utime,after.ru_utime) + elapsed(before.ru_stime,after.ru_stime);
}

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
    cov.printToStream(std::cout);

    int nsample(1000000);
    struct rusage t1,t2,t3;
    
    {
        getrusage(RUSAGE_SELF,&t1);
        boost::shared_array<double> residuals = cov.sample(nsample);
        getrusage(RUSAGE_SELF,&t2);
        lk::CovarianceAccumulator accum(size);
        for(int row = 0; row < nsample; ++row) {
            accum.accumulate(&residuals[size*row]);
        }
        lk::CovarianceMatrixCPtr cptr = accum.getCovariance();
        getrusage(RUSAGE_SELF,&t3);    
        cptr->printToStream(std::cout);

        int nelem = nsample*size;
        double t12 = 1e3*elapsed(t1,t2)/nelem, t23 = 1e3*elapsed(t2,t3)/nelem;
        std::cout << "sample = " << t12 << " ns/elem, accumulate = " << t23
            << " ns/elem, accumulate/sample = " << t23/t12 << std::endl;
        std::cout << cov.getMemoryState() << std::endl;
    }

    // Validate single samples
    {
        std::vector<double> delta(size);
        lk::CovarianceAccumulator accum(size);
        for(int k = 0; k < nsample; ++k) {
            cov.sample(delta);
            accum.accumulate(delta);
        }
        accum.getCovariance()->printToStream(std::cout);
    }

    // Benchmark single samples
    int ntrial = 10000;
    {
        std::vector<double> delta(size);
        boost::shared_array<double> delta2;
        for(nsample = 1; nsample <= 50; ++nsample) {
            getrusage(RUSAGE_SELF,&t1);
            for(int trial = 0; trial < ntrial; ++trial) {
                for(int k = 0; k < nsample; ++k) {
                    cov.sample(delta);
                }
            }
            getrusage(RUSAGE_SELF,&t2);
            for(int trial = 0; trial < ntrial; ++trial) {
                delta2 = cov.sample(nsample);
            }
            getrusage(RUSAGE_SELF,&t3);
            int ntot(nsample*ntrial);
            std::cout << "nsample = " << nsample << ": " << 1e3*elapsed(t1,t2)/ntot << ", "
                << 1e3*elapsed(t2,t3)/ntot << std::endl;
        }
    }
}
