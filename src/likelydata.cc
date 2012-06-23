// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the CovarianceMatrix class.

#include "likely/likely.h"

#include <iostream>
#include <cassert>
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
    std::vector<double> bins(4);
    bins[0] = 0; bins[1] = 0.25; bins[2] = 0.35; bins[3] = 1;
    lk::AbsBinningCPtr
        axis1(new lk::UniformBinning(0.,1.,3)),
        axis2(new lk::UniformSampling(0.,1.,3)),
        axis3(new lk::NonUniformBinning(bins));

    lk::BinnedData data(axis1,axis2,axis3);
    int nAxes(data.getNAxes()), nBins(data.getNBinsTotal());
    std::cout << "naxes = " << nAxes << ", nbins = " << nBins << std::endl;
    std::vector<int> idx(nAxes);
    std::vector<double> centers(nAxes), widths(nAxes);
    for(int index = 0; index < data.getNBinsTotal(); ++index) {
        std::cout << "[" << index << "] =>";
        data.getBinIndices(index,idx);
        assert(data.getIndex(idx) == index);
        for(int k = 0; k < nAxes; ++k) std::cout << ' ' << idx[k];
        data.getBinCenters(index,centers);
        assert(data.getIndex(centers) == index);
        for(int k = 0; k < nAxes; ++k) std::cout << ' ' << centers[k];
        data.getBinWidths(index,widths);
        for(int k = 0; k < nAxes; ++k) std::cout << ' ' << widths[k];
        std::cout << std::endl;
        assert(data.hasData(index) == false);
        data.setData(index,index);
    }
    std::cout << "size = " << data.getMemoryUsage() << std::endl;
    data.compress();
    std::cout << "compressed size = " << data.getMemoryUsage() << std::endl;
    lk::BinnedData copy = data;
    std::cout << "copy size = " << copy.getMemoryUsage() << std::endl;
    assert(copy.isCongruent(data));
    copy += data;
    
    lk::BinnedData::IndexIterator ptr = data.begin();
    for(lk::BinnedData::IndexIterator iter = data.begin(); iter != data.end(); ++iter) {
        std::cout << "[" << *iter << "] = " << data.getData(*iter) << std::endl;
    }
    
    data.finalize();
    try {
        data.setCovariance(0,0,0);
    }
    catch(likely::RuntimeError const &e) {
        std::cout << e.what() << std::endl;
    }
    
    // Test subset combinatorics.
    {
        int n(5),m(2),seqno(0);
        std::vector<int> subset(m);
        while(lk::getSubset(n,seqno,subset)) {
            std::cout << '[' << seqno << "] ";
            for(int k = 0; k < m; ++k) std::cout << subset[k] << ' ';
            std::cout << std::endl;
            seqno++;
        }
    }
    
    // Test unweighted vs weighted data combinations.
    {
        int nbins(3);
        lk::CovarianceMatrixPtr C(lk::createDiagonalCovariance(nbins,1));
        lk::AbsBinningCPtr bins(new lk::UniformBinning(0.,1.,nbins));
        lk::BinnedData d1(bins), d2(bins), d3(bins), c1(bins), c2(bins), c3(bins);
        for(int k = 0; k < nbins; ++k) {
            d1.setData(k,1.); c1.setData(k,1.);
            d2.setData(k,2.); c2.setData(k,2.);
            d3.setData(k,3.); c3.setData(k,3.);
        }
        c1.setCovarianceMatrix(C);
        c2.setCovarianceMatrix(C);
        c3.setCovarianceMatrix(C);

        lk::BinnedData d123(bins), c123(bins);
        d123 += d1;
        c123 += c1;
        d123.printToStream(std::cout);
        c123.printToStream(std::cout);
        d123 += d2;
        c123 += c2;
        d123.printToStream(std::cout);
        c123.printToStream(std::cout);
    }
    
    // Test decorrelated errors
    {
        lk::RandomPtr random(new lk::Random());
        random->setSeed(12345);
        int nbins(5);
        // Generate a random covariance matrix
        lk::CovarianceMatrixPtr C(lk::generateRandomCovariance(nbins,1,random));
        // Initialize an empty dataset.
        lk::AbsBinningCPtr binning(new lk::UniformBinning(0,1,nbins));
        lk::BinnedData data(binning);
        // Generate random prediction and data vectors.
        std::vector<double> pred,noise;
        C->sample(noise,random);
        for(int index = 0; index < nbins; ++index) {
            double truth = random->getUniform();
            pred.push_back(truth);
            data.setData(index, truth + noise[index]);
        }
        // Comment out this line to test without a covariance matrix.
        data.setCovarianceMatrix(C);
        // Calculate the chi-square with the full covariance.
        double chi2 = data.chiSquare(pred);
        // Calculate with decorrelated errors.
        std::vector<double> dwgt;
        data.getDecorrelatedWeights(pred,dwgt);
        double chi2d(0);
        for(int index = 0; index < nbins; ++index) {
            chi2d += noise[index]*noise[index]*dwgt[index];
        }
        std::cout << "chi2 = " << chi2 << ", chi2d = " << chi2d << std::endl;
    }
    
    return 0;
}
