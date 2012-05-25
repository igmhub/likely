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
    
    return 0;
}
