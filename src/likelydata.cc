// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstrates and tests the BinnnedGrid and BinnedData classes.

#include "likely/CovarianceMatrix.h"
#include "likely/CovarianceAccumulator.h"
#include "likely/UniformBinning.h"
#include "likely/UniformSampling.h"
#include "likely/NonUniformBinning.h"
#include "likely/BinnedGrid.h"
#include "likely/BinnedData.h"
#include "likely/BinnedDataResampler.h"
#include "likely/Random.h"
#include "likely/RuntimeError.h"

#include <iostream>
#include <cassert>
#include <cmath>

namespace lk = likely;

bool accumulationMessage(lk::CovarianceAccumulatorCPtr accumulator) {
    std::cout << "accumulated " << accumulator->count() << " samples." << std::endl;
    return true;
}

int main(int argc, char **argv) {
    std::vector<double> bins(4);
    bins[0] = 0; bins[1] = 0.25; bins[2] = 0.35; bins[3] = 1;
    lk::AbsBinningCPtr
        axis1(new lk::UniformBinning(0.,1.,3)),
        axis2(new lk::UniformSampling(0.,1.,3)),
        axis3(new lk::NonUniformBinning(bins));

    lk::BinnedGrid grid(axis1,axis2,axis3);
    lk::BinnedData data(grid);
    int nAxes(grid.getNAxes()), nBins(grid.getNBinsTotal());
    std::cout << "naxes = " << nAxes << ", nbins = " << nBins << std::endl;
    std::vector<int> idx(nAxes);
    std::vector<double> centers(nAxes), widths(nAxes);
    for(lk::BinnedGrid::Iterator iter = grid.begin(); iter != grid.end(); ++iter) {
        int index = *iter;
        std::cout << "[" << index << "] =>";
        grid.getBinIndices(index,idx);
        assert(grid.getIndex(idx) == index);
        for(int k = 0; k < nAxes; ++k) std::cout << ' ' << idx[k];
        grid.getBinCenters(index,centers);
        assert(grid.getIndex(centers) == index);
        for(int k = 0; k < nAxes; ++k) std::cout << ' ' << centers[k];
        grid.getBinWidths(index,widths);
        for(int k = 0; k < nAxes; ++k) std::cout << ' ' << widths[k];
        std::cout << std::endl;
        assert(data.hasData(index) == false);
        data.setData(index,index);
    }
    std::cout << "   initial: " << data.getMemoryState() << std::endl;
    data.printToStream(std::cout);

    data.getData(0,true);
    std::cout << "  weighted: " << data.getMemoryState() << std::endl;
    data.printToStream(std::cout);

    data.getData(0,false);  
    std::cout << "unweighted: " << data.getMemoryState() << std::endl;
    data.printToStream(std::cout);

    data.compress();
    std::cout << "compressed: " << data.getMemoryState() << std::endl;
    data.printToStream(std::cout);

    lk::BinnedData copy = data;
    std::cout << "    copied: " << copy.getMemoryState() << std::endl;
    assert(copy.isCongruent(data));
    copy.printToStream(std::cout);

    copy += data;
    std::cout << "     added: " << copy.getMemoryState() << std::endl;
    copy.printToStream(std::cout);
    
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
        lk::BinnedGrid grid(bins);
        lk::BinnedData d1(grid), d2(grid), d3(grid), c1(grid), c2(grid), c3(grid);
        for(int k = 0; k < nbins; ++k) {
            d1.setData(k,1.); c1.setData(k,1.);
            d2.setData(k,2.); c2.setData(k,2.);
            d3.setData(k,3.); c3.setData(k,3.);
        }
        c1.setCovarianceMatrix(C);
        c2.setCovarianceMatrix(C);
        c3.setCovarianceMatrix(C);

        lk::BinnedData d123(grid), c123(grid);
        double wgt(1.5);
        d123.add(d1,2*wgt);
        c123.add(c1,2*wgt);
        d123.printToStream(std::cout);
        c123.printToStream(std::cout);
        d123.add(d2,wgt);
        c123.add(c2,wgt);
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
        lk::BinnedGrid grid(binning);
        lk::BinnedData data(grid);
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
    
    // Test bootstrap estimated covariance for identically distributed observations.
    {
        std::cout << "== Bootstrap covariance test 1:" << std::endl;
        // Create a prototype dataset.
        int nbins(2);
        lk::AbsBinningCPtr binning(new lk::UniformBinning(0.,1.,nbins));
        lk::BinnedGrid grid(binning);
        lk::BinnedData prototype(grid);
        prototype.setData(0,0);
        prototype.setData(1,1);
        // Define a covariance matrix.
        lk::CovarianceMatrixPtr cov(new lk::CovarianceMatrix(nbins));
        (*cov).setCovariance(0,0,1).setCovariance(0,1,-0.5).setCovariance(1,1,2);
        std::cout << "-- ensemble sample covariance:" << std::endl;
        cov->printToStream(std::cout);
        prototype.setCovarianceMatrix(cov);
        // Generate realizations of this covariance matrix
        int nobs(1000);
        lk::RandomPtr random1(new lk::Random()), random2(new lk::Random());
        lk::BinnedDataResampler resamplerMatrix(false,random1), resamplerScalar(true,random2);
        lk::CovarianceAccumulator accumulator(nbins);
        for(int obs = 0; obs < nobs; ++obs) {
            lk::BinnedDataCPtr data = prototype.sample();
            //data->printToStream(std::cout);
            resamplerMatrix.addObservation(data);
            resamplerScalar.addObservation(data);
            accumulator.accumulate(data);
        }
        std::cout << "-- bootstrap test:" << std::endl;
        lk::BinnedDataCPtr bs;
        bool fixCovariance(false);
        bs = resamplerMatrix.bootstrap(0,fixCovariance);
        bs->printToStream(std::cout);
        bs->getCovarianceMatrix()->printToStream(std::cout);
        bs = resamplerScalar.bootstrap(0,fixCovariance);
        bs->printToStream(std::cout);
        bs->getCovarianceMatrix()->printToStream(std::cout);
        // Dump the combined data.
        std::cout << "-- combined data:" << std::endl;
        lk::BinnedDataCPtr combined;
        combined = resamplerMatrix.combined();
        combined->printToStream(std::cout);
        combined->getCovarianceMatrix()->printToStream(std::cout);
        combined = resamplerScalar.combined();
        combined->printToStream(std::cout);
        combined->getCovarianceMatrix()->printToStream(std::cout);
        // Calculate the covariance of the samples actually generated.
        std::cout << "-- calculated sample covariance:" << std::endl;
        accumulator.getCovariance()->printToStream(std::cout);
        // Estimate the covariance of the observations with bootstrap.
        std::cout << "-- bootstrap covariance estimates:" << std::endl;
        lk::CovarianceMatrixPtr bsCov;
        bsCov = resamplerMatrix.estimateCombinedCovariance(10000)->getCovariance();
        bsCov->applyScaleFactor(nobs);
        bsCov->printToStream(std::cout);
        bsCov = resamplerScalar.estimateCombinedCovariance(10000)->getCovariance();
        bsCov->applyScaleFactor(nobs);
        bsCov->printToStream(std::cout);
    }

    // Test bootstrap estimated covariance for non-identically distributed observations.
    {
        std::cout << "== Bootstrap covariance test 2:" << std::endl;
        // Create two prototype datasets with the same binning and contents (MC truth)
        int nbins(2);
        lk::AbsBinningCPtr binning(new lk::UniformBinning(0.,1.,nbins));
        lk::BinnedGrid grid(binning);
        lk::BinnedData prototype1(grid),prototype2(grid);
        prototype1.setData(0,0);
        prototype1.setData(1,+1);
        prototype2.setData(0,0);
        prototype2.setData(1,+1);
        // Define covariance matrices for each subsample.
        lk::CovarianceMatrixPtr cov1(new lk::CovarianceMatrix(nbins)),cov2(new lk::CovarianceMatrix(nbins));
        (*cov1).setCovariance(0,0,1).setCovariance(0,1,-0.9).setCovariance(1,1,2);
        (*cov2).setCovariance(0,0,1).setCovariance(0,1,-0.9).setCovariance(1,1,2);
        cov2->applyScaleFactor(3);
        //(*cov2).setCovariance(0,0,2).setCovariance(0,1,+0.1).setCovariance(1,1,1);
        std::cout << "-- ensemble sample covariances:" << std::endl;
        cov1->printToStream(std::cout);
        cov2->printToStream(std::cout);
        prototype1.setCovarianceMatrix(cov1);
        prototype2.setCovarianceMatrix(cov2);
        std::cout << "weights: " << std::exp(-cov1->getLogDeterminant()/nbins) << ','
            << std::exp(-cov2->getLogDeterminant()/nbins) << std::endl;
        // Define the estimated covariance matrices we will use below, which are obtained by
        // scaling the true covariances.
        double scale(1);
        lk::CovarianceMatrixPtr cov1e(new lk::CovarianceMatrix(*cov1)), cov2e(new lk::CovarianceMatrix(*cov2));
        cov1e->applyScaleFactor(scale);
        cov2e->applyScaleFactor(scale);
        cov1e->printToStream(std::cout);
        cov2e->printToStream(std::cout);
        // Generate realizations of each covariance matrix.
        int n1(400),n2(600);
        lk::RandomPtr random1(new lk::Random()), random2(new lk::Random());
        lk::BinnedDataResampler resamplerMatrix(false,random1), resamplerScalar(true,random2);
        for(int obs = 0; obs < n1; ++obs) {
            lk::BinnedDataPtr data1 = prototype1.sample();
            data1->unweightData();
            data1->setCovarianceMatrix(cov1e);
            resamplerMatrix.addObservation(data1);
            resamplerScalar.addObservation(data1);
        }
        for(int obs = 0; obs < n2; ++obs) {
            lk::BinnedDataPtr data2 = prototype2.sample();
            data2->unweightData();
            data2->setCovarianceMatrix(cov2e);
            resamplerMatrix.addObservation(data2);
            resamplerScalar.addObservation(data2);
        }
        std::cout << "-- bootstrap test:" << std::endl;
        lk::BinnedDataCPtr bs;
        bool fixCovariance(false);
        bs = resamplerMatrix.bootstrap(0,fixCovariance);
        bs->printToStream(std::cout);
        bs->getCovarianceMatrix()->printToStream(std::cout);
        bs = resamplerScalar.bootstrap(0,fixCovariance);
        bs->printToStream(std::cout);
        bs->getCovarianceMatrix()->printToStream(std::cout);
        std::cout << "-- combined data:" << std::endl;
        resamplerMatrix.combined()->printToStream(std::cout);
        resamplerScalar.combined()->printToStream(std::cout);
        // Calculate the ensemble covariance using Cinv12 = n1*Cinv1 + n2*Cinv2
        std::cout << "-- calculated sample covariance:" << std::endl;
        lk::CovarianceMatrix cov12(*cov1);
        cov12.applyScaleFactor(1./n1);
        cov12.addInverse(*cov2,n2);
        cov12.printToStream(std::cout);
        // Estimate the covariance of the observations with bootstrap.
        int ntrials(1000);
        std::cout << "-- bootstrap covariance estimates:" << std::endl;
        lk::BinnedDataResampler::AccumulationCallback callback(accumulationMessage);
        lk::CovarianceAccumulatorPtr accum;
        lk::CovarianceMatrixPtr bsCov;
        accum = resamplerMatrix.estimateCombinedCovariance(ntrials,callback,5000);
        accum->dump(std::cout);
        bsCov = accum->getCovariance();
        bsCov->printToStream(std::cout);
        accum = resamplerScalar.estimateCombinedCovariance(ntrials,callback,5000);
        accum->dump(std::cout);
        bsCov = accum->getCovariance();
        bsCov->printToStream(std::cout);
    }
    
    return 0;
}
