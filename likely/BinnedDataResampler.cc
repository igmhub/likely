// Created 29-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/BinnedDataResampler.h"
#include "likely/RuntimeError.h"
#include "likely/BinnedData.h"
#include "likely/CovarianceMatrix.h"

#include <algorithm>

namespace local = likely;

local::BinnedDataResampler::BinnedDataResampler(int randomSeed)
{
    setSeed(randomSeed);
}

local::BinnedDataResampler::~BinnedDataResampler() { }

void local::BinnedDataResampler::addObservation(BinnedDataCPtr observation) {
    if(std::find(_observations.begin(), _observations.end(), observation) != _observations.end()) {
        throw RuntimeError("BinnedDataResampler::addObservation: cannot add duplicate.");
    }
    if(getNObservations() > 0 && !_observations[0]->isCongruent(*observation)) {
        throw RuntimeError("BinnedDataResampler::addObservation: new observation is incongruent.");
    }
    _observations.push_back(observation);
}

local::BinnedDataPtr local::BinnedDataResampler::combined() const {
    BinnedDataPtr all;
    int size(_observations.size());
    // Return an unassigned shared pointer if we don't have any observations yet.
    if(0 == size) return all;
    // Create an empty dataset with the right axis binning.
    all.reset(_observations[0]->clone(true));
    for(int obsIndex = 0; obsIndex < size; ++obsIndex) {
        *all += *_observations[obsIndex];
    }
    return all;
}

// This should not be random. Instead: jackknife(int ndrop, int &seqno) updates seqno to
// iterate over all ways to drop ndrop observations and returns an empty ptr after the last one.
local::BinnedDataPtr local::BinnedDataResampler::jackknife(int size) const {
    if(size <= 0 || size > _observations.size()) {
        throw RuntimeError("BinnedDataResampler::jackknife: invalid size.");
    }
    // Do we need to (re)initialize our shuffle vector?
    if(_shuffle.size() != _observations.size()) {
        _shuffle.resize(0);
        _shuffle.reserve(_observations.size());
        for(int obsIndex = 0; obsIndex < _observations.size(); ++obsIndex) {
            _shuffle.push_back(obsIndex);
        }
    }
    // Do a partial shuffling so that the first size elements index our jackknife sample.
    _random.partialShuffle(_shuffle,size);
    // Create an empty dataset with the right axis binning.
    BinnedDataPtr resample(_observations[0]->clone(true));
    // Add each observation from the generated sample.
    for(int obsIndex = 0; obsIndex < size; ++obsIndex) {
        *resample += *_observations[_shuffle[obsIndex]];
    }
    return resample;
}

local::BinnedDataPtr local::BinnedDataResampler::bootstrap(int size, bool fixCovariance) const {
    if(size <= 0) {
        throw RuntimeError("BinnedDataResampler::bootstrap: invalid size.");
    }
    if(0 == getNObservations()) return BinnedDataPtr();
    // Do we need to (re)initialize our counts vector?
    if(_counts.size() != _observations.size()) {
        _counts.resize(_observations.size(),0);
    }
    // Generate a random sample with replacement.
    _random.sampleWithReplacement(_counts,size);
    // Create an empty dataset with the right axis binning.
    BinnedDataPtr resample(_observations[0]->clone(true));
    // Initialize matrix needed to fix final covariance.
    likely::CovarianceMatrixPtr D;
    int nbins = _observations[0]->getNBinsWithData();
    if(fixCovariance) D.reset(new likely::CovarianceMatrix(nbins));
    // Loop over observations, adding each one the appropriate number of times.
    bool duplicatesFound(false);
    for(int obsIndex = 0; obsIndex < _observations.size(); ++obsIndex) {
        int count(_counts[obsIndex]);
        if(0 == count) continue;
        BinnedDataCPtr observation = _observations[obsIndex];
        resample->add(*observation,count);
        if(fixCovariance) D->addInverse(*(observation->getCovarianceMatrix()),count*count);
    }
    if(fixCovariance) resample->transformCovariance(D);
    return resample;
}
