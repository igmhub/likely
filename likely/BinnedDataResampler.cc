// Created 29-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/BinnedDataResampler.h"
#include "likely/RuntimeError.h"
#include "likely/BinnedData.h"
#include "likely/Random.h"
#include "likely/CovarianceMatrix.h"
#include "likely/CovarianceAccumulator.h"

#include "boost/math/special_functions/binomial.hpp"

#include <algorithm>
#include <iostream>

namespace local = likely;

local::BinnedDataResampler::BinnedDataResampler(bool useScalarWeights, RandomPtr random)
: _useScalarWeights(useScalarWeights), _random(random), _combinedScalarWeight(0)
{
    if(!_random) _random = Random::instance();
}

local::BinnedDataResampler::~BinnedDataResampler() { }

int local::BinnedDataResampler::addObservation(BinnedDataCPtr observation, int reuseCovIndex) {
    std::cout << "  add-in: " << observation->getMemoryState() << std::endl;
    // Check that this new observation is congruent with what we have so far. Ignore covariance
    // in the congruence test if we will be reusing a previous covariance matrix.
    if(getNObservations() > 0) {
        if(!_combined->isCongruent(*observation,false,reuseCovIndex >= 0)) {
            throw RuntimeError("BinnedDataResampler::addObservation: new observation is incongruent.");
        }
    }
    else {
        bool binningOnly(true);
        _combined.reset(observation->clone(binningOnly));
    }
    int newIndex = _observations.size();
    // Make a copy of this observation that we will keep.
    BinnedDataPtr copy(observation->clone());
    // Reuse the covariance of a previously added dataset?
    if(reuseCovIndex >= 0) {
        if(reuseCovIndex >= newIndex) {
            throw RuntimeError("BinnedDataResampler::addObservation: invalid reuseCovIndex.");
        }
        BinnedDataCPtr reuseData = _observations[reuseCovIndex];
        // This operation will need an input covariance if the dataset was previously weighted.
        copy->setWeighted(false);
        if(_useScalarWeights) {
            // We already dropped the previously added dataset's covariance, so we need to
            // temporarily reconstruct something here that we can add to our combined dataset.
            double weight = reuseData->getScalarWeight();
            CovarianceMatrixPtr reconstructed(new CovarianceMatrix(*_combined->getCovarianceMatrix()));
            reconstructed->applyScaleFactor(_combinedScalarWeight/weight);
            copy->setCovarianceMatrix(reconstructed);
        }
        else {
            if(!reuseData->hasCovariance()) {
                throw RuntimeError("BinnedDataResampler::addObservation: no covariance to reuse.");
            }
            // Reuse the previously added dataset's covariance matrix.
            copy->shareCovarianceMatrix(*reuseData);
        }
    }
    else {
        // Make a copy of the covariance matrix so changes to the input observation's
        // covariance will not affect us (this is harmless if the copy doesn't actually
        // have any covariance) and so that the operations below do not alter the input
        // observation at all.
        copy->cloneCovariance();
    }
    std::cout << "add-tmp1: " << copy->getMemoryState() << std::endl;
    // Add this copy to our combined dataset.
    *_combined += *copy;
    std::cout << "add-tmp2: " << copy->getMemoryState() << std::endl;
    std::cout << "add-comb: " << _combined->getMemoryState() << std::endl;
    if(_useScalarWeights) {
        // Replace Cinv with the scalar weight |Cinv|^(1/n)
        double weight = copy->getScalarWeight();
        std::cout << "add-tmp3: " << copy->getMemoryState() << std::endl;
        copy->dropCovariance(weight);
        _combinedScalarWeight += weight;
    }
    // Compress the copy before we save it (none of our resampling methods should uncompress it)
    copy->compress();
    // Remember this (copied) observation
    std::cout << " add-out: " << observation->getMemoryState() << std::endl;
    std::cout << "add-copy: " << copy->getMemoryState() << std::endl;
    _observations.push_back(copy);
    return newIndex;
}

local::BinnedDataPtr local::BinnedDataResampler::combined() const {
    BinnedDataPtr all(_combined->clone());
    return all;
}

bool local::getSubset(int n, unsigned long seqno, std::vector<int> &subset) {
    if(n <= 0) {
        throw RuntimeError("BinnedDataResampler::getSubset: invalid n.");
    }
    int m = subset.size();
    if(0 == m) {
        throw RuntimeError("BinnedDataResampler::getSubset: subset is empty.");
    }
    if(seqno < 0) {
        throw RuntimeError("BinnedDataResampler::getSubset: expected seqno >= 0.");
    }
    // See http://en.wikipedia.org/wiki/Combinatorial_number_system
    for(int k = m; k > 0; --k) {
        int next = k-1;
        double last,nCk = 0;
        while(nCk <= seqno) {
            last = nCk;
            nCk = boost::math::binomial_coefficient<double>(++next,k);
            if(next > n) return false;
        }
        subset[k-1] = next-1;
        seqno -= last;
    }
    return true;
}

void local::BinnedDataResampler::_addCovariance(BinnedDataPtr sample) const {
    if(!_useScalarWeights) return;
    sample->setWeighted(false);
    double weight = sample->getScalarWeight();
    CovarianceMatrixPtr cov(new CovarianceMatrix(*_combined->getCovarianceMatrix()));
    cov->applyScaleFactor(_combinedScalarWeight/weight);
    sample->setCovarianceMatrix(cov);
}

local::BinnedDataCPtr local::BinnedDataResampler::getObservation(int index) const {
    if(index < 0 || index >= getNObservations()) {
        throw RuntimeError("BinnedDataResampler::getObservation: index of our range.");
    }
    return _observations[index];
}

local::BinnedDataPtr local::BinnedDataResampler::getObservationCopy(int index, bool addCovariance) const {
    if(index < 0 || index >= getNObservations()) {
        throw RuntimeError("BinnedDataResampler::getObservation: index of our range.");
    }
    BinnedDataPtr copy(_observations[index]->clone());
    if(addCovariance) _addCovariance(copy);
    return copy;
}

local::BinnedDataPtr local::BinnedDataResampler::jackknife(int ndrop, unsigned long seqno,
bool addCovariance) const {
    int nobs(_observations.size());
    if(ndrop < 0 || ndrop >= nobs) {
        throw RuntimeError("BinnedDataResampler::jackknife: invalid ndrop.");
    }
    // Fill our _subset vector with the subset indices corresponding to seqno.
    int nkeep = nobs - ndrop;
    _subset.resize(nkeep);
    if(!getSubset(nobs,seqno,_subset)) return BinnedDataPtr();
    // Create an empty dataset with the right axis binning.
    BinnedDataPtr resample(_observations[0]->clone(true));
    // Add each observation from the generated sample.
    for(int obsIndex = 0; obsIndex < nkeep; ++obsIndex) {
        *resample += *_observations[_subset[obsIndex]];
    }
    if(addCovariance) _addCovariance(resample);
    return resample;
}

local::BinnedDataPtr local::BinnedDataResampler::bootstrap(int size, bool fixCovariance,
bool addCovariance) const {
    if(size < 0) {
        throw RuntimeError("BinnedDataResampler::bootstrap: invalid size.");
    }
    if(0 == size) size = getNObservations();
    if(0 == getNObservations()) return BinnedDataPtr();
    // Do we need to (re)initialize our counts vector?
    if(_counts.size() != _observations.size()) {
        _counts.resize(_observations.size(),0);
    }
    // Generate a random sample with replacement.
    _random->sampleWithReplacement(_counts,size);
    // Create an empty dataset with the right axis binning.
    BinnedDataPtr resample(_observations[0]->clone(true));
    // We cannot fix a non-existent covariance.
    if(!_observations[0]->hasCovariance() || _useScalarWeights) fixCovariance = false;
    // Initialize matrix needed to fix final covariance.
    likely::CovarianceMatrixPtr D;
    int nbins = _observations[0]->getNBinsWithData();
    if(fixCovariance) D.reset(new likely::CovarianceMatrix(nbins));
    // Loop over observations, adding each one the appropriate number of times.
    bool duplicatesFound(false);
    for(int obsIndex = 0; obsIndex < _observations.size(); ++obsIndex) {
        int count(_counts[obsIndex]);
        if(0 == count) continue;
        if(count > 1) duplicatesFound = true;
        BinnedDataCPtr observation = _observations[obsIndex];
        resample->add(*observation,count);
        if(fixCovariance) D->addInverse(*(observation->getCovarianceMatrix()),count*count);
    }
    // We can skip this relatively expensive operation if all counts are 0,1.
    if(duplicatesFound && fixCovariance) resample->transformCovariance(D);
    if(addCovariance) _addCovariance(resample);
    return resample;
}

local::CovarianceMatrixPtr
local::BinnedDataResampler::estimateCombinedCovariance(int nSamples, int messageInterval) const {
    if(nSamples <= 0) {
        throw RuntimeError("BinnedDataResampler::estimateCombinedCovariance: expected nSamples > 0.");
    }
    if(0 == getNObservations()) return CovarianceMatrixPtr();
    CovarianceAccumulator accumulator(_observations[0]->getNBinsWithData());
    bool fixCovariance(false),addCovariance(false);
    for(int sample = 0; sample < nSamples; ++sample) {
        if(messageInterval > 0 && sample > 0 && sample % messageInterval == 0) {
            std::cout << "generated " << sample << " bootstrap samples." << std::endl;
        }
        BinnedDataPtr data = bootstrap(0,fixCovariance,addCovariance);
        accumulator.accumulate(data);
    }
    try {
        return accumulator.getCovariance();
    }
    catch(RuntimeError const &e) {
        throw RuntimeError("BinnedDataResampler::estimateCombinedCovariance: failed - try more samples?");
    }
}
