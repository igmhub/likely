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
: _useScalarWeights(useScalarWeights), _random(random), _combinedWeight(0)
{
    if(!_random) _random = Random::instance();
}

local::BinnedDataResampler::~BinnedDataResampler() { }

void local::BinnedDataResampler::addObservation(BinnedDataCPtr observation) {
    BinnedDataCPtr keeper = observation;
    if(_useScalarWeights) {
        // Keep a copy of the observation with its covariance replaced by a scalar weight
        // equal to |Cinv|^(1/n) where n is the covariance matrix size.
        BinnedDataPtr copy(observation->clone());
        double weight = std::exp(-copy->getCovarianceMatrix()->getLogDeterminant()/copy->getNBinsWithData());
        copy->dropCovariance(weight);
        _combinedWeight += weight;
        keeper = copy;
    }
    if(getNObservations() > 0) {
        if(!_observations[0]->isCongruent(*keeper)) {
            throw RuntimeError("BinnedDataResampler::addObservation: new observation is incongruent.");
        }
        if(_useScalarWeights) {
            _combinedCovariance->addInverse(*observation->getCovarianceMatrix());
        }
    }
    else {
        if(_useScalarWeights) {
            _combinedCovariance.reset(new CovarianceMatrix(*observation->getCovarianceMatrix()));
        }
    }
    _observations.push_back(keeper);
}

local::BinnedDataPtr local::BinnedDataResampler::combined() const {
    // The reason we don't build and cache the combined dataset in addObservation is
    // that we have no way to guarantee that the individual observations won't be
    // changed after they are added. Instead, we build the combination each time we
    // are called and leave it up to the user to cache the result when they know that
    // nothing has changed.
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

local::BinnedDataCPtr local::BinnedDataResampler::getObservation(int index) const {
    if(index < 0 || index >= getNObservations()) {
        throw RuntimeError("BinnedDataResampler::getObservation: index of our range.");
    }
    return _observations[index];
}

local::BinnedDataPtr local::BinnedDataResampler::getObservationCopy(int index) const {
    if(index < 0 || index >= getNObservations()) {
        throw RuntimeError("BinnedDataResampler::getObservation: index of our range.");
    }
    BinnedDataPtr copy(_observations[index]->clone());
    return copy;
}

local::BinnedDataPtr local::BinnedDataResampler::jackknife(int ndrop, unsigned long seqno) const {
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
    return resample;
}

local::BinnedDataPtr local::BinnedDataResampler::bootstrap(int size, bool fixCovariance,
bool scalarWeights) const {
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
    if(!_observations[0]->hasCovariance() || scalarWeights) fixCovariance = false;
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
        if(scalarWeights) {
            double weight = count*std::exp(-observation->getCovarianceMatrix()->getLogDeterminant()/nbins);
            observation->setWeighted(false);
            BinnedDataPtr copy(observation->clone());
            copy->dropCovariance();
            resample->add(*copy,weight);
        }
        else {
            resample->add(*observation,count);
        }
        if(fixCovariance) D->addInverse(*(observation->getCovarianceMatrix()),count*count);
    }
    // We can skip this relatively expensive operation if all counts are 0,1.
    if(duplicatesFound && fixCovariance) resample->transformCovariance(D);
    return resample;
}

local::CovarianceMatrixPtr
local::BinnedDataResampler::estimateCombinedCovariance(int nSamples, int messageInterval,
bool scalarWeights) const {
    if(nSamples <= 0) {
        throw RuntimeError("BinnedDataResampler::estimateCombinedCovariance: expected nSamples > 0.");
    }
    if(0 == getNObservations()) return CovarianceMatrixPtr();
    CovarianceAccumulator accumulator(_observations[0]->getNBinsWithData());
    bool fixCovariance(false);
    for(int sample = 0; sample < nSamples; ++sample) {
        if(messageInterval > 0 && sample > 0 && sample % messageInterval == 0) {
            std::cout << "generated " << sample << " bootstrap samples." << std::endl;
        }
        BinnedDataPtr data = bootstrap(0,fixCovariance,scalarWeights);
        accumulator.accumulate(data);
    }
    try {
        return accumulator.getCovariance();
    }
    catch(RuntimeError const &e) {
        throw RuntimeError("BinnedDataResampler::estimateCombinedCovariance: failed - try more samples?");
    }
}
