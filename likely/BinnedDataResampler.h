// Created 29-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_DATA_RESAMPLER
#define LIKELY_BINNED_DATA_RESAMPLER

#include "likely/types.h"
#include "likely/Random.h"

#include <vector>

namespace likely {
	class BinnedDataResampler {
	// Collects and resamples a set of congruent BinnedData observations using jackknife
	// and bootstrap techniques.
	public:
	    // Creates a new resampler using the specified random seed.
		BinnedDataResampler(int randomSeed);
		virtual ~BinnedDataResampler();
		// Sets the random seed value to use for subsequent random resampling.
        void setSeed(int seedValue);
		// Adds a new observation for resampling. Throws a RuntimeError if this observation
		// has already been added or is not congruent with existing observations.
        void addObservation(BinnedDataCPtr observation);
        // Returns the number of observations available for resampling.
        int getNObservations() const;
        // Returns a shared pointer to a new BinnedData that represents all observations combined.
        // Each call to this method builds a new combined dataset so save the result unless
        // you actually want many independent copies.
        BinnedDataPtr combined() const;
        // Returns a shared pointer to a new BinnedData that represents a jackknife resampling
        // of our observations of the specified size.
        BinnedDataPtr jackknife(int size) const;
        // Returns a shared pointer to a new BinnedData that represents a bootstrap resampling
        // of our observations of the specified size.
        BinnedDataPtr bootstrap(int size, bool accurateWeights = true) const;
	private:
        mutable Random _random;
        std::vector<BinnedDataCPtr> _observations;
        mutable std::vector<int> _shuffle, _counts;
	}; // BinnedDataResampler
	
    inline void BinnedDataResampler::setSeed(int seedValue) { _random.setSeed(seedValue); }
    inline int BinnedDataResampler::getNObservations() const { return _observations.size(); }

} // likely

#endif // LIKELY_BINNED_DATA_RESAMPLER
