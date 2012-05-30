// Created 29-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_DATA_RESAMPLER
#define LIKELY_BINNED_DATA_RESAMPLER

#include "likely/types.h"
#include "likely/Random.h"

#include <vector>

namespace likely {
	class BinnedDataResampler {
	// Collects and resamples a set of congruent BinnedData observations using jackknife
	// and bootstrap techniques. Resampler works with subclasses X of BinnedData as long
	// as they implement their own X::clone() method. Subclasses will need to use
	// boost::dynamic_pointer_cast<...>, e.g.
	//
	//  typedef boost::shared_ptr<const X> XCPtr;
	//  typedef boost::shared_ptr<X> XPtr;
	//  XCPtr x1,x2;
	//  likely::BinnedResampler resampler;
	//  resampler.add(boost::dynamic_pointer_cast<likely::BinnedData>(x1));
	//  resampler.add(boost::dynamic_pointer_cast<likely::BinnedData>(x2));
	//  XPtr bs = boost::dynamic_pointer_cast<X>(resampler.bootstrap(100));
	//
	public:
	    // Creates a new resampler using the specified random seed.
		BinnedDataResampler(int randomSeed);
		virtual ~BinnedDataResampler();
		// Sets the random seed value to use for subsequent random resampling.
        void setSeed(int seedValue);
		// Adds a new observation for resampling. Throws a RuntimeError if this observation
		// has already been added or is not congruent with existing observations. Calls to
		// this method can be interspersed with calls to resampling methods below.
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
