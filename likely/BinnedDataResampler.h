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
        // of our observations of the specified size. The fixCovariance option requests that
        // the final covariance matrix be corrected for double counting of identical observations.
        // However, this is a relatively slow operation for large datasets, involving a triple
        // matrix product, so you might want to use fixCovariance = false if you don't need
        // an accurate covariance matrix. Chi-square values calculated with fixCovariance = false
        // will be roughly twice as large as the correct values obtained with fixCovariance = true.
        BinnedDataPtr bootstrap(int size, bool fixCovariance = true) const;
	private:
        mutable Random _random;
        std::vector<BinnedDataCPtr> _observations;
        mutable std::vector<int> _shuffle, _counts;
	}; // BinnedDataResampler
	
    inline void BinnedDataResampler::setSeed(int seedValue) { _random.setSeed(seedValue); }
    inline int BinnedDataResampler::getNObservations() const { return _observations.size(); }
    
    // Fills the integer vector provided with a subset of [0,1,...,n-1] of length m=subset.size().
    // The value of seqno determines which subset is selected and values of seqno from zero
    // up to nCm-1 represent a complete enumeration of all possible m-subsets.
    // Returns true for a valid seqno >= 0 (in which case the subset elements will be in
    // increasing order) or false if seqno >= nCm (in which case the input subset is unchanged).
    bool getSubset(int n, unsigned long seqno, std::vector<int> &subset);

} // likely

#endif // LIKELY_BINNED_DATA_RESAMPLER
