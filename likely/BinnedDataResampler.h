// Created 29-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_DATA_RESAMPLER
#define LIKELY_BINNED_DATA_RESAMPLER

#include "likely/types.h"

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
	    // Creates a new resampler using the random generator provided, or else the default
	    // Random::instance().
		BinnedDataResampler(RandomPtr random = RandomPtr());
		virtual ~BinnedDataResampler();
		// Adds a new observation for resampling. Throws a RuntimeError if this observation
		// has already been added or is not congruent with existing observations. Calls to
		// this method can be interspersed with calls to resampling methods below.
        void addObservation(BinnedDataCPtr observation);
        // Returns the number of observations available for resampling.
        int getNObservations() const;
        // Returns a shared pointer to the specified (readonly) observation.
        BinnedDataCPtr getObservation(int index) const;
        // Returns a shared pointer to a (modifiable) copy of the specified observation.
        BinnedDataPtr getObservationCopy(int index) const;
        // Returns a shared pointer to a new BinnedData that combines all observations added
        // so far. Each call to this method builds a new combined dataset so save the result
        // unless you actually want many independent copies.
        BinnedDataPtr combined() const;
        // Returns a shared pointer to a new BinnedData that represents a jackknife resampling
        // of our observations with the specified number of observations dropped. The specific
        // resampling is determined by the value of seqno. Use the following to generate the
        // full set of jackknife samples:
        //
	    //  likely::BinnedResampler resampler;
	    //  BinnedDataPtr sample;
        //  unsigned long seqno(0);
        //  while(sample = resampler.jackknife(ndrop,seqno++)) {
        //    ...
        //  }
        // Note that the number of jackknife samples generated this way gets large quickly
        // as ndrop increases. There is no requirement that seqno increase by one for successive
        // calls, so jackknifing can easily be parallelized in various ways.
        BinnedDataPtr jackknife(int ndrop, unsigned long seqno) const;
        // Returns a shared pointer to a new BinnedData that represents a bootstrap resampling
        // of our observations of the specified size, which defaults to the number of observations
        // when zero. The fixCovariance option requests that the final covariance matrix be corrected
        // for double counting of identical observations. However, this is a relatively slow operation
        // for large datasets, involving a triple matrix product, so you might want to use
        // fixCovariance = false if you don't need an accurate covariance matrix. Chi-square values
        // calculated with fixCovariance = false will be roughly twice as large as the correct values
        // obtained with fixCovariance = true. With scalarWeights = true, observations are combined
        // using scalar weights det(Cinv)^(1/n) instead of matrix weights Cinv. This method is faster
        // but only valid when all observations are identically distributed up to a scale factor, i.e.
        // each observation's true covariance is proportional to some global covariance. The data
        // returned with this method will have no covariance matrix associated with it.
        BinnedDataPtr bootstrap(int size = 0, bool fixCovariance = true, bool scalarWeights = false) const;
        // Returns a shared pointer to a new CovarianceMatrix that estimates the covariance of
        // our combined observations using the specified number of bootstrap samples with a
        // CovarianceAccumulator. Throws a RuntimeError if the estimated covariance is not positive
        // definite, which can usually be fixed with more samples. See the description of the bootstrap()
        // method for details on the scalarWeights option.
        CovarianceMatrixPtr estimateCombinedCovariance(int nSamples, int messageInterval = 0,
            bool scalarWeights = false) const;
	private:
        mutable RandomPtr _random;
        std::vector<BinnedDataCPtr> _observations;
        mutable std::vector<int> _subset, _counts;
	}; // BinnedDataResampler
	
    inline int BinnedDataResampler::getNObservations() const { return _observations.size(); }
    
    // Fills the integer vector provided with a subset of [0,1,...,n-1] of length m=subset.size().
    // The value of seqno determines which subset is selected and values of seqno from zero
    // up to nCm-1 represent a complete enumeration of all possible m-subsets.
    // Returns true for a valid seqno >= 0 (in which case the subset elements will be in
    // increasing order) or false if seqno >= nCm (in which case the input subset is unchanged).
    bool getSubset(int n, unsigned long seqno, std::vector<int> &subset);

} // likely

#endif // LIKELY_BINNED_DATA_RESAMPLER
