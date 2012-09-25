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
		BinnedDataResampler(bool useScalarWeights = false, RandomPtr random = RandomPtr());
		virtual ~BinnedDataResampler();
		// Adds a copy of the specified observation. Throws a RuntimeError if this observation
		// is not congruent with existing observations. You are allowed to add the same
		// observation several times, but you normally don't want to do this. Calls to
		// this method can be interspersed with calls to resampling methods below. Since
		// observations are copied when you add them, subsequent changes to the observation
		// will have not effect on the resampler. Also, there is no need to compress observations
		// before adding them (adding compressed observations can be slower, but will not uncompress
		// the input observation). Returns the index of the added observation. Passing the index
		// of a previously added observation as reuseCovIndex will result in the two observations
		// sharing the same covariance matrix object (and the input observation does not need to
		// have any covariance matrix as long as it is unweighted). If useScalarWeights is false,
		// the re-using covariances should give identical results but using less memory. However,
		// if useScalarWeights is true, then results are only identical in the limit that all
		// covariances are proportional (since we assume that the re-used covariance is proportional
		// to the combined covariances seen so far).
        int addObservation(BinnedDataCPtr observation, int reuseCovIndex = -1);
        // Returns the number of observations available for resampling.
        int getNObservations() const;
        // Returns a shared pointer to the specified (readonly) observation.
        BinnedDataCPtr getObservation(int index) const;
        // Returns a shared pointer to a (modifiable) copy of the specified observation.
        BinnedDataPtr getObservationCopy(int index, bool addCovariance = true) const;
        // Returns a shared pointer to a new BinnedData that combines all observations added
        // so far. Each call to this method builds a new combined dataset so save the result
        // if you don't want independent copies (and you know that the observations have not
        // changed since the last combination).
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
        BinnedDataPtr jackknife(int ndrop, unsigned long seqno, bool addCovariance = true) const;
        // Returns a shared pointer to a new BinnedData that represents a bootstrap resampling
        // of our observations of the specified size, which defaults to the number of observations
        // when zero. The fixCovariance option requests that the final covariance matrix be corrected
        // for double counting of identical observations. However, this is a relatively slow operation
        // for large datasets, involving a triple matrix product, so you might want to use
        // fixCovariance = false if you don't need an accurate covariance matrix. Chi-square values
        // calculated with fixCovariance = false will be roughly twice as large as the correct values
        // obtained with fixCovariance = true.
        BinnedDataPtr bootstrap(int size = 0, bool fixCovariance = true, bool addCovariance = true) const;
        // Returns a shared pointer to a new CovarianceMatrix that estimates the covariance of
        // our combined observations using the specified number of bootstrap samples with a
        // CovarianceAccumulator. Throws a RuntimeError if the estimated covariance is not positive
        // definite, which can usually be fixed with more samples.
        CovarianceMatrixPtr estimateCombinedCovariance(int nSamples, int messageInterval = 0) const;
	private:
	    // Adds a covariance matrix to a resampling built with scalar weights. The matrix will be
	    // a copy of our combined covariance scaled by the ratio of our _combinedScalarWeight to
	    // the sample's scalar weight.
        void _addCovariance(BinnedDataPtr sample) const;
        bool _useScalarWeights;
        mutable RandomPtr _random;
        std::vector<BinnedDataCPtr> _observations;
        double _combinedScalarWeight;
        BinnedDataPtr _combined;
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
