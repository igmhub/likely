// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_NON_UNIFORM_SAMPLING
#define LIKELY_NON_UNIFORM_SAMPLING

#include "likely/AbsBinning.h"

#include <vector>

namespace likely {
	// Represents a non-uniform 1D sampling of a finite interval. A sampling is a binning with zero-width
	// bins, and represents point samples of an interval rather than bin integrals.
	class NonUniformSampling : public AbsBinning {
	public:
	    // Creates a new sampling using the specified list of sample points, which must be in increasing
	    // order. Throws a BinningError unless maxValue >= minValue.
		explicit NonUniformSampling(std::vector<double> const &samplePoints);
		virtual ~NonUniformSampling();
        // Returns the total number of bins.
        virtual int getNBins() const;
        // Returns the lower bound of the specified bin. Throws a BinningError if index is out of range.
        virtual double getBinLowEdge(int index) const;
        // Returns the upper bound of the specified bin. Throws a BinningError if index is out of range.
        virtual double getBinHighEdge(int index) const;
        // Returns the full width (hi-lo) of the specified bin, which might be zero if this bin represents
        // a point sample rather than an integral over some interval. Throws a BinningError if index
        // is out of range.
        virtual double getBinWidth(int index) const;
        // Returns the midpoint value (lo+hi)/2 of the specified bin. Throws a BinningError if index is
        // out of range.
        virtual double getBinCenter(int index) const;
	private:
        std::vector<double> _samplePoints;
	}; // NonUniformSampling
} // likely

#endif // LIKELY_NON_UNIFORM_SAMPLING
