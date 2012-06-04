// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_UNIFORM_SAMPLING
#define LIKELY_UNIFORM_SAMPLING

#include "likely/AbsBinning.h"

namespace likely {
	// Represents a uniform 1D sampling of a finite interval. A sampling is a binning with zero-width
	// bins, and represents point samples of an interval rather than bin integrals.
	class UniformSampling : public AbsBinning {
	public:
	    // Creates a new uniform sampling for the interval [minValue,maxValue] using the specified
	    // number of samples. Throws a BinningError unless maxValue > minValue and nSamples > 1,
	    // or else maxValue==minValue and nSamples==1. See getBinIndex for a description of ftol.
		UniformSampling(double minValue, double maxValue, int nSamples, double ftol = 1e-4);
		virtual ~UniformSampling();
        // Returns the bin index [0,nBins-1] corresponding to the specified value, or throws a
        // BinningError if value is not within ftol*spacing of a sample point, where spacing is
        // the fixed distance between samples returned by getBinWidth().
        virtual int getBinIndex(double value) const;
        // Returns the total number of bins, which is equal to the number of samples.
        virtual int getNBins() const;
        // Returns the lower bound of the specified bin. Throws a BinningError if index is out of range.
        virtual double getBinLowEdge(int index) const;
        // Returns the upper bound of the specified bin. Throws a BinningError if index is out of range.
        virtual double getBinHighEdge(int index) const;
        // Returns the full width (hi-lo) of the specified bin, which is zero for a sampling.
        // Throws a BinningError if index is out of range.
        virtual double getBinWidth(int index) const;
        // Returns the midpoint value (lo+hi)/2 of the specified bin. Throws a BinningError if index is
        // out of range.
        virtual double getBinCenter(int index) const;
	private:
        double _minValue, _maxValue, _sampleSpacing, _ftol;
        int _nSamples;
	}; // UniformSampling
} // likely

#endif // LIKELY_UNIFORM_SAMPLING
