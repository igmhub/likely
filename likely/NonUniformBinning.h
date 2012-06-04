// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_NON_UNIFORM_BINNING
#define LIKELY_NON_UNIFORM_BINNING

#include "likely/AbsBinning.h"

#include <vector>

namespace likely {
    // Represents a non-uniform 1D binning, without gaps, of a finite interval.
	class NonUniformBinning : public AbsBinning {
	public:
	    // Creates a new binning using the specified list of bin edges, which must be in increasing
	    // order (but zero-width bins are allowed). The number of bins will be equal to the
	    // number of input vector elements minus one.
		explicit NonUniformBinning(std::vector<double> const &binEdges);
		virtual ~NonUniformBinning();
        // Returns the bin index [0,nBins-1] corresponding to the specified value, or throws a
        // BinningError if value does not fall in any bin.
        virtual int getBinIndex(double value) const;
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
        std::vector<double> _binEdges;
	}; // NonUniformBinning
} // likely

#endif // LIKELY_NON_UNIFORM_BINNING
