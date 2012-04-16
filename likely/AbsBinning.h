// Created 14-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_BINNING
#define LIKELY_ABS_BINNING

#include <iosfwd>

namespace likely {
	class AbsBinning {
	// Represents an abstract 1D binning of a real-valued independent variable. Each bin is associated
	// with an integer non-negative index and represented by low- and high-edge values. Bins are required
	// to be non-overlapping but gaps between bins are permitted. Increasing bin index corresponds to
	// increasing low and high edge values.
	public:
		AbsBinning();
		virtual ~AbsBinning();
        // Returns the bin index [0,nBins-1] corresponding to the specified value, or throws a
        // BinningError if value does not fall in any bin.
        virtual int getBinIndex(double value) const = 0;
        // Returns the total number of bins.
        virtual int getNBins() const = 0;
        // Returns the lower bound of the specified bin. Throws a BinningError if index is out of range.
        virtual double getBinLowEdge(int index) const = 0;
        // Returns the upper bound of the specified bin. Throws a BinningError if index is out of range.
        virtual double getBinHighEdge(int index) const = 0;
        // Returns the full width (hi-lo) of the specified bin, which might be zero if this bin represents
        // a point sample rather than an integral over some interval. Throws a BinningError if index
        // is out of range.
        virtual double getBinWidth(int index) const;
        // Returns the midpoint value (lo+hi)/2 of the specified bin. Throws a BinningError if index is
        // out of range.
        virtual double getBinCenter(int index) const;
        // Dumps this binning to the specified output stream in a standard one-line format.
        void dump(std::ostream &os) const;
	private:
	}; // AbsBinning
} // likely

#endif // LIKELY_ABS_BINNING
