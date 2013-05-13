// Created 14-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_BINNING
#define LIKELY_ABS_BINNING

#include "likely/types.h"

#include <string>
#include <iosfwd>

namespace likely {
	// Represents an abstract 1D binning of a real-valued independent variable. Each bin is associated
	// with an integer non-negative index and represented by low- and high-edge values. Bins are required
	// to be non-overlapping but gaps between bins are permitted. Increasing bin index corresponds to
	// increasing low and high edge values.
	class AbsBinning {
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
        // Prints this binning to the specified output stream in a format compatible with createBinning.
        virtual void printToStream(std::ostream &os) const = 0;
        // Returns true if the specified bin index is within range. Otherwise, throws a BinningError
        // if a non-empty error message format is provided, with the first %d replaced by the index,
        // or else returns false if no error message format is specified.
        bool isValidBinIndex(int index, std::string const &errorFormat = std::string("")) const;
	private:
	}; // AbsBinning
	
	// Creates a new binning object using the specification string provided, or throws a RuntimeError.
	// The supported specification string formats are:
	// - "[x1,x2,...,xn]" => NonUniformBinning where xi are the bin edges (so n-1 bins)
	// - "{x1,x2,...,xn}" => NonUniformSampling (or UniformSampling for only 2 points)
	// - "[lo:hi]*n" => UniformBinning covering [lo,hi] with n bins of size (lo-hi)/n
	// - "{lo:hi}*n" => UniformSampling covering {lo,hi} with n samples separated by (hi-lo)/(n-1)
	// Use AbsBinning::printToStream to write any binning object in the appropriate format.
    AbsBinningCPtr createBinning(std::string const &binningSpec);
    
} // likely

#endif // LIKELY_ABS_BINNING
