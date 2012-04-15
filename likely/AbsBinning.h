// Created 14-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_BINNING
#define LIKELY_ABS_BINNING

#include <iosfwd>

namespace likely {
	class AbsBinning {
	public:
		AbsBinning();
		virtual ~AbsBinning();
        // Returns the bin index [0,nBins-1] or else -1.
        virtual int getBinIndex(double value) const = 0;
        // Returns the total number of bins.
        virtual int getNBins() const = 0;
        // Returns the full width of the specified bin.
        virtual double getBinSize(int index) const = 0;
        // Returns the lower bound of the specified bin. Use index=nbins for the upper bound of the last bin.
        virtual double getBinLowEdge(int index) const = 0;
        // Returns the midpoint value of the specified bin.
        virtual double getBinCenter(int index) const;
        // Dumps this binning to the specified output stream in a standard format.
        void dump(std::ostream &os) const;
	private:
	}; // AbsBinning
} // likely

#endif // LIKELY_ABS_BINNING
