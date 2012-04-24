// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_DATA
#define LIKELY_BINNED_DATA

#include "likely/types.h"

#include <vector>

namespace likely {
    // Represents data that is binned independently along one or more axes.
	class BinnedData {
	public:
	    // Creates a new dataset with an arbitrary number of axes.
        explicit BinnedData(std::vector<AbsBinningCPtr> axes);
        // Creates a new dataset with a single axis.
		explicit BinnedData(AbsBinningCPtr axis1);
        // Creates a new dataset with two axes.
        BinnedData(AbsBinningCPtr axis1, AbsBinningCPtr axis2);
        // Creates a new dataset with three axes.
        BinnedData(AbsBinningCPtr axis1, AbsBinningCPtr axis2, AbsBinningCPtr axis3);
		virtual ~BinnedData();
		// Returns the number of axes used to bin this data.
        int getNAxes() const;
        // Returns the total number of bins covering the rectangular volume of our axes.
        int getNBinsTotal() const;
        // Returns the number of bins with data, which is never more than getNumBinsTotal().
        int getNBinsWithData() const;
	private:
        std::vector<AbsBinningCPtr> _axisBinning;
        int _nbins, _ndata;
        void _initialize();
	}; // BinnedData
	
    inline int BinnedData::getNAxes() const { return _axisBinning.size(); }
    inline int BinnedData::getNBinsTotal() const { return _nbins; }
    inline int BinnedData::getNBinsWithData() const { return _ndata; }

} // likely

#endif // LIKELY_BINNED_DATA
