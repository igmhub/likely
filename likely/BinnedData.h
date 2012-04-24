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
        int getDimension() const;
	private:
        std::vector<AbsBinningCPtr> _axisBinning;
	}; // BinnedData
	
    inline int BinnedData::getDimension() const { return _axisBinning.size(); }

} // likely

#endif // LIKELY_BINNED_DATA
