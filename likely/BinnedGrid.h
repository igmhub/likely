// Created 10-May-2013 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_GRID
#define LIKELY_BINNED_GRID

#include "likely/types.h"

#include <vector>

namespace likely {
    // Represents a multidimensional grid defined by an AbsBinning along each axis, and
    // defines a mapping between a global integer index and the individual axis indices.
	class BinnedGrid {
	public:
	    // Creates a new grid with the specified binning for each axis.
        explicit BinnedGrid(std::vector<AbsBinningCPtr> axes);
		virtual ~BinnedGrid();
		// Returns the number of axes for this grid.
        int getNAxes() const;
        // Returns the total number of bins covering the rectangular volume of this grid.
        int getNBinsTotal() const;
        // Returns the global index corresponding to the specified bin index values along
        // each axis. The global index is defined as (i0*n1+i1)*n2+…) where ik, nk are
        // the bin index and number of bins for axis k, respectively. The global
        // index will always be >= 0 and < getNBinsTotal().
        int getIndex(std::vector<int> const &binIndices) const;
        // Returns the global index corresponding to the specified coordinate values along
        // each axis.
        int getIndex(std::vector<double> const &values) const;
        // Fills the vector provided with the bin index values along each axis for the specified
        // global index.
        void getBinIndices(int index, std::vector<int> &binIndices) const;
        // Fills the vector provided with the bin centers along each axis for the specified
        // global index.
        void getBinCenters(int index, std::vector<double> &binCenters) const;
        // Fills the vector provided with the full bin widths along each axis for the specified
        // global index.
        void getBinWidths(int index, std::vector<double> &binWidths) const;
	private:
        // Throws a RuntimeError unless the specified global index is valid.
        void _checkIndex(int index) const;
        int _nbins;
        std::vector<AbsBinningCPtr> _axisBinning;
	}; // BinnedGrid

    inline int BinnedGrid::getNAxes() const { return _axisBinning.size(); }
    inline int BinnedGrid::getNBinsTotal() const { return _nbins; }

} // likely

#endif // LIKELY_BINNED_GRID
