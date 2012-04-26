// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_DATA
#define LIKELY_BINNED_DATA

#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <vector>

namespace likely {
    // Represents data that is binned independently along one or more axes. Not all possible
    // bins within the rectangular volumed defined by the binning are assumed to be filled.
    // The binned data may have an associated covariance matrix. Most of this class' methods
    // refer to multidimensional bins with the single integer "global index" defined by
    // the getIndex method.
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
        // The hasData method defines exactly what constitutes a bin with data.
        int getNBinsWithData() const;
        
        // Returns the global index corresponding to the specified bin index values along each
        // axis. The global index is defined as i0 + n0*(i1 + n1*(i2 + n2*(...))) where
        // ik, nk are the bin index and number of bins for axis k, respectively. The global
        // index will always be >= 0 and < getNBinsTotal().
        int getIndex(std::vector<int> const &binIndices) const;
        // Returns the global index corresponding to the specified coordinate values along each
        // axis.
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
        
        // Returns true if the bin corresponding to the specified global index has data, or
        // else returns false. Note that a data whose contents is zero is not considered empty.
        // An empty bin is one that has never had any value assigned to it.
        bool hasData(int index) const;
        // Returns the data associated with the specified global index or else throws a
        // RuntimeError if this bin does not contain any data.
        double getData(int index) const;
        // Sets the data value for the bin associated with the specified global index. After
        // calling this method successfully, hasData(index) will be true. Note that once this
        // object has started filling a covariance matrix, i.e., hasCovariance() == true, then
        // no new bins can be filled (but already filled bins can have their values changed.)
        void setData(int index, double value);
        // Adds the specified offset to the value for the specified bin, which must already
        // have data associated with it.
        void addData(int index, double offset);

        // Returns true if covariance data is available.
        bool hasCovariance() const;
        // Returns the (inverse) covariance matrix element for the specified pair of global
        // indices, or throws a RuntimeError if either of the corresponding bins has no data,
        // or if no covariance has been specified for this data.
        double getCovariance(int index1, int index2) const;
        double getInverseCovariance(int index1, int index2) const;
        // Sets the (inverse) covariance matrix element for the specified pair of global
        // indices, or throws a RuntimeError if either of the corresponding bins has no data.
        // After the first call to one of these methods, hasCovariance() == true and no
        // further calls to setData are allowed for bins that do not already have data.
        void setCovariance(int index1, int index2, double value);
        void setInverseCovariance(int index1, int index2, double value);

        // Requests that this object be compressed to reduce its memory usage,
        // if possible. Returns immediately if we are already compressed. Any compression
        // is lossless. Any subsequent reading or writing of covariance matrix elements
        // will automatically trigger a decompression. Return value
        // indicates if any compression was actually performed.
        bool compress() const;
        // Returns true if this covariance matrix is currently compressed.
        bool isCompressed() const;
        // Returns the memory usage of this object. Does not include the memory used by
        // the binning objects whose pointers are passed to our constructor.
        std::size_t getMemoryUsage() const;

	private:
        int _nbins, _ndata;
        std::vector<AbsBinningCPtr> _axisBinning;
        enum { EMPTY_BIN = -1 };
        std::vector<int> _offset, _index;
        std::vector<double> _data;
        boost::shared_ptr<CovarianceMatrix> _covariance;
        void _initialize();
        void _checkIndex(int index) const;
	}; // BinnedData
	
    inline int BinnedData::getNAxes() const { return _axisBinning.size(); }
    inline int BinnedData::getNBinsTotal() const { return _nbins; }
    inline int BinnedData::getNBinsWithData() const { return _ndata; }
    inline bool BinnedData::hasCovariance() const { return _covariance.get() != 0; }

} // likely

#endif // LIKELY_BINNED_DATA
