// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNED_DATA
#define LIKELY_BINNED_DATA

#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <vector>
#include <set>

namespace likely {
    // Represents data that is binned independently along one or more axes. Not all possible
    // bins within the rectangular volumed defined by the binning are assumed to be filled.
    // The binned data may have an associated covariance matrix. Most of this class' methods
    // refer to multidimensional bins with the single integer "global index" defined by
    // the getIndex method. BinnedData may be copied (using the default copy constructor or
    // the provided assignment operator), in which case the copy will have smart pointers
    // the the same underlying binning objects and covariance matrix (if one is present in
    // the original).
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
		
		// Shallow copying is supported via the default copy constructor, which makes copies of
		// our data, but just adds smart pointer references to the original object's binning
		// objects and covariance matrix (if any). Any covariance matrix with more than one
		// reference count is being shared and cannot be modified via the newly created object
		// or the original. Use the isCovarianceModifiable() method to test for this condition.
		// See the cloneCovariance() method if you actually want each object to have separate
		// and modifiable covariance matrices.
		
		// Clone method performs the same shallow copy but is virtual to allow polymorphic
		// copying. Set onlyBinning = true to only clone our binning specifications and
		// create an empty dataset. Subclasses X must override this method, normally with:
		//
		//   return binningOnly ? new X(getAxisBinning()) : new X(*this)
		//
		// This means that X::X(std::vector<AbsBinningCPtr> axes) and X::X(X const &other)
		// must also be valid (but the default copy ctor is usually ok).
        virtual BinnedData *clone(bool binningOnly = false) const;
		
		// Assignment operator supports the same shallow copy semantics.
        BinnedData& operator=(BinnedData other);
        friend void swap(BinnedData& a, BinnedData& b);

        // Adds another congruent binned dataset to our dataset (with weight 1).
        // Equivalent to add(other).
        BinnedData& operator+=(BinnedData const &other);
        // Adds another congruent binned dataset to our dataset with an arbitrary weight.
        // Note that using weights different from 1 will generally produce an incorrect
        // covariance matrix. For some common cases of correctly weighted combinations,
        // use a BinnedDataResampler.
        virtual BinnedData& add(BinnedData const &other, double weight = 1);
        // Tests if another binned dataset is congruent with ours. Congruence requires
        // identical binning specifications along each axis and (unless onlyBinning = true)
        // that the same bins be occupied and that both datasets either have or do not
        // have covariance matrices.
        virtual bool isCongruent(BinnedData const &other, bool onlyBinning = false) const;
		
		// Returns the number of axes used to bin this data.
        int getNAxes() const;
        // Returns the total number of bins covering the rectangular volume of our axes.
        int getNBinsTotal() const;
        // Returns the number of bins with data, which is never more than getNumBinsTotal().
        // The hasData method defines exactly what constitutes a bin with data.
        int getNBinsWithData() const;
        // Returns a vector of shared pointers to our axis specification objects.
        std::vector<AbsBinningCPtr> getAxisBinning() const;
        
        // Returns the global index corresponding to the specified bin index values along
        // each axis. The global index is defined as (i0*n1+i1)*n2+…) where ik, nk are
        // the bin index and number of bins for axis k, respectively. The global
        // index will always be >= 0 and < getNBinsTotal().
        int getIndex(std::vector<int> const &binIndices) const;
        // Returns the global index corresponding to the specified coordinate values along
        // each axis.
        int getIndex(std::vector<double> const &values) const;
        
        // Returns iterators pointing to the first and last global indices for bins with data.
        // Iteration order is defined by the order of setData(...) calls, and not by the global
        // index value.
        typedef std::vector<int>::const_iterator IndexIterator;
        IndexIterator begin() const;
        IndexIterator end() const;

        // Returns the global index corresponding to the specified offset, where offset = 0
        // is the first data value loaded by setData, offset = 1 is the next data value, etc.
        // This method is useful for matching up BinnedData entries with an ordered list of
        // values that was used to create it. Throws a RuntimeError if offset is out of range.
        int getIndexAtOffset(int offset) const;
        // Returns the offset corresponding to a global index, or throws a RuntimeError for
        // an index with no associated data or out of range.
        int getOffsetForIndex(int index) const;
        
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
        // RuntimeError if this bin does not contain any data. If weighted is true, then
        // the value returned is (Cinv.data)[index] instead of data[index]. Be aware that
        // going back and forth between weighted and unweighted access requires potentially
        // expensive covariance matrix operations. If this data has no covariance, then
        // weighted and unweighted values are equivalent.
        double getData(int index, bool weighted = false) const;
        // Sets the data value for the bin associated with the specified global index. After
        // calling this method successfully, hasData(index) will be true. Note that once this
        // object has started filling a covariance matrix, i.e., hasCovariance() == true, then
        // no new bins can be filled (but already filled bins can have their values changed.)
        // If weighted is true, then the value set is (Cinv.data)[index] instead of data[index].
        // Be aware that going back and forth between weighted and unweighted access requires
        // potentially expensive covariance matrix operations. If this data has no covariance, then
        // weighted and unweighted values are equivalent.
        void setData(int index, double value, bool weighted = false);
        // Adds the specified offset to the value for the specified bin, which must already
        // have data associated with it. If weighted is true, then the value being updated is
        // (Cinv.data)[index] instead of data[index]. Be aware that going back and forth between
        // weighted and unweighted access requires potentially expensive covariance matrix operations.
        // If this data has no covariance, then weighted and unweighted values are equivalent.
        void addData(int index, double offset, bool weighted = false);
        // Returns true if our internal data representation is currently weighted, i.e.,
        // stored as Cinv.data rather than data. Changes to our internal representation are
        // triggered automatically, so this method simply allows these changes to be tracked.
        bool isDataWeighted() const;

        // Returns true if covariance data is available.
        bool hasCovariance() const;
        // Returns true if covariance data can be modified. Covariance data
        // might be present but not modifiable if we are sharing our matrix (via a smart pointer)
        // with some other consumer. This will generally be the case when a BinnedData object
        // has been copied (via the copy ctor or assignment operator). Note that our modifiable
        // state depends on the creation and destruction of smart pointers in other objects, so
        // might change even when nothing changes internally within this object itself.
        bool isCovarianceModifiable() const;
        // Replaces our covariance, if any, with a clone of itself. Since a covariance matrix
        // cannot be modified via this class when its reference count is greater than one, this
        // method allows an object created via the copy constructor or assignment operator to
        // modify its covariance matrix, with a corresponding increase in memory usage.
        void cloneCovariance();
        // Drops any covariance matrix.
        void dropCovariance();
        // Returns the (inverse) covariance matrix element for the specified pair of global
        // indices, or throws a RuntimeError if either of the corresponding bins has no data,
        // or if no covariance has been specified for this data. Be aware that going back and
        // forth between Covariance and InverseCovariance operations requires potentially
        // expensive matrix operations.
        double getCovariance(int index1, int index2) const;
        double getInverseCovariance(int index1, int index2) const;
        // Sets the (inverse) covariance matrix element for the specified pair of global
        // indices, or throws a RuntimeError if either of the corresponding bins has no data.
        // After the first call to one of these methods, hasCovariance() == true and no
        // further calls to setData are allowed for bins that do not already have data.
        // Throws a RuntimeError if we have an unmodifiable covariance matrix. Be aware that
        // going back and forth between Covariance and InverseCovariance operations requires
        // potentially expensive matrix operations. Note that changing the (inverse) covariance
        // with these methods does not directly change the contents of our data vector, but
        // it does change the meaning of weighted data. For example, if isDataWeighted() is true,
        // then setCovariance() changes the subsequent result of getData(...,weighted=false) but
        // not of getData(...,weighted=true).
        void setCovariance(int index1, int index2, double value);
        void setInverseCovariance(int index1, int index2, double value);
        // Returns a const shared pointer to our covariance matrix, if any.
        CovarianceMatrixCPtr getCovarianceMatrix() const;
        // Transforms our covariance matrix C by replacing it with C.Dinv.C. On return, D
        // contains our original covariance matrix.
        void transformCovariance(CovarianceMatrixPtr D);

        // Calculates the chi-square = (data-pred).Cinv.(data-pred) for the specified
        // vector of predicted data, or throws a RuntimeError. The predicted data vector
        // must use the same index sequence as our index iterator. (The copy by value
        // used here is an optimization, not a mistake.)
        double chiSquare(std::vector<double> pred) const;
        
        // Prunes our data to the subset of bins listed by their global index in the
        // specified keep vector. Throws a RuntimeError if any indices are out of range.
        // Pruning is done in place and does not require any new memory allocation.
        // The axis binning is unchanged by pruning. A covariance matrix, if present,
        // will be changed. If an existing covariance matrix is not modifiable, it will
        // be cloned before pruning.
        void prune(std::set<int> const &keep);

        // Finalizes this object by preventing any further changes to our "shape", as
        // implemented by isCongruent(). Specifically, an existing covariance cannot be
        // dropped, a new covariance cannot be created, and previously unused bins
        // cannot be used. A finalized object cannot be pruned but it can be compressed.
        // The index iteration sequence of a finalized object is guaranteed never to change.
        virtual void finalize();
        // Retrurns true if this object has been finalized, or else false.
        bool isFinalized() const;
        
        // Requests that this object be compressed to reduce its memory usage,
        // if possible. Returns immediately if we are already compressed. Any compression
        // is lossless. Any subsequent reading or writing of covariance matrix elements
        // will automatically trigger a decompression. Return value indicates if any
        // compression was actually performed. Compression is considered a logically-const
        // operation since it is lossless, and an unmodifiable covariance matrix can still
        // be compressed. If weighted is true, then the data vector is converted to weighted
        // form, Cinv.data, (if a covariance is available) before compressing the covariance.
        // This allows this dataset to be added to other datasets while compressed but means
        // that calling getData() will trigger an automatic decompression of the covariance.
        bool compress(bool weighted = true) const;
        // Returns true if this covariance matrix is currently compressed. Note that uncompression
        // happens automatically, on demand, so there is no guarantee that a compressed object
        // will remain compressed.
        bool isCompressed() const;
        // Returns the memory usage of this object. Does not include the memory used by
        // the binning objects whose pointers are passed to our constructor. Use includeCovariance
        // to specify if the memory usage of any covariance matrix should be included in the result.
        std::size_t getMemoryUsage(bool includeCovariance = true) const;

	private:
        int _nbins;
        std::vector<AbsBinningCPtr> _axisBinning;
        enum { EMPTY_BIN = -1 };
        std::vector<int> _offset, _index;
        // Our data vector which might be weighted.
        mutable std::vector<double> _data;
        CovarianceMatrixPtr _covariance;
        // Is our _data vector weighted by _Cinv?
        mutable bool _weighted;
        // Have we been finalized?
        bool _finalized;
        // Changes whether our _data vector is weighted by _Cinv by multiplying
        // by Cinv or C, as needed.
        void _setWeighted(bool weighted) const;
        // Initializes a new object.
        void _initialize();
        // Throws a RuntimeError unless the specified global index is valid.
        void _checkIndex(int index) const;
	}; // BinnedData
	
    void swap(BinnedData& a, BinnedData& b);
	
    inline int BinnedData::getNAxes() const { return _axisBinning.size(); }
    inline int BinnedData::getNBinsTotal() const { return _nbins; }
    inline int BinnedData::getNBinsWithData() const { return _index.size(); }
    inline std::vector<AbsBinningCPtr> BinnedData::getAxisBinning() const { return _axisBinning; }
    inline bool BinnedData::hasCovariance() const { return _covariance.get() != 0; }
    inline bool BinnedData::isDataWeighted() const { return _weighted; }
    inline CovarianceMatrixCPtr BinnedData::getCovarianceMatrix() const { return _covariance; }
    inline bool BinnedData::isCovarianceModifiable() const {
        return 0 == _covariance.get() || _covariance.unique();
    }
    inline BinnedData::IndexIterator BinnedData::begin() const { return _index.begin(); }
    inline BinnedData::IndexIterator BinnedData::end() const { return _index.end(); }
    inline bool BinnedData::isFinalized() const { return _finalized; }
    inline BinnedData& BinnedData::operator+=(BinnedData const& other) { return add(other); }

} // likely

#endif // LIKELY_BINNED_DATA