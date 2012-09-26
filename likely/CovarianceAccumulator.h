// Created 20-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_COVARIANCE_ACCUMULATOR
#define LIKELY_COVARIANCE_ACCUMULATOR

#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <vector>
#include <iosfwd>

namespace likely {
    // Accumulates statistics to estimate the covariance of a dataset.
	class CovarianceAccumulator {
	public:
	    // Create a new accumulator for vectors of the specified size.
		explicit CovarianceAccumulator(int size);
		virtual ~CovarianceAccumulator();
		// Accumulate a single vector using the specified weight.
        void accumulate(std::vector<double> const &vector, double wgt = 1);
        void accumulate(double const *vector, double wgt = 1);
        // Accumulate the data vector of a BinnedData object.
        void accumulate(BinnedDataCPtr data, double wgt = 1);
        // Returns the number of vectors accumulated so far.
        int count() const;
        // Return the estimated covariance matrix of all vectors accumulated so far.
        // The result may not be positive definite (and this is not checked here)
        // but this can usually be fixed by accumulating more samples.
        CovarianceMatrixPtr getCovariance() const;
        // Dumps our internal state to the specified output stream in the following format:
        //
        //   size
        //   count
        //   sum-of-weights  (equal to count when wgt=1)
        //   col weighted-mean[col]  (0 <= col < size)
        //   ...
        //   row col weighted-covariance[row,col]  (0 <= col < size and 0 <= row <= col)
        //   ...
        //
        // The total number of lines is 3 + size*(size+3)/2. Floating point values are
        // written using the full internal precision.
        void dump(std::ostream &os) const;
	private:
        int _size;
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
	}; // CovarianceAccumulator
} // likely

#endif // LIKELY_COVARIANCE_ACCUMULATOR
