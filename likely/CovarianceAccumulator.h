// Created 20-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_COVARIANCE_ACCUMULATOR
#define LIKELY_COVARIANCE_ACCUMULATOR

#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <vector>

namespace likely {
    // Accumulates statistics to estimate the covariance of a dataset.
	class CovarianceAccumulator {
	public:
	    // Create a new accumulator for vectors of the specified size.
		explicit CovarianceAccumulator(int size);
		virtual ~CovarianceAccumulator();
		// Accumulate a single vector.
        void accumulate(std::vector<double> const &vector);
        void accumulate(double const *vector);
        // Return the estimated covariance matrix of all vectors accumulated so far.
        CovarianceMatrixPtr getCovariance() const;
	private:
        int _size;
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
	}; // CovarianceAccumulator
} // likely

#endif // LIKELY_COVARIANCE_ACCUMULATOR
