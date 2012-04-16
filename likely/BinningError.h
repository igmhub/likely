// Created 16-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_BINNING_ERROR
#define LIKELY_BINNING_ERROR

#include "likely/RuntimeError.h"

namespace likely {
	class BinningError : public RuntimeError {
	public:
		BinningError(std::string const &reason);
		virtual ~BinningError() throw ();
	private:
	}; // BinningError

	inline BinningError::BinningError(std::string const &reason)
	: RuntimeError(reason) { }

    inline BinningError::~BinningError() throw () { }

} // likely

#endif // LIKELY_BINNING_ERROR
