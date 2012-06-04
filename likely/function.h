// Created 11-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FUNCTION
#define LIKELY_FUNCTION

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

namespace likely {

	// Creates and returns a shared pointer to a generic function object that wraps a
	// shared pointer pimpl to an implementation function object of class P. The
	// returned shared pointer creates a new reference to the input shared pointer so that
	// the input object is guaranteed to stay alive as long as the returned object does.

	typedef boost::function<double (double)> GenericFunction;
    typedef boost::shared_ptr<GenericFunction> GenericFunctionPtr;
    template <class P> GenericFunctionPtr createFunctionPtr(boost::shared_ptr<P> pimpl);

} // likely

#endif // LIKELY_FUNCTION
