// Created 28-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_ENGINE
#define LIKELY_ABS_ENGINE

#include "boost/utility.hpp"

namespace likely {
	class AbsEngine : public boost::noncopyable {
	public:
		AbsEngine();
		virtual ~AbsEngine();
	private:
	}; // AbsEngine
} // likely

#endif // LIKELY_ABS_ENGINE
