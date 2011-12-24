// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TRI_CUBIC_INTERPOLATOR
#define LIKELY_TRI_CUBIC_INTERPOLATOR

#include "boost/smart_ptr.hpp"

namespace likely {
	class TriCubicInterpolator {
	// Performs tri-cubic interpolation within a 3D periodic grid.
	public:
        typedef boost::shared_array<double> DataCube;
		TriCubicInterpolator(DataCube data, int n1, int n2 = 0, int n3 = 0);
		virtual ~TriCubicInterpolator();
	private:
        DataCube _data;
        int _n1, _n2, _n3;
	}; // TriCubicInterpolator
} // likely

#endif // LIKELY_TRI_CUBIC_INTERPOLATOR
