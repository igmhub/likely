// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TRI_CUBIC_INTERPOLATOR
#define LIKELY_TRI_CUBIC_INTERPOLATOR

#include "boost/smart_ptr.hpp"

namespace likely {
	class TriCubicInterpolator {
	// Performs tri-cubic interpolation within a 3D periodic grid.
	public:
        typedef boost::shared_array<double> DataCube;
        // Initializes an interpolator using the specified datacube of length n1*n2*n3 where
        // data is ordered first along the n1 axis [0,0,0], [1,0,0], ..., [n1-1,0,0], [0,1,0], ...
        // If n2 and n3 are both omitted, then n1=n2=n3 is assumed. Data is assumed to be
        // equally spaced and periodic along each axis.
		TriCubicInterpolator(DataCube data, int n1, int n2 = 0, int n3 = 0);
		virtual ~TriCubicInterpolator();
        // Returns the interpolated data value for the specified x,y,z point in the unit cube.
        double operator()(double x, double y, double z) const;
	private:
	    // Returns the unrolled 1D index corresponding to [i1,i2,i3] after mapping to each ik into [0,nk).
	    // Assumes that i1 increases fastest in the 1D array.
        int _index(int i1, int i2, int i3) const;
        DataCube _data;
        int _n1, _n2, _n3;
        mutable int _i1, _i2, _i3;
        mutable double _coefs[64];
        mutable bool _initialized;
        static double _C[64][64];
	}; // TriCubicInterpolator
	
	inline int TriCubicInterpolator::_index(int i1, int i2, int i3) const {
        if((i1 %= _n1) < 0) i1 += _n1;
        if((i2 %= _n2) < 0) i2 += _n2;
        if((i3 %= _n3) < 0) i3 += _n3;
        return i1 + _n1*(i2 + _n2*i3);
	}

} // likely

#endif // LIKELY_TRI_CUBIC_INTERPOLATOR
