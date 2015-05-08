// Created 29-Aug-2012 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef LIKELY_BI_CUBIC_INTERPOLATOR
#define LIKELY_BI_CUBIC_INTERPOLATOR

#include "boost/smart_ptr.hpp"

namespace likely {
	class BiCubicInterpolator {
	// Performs bi-cubic interpolation within a 2D periodic grid.
	public:
        typedef boost::shared_array<double> DataPlane;
        // Initializes an interpolator using the specified dataplane of length n1*n2 where
        // data is ordered first along the n1 axis [0,0], [1,0], ..., [n1-1,0], [0,1], ...
        // If n2 is omitted, then n1=n2 is assumed. Data is assumed to be equally spaced
        // and periodic along each axis, with the coordinate origin (x0,y0) at grid index [0,0].
		BiCubicInterpolator(DataPlane data, double spacing, int n1, int n2 = 0, double x0 = 0, double y0 = 0);
		virtual ~BiCubicInterpolator();
        // Returns the interpolated data value for the specified x,y point. If the point lies
        // outside the box [0,n1*spacing) x [0,n2*spacing), it will be folded
        // back assuming periodicity along each axis.
        double operator()(double x, double y) const;
        // Returns the grid parameters.
        double getSpacing() const;
        int getN1() const;
        int getN2() const;
	private:
	    // Returns the unrolled 1D index corresponding to [i1,i2] after mapping to each ik into [0,nk).
	    // Assumes that i1 increases fastest in the 1D array.
        int _index(int i1, int i2) const;
        DataPlane _data;
        double _spacing, _x0, _y0;
        int _n1, _n2;
        mutable int _i1, _i2;
        mutable double _coefs[16];
        mutable bool _initialized;
        static int _C[16][16];
	}; // BiCubicInterpolator

    inline double BiCubicInterpolator::getSpacing() const { return _spacing; }
    inline int BiCubicInterpolator::getN1() const { return _n1; }
    inline int BiCubicInterpolator::getN2() const { return _n2; }

	inline int BiCubicInterpolator::_index(int i1, int i2) const {
        if((i1 %= _n1) < 0) i1 += _n1;
        if((i2 %= _n2) < 0) i2 += _n2;
        return i1 + _n1*i2;
	}

} // likely

#endif // LIKELY_BI_CUBIC_INTERPOLATOR