// Created 29-Aug-2012 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#ifndef LIKELY_BI_CUBIC_INTERPOLATOR
#define LIKELY_BI_CUBIC_INTERPOLATOR

#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <vector>

namespace likely {
	class BiCubicInterpolator {
	// Performs bi-cubic interpolation within a 2D periodic grid.
	public:
        typedef boost::shared_array<double> DataPlane;
        // Initializes an interpolator using the specified dataplane of length nx*ny where
        // data is ordered first along the nx axis [0,0], [1,0], ..., [nx-1,0], [0,1], ...
        // If ny is omitted, then nx=ny is assumed. If yspacing is omitted, then 
        // xspacing=yspacing is assumed. Data is assumed to be periodic along each axis,
        // with the coordinate origin (x0,y0) at grid index [0,0].
		BiCubicInterpolator(DataPlane data, double xspacing, int nx, int ny = 0,
		double yspacing = 0, double x0 = 0, double y0 = 0);
		virtual ~BiCubicInterpolator();
        // Returns the interpolated data value for the specified x,y point. If the point lies
        // outside the box [0,nx*xspacing) x [0,ny*yspacing), it will be folded
        // back assuming periodicity along each axis.
        double operator()(double x, double y) const;
        // Returns the grid parameters.
        double getXSpacing() const;
        double getYSpacing() const;
        int getNX() const;
        int getNY() const;
        double getX0() const;
        double getY0() const;
	private:
	    // Returns the unrolled 1D index corresponding to [i1,i2] after mapping to each ik into [0,nk).
	    // Assumes that i1 increases fastest in the 1D array.
        int _index(int i1, int i2) const;
        DataPlane _data;
        double _xspacing, _yspacing, _x0, _y0;
        int _nx, _ny;
        mutable int _i1, _i2;
        mutable double _coefs[16];
        mutable bool _initialized;
        static int _C[16][16];
    }; // BiCubicInterpolator

    inline double BiCubicInterpolator::getXSpacing() const { return _xspacing; }
    inline double BiCubicInterpolator::getYSpacing() const { return _yspacing; }
    inline int BiCubicInterpolator::getNX() const { return _nx; }
    inline int BiCubicInterpolator::getNY() const { return _ny; }
    inline double BiCubicInterpolator::getX0() const { return _x0; }
    inline double BiCubicInterpolator::getY0() const { return _y0; }

	inline int BiCubicInterpolator::_index(int i1, int i2) const {
        if((i1 %= _nx) < 0) i1 += _nx;
        if((i2 %= _ny) < 0) i2 += _ny;
        return i1 + _nx*i2;
	}
	
    // Returns a smart pointer to a bicubic interpolator based on control points read
    // from the specified file name.
    BiCubicInterpolatorPtr createBiCubicInterpolator(std::string const &filename);

} // likely

#endif // LIKELY_BI_CUBIC_INTERPOLATOR