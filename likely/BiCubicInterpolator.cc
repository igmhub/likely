// Created 29-Aug-2012 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>

#include "likely/BiCubicInterpolator.h"
#include "likely/RuntimeError.h"

#include <cmath>

namespace local = likely;

local::BiCubicInterpolator::BiCubicInterpolator(DataPlane data, double xspacing, int nx, int ny,
    double yspacing, double x0, double y0)
: _data(data), _xspacing(xspacing), _nx(nx), _ny(ny), _yspacing(yspacing), _x0(x0), _y0(y0),
_initialized(false)
{
    if(_ny == 0) {
        _ny = _nx;
    }
    if(_yspacing == 0) {
        _yspacing = _xspacing;
    }
    if(_nx <= 0 || _ny <= 0) throw RuntimeError("Bad dataplane dimensions.");
    if(_xspacing <= 0 || _yspacing <= 0) throw RuntimeError("Bad dataplane grid spacing.");
}

local::BiCubicInterpolator::~BiCubicInterpolator() { }

double local::BiCubicInterpolator::operator()(double x, double y) const {
    // Code here is based on:
    // https://svn.blender.org/svnroot/bf-blender/branches/volume25/source/blender/blenlib/intern/voxel.c
    
    // Map x,y to a point dx,dy in the plane [0,nx) x [0,ny)
    double dx(std::fmod((x-_x0)/_xspacing,_nx)), dy(std::fmod((y-_y0)/_yspacing,_ny));
    if(dx < 0) dx += _nx;
    if(dy < 0) dy += _ny;
    // Calculate the corresponding lower-bound grid indices.
    int xi = (int)std::floor(dx);
    int yi = (int)std::floor(dy);
    // Check if we can re-use coefficients from the last interpolation.
    if(!_initialized || xi != _i1 || yi != _i2) {
        // Extract the local vocal values and calculate partial derivatives.
		double x[16] = {
		    // values of f(x,y) at each corner.
		    _data[_index(xi,yi)],
		    _data[_index(xi+1,yi)],
		    _data[_index(xi,yi+1)],
		    _data[_index(xi+1,yi+1)],
            // values of df/dx at each corner.
		    0.5*(_data[_index(xi+1,yi)]-_data[_index(xi-1,yi)]),
		    0.5*(_data[_index(xi+2,yi)]-_data[_index(xi,yi)]),
			0.5*(_data[_index(xi+1,yi+1)]-_data[_index(xi-1,yi+1)]),
			0.5*(_data[_index(xi+2,yi+1)]-_data[_index(xi,yi+1)]),
            // values of df/dy at each corner.
		    0.5*(_data[_index(xi,yi+1)]-_data[_index(xi,yi-1)]),
		    0.5*(_data[_index(xi+1,yi+1)]-_data[_index(xi+1,yi-1)]),
			0.5*(_data[_index(xi,yi+2)]-_data[_index(xi,yi)]),
			0.5*(_data[_index(xi+1,yi+2)]-_data[_index(xi+1,yi)]),
            // values of d2f/dxdy at each corner.
		    0.25*(_data[_index(xi+1,yi+1)]-_data[_index(xi-1,yi+1)]-_data[_index(xi+1,yi-1)]+_data[_index(xi-1,yi-1)]),
			0.25*(_data[_index(xi+2,yi+1)]-_data[_index(xi,yi+1)]-_data[_index(xi+2,yi-1)]+_data[_index(xi,yi-1)]),
			0.25*(_data[_index(xi+1,yi+2)]-_data[_index(xi-1,yi+2)]-_data[_index(xi+1,yi)]+_data[_index(xi-1,yi)]),
			0.25*(_data[_index(xi+2,yi+2)]-_data[_index(xi,yi+2)]-_data[_index(xi+2,yi)]+_data[_index(xi,yi)])
		};
		// Convert pixel values and partial derivatives to interpolation coefficients.
    	for (int i=0;i<16;++i) {
    		_coefs[i] = 0.0;
    		for (int j=0;j<16;++j) {
    			_coefs[i] += _C[i][j]*x[j];
    		}
    	}
        // Remember this voxel for next time.
        _i1 = xi;
        _i2 = yi;
        _initialized = true;
    }
    // Evaluate the interpolation within this grid voxel.
    dx -= xi;
    dy -= yi;
    int ijkn(0);
    double dypow(1);
    double result(0);
    for(int j = 0; j < 4; ++j) {
        result += dypow*
            (_coefs[ijkn] + dx*(_coefs[ijkn+1] + dx*(_coefs[ijkn+2] + dx*_coefs[ijkn+3])));
        ijkn += 4;
        dypow *= dy;
    }
    return result;
}

int local::BiCubicInterpolator::_C[16][16] = {
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
    {-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0},
    { 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1},
    {-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1},
    { 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
    {-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1},
    { 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1}
};
