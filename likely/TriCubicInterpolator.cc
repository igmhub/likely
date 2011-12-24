// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/TriCubicInterpolator.h"
#include "likely/RuntimeError.h"

#include <cmath>

namespace local = likely;

local::TriCubicInterpolator::TriCubicInterpolator(DataCube data, int n1, int n2, int n3)
: _data(data), _n1(n1), _n2(n2), _n3(n3), _initialized(false)
{
    if(_n2 == 0 && _n3 == 0) {
        _n3 = _n2 = _n1;
    }
    if(_n1 <= 0 || _n2 <= 0 || _n3 <= 0) throw RuntimeError("Bad datacube dimensions.");
}

local::TriCubicInterpolator::~TriCubicInterpolator() { }

double local::TriCubicInterpolator::operator()(double x, double y, double z) const {
    // Map x,y,z to a point dx,dy,dz in the unit cube [0,1)^3
    double dx(std::fmod(x,1)), dy(std::fmod(y,1)), dz(std::fmod(z,1));
    if(dx < 0) dx += 1;
    if(dy < 0) dy += 1;
    if(dz < 0) dz += 1;
    // Calculate the corresponding lower-bound grid indices.
    int xi = (int)std::floor(dx*_n1);
    int yi = (int)std::floor(dy*_n2);
    int zi = (int)std::floor(dz*_n3);
    // Check if we can re-use coefficients from the last interpolation.
    if(!_initialized || xi != _i1 || yi != _i2 || zi != _i3) {
		double fval[8] = {
		    _data[_index(xi,yi,zi)],_data[_index(xi+1,yi,zi)],_data[_index(xi,yi+1,zi)],
		    _data[_index(xi+1,yi+1,zi)],_data[_index(xi,yi,zi+1)],_data[_index(xi+1,yi,zi+1)],
		    _data[_index(xi,yi+1,zi+1)],_data[_index(xi+1,yi+1,zi+1)]
		};
		double dfdxval[8] = {
		    0.5*(_data[_index(xi+1,yi,zi)]-_data[_index(xi-1,yi,zi)]),
		    0.5*(_data[_index(xi+2,yi,zi)]-_data[_index(xi,yi,zi)]),
			0.5*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi-1,yi+1,zi)]),
			0.5*(_data[_index(xi+2,yi+1,zi)]-_data[_index(xi,yi+1,zi)]),
			0.5*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi-1,yi,zi+1)]),
			0.5*(_data[_index(xi+2,yi,zi+1)]-_data[_index(xi,yi,zi+1)]),
			0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]),
			0.5*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)])
		};
		double dfdyval[8] = {
		    0.5*(_data[_index(xi,yi+1,zi)]-_data[_index(xi,yi-1,zi)]),
		    0.5*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi+1,yi-1,zi)]),
			0.5*(_data[_index(xi,yi+2,zi)]-_data[_index(xi,yi,zi)]),
			0.5*(_data[_index(xi+1,yi+2,zi)]-_data[_index(xi+1,yi,zi)]),
			0.5*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi-1,zi+1)]),
			0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]),
			0.5*(_data[_index(xi,yi+2,zi+1)]-_data[_index(xi,yi,zi+1)]),
			0.5*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)])
		};
		double dfdzval[8] = {
		    0.5*(_data[_index(xi,yi,zi+1)]-_data[_index(xi,yi,zi-1)]),
		    0.5*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi+1,yi,zi-1)]),
			0.5*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi+1,zi-1)]),
			0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]),
			0.5*(_data[_index(xi,yi,zi+2)]-_data[_index(xi,yi,zi)]),
			0.5*(_data[_index(xi+1,yi,zi+2)]-_data[_index(xi+1,yi,zi)]),
			0.5*(_data[_index(xi,yi+1,zi+2)]-_data[_index(xi,yi+1,zi)]),
			0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi+1,yi+1,zi)])
		};
		double d2fdxdyval[8] = {
		    0.25*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi-1,yi+1,zi)]-_data[_index(xi+1,yi-1,zi)]+_data[_index(xi-1,yi-1,zi)]),
			0.25*(_data[_index(xi+2,yi+1,zi)]-_data[_index(xi,yi+1,zi)]-_data[_index(xi+2,yi-1,zi)]+_data[_index(xi,yi-1,zi)]),
			0.25*(_data[_index(xi+1,yi+2,zi)]-_data[_index(xi-1,yi+2,zi)]-_data[_index(xi+1,yi,zi)]+_data[_index(xi-1,yi,zi)]),
			0.25*(_data[_index(xi+2,yi+2,zi)]-_data[_index(xi,yi+2,zi)]-_data[_index(xi+2,yi,zi)]+_data[_index(xi,yi,zi)]),
			0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]+_data[_index(xi-1,yi-1,zi+1)]),
			0.25*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi-1,zi+1)]+_data[_index(xi,yi-1,zi+1)]),
			0.25*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi-1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]+_data[_index(xi-1,yi,zi+1)]),
			0.25*(_data[_index(xi+2,yi+2,zi+1)]-_data[_index(xi,yi+2,zi+1)]-_data[_index(xi+2,yi,zi+1)]+_data[_index(xi,yi,zi+1)])
        };
		double d2fdxdzval[8] = {
		    0.25*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi-1,yi,zi+1)]-_data[_index(xi+1,yi,zi-1)]+_data[_index(xi-1,yi,zi-1)]),
			0.25*(_data[_index(xi+2,yi,zi+1)]-_data[_index(xi,yi,zi+1)]-_data[_index(xi+2,yi,zi-1)]+_data[_index(xi,yi,zi-1)]),
			0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi-1,yi+1,zi-1)]),
			0.25*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi+1,zi-1)]+_data[_index(xi,yi+1,zi-1)]),
			0.25*(_data[_index(xi+1,yi,zi+2)]-_data[_index(xi-1,yi,zi+2)]-_data[_index(xi+1,yi,zi)]+_data[_index(xi-1,yi,zi)]),
			0.25*(_data[_index(xi+2,yi,zi+2)]-_data[_index(xi,yi,zi+2)]-_data[_index(xi+2,yi,zi)]+_data[_index(xi,yi,zi)]),
			0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi-1,yi+1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi-1,yi+1,zi)]),
			0.25*(_data[_index(xi+2,yi+1,zi+2)]-_data[_index(xi,yi+1,zi+2)]-_data[_index(xi+2,yi+1,zi)]+_data[_index(xi,yi+1,zi)])
		};
		double d2fdydzval[8] = {
		    0.25*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi-1,zi+1)]-_data[_index(xi,yi+1,zi-1)]+_data[_index(xi,yi-1,zi-1)]),
			0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi+1,yi-1,zi-1)]),
			0.25*(_data[_index(xi,yi+2,zi+1)]-_data[_index(xi,yi,zi+1)]-_data[_index(xi,yi+2,zi-1)]+_data[_index(xi,yi,zi-1)]),
			0.25*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]-_data[_index(xi+1,yi+2,zi-1)]+_data[_index(xi+1,yi,zi-1)]),
			0.25*(_data[_index(xi,yi+1,zi+2)]-_data[_index(xi,yi-1,zi+2)]-_data[_index(xi,yi+1,zi)]+_data[_index(xi,yi-1,zi)]),
			0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi+1,yi-1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi+1,yi-1,zi)]),
			0.25*(_data[_index(xi,yi+2,zi+2)]-_data[_index(xi,yi,zi+2)]-_data[_index(xi,yi+2,zi)]+_data[_index(xi,yi,zi)]),
			0.25*(_data[_index(xi+1,yi+2,zi+2)]-_data[_index(xi+1,yi,zi+2)]-_data[_index(xi+1,yi+2,zi)]+_data[_index(xi+1,yi,zi)])
		};
		double d3fdxdydzval[8] = {
		    0.125f*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]+_data[_index(xi-1,yi-1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi-1,yi+1,zi-1)]+_data[_index(xi+1,yi-1,zi-1)]-_data[_index(xi-1,yi-1,zi-1)]),
			0.125f*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi-1,zi+1)]+_data[_index(xi,yi-1,zi+1)]-_data[_index(xi+2,yi+1,zi-1)]+_data[_index(xi,yi+1,zi-1)]+_data[_index(xi+2,yi-1,zi-1)]-_data[_index(xi,yi-1,zi-1)]),
			0.125f*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi-1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]+_data[_index(xi-1,yi,zi+1)]-_data[_index(xi+1,yi+2,zi-1)]+_data[_index(xi-1,yi+2,zi-1)]+_data[_index(xi+1,yi,zi-1)]-_data[_index(xi-1,yi,zi-1)]),
			0.125f*(_data[_index(xi+2,yi+2,zi+1)]-_data[_index(xi,yi+2,zi+1)]-_data[_index(xi+2,yi,zi+1)]+_data[_index(xi,yi,zi+1)]-_data[_index(xi+2,yi+2,zi-1)]+_data[_index(xi,yi+2,zi-1)]+_data[_index(xi+2,yi,zi-1)]-_data[_index(xi,yi,zi-1)]),
			0.125f*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi-1,yi+1,zi+2)]-_data[_index(xi+1,yi-1,zi+2)]+_data[_index(xi-1,yi-1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi-1,yi+1,zi)]+_data[_index(xi+1,yi-1,zi)]-_data[_index(xi-1,yi-1,zi)]),
			0.125f*(_data[_index(xi+2,yi+1,zi+2)]-_data[_index(xi,yi+1,zi+2)]-_data[_index(xi+2,yi-1,zi+2)]+_data[_index(xi,yi-1,zi+2)]-_data[_index(xi+2,yi+1,zi)]+_data[_index(xi,yi+1,zi)]+_data[_index(xi+2,yi-1,zi)]-_data[_index(xi,yi-1,zi)]),
			0.125f*(_data[_index(xi+1,yi+2,zi+2)]-_data[_index(xi-1,yi+2,zi+2)]-_data[_index(xi+1,yi,zi+2)]+_data[_index(xi-1,yi,zi+2)]-_data[_index(xi+1,yi+2,zi)]+_data[_index(xi-1,yi+2,zi)]+_data[_index(xi+1,yi,zi)]-_data[_index(xi-1,yi,zi)]),
			0.125f*(_data[_index(xi+2,yi+2,zi+2)]-_data[_index(xi,yi+2,zi+2)]-_data[_index(xi+2,yi,zi+2)]+_data[_index(xi,yi,zi+2)]-_data[_index(xi+2,yi+2,zi)]+_data[_index(xi,yi+2,zi)]+_data[_index(xi+2,yi,zi)]-_data[_index(xi,yi,zi)])
		};
        // Remember this voxel for next time.
        _i1 = xi;
        _i2 = yi;
        _i3 = zi;
        _initialized = true;
    }
    double result(0);
    for (int i=0;i<4;++i) {
		for (int j=0;j<4;++j) {
			for (int k=0;k<4;++k) {
				result += _coefs[i+4*j+16*k]*pow(x,i)*pow(y,j)*pow(z,k);
			}
		}
	}
    return result;
}
