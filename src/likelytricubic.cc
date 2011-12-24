// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A tri-cubic 3D interpolation algorithm test program.

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"

#include <iostream>
#include <cmath>

namespace lk = likely;
namespace po = boost::program_options;

class Periodic3D {
public:
    Periodic3D(double delta, int nx, int ny, int nz)
    : _lx(delta*nx), _ly(delta*ny), _lz(delta*nz), _twopi(4*atan2(1,0))
    { }
    double operator()(double x, double y, double z) const {
        return std::sin(_twopi*(x/_lx-0.3))*std::sin(2*_twopi*(y/_ly+0.3))*std::sin(3*_twopi*(z/_lz));
    }
private:
    double _lx, _ly, _lz,_twopi;
    int _nx, _ny, _nz;
};

int main(int argc, char **argv) {

    // Configure command-line option processing
    int nx,ny,nz,ntrial;
    double spacing;
    po::options_description cli("Tri-cubic interpolation test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("nx", po::value<int>(&nx)->default_value(20),
            "Number of subdivisions along grid x axis.")
        ("ny", po::value<int>(&ny)->default_value(10),
            "Number of subdivisions along grid y axis.")
        ("nz", po::value<int>(&nz)->default_value(15),
            "Number of subdivisions along grid z axis.")
        ("spacing", po::value<double>(&spacing)->default_value(0.1),
            "Spacing between grid points")
        ("ntrial", po::value<int>(&ntrial)->default_value(1000),
            "Number of random points for testing the interpolator.")
        ;

    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose"));

    if(nx <= 0 || ny <= 0 || nz <= 0) {
        std::cerr << "Bad dimensions nx,ny,nz." << std::endl;
        return -1;
    }
    
    try {
        // Create a periodic function for testing.
        Periodic3D f(spacing,nx,ny,nz);
        // Create a new datacube with the requested size that samples the testing function.
        int n(nx*ny*nz);
        lk::TriCubicInterpolator::DataCube data(new double[n]);
        for(int ix = 0; ix < nx; ++ix) {
            for(int iy = 0; iy < ny; ++iy) {
                for(int iz = 0; iz < nz; ++iz) {
                    data[ix + nx*(iy + ny*iz)] = f(ix*spacing,iy*spacing,iz*spacing);
                }
            }
        }
        // Interpolate in this datacube at random points.
        lk::TriCubicInterpolator interpolator(data,spacing,nx,ny,nz);
        lk::Random &random(lk::Random::instance());
        random.setSeed(1234);
        lk::WeightedAccumulator stats;
        double lx(nx*spacing), ly(ny*spacing), lz(nz*spacing);
        for(int trial = 0; trial < ntrial; ++trial) {
            double x = random.getUniform()*lx;
            double y = random.getUniform()*ly;
            double z = random.getUniform()*lz;
            double error = interpolator(x,y,z) - f(x,y,z);
            stats.accumulate(error);
        }
        std::cout << "mean error = " << stats.mean() << " (should be nearly zero)" << std::endl;
        std::cout << "RMS error = " << stats.error() << " (the smaller, the better)" << std::endl;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
