// Created 31-Aug-2012 by Michael Blomqvist (University of California, Irvine) <cblomqvi@uci.edu>
// A bi-cubic 2D interpolation algorithm test program.

#include "likely/Random.h"
#include "likely/BiCubicInterpolator.h"
#include "likely/WeightedAccumulator.h"
#include "likely/RuntimeError.h"

#include "boost/program_options.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"

#include <iostream>
#include <cmath>

namespace lk = likely;
namespace po = boost::program_options;

class Periodic2D {
public:
    Periodic2D(double delta, int nx, int ny)
    : _lx(delta*nx), _ly(delta*ny), _twopi(4*atan2(1,0))
    { }
    double operator()(double x, double y) const {
        return std::sin(_twopi*(x/_lx-0.3))*std::sin(2*_twopi*(y/_ly+0.3));
    }
private:
    double _lx, _ly, _twopi;
    int _nx, _ny;
};

int main(int argc, char **argv) {

    // Configure command-line option processing
    int nx,ny,ntrial;
    double spacing;
    po::options_description cli("Bi-cubic interpolation test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("nx", po::value<int>(&nx)->default_value(20),
            "Number of subdivisions along grid x axis.")
        ("ny", po::value<int>(&ny)->default_value(10),
            "Number of subdivisions along grid y axis.")
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

    if(nx <= 0 || ny <= 0) {
        std::cerr << "Bad dimensions nx,ny." << std::endl;
        return -1;
    }
    
    try {
        // Create a periodic function for testing.
        Periodic2D f(spacing,nx,ny);
        // Create a new dataplane with the requested size that samples the testing function.
        int n(nx*ny);
        lk::BiCubicInterpolator::DataPlane data(new double[n]);
        for(int ix = 0; ix < nx; ++ix) {
            for(int iy = 0; iy < ny; ++iy) {
                data[ix + nx*iy] = f(ix*spacing,iy*spacing);
            }
        }
        // Interpolate in this dataplane at random points.
        lk::BiCubicInterpolator interpolator(data,spacing,nx,ny);
        lk::RandomPtr random = lk::Random::instance();
        random->setSeed(1234);
        lk::WeightedAccumulator stats;
        double lx(nx*spacing), ly(ny*spacing);
        for(int trial = 0; trial < ntrial; ++trial) {
            double x = (random->getUniform()-0.25)*4*lx;
            double y = (random->getUniform()-0.25)*4*ly;
            double error = interpolator(x,y) - f(x,y);
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
