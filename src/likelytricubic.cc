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

int main(int argc, char **argv) {

    // Configure command-line option processing
    int nx,ny,nz;
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
        // Create a new datacube with the requested size.
        int n(nx*ny*nz);
        lk::TriCubicInterpolator::DataCube data(new double[n]);
        double twopi(4*atan2(1,0));
        for(int ix = 0; ix < nx; ++ix) {
            for(int iy = 0; iy < ny; ++iy) {
                for(int iz = 0; iz < nz; ++iz) {
                    data[ix + nx*(iy + ny*iz)] = 1;
                        //std::sin(twopi*(ix+3)/nx)*std::sin(2*twopi*(iy-3)/ny)*std::sin(3*twopi*iz/nz);
                }
            }
        }        
        // Interpolate in this datacube.
        lk::TriCubicInterpolator interpolator(data,spacing,nx,ny,nz);
        std::cout << "interpolated = " << interpolator(1,2,3) << std::endl;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
