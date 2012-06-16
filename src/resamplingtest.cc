// Created 16-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying resampling methods.

#include "likely/likely.h"

#include "boost/program_options.hpp"

#include <iostream>

namespace lk = likely;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int ndim,seed;
    po::options_description cli("Resampling methods test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("ndim", po::value<int>(&ndim)->default_value(3),
            "Number of dimensions for binning data.")
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

    try {
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
