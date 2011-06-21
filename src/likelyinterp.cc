// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// An interpolation algorithm test program.

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>

namespace lk = likely;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string inputName,algorithms;
    int numSamples;
    po::options_description cli("Interpolation test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&inputName)->default_value(""),
            "Input filename to read interpolation control points from.")
        ("algorithms,a", po::value<std::string>(&algorithms)->default_value("linear"),
            "Comma-separated list of algorithms to use.")
        ("nsamples,n", po::value<int>(&numSamples)->default_value(100),
            "Number of interpolated samples to calculate and print.")
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
