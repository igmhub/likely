// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying likelihood analysis methods.

#include "likely/likely.h"

#include "boost/program_options.hpp"

#include <iostream>

namespace test = likely::test;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Likelihood analysis test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
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

    test::TestLikelihood likelyfn(3,1,-0.75,0.25);
    std::vector<double> params(3,0);
    params[1] = 1;
    std::cout << likelyfn(params) << std::endl;

    return 0;
}
