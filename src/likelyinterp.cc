// Created 21-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// An interpolation algorithm test program.

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace lk = likely;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string inputName,algorithms,formatString;
    int numSamples;
    po::options_description cli("Interpolation test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&inputName)->default_value(""),
            "Input filename to read interpolation control points from.")
        ("ignore-extra", "Ignores extra input beyond x y columns.")
        ("algorithms,a", po::value<std::string>(&algorithms)->default_value("linear"),
            "Comma-separated list of algorithms to use.")
        ("nsamples,n", po::value<int>(&numSamples)->default_value(100),
            "Number of interpolated samples to calculate and print.")
        ("format", po::value<std::string>(&formatString)->default_value("%.5f"),
            "Printf format to use for x and y values.")
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
    bool verbose(vm.count("verbose")), ignoreExtra(vm.count("ignore-extra"));
    
    if(numSamples <= 0) {
        std::cerr << "nsamples must be > 0." << std::endl;
        return -2;
    }
    
    try {
        // Open the requested input file or use stdin.
        std::istream *input(&std::cin);
        boost::scoped_ptr<std::istream> closeme;
        if(0 < inputName.length()) {
            closeme.reset(new std::ifstream(inputName.c_str()));
            input = closeme.get();
        }
        // Read two columns of x,y control points from the input.
        std::vector<std::vector<double> > columns(2);
        lk::readVectors(*input, columns, ignoreExtra);
        closeme.reset();
        if(verbose) {
            std::cout << "Read " << columns[0].size() << " control points." << std::endl;
        }
        // Loop over the requested interpolation algorithms.
        std::vector<lk::InterpolatorPtr> interpolators;
        std::list<std::string> tokens;
        boost::split(tokens, algorithms, boost::is_any_of(","));
        BOOST_FOREACH(std::string const &token, tokens) {
            interpolators.push_back(
                lk::InterpolatorPtr(new lk::Interpolator(columns[0],columns[1],token)));
            if(verbose) {
                std::cout << "Created '" << token << "' interpolator." << std::endl;
            }
        }
        // Prepare the output formatter.
        boost::format fmt(formatString);
        // Loop over the x range with the requested number of steps.
        double xlo(columns[0].front()), xhi(columns[0].back()), dx((xhi-xlo)/numSamples);
        for(int i = 0; i <= numSamples; i++) {
            double x = xlo + i*dx;
            std::cout << fmt % x;
            // Print the results of each interpolator on a single line.
            BOOST_FOREACH(lk::InterpolatorPtr interpolator, interpolators) {
                double y = (*interpolator)(x);
                std::cout << ' ' << fmt % y;
            }
            std::cout << std::endl;
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
