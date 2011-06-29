// Created 25-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A weighted sum accumulator test program.

#include "likely/likely.h"

#include "boost/program_options.hpp"

#include <string>
#include <iostream>
#include <fstream>

namespace lk = likely;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string inputName;
    po::options_description cli("Weighted sum accumulator test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&inputName)->default_value(""),
            "Input filename to read sample x wgt values from.")
        ("ignore-extra", "Ignores extra input beyond x wgt columns.")
        ("sigmas", "Second column gives sigmas instead of weights.")
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
    bool verbose(vm.count("verbose")), sigmas(vm.count("sigmas")),
        ignoreExtra(vm.count("ignore-extra"));;
    
    try {
        // Open the requested input file or use stdin.
        std::istream *input(&std::cin);
        boost::scoped_ptr<std::istream> closeme;
        if(0 < inputName.length()) {
            closeme.reset(new std::ifstream(inputName.c_str()));
            input = closeme.get();
        }
        // Read two columns of x,wgt values from the input.
        std::vector<std::vector<double> > columns(2);
        lk::readVectors(*input, columns, ignoreExtra);
        closeme.reset();
        int nSamples(columns[0].size());
        if(verbose) {
            std::cout << "Read " << nSamples << " sample values." << std::endl;
        }
        // Prepare some weighted accumulators
        lk::WeightedAccumulator all, even, odd;
        // Loop over the input samples.
        for(int sample = 0; sample < nSamples; ++sample) {
            double wgt(columns[1][sample]);
            if(wgt <= 0) continue;
            if(sigmas) wgt = 1/(wgt*wgt);
            double val(columns[0][sample]);
            all.accumulate(val, wgt);
            if(sample % 2) {
                odd.accumulate(val, wgt);
            }
            else {
                even.accumulate(val, wgt);
            }
        }
        // Print the results.
        std::cout << "Accumulated " << all.count()
            << " samples with mean " << all.mean()
            << ", sum " << all.sum()
            << ", sqrt(variance) " << all.error()
            << ", sum(weights) " << all.sumOfWeights() << std::endl;
        if(sigmas) {
            std::cout << "Error on the mean is " << all.errorOnMean() << std::endl;
        }
        // Compare the results of combining the even/odd subsamples
        lk::WeightedCombiner combined;
        combined.combine(even);
        combined.combine(odd);
        std::cout << "Even + odd: " << combined.count()
            << " samples with mean " << combined.mean()
            << ", sum " << combined.sum()
            << ", sqrt(variance) " << combined.error()
            << ", sum(weights) " << combined.sumOfWeights() << std::endl;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
