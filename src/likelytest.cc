// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying likelihood analysis methods.

#include "likely/likely.h"

#include "Minuit2/MnPrint.h"

#include "gsl/gsl_multimin.h"

#include "boost/program_options.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_real.hpp"
#include "boost/random/variate_generator.hpp"

#include <iostream>
#include <cmath>

namespace lk = likely;
namespace test = likely::test;
namespace po = boost::program_options;
namespace mn = ROOT::Minuit2;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int npar,ntrial,seed;
    double rho,alpha;
    po::options_description cli("Likelihood analysis test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("trace", "Traces all calls to the likelihood function.")
        ("ntrial", po::value<int>(&ntrial)->default_value(100),
            "Number of minimization trials to perform.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed for generating initial parameter values.")
        ("npar", po::value<int>(&npar)->default_value(3),
            "Number of floating parameters to use.")
        ("rho", po::value<double>(&alpha)->default_value(0),
            "Likelihood linear correlation parameter.")
        ("alpha", po::value<double>(&alpha)->default_value(0),
            "Likelihood non-linearity parameter.")
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
    bool verbose(vm.count("verbose")), trace(vm.count("trace"));

    if(npar <= 0) {
        std::cerr << "Number of parameters (npar) must be > 0." << std::endl;
        return -2;
    }

    // Initialize a uniform random number generator.
    boost::mt19937 gen(seed);
    boost::uniform_real<> flat(-1,+1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random(gen,flat);
    
    // Print out column headings for our output below.
    std::cout << "norm gsl_simplex2 gsl_simplex2rand" << std::endl;

    try {
        // Create a likelihood function using the command-line parameters.
        test::TestLikelihood testfn(npar,1,rho,alpha);
        if(trace) testfn.setTrace(true);
        double trueMinVal(testfn.getMinimum());
        lk::FunctionMinimumPtr fmin;
            
        // Loop over minimization trials.
        for(int trial = 0; trial < ntrial; ++trial) {
            // Choose random initial parameter values within [-1,+1]
            lk::Parameters initial(npar);
            double norm(0);
            for(int par = 0; par < npar; ++par) {
                double value = random();
                initial[par] = value;
                norm += value*value;
            }
            norm = std::sqrt(norm);
            std::cout << norm;
            // Use fixed initial error estimates.
            lk::Parameters errors(npar,1);
            // Use methods that do not use the function gradient.
            fmin = lk::findMinimum(testfn,initial,errors,"gsl::simplex2");
            std::cout << ' ' << testfn.getCount() << ' ' << fmin->getMinValue()-trueMinVal;
            testfn.resetCount();
            fmin = lk::findMinimum(testfn,initial,errors,"gsl::simplex2rand");
            std::cout << ' ' << testfn.getCount() << ' ' << fmin->getMinValue()-trueMinVal;
            testfn.resetCount();
            std::cout << std::endl;
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
