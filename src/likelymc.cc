// Created 9-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A Markov-chain Monte Carlo test program.

#include "config.h"

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_on_sphere.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/format.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <iostream>

namespace lk = likely;
namespace test = likely::test;
namespace po = boost::program_options;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int npar,ntrial,seed;
    double rho,alpha,radius;
    po::options_description cli("Markov-chain Monte Carlo test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("ntrial", po::value<int>(&ntrial)->default_value(100),
            "Number of minimization trials to perform.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed for generating initial parameter values.")
        ("npar", po::value<int>(&npar)->default_value(3),
            "Number of floating parameters to use.")
        ("rho", po::value<double>(&rho)->default_value(0),
            "NLL correlation coefficient in the range (-1,+1).")
        ("alpha", po::value<double>(&alpha)->default_value(0),
            "Size of NLL non-parabolic effects.")
        ("radius", po::value<double>(&radius)->default_value(2),
            "Radius of the initial-parameter value sphere.")
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

    if(npar <= 0) {
        std::cerr << "Number of parameters (npar) must be > 0." << std::endl;
        return -2;
    }

    // Initialize the random number generator.
    lk::Random &random(lk::Random::instance());
    random.setSeed(seed);
    
    // Print out column headings for the output we generate below.
    //std::cout << "method ncall ngrad accuracy" << std::endl;

    try {
        // Create a likelihood function using the command-line parameters.
        test::TestLikelihood tester(npar,1,rho,alpha);        

        // Since our function object has internal state (call counters), we
        // want to ensure that the original object is passed around and
        // never copied. Here are two ways to do this:
        lk::FunctionPtr f(new lk::Function(boost::ref(tester)));

        // Use boost::bind for the gradient calculator.
        lk::GradientCalculatorPtr gc(new lk::GradientCalculator(
            boost::bind(&test::TestLikelihood::evaluateGradient,&tester,_1,_2)));

        // Loop over minimization trials.
        for(int trial = 0; trial < ntrial; ++trial) {
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
