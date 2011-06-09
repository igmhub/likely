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

void saveSample(lk::Parameters const &params, double fVal, bool accepted) {
    boost::format real(" %.5f");
    std::cout << (accepted ? 1 : 0) << real % fVal;
    for(int i = 0; i < params.size(); ++i) std::cout << real % params[i];
    std::cout << std::endl;
}

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int npar,ncycle,nstep,seed;
    double rho,alpha,initial;
    po::options_description cli("Markov-chain Monte Carlo test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed for generating initial parameter values.")
        ("ncycle", po::value<int>(&ncycle)->default_value(10),
            "Number of MCMC cycles to perform.")
        ("nstep", po::value<int>(&nstep)->default_value(10),
            "Number of MCMC steps to take per parameter in each cycle.")
        ("npar", po::value<int>(&npar)->default_value(3),
            "Number of floating parameters to use.")
        ("rho", po::value<double>(&rho)->default_value(0),
            "NLL correlation coefficient in the range (-1,+1).")
        ("alpha", po::value<double>(&alpha)->default_value(0),
            "Size of NLL non-parabolic effects.")
        ("initial", po::value<double>(&initial)->default_value(2),
            "Initial value of all parameters.")
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
        lk::FunctionPtr f(new lk::Function(boost::ref(tester)));

        // Set the initial parameters and errors.
        lk::Parameters params(npar,initial);
        lk::PackedCovariance errors(npar,1);
        
        // Set the initial function minimum to use.
        double fval((*f)(params));
        lk::FunctionMinimum fmin(fval,params,errors,true);
        
        // Create an MCMC engine to use.
        lk::MarkovChainEngine mcmc(f,npar);

        // Loop over minimization trials.
        for(int cycle = 0; cycle < ncycle; ++cycle) {
            fval = mcmc.generate(fmin,params,fval,saveSample,nstep);
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
