// Created 9-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A Markov-chain Monte Carlo test program.

#include "config.h"

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <iostream>
#include <fstream>

namespace lk = likely;
namespace test = likely::test;
namespace po = boost::program_options;

std::ofstream cycleOut;

void saveSample(lk::Parameters const &params, double fVal, bool accepted) {
    boost::format real(" %.5f");
    cycleOut << (accepted ? 1 : 0) << real % fVal;
    for(int i = 0; i < params.size(); ++i) cycleOut << real % params[i];
    cycleOut << std::endl;
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
    
    try {
        // Create a likelihood function using the command-line parameters.
        test::TestLikelihood tester(npar,1,rho,alpha);
        lk::FunctionPtr f(new lk::Function(boost::ref(tester)));

        // Set the initial parameters and errors.
        lk::Parameters params(npar,initial);
        lk::PackedCovariance errors(npar,1);
        
        // Set the initial function minimum to use.
        lk::FunctionMinimumPtr fmin(new lk::FunctionMinimum(
            (*f)(params),params,errors,true));
        
        // Create an MCMC engine to use.
        lk::MarkovChainEngine mcmc(f,npar);

        // Loop over MCMC cycles.
        boost::format cycleOutName("cycle-%d.dat"),valueFmt(" %.5f");
        for(int cycle = 0; cycle < ncycle; ++cycle) {
            // Open a file to save this cycle's steps to.
            cycleOut.open(boost::str(cycleOutName % cycle).c_str());
            // Run the MCMC generator.
            int naccept = mcmc.generate(fmin,saveSample,nstep);
            cycleOut.close();
            // Print a summary of this cycle.
            std::cout << boost::format("cycle %d accepted %d / %d\n")
                % cycle % naccept % nstep;
            lk::Parameters where(fmin->getParameters());
            std::cout << "where = {" << valueFmt % where[0];
            for(int i = 1; i < npar; ++i) std::cout << ',' << valueFmt % where[i];
            std::cout << " };" << std::endl;
            lk::PackedCovariancePtr covar(fmin->getCovariance());
            std::cout << "covar = {" << std::endl;
            for(int i = 0; i < npar; ++i) {
                int index = i*(i+1)/2;
                std::cout << "  {" << valueFmt % (*covar)[index];
                for(int j = 1; j < npar; ++j) {
                    index = (i <= j) ? i+j*(j+1)/2 : j+i*(i+1)/2;
                    std::cout << ',' << valueFmt % (*covar)[index];
                }
                std::cout << " }" << (i == npar-1 ? ' ':',') << std::endl;
            }
            std::cout << "};" << std::endl;
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
