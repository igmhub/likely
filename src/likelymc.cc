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

std::ofstream cycleOut,resultsOut;

void saveSample(lk::Parameters const &params, double fVal, bool accepted) {
    boost::format real(" %.5f");
    cycleOut << (accepted ? 1 : 0) << real % fVal;
    for(int i = 0; i < params.size(); ++i) cycleOut << real % params[i];
    cycleOut << std::endl;
}

void printSummary(lk::FunctionMinimumPtr fmin, int index, std::string const &tag) {
    boost::format valueFmt(" %.5f");
    lk::Parameters where(fmin->getParameters());
    int npar(where.size());
    resultsOut << boost::format("min%d[\"%s\"] = { %.5f") % index % tag % where[0];
    for(int i = 1; i < npar; ++i) resultsOut << ',' << valueFmt % where[i];
    resultsOut << " };" << std::endl;
    lk::PackedCovariancePtr covar(fmin->getCovariance());
    resultsOut << boost::format("covar%d[\"%s\"] = {\n") % index % tag;
    for(int i = 0; i < npar; ++i) {
        int index = i*(i+1)/2;
        resultsOut << "  {" << valueFmt % (*covar)[index];
        for(int j = 1; j < npar; ++j) {
            index = (i <= j) ? i+j*(j+1)/2 : j+i*(i+1)/2;
            resultsOut << ',' << valueFmt % (*covar)[index];
        }
        resultsOut << " }" << (i == npar-1 ? ' ':',') << std::endl;
    }
    resultsOut << "};" << std::endl;
}

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string tag;
    int npar,ncycle,naccept,maxtrials,seed;
    double rho,alpha,initial;
    po::options_description cli("Markov-chain Monte Carlo test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed for generating initial parameter values.")
        ("tag", po::value<std::string>(&tag)->default_value("mc"),
            "Tag to identify the program output.")
        ("ncycle", po::value<int>(&ncycle)->default_value(4),
            "Number of MCMC cycles to perform.")
        ("naccept", po::value<int>(&naccept)->default_value(1000),
            "Number of MCMC trial steps to accept in each cycle.")
        ("maxtrials", po::value<int>(&maxtrials)->default_value(0),
            "Maximum number of MCMC trials to generate in each cycle (0=unlimited).")
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
        // Open the summary output file.
        std::string resultsName = tag + "/results.m";
        resultsOut.open(resultsName.c_str());

        // Create a likelihood function using the command-line parameters.
        test::TestLikelihood tester(npar,1,rho,alpha);
        lk::FunctionPtr f(new lk::Function(boost::ref(tester)));

        // Set the initial parameters and errors.
        lk::Parameters params(npar,initial);
        lk::PackedCovariance errors(npar,1);
        
        // Set the initial function minimum to use.
        lk::FunctionMinimumPtr fmin(new lk::FunctionMinimum(
            (*f)(params),params,errors,true));
        printSummary(fmin,0,tag);
        
        // Create an MCMC engine to use.
        lk::MarkovChainEngine mcmc(f,lk::GradientCalculatorPtr(),npar,"mc::saunter");
        
        // Loop over MCMC cycles.
        boost::format cycleOutName("%s/cycle-%d.dat"),valueFmt(" %.5f");
        for(int cycle = 0; cycle < ncycle; ++cycle) {
            // Open a file to save this cycle's steps to.
            cycleOut.open(boost::str(cycleOutName % tag % cycle).c_str());
            // Run the MCMC generator.
            int nsample = mcmc.generate(fmin,naccept,maxtrials,saveSample);
            cycleOut.close();
            // Print a summary of this cycle.
            std::cout << boost::format("cycle %d accepted %d / %d\n")
                % cycle % naccept % nsample;
            printSummary(fmin,cycle+1,tag);
        }
        resultsOut << boost::format("trueNLL[\"%s\"] = nonlinearNLL[1,%f,%f,%d];\n")
            % tag % rho % alpha % npar;
        resultsOut.close();
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
