// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying likelihood analysis methods.

#include "likely/likely.h"

#include "Minuit2/MnPrint.h"

#include "gsl/gsl_multimin.h"

#include "boost/program_options.hpp"

#include <iostream>

namespace lk = likely;
namespace test = likely::test;
namespace po = boost::program_options;
namespace mn = ROOT::Minuit2;

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int npar;
    double rho,alpha;
    po::options_description cli("Likelihood analysis test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("npar,n", po::value<int>(&npar)->default_value(3),
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
    bool verbose(vm.count("verbose"));

    if(npar <= 0) {
        std::cerr << "Number of parameters (npar) must be > 0." << std::endl;
        return -2;
    }

    try {
        // Create the likelihood function to use.
        test::TestLikelihood testfn(npar,1,rho,alpha);
        testfn.setTrace(true);
    
        // Specify the initial parameter values and error estimates.
        std::vector<double> initial(npar,1),errors(npar,1);
        std::cout << "f0 = " << testfn(initial) << std::endl;
    
        //lk::AbsMinimizerPtr minimizer(new lk::GslMinimizer(testfn,npar));
        //lk::Parameters final(minimizer->minimize(initial,errors));

        lk::FunctionMinimum fmin = lk::findMinimum(testfn,initial,errors,"gsl::simplex2");
        //fmin.print(std::cout);
        //fmin.estimateCovariance("mn2:hesse");
        //fmin.estimateError();

        //lk::MinuitEngine minuit(testfn,npar);
        //mn::FunctionMinimum mfit = minuit.simplex(initial,errors);
        //mn::FunctionMinimum mfit = minuit.variableMetric(initial,errors);
        //std::cout << mfit;

        //lk::GslEngine gsl(testfn,npar);
        //gsl.minimize(gsl_multimin_fminimizer_nmsimplex2,initial,errors);
    
        //lk::MarkovChainEngine mc(testfn,npar);
        //lk::Parameters sample = mc.advance(initial,errors,10);
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
