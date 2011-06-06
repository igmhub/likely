// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying likelihood analysis methods.

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

#include <cassert> //!!

namespace lk = likely;
namespace test = likely::test;
namespace po = boost::program_options;

void useMethod(int methodId, std::string const &methodName,
lk::FunctionPtr f, lk::Parameters const &initial,
lk::Parameters const &errors, double prec) {
    boost::format fmt("%d %.4f %.4f %.4f\n");
    lk::FunctionMinimumPtr fmin = lk::findMinimum(f,initial,errors,methodName,prec);
    std::cout << fmt % methodId % std::log10(lk::lastMinEvalCount)
        % (lk::lastMinGradCount ? std::log10(lk::lastMinGradCount) : 0.)
        % -std::log10(fmin->getMinValue());
}

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int npar,ntrial,seed;
    double rho,alpha,radius;
    po::options_description cli("Likelihood analysis test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("eval", "Prints the NLL and its gradient at a fixed point.")
        ("trace", "Traces all calls to the NLL function.")
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
    bool verbose(vm.count("verbose")), trace(vm.count("trace")),eval(vm.count("eval"));

    if(npar <= 0) {
        std::cerr << "Number of parameters (npar) must be > 0." << std::endl;
        return -2;
    }

    // Initialize a uniform random number generator.
    boost::mt19937 gen(seed);
    boost::uniform_on_sphere<> spherical(npar);
    boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<> >
        randomOnSphere(gen,spherical);
    
    // Print out column headings for the output we generate below.
    std::cout << "method ncall ngrad accuracy" << std::endl;

    try {
        // Create a likelihood function using the command-line parameters.
        test::TestLikelihood tester(npar,1,rho,alpha);        
        if(eval) {
            tester.setTrace(true);
            lk::Parameters evalAt(npar,0);
            evalAt[0] = radius;
            tester.evaluate(evalAt);
            lk::Gradient grad(npar);
            tester.evaluateGradient(evalAt,grad);
        }
        tester.setTrace(trace);

        // Since our function object has internal state (call counters), we
        // want to ensure that the original object is passed around and
        // never copied. There are two ways to do this:
        lk::FunctionPtr f;
        if(true) {
            // First, using boost::ref with the tester operator() method
            f.reset(new lk::Function(boost::ref(tester)));
        }
        else {
            // Second, bind the object's evaluate() method with a pointer to tester
            f.reset(new lk::Function(
                boost::bind(&test::TestLikelihood::evaluate,&tester,_1)));
        }

        // Specify the different precision values to use for each trial.
        std::vector<double> precision;
        precision.push_back(1e-1);
        precision.push_back(1e-2);
        precision.push_back(1e-3);
        precision.push_back(1e-4);

        // Use fixed initial error estimates.
        lk::Parameters errors(npar,1);

        // Loop over minimization trials.
        for(int trial = 0; trial < ntrial; ++trial) {
            // Choose a random point on a sphere (in npar dimensions)
            // for the initial parameter values.
            lk::Parameters initial(randomOnSphere());
            for(int i = 0; i < initial.size(); ++i) initial[i] *= radius;
           // Loop over precision goals.
            for(int precIndex = 0; precIndex < precision.size(); ++precIndex) {
                double precValue(precision[precIndex]);
                // Use methods that do not use the function gradient.
#ifdef HAVE_LIBGSL
                useMethod(1,"gsl::nmsimplex2",f,initial,errors,precValue);
                useMethod(2,"gsl::nmsimplex2rand",f,initial,errors,precValue);
#endif
#ifdef HAVE_LIBMINUIT2
                useMethod(3,"mn::simplex",f,initial,errors,precValue);
                useMethod(4,"mn::vmetric",f,initial,errors,precValue);
                useMethod(5,"mn::vmetric_fast",f,initial,errors,precValue);
#endif
                // Use methods that require a gradient calculator.
#ifdef HAVE_LIBGSL
                //useMethod(6,"gsl::conjugate_fr",f,gc,initial,errors,precValue);
#endif
            }
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
