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

namespace lk = likely;
namespace test = likely::test;
namespace po = boost::program_options;

void useMethod(int methodId, std::string const &methodName,
lk::FunctionPtr f, lk::GradientCalculatorPtr gc, lk::Parameters const &initial,
lk::Parameters const &errors, double prec) {
    lk::FunctionMinimumPtr fmin;
    long maxIterations(1000000); //1e6
    if(gc) {
        // Use an algorithm that requires a gradient calculator.
        fmin = lk::findMinimum(f,gc,initial,errors,methodName,prec,maxIterations);
    }
    else {
        fmin = lk::findMinimum(f,initial,errors,methodName,prec,maxIterations);
    }
    double error(0);
    if(fmin->haveCovariance()) error = fmin->getErrors()[0];
    boost::format fmt("%d %.4f %.4f %.4f %.4f\n");
    std::cout << fmt % methodId % std::log10((double)lk::lastMinEvalCount)
        % (lk::lastMinGradCount ? std::log10((double)lk::lastMinGradCount) : 0.)
        % (fmin->getMinValue() > 0 ? -std::log10(fmin->getMinValue()) : 0) % error;
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
        ("ntrial", po::value<int>(&ntrial)->default_value(20),
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
    lk::Random &random(lk::Random::instance());
    random.setSeed(seed);
    boost::uniform_on_sphere<> spherical(npar);
    boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<> >
        randomOnSphere(random.getGenerator(),spherical);
    
    // Print out column headings for the output we generate below.
    std::cout << "method ncall ngrad accuracy error" << std::endl;

    try {
        // Create a likelihood function using the command-line parameters.
        test::TestLikelihood tester(npar,1,rho,alpha);
        
        // Trace a single function and gradient evaluation, if requested.
        if(eval) {
            tester.setTrace(true);
            lk::Parameters evalAt(npar,radius);
            tester.evaluate(evalAt);
            lk::Gradient grad(npar);
            tester.evaluateGradient(evalAt,grad);
        }
        tester.setTrace(trace);

        // Since our function object has internal state (call counters), we
        // want to ensure that the original object is passed around and
        // never copied. Here are two ways to do this:
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
        // Use boost::bind for the gradient calculator.
        lk::GradientCalculatorPtr gc(new lk::GradientCalculator(
            boost::bind(&test::TestLikelihood::evaluateGradient,&tester,_1,_2)));
        // Initialize a null gradient calculator for methods that don't need one.
        lk::GradientCalculatorPtr noGC;

        // Specify the different precision values to use for each trial.
        std::vector<double> precision;
        precision.push_back(1e-1);
        precision.push_back(1e-2);
        precision.push_back(1e-3);
        precision.push_back(1e-4);

        // Use fixed initial error estimates of 1.3 for each parameter, which
        // overestimate the true error, which is exactly 1 for alpha = 0.
        lk::Parameters errors(npar,1.3);

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
                useMethod(1,"gsl::nmsimplex2",f,noGC,initial,errors,precValue);
                useMethod(2,"gsl::nmsimplex2rand",f,noGC,initial,errors,precValue);
#endif
#ifdef HAVE_LIBMINUIT2
                useMethod(3,"mn2::simplex",f,noGC,initial,errors,precValue);
                useMethod(4,"mn2::vmetric",f,noGC,initial,errors,precValue);
                useMethod(5,"mn2::vmetric_fast",f,noGC,initial,errors,precValue);
#endif
                useMethod(6,"mc::saunter",f,noGC,initial,errors,precValue);
                useMethod(7,"mc::stroll",f,noGC,initial,errors,precValue);
                // Use methods that require a gradient calculator.
#ifdef HAVE_LIBGSL
                useMethod(11,"gsl::conjugate_fr",f,gc,initial,errors,precValue);
                useMethod(12,"gsl::conjugate_pr",f,gc,initial,errors,precValue);
                useMethod(13,"gsl::vector_bfgs2",f,gc,initial,errors,precValue);
                useMethod(14,"gsl::steepest_descent",f,gc,initial,errors,precValue);
#endif
#ifdef HAVE_LIBMINUIT2
                useMethod(15,"mn2::vmetric_grad",f,gc,initial,errors,precValue);
                useMethod(16,"mn2::vmetric_grad_fast",f,gc,initial,errors,precValue);
#endif
            }
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
