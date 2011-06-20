// Created 20-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

// Uses the likely package to minimize a C++ likelihood class that provides
// function gradient info.

#include "config.h"
#include "likely/likely.h"

#include "boost/bind.hpp"

#include <iostream>

namespace lk = likely;

class Parabola {
public:
    Parabola(double constant) : _constant(constant) { }
    typedef std::vector<double> Parameters;
    // Returns constant+|p|^2 for an arbitrary-sized input parameter vector.
    double operator()(Parameters const &parameters) const {
        double result = _constant;
        for(int pIndex = 0; pIndex < parameters.size(); ++pIndex) {
            result += parameters[pIndex]*parameters[pIndex];
        }
        return result;        
    }
    // Saves the components of 2*p in the gradient vector provided.
    void calculateGradient(Parameters const &parameters, Parameters &gradients) const {
        for(int pIndex = 0; pIndex < parameters.size(); ++pIndex) {
            gradients[pIndex] = 2*parameters[pIndex];
        }
    }
private:
    double _constant;
};

int main(int argc, char **argv) {

    // Specify the number of parameters to use.
    int npar = 4;

    // Specify our initial guess at the minimum location and errors.
    lk::Parameters initial(npar,1), errors(npar,1);

    // Results for all algorithms are returned as a FunctionMinimum object.
    lk::FunctionMinimumPtr fmin;

    // Create an instance of the Parabola class above.
    Parabola parabola(1.23);

    // Adapt the parabola object's operator() and calculateGradient() to
    // what findMinimum(...) expects.
    lk::FunctionPtr fptr(new lk::Function(parabola));
    lk::GradientCalculatorPtr gcptr(new lk::GradientCalculator(
        boost::bind(&Parabola::calculateGradient,&parabola,_1,_2)));

#ifdef HAVE_LIBGSL
    // Estimate the parabola minimum using a GSL algorithm.
    fmin = lk::findMinimum(fptr,gcptr,initial,errors,"gsl::conjugate_pr");
    std::cout << "=== GSL" << std::endl;
    fmin->printToStream(std::cout);
#endif

#ifdef HAVE_LIBMINUIT2
    // Estimate the parabola minimum and parameter errors using a Minuit algorithm.
    fmin = lk::findMinimum(fptr,gcptr,initial,errors,"mn2::vmetric_grad");
    std::cout << "=== Minuit2" << std::endl;
    fmin->printToStream(std::cout);
#endif
}