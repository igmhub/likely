// Created 20-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

// Uses the likely package to minimize a global C-style function that does not provide
// any function gradient info.

#include "config.h"
#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

// Returns |p|^2 for an arbitrary-sized input parameter vector.

double parabola(std::vector<double> const &parameters) {
    double result(0);
    for(int pIndex = 0; pIndex < parameters.size(); ++pIndex) {
        result += parameters[pIndex]*parameters[pIndex];
    }
    return result;
}

int main(int argc, char **argv) {

    // Specify the number of parameters to use.
    int npar = 4;

    // Specify our initial guess at the minimum location and errors.
    lk::Parameters initial(npar,1), errors(npar,1);

    // Results for all algorithms are returned as a FunctionMinimum object.
    lk::FunctionMinimumPtr fmin;

    // Adapt the parabola global function above to what findMinimum(...) expects.
    lk::FunctionPtr fptr(new lk::Function(parabola));

#ifdef HAVE_LIBGSL
    // Estimate the parabola minimum using a GSL algorithm.
    fmin = lk::findMinimum(fptr,initial,errors,"gsl::nmsimplex2");
    std::cout << "=== GSL" << std::endl;
    fmin->printToStream(std::cout);
#endif

#ifdef HAVE_LIBMINUIT2
    // Estimate the parabola minimum and parameter errors using a Minuit algorithm.
    fmin = lk::findMinimum(fptr,initial,errors,"mn2::vmetric");
    std::cout << "=== Minuit2" << std::endl;
    fmin->printToStream(std::cout);
#endif

    // Estimate the parabola minimum and parameter errors using a Markov chain MC
    fmin = lk::findMinimum(fptr,initial,errors,"mc::stroll");
    std::cout << "=== Markov chain MC" << std::endl;
    fmin->printToStream(std::cout);
}