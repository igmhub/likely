// Created 7-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates and tests the Integrator class.

#include "likely/likely.h"

#include <cmath>
#include <iostream>

namespace lk = likely;

double integrand1(double x) { return std::log(x)/(x*x); }

double integrand2(double x) { return std::log(x)/std::sqrt(x); }

int main(int argc, char **argv) {

    // Create function wrappers.
    lk::Integrator::IntegrandPtr
        int1(new lk::Integrator::Integrand(integrand1)),
        int2(new lk::Integrator::Integrand(integrand2));

    // Create integrators.
    double epsAbs(1e-7), epsRel(0);
    lk::Integrator integrator1(int1,epsAbs,epsRel), integrator2(int2,epsAbs,epsRel);

    double exact1 = (1 - std::log(2))/2;
    double result1 = integrator1.integrateSmooth(1,2);
    std::cout << "I1[1,2]: " << exact1 << " (exact) - " << result1 << " +/- "
        << integrator1.getAbsError() << " (estimate) = "
        << exact1 - result1 << std::endl;

    double exact2 = (2 - std::log(3))/3;
    double result2 = integrator1.integrateSmooth(1,3);
    std::cout << "I1[2,3]: " << exact2 << " (exact) - " << result2 << " +/- "
        << integrator1.getAbsError() << " (estimate) = "
        << exact2 - result2 << std::endl;
}
