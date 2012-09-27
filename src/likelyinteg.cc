// Created 7-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates and tests the Integrator class.

#include "likely/Integrator.h"
#include "likely/RuntimeError.h"

#include <cmath>
#include <iostream>

namespace lk = likely;

double integrand1(double x) { return std::log(x)/(x*x); }

double integrand2(double x) { return std::log(x)/std::sqrt(x); }

double integrand3(double x) { return 1./0.; }

int main(int argc, char **argv) {

    // Create function wrappers.
    lk::Integrator::IntegrandPtr
        int1(new lk::Integrator::Integrand(integrand1)),
        int2(new lk::Integrator::Integrand(integrand2)),
        int3(new lk::Integrator::Integrand(integrand3));

    // Create integrators.
    double epsAbs(1e-7), epsRel(0);
    lk::Integrator integrator1(int1,epsAbs,epsRel), integrator2(int2,epsAbs,epsRel),
         integrator3(int3,epsAbs,epsRel);

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
        
    double exact3 = 1;
    double result3 = integrator1.integrateUp(1);
    std::cout << "I1[1,+inf): " << exact3 << " (exact) - " << result3 << " +/- "
        << integrator1.getAbsError() << " (estimate) = "
        << exact3 - result3 << std::endl;

    double exact4 = -4;
    double result4 = integrator2.integrateSingular(0,1);
    std::cout << "I2[0,1]: " << exact4 << " (exact) - " << result4 << " +/- "
        << integrator2.getAbsError() << " (estimate) = "
        << exact4 - result4 << std::endl;
    
    double exact5 = -4;
    double result5 = integrator2.integrateRobust(0,1);
    std::cout << "I2[0,1]: " << exact5 << " (exact) - " << result5 << " +/- "
        << integrator2.getAbsError() << " (estimate) = "
        << exact5 - result5 << std::endl;
    
    double exact6 = -0.060498064358031995849;
    double result6 = integrator1.integrateOsc(1,15,3);
    std::cout << "I1[1,15]*sin(3x): " << exact6 << " (exact) - " << result6 << " +/- "
        << integrator2.getAbsError() << " (estimate) = "
        << exact6 - result6 << std::endl;
    
    try {
        integrator3.integrateSingular(0,1);
    }
    catch(lk::RuntimeError const &err) {
        std::cout << "Integral of 1/0 failed as expected." << std::endl;
    }
}
