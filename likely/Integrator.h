// Created 07-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_INTEGRATOR
#define LIKELY_INTEGRATOR

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

#include <stack>

namespace likely {
    // Implements one-dimensional numerical integration algorithms.
	class Integrator {
	public:
        typedef boost::function<double (double)> Integrand;
        typedef boost::shared_ptr<Integrand> IntegrandPtr;
        // Creates a new integrator of the specified integrand.
		Integrator(IntegrandPtr integrand, double epsAbs, double epsRel);
		virtual ~Integrator();
		// Returns the integral over an interval [a,b] where the integrand is smooth
		// and non-singular. Updates the value returned by getAbsError(). Uses GSL QAG.
        double integrateSmooth(double a, double b);
		// Returns the integral over an interval [a,b] where the integrand is singular
		// at the endpoints and/or interior points. Updates the value returned by
		// getAbsError(). Uses GSL QAGS.
        double integrateSingular(double a, double b);
        // Returns the integral over an interval [a,b] using a robust but slower method
        // that can handle singularities and inf,nan values. Updates the value returned
        // by getAbsError(). Uses GSL CQUAD (added in GSL version 1.15)
        double integrateRobust(double a, double b);
        // Returns the integral from [a,+infinity). Updates the value returned by
        // getAbsError(). Uses GSL QAGIU.
        double integrateUp(double a);
        // Returns the integral from (-infinity,b]. Updates the value returned by
        // getAbsError(). Uses GSL QAGIL.
        double integrateDown(double b);
        // Returns the integral from (-infinity,+infinity). Updates the value return
        // by getAbsError(). Uses GSL QAGI.
        double integrateAll();
        // Returns the integral of integrand(x)*osc(omega*x) where osc = sin or cos.
        // Updates the value return by getAbsError(). Uses GSL QAWO.
        double integrateOsc(double a, double b, double omega, bool useSin = true);
        // Returns the estimated absolute error from the last integration or zero if
        // no integrations have been performed yet.
        double getAbsError() const;
	private:
        IntegrandPtr _integrand;
        double _epsAbs, _epsRel, _absError;
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
        // Global C-style callbacks that evaluate the top engine on the stack.
        static double _evaluate(double x, void *params);
        // Maintains an integrand stack that the global callback uses to determine which
        // function to invoke when called.
        static std::stack<const Integrator*> &_getStack();
	}; // Integrator
	
    inline double Integrator::getAbsError() const { return _absError; }

} // likely

#endif // LIKELY_INTEGRATOR
