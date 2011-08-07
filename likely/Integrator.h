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
		// and non-singular. Updates the value returned by getAbsError().
        double integrateSmooth(double a, double b);
		// Returns the integral over an interval [a,b] where the integrand is singular
		// at the endpoints and/or interior points. Updates the value returned by
		// getAbsError().
        double integrateSingular(double a, double b);
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
