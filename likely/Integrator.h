// Created 07-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_INTEGRATOR
#define LIKELY_INTEGRATOR

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

namespace likely {
    // Implements one-dimensional numerical integration algorithms.
	class Integrator {
	public:
        typedef boost::function<double (double)> Integrand;
        typedef boost::shared_ptr<Integrand> IntegrandPtr;
        // Creates a new integrator of the specified integrand.
		Integrator(IntegrandPtr integrand);
		virtual ~Integrator();
	private:
        IntegrandPtr _integrand;
        class Implementation;
        boost::scoped_ptr<Implementation> _pimpl;
	}; // Integrator
} // likely

#endif // LIKELY_INTEGRATOR
