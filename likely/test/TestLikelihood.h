// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TEST_TEST_LIKELIHOOD
#define LIKELY_TEST_TEST_LIKELIHOOD

#include "likely/types.h"

#include "boost/utility.hpp"

#include <vector>

namespace likely {
namespace test {
	class TestLikelihood : public boost::noncopyable {
	public:
	    // Creates a test likelihood function with the specified number of parameters.
	    // The common variance of all parameters is sigma^2.
	    // Use -1 < rho < +1 to specify a common correlation coefficient between all
	    // parameters. Use alpha != 0 to introduce a non-linear mapping of the input
	    // parameters to the internal Gaussian parameters.
	    // The function has a single global minimum at the origin where the function
	    // value is zero.
		TestLikelihood(int npar, double sigma = 1, double rho = 0, double alpha = 0);
		virtual ~TestLikelihood();
		// Evaluates this likelihood function at the specified parameter values.
        double operator()(Parameters const &params) const;
        // Turns evaluation tracing on/off.
        void setTrace(bool value);
        // Returns the number of function evaluations since this object was created
        // or resetCount() was called.
        long getCount() const;
        // Resets the evaluation counter.
        void resetCount();
	private:
        int _npar;
        mutable long _count;
        double _sigma, _rho, _alpha;
        double _inverseDiagonal, _inverseOffDiagonal;
        bool _trace;
	}; // TestLikelihood
	
    inline void TestLikelihood::setTrace(bool value) { _trace = value; }
    inline long TestLikelihood::getCount() const { return _count; }
    inline void TestLikelihood::resetCount() { _count = 0; }

}} // likely::test

#endif // LIKELY_TEST_TEST_LIKELIHOOD
