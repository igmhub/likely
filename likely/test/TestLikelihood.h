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
	    // Creates a negative log-likelihood (NLL) function with the specified number
	    // of parameters, all having the common variance sigma^2.
	    // Use -1 < rho < +1 to specify a common correlation coefficient between all
	    // parameters. Use alpha != 0 to introduce a non-linear mapping of the input
	    // parameters to the internal Gaussian parameters.
	    // The function has a single global minimum at the origin where the function
	    // value is zero.
		TestLikelihood(int npar, double sigma = 1, double rho = 0, double alpha = 0);
		virtual ~TestLikelihood();
		// Evaluates this NLL function at the specified parameter values.
        double evaluate(Parameters const &params) const;
        double operator()(Parameters const &params) const;
        // Evaluates this NLL function's gradients at the specified parameter values.
        void evaluateGradient(Parameters const &params, Gradient &grad) const;
        // Turns evaluation tracing on/off.
        void setTrace(bool value);
	private:
        int _npar;
        double _sigma, _rho, _alpha, _delta;
        double _inverseDiagonal, _inverseOffDiagonal;
        bool _trace;
        // Applies the transformation y_i(x) = x_i - alpha*(|x|^2 - x_i^2)
        Parameters _transform(Parameters const &x) const;
	}; // TestLikelihood
	
    inline void TestLikelihood::setTrace(bool value) { _trace = value; }

    inline double TestLikelihood::operator()(Parameters const &params) const {
        return evaluate(params);
    }

}} // likely::test

#endif // LIKELY_TEST_TEST_LIKELIHOOD
