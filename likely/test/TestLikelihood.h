// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TEST_TEST_LIKELIHOOD
#define LIKELY_TEST_TEST_LIKELIHOOD

#include "likely/types.h"

#include <vector>

namespace likely {
namespace test {
	class TestLikelihood {
	public:
	    // Creates a test likelihood function with the specified number of parameters.
	    // The common variance of all parameters is sigma^2.
	    // Use -1 < rho < +1 to specify a common correlation coefficient between all
	    // parameters. Use alpha != 0 to introduce a non-linear mapping of the input
	    // parameters to the internal Gaussian parameters.
		TestLikelihood(int npar, double sigma = 1, double rho = 0, double alpha = 0);
		virtual ~TestLikelihood();
        double operator()(std::vector<double> const &params) const;
	private:
        int _npar;
        double _sigma, _rho, _alpha;
        double _inverseDiagonal, _inverseOffDiagonal, _norm;
	}; // TestLikelihood
}} // likely::test

#endif // LIKELY_TEST_TEST_LIKELIHOOD
