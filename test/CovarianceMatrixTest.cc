// Created 6-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// CovarianceMatrix class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

struct CovarianceMatrixFixture
{
    CovarianceMatrixFixture() {
    	size = 3;
    	cov.reset(new lk::CovarianceMatrix(size));
		for(int k = 0; k < size; ++k) {
			cov->setCovariance(k,k,k+1);
		}
		cov->setCovariance(0,1,0.1);
		cov->setCovariance(1,2,-0.2);
	}
    ~CovarianceMatrixFixture() {
    }
    int size;
    boost::shared_ptr<lk::CovarianceMatrix> cov;
};


BOOST_FIXTURE_TEST_SUITE( CovarianeMatrix, CovarianceMatrixFixture )

BOOST_AUTO_TEST_CASE( shouldMakeCovarianceMatrixWithCorrectSize ) {
	int nelem(10);
	lk::CovarianceMatrix myCov(std::vector<double>(nelem,1));
	BOOST_CHECK_EQUAL(myCov.getSize(), 4);
}

BOOST_AUTO_TEST_CASE( shouldThrowErrorSettingInvalidMatrixElement ) {
	BOOST_CHECK_THROW(cov->setCovariance(0,0,-2), lk::RuntimeError);
	BOOST_CHECK_THROW(cov->setInverseCovariance(1,1,0), lk::RuntimeError);
}

// =, swap
// getCovariance, getInverseCovariance
// setCovariance, setInverseCovariance
// multiplyByCovariance, multiplyByInverseCovariance
// chisquare
// replaceWithTripleProduct
// addInverse
// sample
// prune
// compress
// choleskyDecompose
// invertCholesky

BOOST_AUTO_TEST_CASE( shouldConvertMatrixIndexToPackedArrayIndex ) {
	BOOST_CHECK_EQUAL(lk::symmetricMatrixIndex(1,2,3), 4);
	BOOST_CHECK_EQUAL(lk::symmetricMatrixIndex(2,1,3), 4);
}

BOOST_AUTO_TEST_CASE( shouldCalculateSizeOfMatrixFromNumberOfElements ) {
	BOOST_CHECK_EQUAL(lk::symmetricMatrixSize(15), 5);
}

BOOST_AUTO_TEST_CASE( shouldMultiplyVectorByCovaraianceMatrix ) {
	std::vector<double> vec(size);
	vec[0] = 1; vec[1] = 2; vec[2] = 3;
	cov->multiplyByCovariance(vec);
	BOOST_CHECK_CLOSE(vec[0], 1.2, 1e-6);
	BOOST_CHECK_CLOSE(vec[1], 3.5, 1e-6);
	BOOST_CHECK_CLOSE(vec[2], 8.6, 1e-6);
	
}

BOOST_AUTO_TEST_SUITE_END()