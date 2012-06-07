// Created 6-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// CovarianceMatrix class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

BOOST_AUTO_TEST_SUITE( CovarianeMatrix )

BOOST_AUTO_TEST_CASE( shouldMakeCovarianceMatrixWithCorrectSize ) {
	int nelem(10);
	lk::CovarianceMatrix cov2(std::vector<double>(nelem,1));
	int cov2Size = cov2.getSize();
	int expected(4);
	BOOST_CHECK_EQUAL(expected, cov2Size);
}

BOOST_AUTO_TEST_CASE( shouldConvertMatrixIndexToPackedArrayIndex ) {
	int expected(4);
	int result(lk::symmetricMatrixIndex(1,2,3));
	BOOST_CHECK_EQUAL(expected, result);
	result = lk::symmetricMatrixIndex(2,1,3);
	BOOST_CHECK_EQUAL(expected, result);
}

BOOST_AUTO_TEST_CASE( shouldCalculateSizeOfMatrixFromNumberOfElements ) {
	int expected(5);
	int result(lk::symmetricMatrixSize(15));
	BOOST_CHECK_EQUAL(expected, result);
}

BOOST_AUTO_TEST_SUITE_END()