#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

BOOST_AUTO_TEST_SUITE( covariance_matrix )

BOOST_AUTO_TEST_CASE( shouldMakeCovarianceMatrixWithCorrectSize ) {
	int nelem(10);
	lk::CovarianceMatrix cov2(std::vector<double>(nelem,1));
	int cov2Size = cov2.getSize();
	int expected(4);
	BOOST_CHECK_EQUAL(expected, cov2Size);
}

BOOST_AUTO_TEST_CASE( shouldConvertMatrixIndexToPackedArrayIndex ) {
	int expected(4);
	BOOST_CHECK_EQUAL(expected, lk::symmetricMatrixIndex(1,2,3));
	BOOST_CHECK_EQUAL(expected, lk::symmetricMatrixIndex(2,1,3));
}

BOOST_AUTO_TEST_CASE( shouldCalculateSizeOfMatrixFromNumberOfElements ) {
	int expected(5);
	BOOST_CHECK_EQUAL(expected, lk::symmetricMatrixSize(15));
}

BOOST_AUTO_TEST_SUITE_END()