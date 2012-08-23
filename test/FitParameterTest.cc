// Created 23-Aug-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// FitParameter class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/assign.hpp>

#include <string>

#include "likely/likely.h"


namespace lk = likely;

struct FitParameterFixture
{
    FitParameterFixture() {
    }
    ~FitParameterFixture() { }
    double value;
    std::vector<double> errors; 
};

BOOST_FIXTURE_TEST_SUITE( FitParameter, FitParameterFixture )

BOOST_AUTO_TEST_CASE( shouldFormatValueWithError ) {
	errors.push_back(.547);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(105.658, errors), "105.7 +/- 0.5");
	errors.push_back(.12);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(105.658, errors), "105.66 +/- 0.55 +/- 0.12");
	errors.push_back(.0097);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(105.658, errors), "105.658 +/- 0.547 +/- 0.120 +/- 0.010");
}

BOOST_AUTO_TEST_CASE( shouldFormatNegativeValueWithError ) {
	errors.push_back(.547);
	errors.push_back(.12);
	errors.push_back(.0097);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(-105.658, errors), "-105.658 +/- 0.547 +/- 0.120 +/- 0.010");
}

BOOST_AUTO_TEST_CASE( shouldInsertCustomSeperator ) {
	errors.push_back(.5);
	errors.push_back(.12);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(105.658, errors, "\\pm"), "105.66 \\pm 0.50 \\pm 0.12");
}

BOOST_AUTO_TEST_CASE( shouldFormatAnyMagnitude ) {
	errors.push_back(500);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(987654.321, errors), "987700 +/- 500");
	errors.push_back(12);
	BOOST_REQUIRE_EQUAL(lk::roundValueWithError(987654.321, errors), "987654 +/- 500 +/- 12");
}



BOOST_AUTO_TEST_SUITE_END() // FitParameter