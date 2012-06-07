// Created 6-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// UniformBinning class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

// TODO: implement common setup here
struct UniformBinningTestFixture
{
    UniformBinningTestFixture() {
    	BOOST_TEST_MESSAGE("Setup UniformBinningTestFixture...");
    }
    
    ~UniformBinningTestFixture() {
    	BOOST_TEST_MESSAGE("Teardown UniformBinningTestFixture...");
    }
        
};

BOOST_FIXTURE_TEST_SUITE( UniformBinningTest, UniformBinningTestFixture )

BOOST_AUTO_TEST_CASE( shouldThrowBinningErrorInConstructorWithInvalidParameters ) {
	// nBins == 0, should throw BinningError
	BOOST_CHECK_THROW(lk::UniformBinning(0.,1.,0), lk::BinningError);
	// nBins < 0, should throw BinningError
	BOOST_CHECK_THROW(lk::UniformBinning(0.,1.,-1), lk::BinningError);
	// maxValue == minValue, should throw BinningError
	BOOST_CHECK_THROW(lk::UniformBinning(2.3,2.3,5), lk::BinningError);
	// maxValue < minValue, should throw BinningError
	BOOST_CHECK_THROW(lk::UniformBinning(10,0.,5), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectBinIndexOrThrowErrorForOutOfRangeIndex ) {
	lk::UniformBinning axis(0.,5.,20);
	int expected(7);
	int result = axis.getBinIndex(1.9);
	BOOST_CHECK_EQUAL(result, expected);
	
	// value out of range, should throw BinningError
	BOOST_CHECK_THROW(axis.getBinIndex(-0.1), lk::BinningError);
	BOOST_CHECK_THROW(axis.getBinIndex(100.), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectNumberOfBins ) {
	lk::UniformBinning axis(0.,5.,20);
	int expected(20);
	int result = axis.getNBins();
	BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE( shouldReturnFalseOrThrowErrorForIndexOutOfRange ) {
	lk::UniformBinning axis(0.,5.,20);
	// index out of range, return false
	BOOST_CHECK(not axis.isValidBinIndex(-1));
	BOOST_CHECK(not axis.isValidBinIndex(20));
	// index out of range, w/ specified errorFormat, should throw BinningError
	BOOST_CHECK_THROW(axis.isValidBinIndex(-1,"invalid bin index: %d"), lk::BinningError);
	BOOST_CHECK_THROW(axis.isValidBinIndex(20,"invalid bin index: %d"), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnLowerBoundOfSpecifiedBin ) {
	lk::UniformBinning axis(0.,5.,20);
	double expectedLow(1.75);
	double resultLow = axis.getBinLowEdge(7);
	BOOST_CHECK_EQUAL(resultLow, expectedLow);
	double expectedHigh(2.);
	double resultHigh = axis.getBinHighEdge(7);
	BOOST_CHECK_EQUAL(resultHigh, expectedHigh);
}

BOOST_AUTO_TEST_CASE( shouldReturnBinWidthOfSpecifiedBin ) {
	lk::UniformBinning axis(0.,5.,20);
	double expected(.25);
	double result = axis.getBinWidth(7);
	BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE( shouldReturnMidpointValueOfSpecifiedBin ){
	lk::UniformBinning axis(0.,5.,20);
	double expected(1.875);
	double result = axis.getBinCenter(7);
	BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END() // UniformBinningTest