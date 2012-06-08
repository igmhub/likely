// Created 6-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// UniformBinning class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

struct UniformBinningFixture
{
    UniformBinningFixture() {
    	axis.reset(new lk::UniformBinning(0.,5.,20));
    }
    ~UniformBinningFixture() { }   
    lk::AbsBinningCPtr axis;
};

BOOST_FIXTURE_TEST_SUITE( UniformBinning, UniformBinningFixture )

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
	BOOST_CHECK_EQUAL(axis->getBinIndex(1.9), 7);
	// value out of range, should throw BinningError
	BOOST_CHECK_THROW(axis->getBinIndex(-0.1), lk::BinningError);
	BOOST_CHECK_THROW(axis->getBinIndex(100.), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectNumberOfBins ) {
	BOOST_CHECK_EQUAL(axis->getNBins(), 20);
}

BOOST_AUTO_TEST_CASE( shouldReturnFalseOrThrowErrorForIndexOutOfRange ) {
	// index out of range, return false
	BOOST_CHECK(!axis->isValidBinIndex(-1));
	BOOST_CHECK(!axis->isValidBinIndex(20));
	// index out of range, w/ specified errorFormat, should throw BinningError
	BOOST_CHECK_THROW(axis->isValidBinIndex(-1,"invalid bin index: %d"), lk::BinningError);
	BOOST_CHECK_THROW(axis->isValidBinIndex(20,"invalid bin index: %d"), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnLowerBoundOfSpecifiedBin ) {
	BOOST_CHECK_EQUAL(axis->getBinLowEdge(7), 1.75);
	BOOST_CHECK_EQUAL(axis->getBinHighEdge(7), 2.);
}

BOOST_AUTO_TEST_CASE( shouldReturnBinWidthOfSpecifiedBin ) {
	BOOST_CHECK_EQUAL(axis->getBinWidth(7), .25);
	BOOST_CHECK_EQUAL(axis->getBinWidth(4), .25);
}

BOOST_AUTO_TEST_CASE( shouldReturnMidpointValueOfSpecifiedBin ){
	BOOST_CHECK_EQUAL(axis->getBinCenter(7), 1.875);
}

BOOST_AUTO_TEST_SUITE_END() // UniformBinningTest