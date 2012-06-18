// Created 6-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// NonUniformBinning class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

struct NonUniformBinningFixture
{
    NonUniformBinningFixture() {
    	std::vector<double> bins(4);
    	bins[0] = 0; bins[1] = 0.25; bins[2] = 0.35; bins[3] = 1;
    	axis.reset(new lk::NonUniformBinning(bins));
    }
    ~NonUniformBinningFixture() { }   
    lk::AbsBinningCPtr axis;
};

BOOST_FIXTURE_TEST_SUITE( NonUniformBinning, NonUniformBinningFixture )

BOOST_AUTO_TEST_CASE( shouldThrowBinningErrorInConstructorWithInvalidParameters ) {
	// nBinEdges < 2, should throw BinningError
	BOOST_CHECK_THROW(axis.reset(new lk::NonUniformBinning(std::vector<double>(1))), lk::BinningError);
	// binEdges are not in increasing order, should throw BinningError
	std::vector<double> bins(3);
	bins[0] = 0; bins[1] = 13.7; bins[2] = 4.54;
	BOOST_CHECK_THROW(axis.reset(new lk::NonUniformBinning(bins)), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectBinIndexOrThrowErrorForOutOfRangeIndex ) {
	BOOST_CHECK_EQUAL(axis->getBinIndex(.75), 2);
	BOOST_CHECK_EQUAL(axis->getBinIndex(.3), 1);
	BOOST_CHECK_EQUAL(axis->getBinIndex(.123456), 0);
	// value out of range, should throw BinningError
	BOOST_CHECK_THROW(axis->getBinIndex(-0.1), lk::BinningError);
	BOOST_CHECK_THROW(axis->getBinIndex(1.1), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectNumberOfBins ) {
	BOOST_CHECK_EQUAL(axis->getNBins(), 3);
}

BOOST_AUTO_TEST_CASE( shouldReturnFalseOrThrowErrorForIndexOutOfRange ) {
	// index out of range, return false
	BOOST_CHECK(!axis->isValidBinIndex(-1));
	BOOST_CHECK(!axis->isValidBinIndex(4));
	// index out of range, w/ specified errorFormat, should throw BinningError
	BOOST_CHECK_THROW(axis->isValidBinIndex(-1,"invalid bin index: %d"), lk::BinningError);
	BOOST_CHECK_THROW(axis->isValidBinIndex(4,"invalid bin index: %d"), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnEdgeOfSpecifiedBin ) {
	BOOST_CHECK_EQUAL(axis->getBinLowEdge(1), .25);
	BOOST_CHECK_EQUAL(axis->getBinHighEdge(1), .35);
	BOOST_CHECK_EQUAL(axis->getBinLowEdge(2), .35);
	BOOST_CHECK_EQUAL(axis->getBinHighEdge(2), 1);
}

BOOST_AUTO_TEST_CASE( shouldReturnBinWidthOfSpecifiedBin ) {
	BOOST_CHECK_EQUAL(axis->getBinWidth(0), .25);
	BOOST_CHECK_EQUAL(axis->getBinWidth(1), .35-.25); //.35-.25 != .1 => [0.099999999999999978 != 0.10000000000000001]
	BOOST_CHECK_EQUAL(axis->getBinWidth(2), .65);
}

BOOST_AUTO_TEST_CASE( shouldReturnMidpointValueOfSpecifiedBin ){
	BOOST_CHECK_EQUAL(axis->getBinCenter(0), .125);
	BOOST_CHECK_EQUAL(axis->getBinCenter(1), .3);
}

BOOST_AUTO_TEST_SUITE_END() // NonUniformBinningTest