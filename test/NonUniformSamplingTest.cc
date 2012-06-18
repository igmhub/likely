// Created 8-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// NonUniformSampling class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

struct NonUniformSamplingFixture
{
    NonUniformSamplingFixture() {
    	std::vector<double> points(7);
    	points[0] = -.2; points[1] = -.05; points[2] = 0.; points[3] = 1.;
    	points[4] = 1.35; points[5] = 5.; points[6] = 10.;
    	axis.reset(new lk::NonUniformSampling(points));
    }
    ~NonUniformSamplingFixture() {
    }  
    lk::AbsBinningCPtr axis;
};

BOOST_FIXTURE_TEST_SUITE( NonUniformSampling, NonUniformSamplingFixture )

BOOST_AUTO_TEST_CASE( shouldThrowBinningErrorInConstructorWithInvalidParameters ) {
	// nSamples < 3, should throw BinningError
	BOOST_CHECK_THROW(axis.reset(new lk::NonUniformSampling(std::vector<double>(2))), lk::BinningError);
	// points are not in increasing order, should throw BinningError
    std::vector<double> points(3);
    points[0] = 0; points[1] = 13.7; points[2] = 4.54;
	BOOST_CHECK_THROW(axis.reset(new lk::NonUniformSampling(points)), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectBinIndexOrThrowError ) {
	// check correct bin
	BOOST_CHECK_EQUAL(axis->getBinIndex(1.35), 4);
	// check precision
	BOOST_CHECK_EQUAL(axis->getBinIndex(1.35020), 4);
	BOOST_CHECK_THROW(axis->getBinIndex(1.35021), lk::BinningError);
	BOOST_CHECK_THROW(axis->getBinIndex(1.34979), lk::BinningError);
	// value out of range, should throw BinningError
	BOOST_CHECK_THROW(axis->getBinIndex(10.1), lk::BinningError);
	BOOST_CHECK_THROW(axis->getBinIndex(-3), lk::BinningError);
	// value out of range but within precision, should be OK
	BOOST_CHECK_EQUAL(axis->getBinIndex(10.0001), 6);
	BOOST_CHECK_EQUAL(axis->getBinIndex(-.20001), 0);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectNumberOfBins ) {
	BOOST_CHECK_EQUAL(axis->getNBins(), 7);
}

BOOST_AUTO_TEST_CASE( shouldReturnFalseOrThrowErrorForIndexOutOfRange ) {
	// index out of range, return false
	BOOST_CHECK(!axis->isValidBinIndex(-1));
	BOOST_CHECK(!axis->isValidBinIndex(7));
	// index out of range, w/ specified errorFormat, should throw BinningError
	BOOST_CHECK_THROW(axis->isValidBinIndex(-1,"invalid bin index: %d"), lk::BinningError);
	BOOST_CHECK_THROW(axis->isValidBinIndex(7,"invalid bin index: %d"), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnZeroAsWidthOfSpecifiedBin ) {
	BOOST_CHECK_EQUAL(axis->getBinWidth(4), 0);
}

BOOST_AUTO_TEST_CASE( shouldReturnValueOfSpecifiedBin ){
	BOOST_CHECK_EQUAL(axis->getBinCenter(4), 1.35);
	BOOST_CHECK_EQUAL(axis->getBinLowEdge(4), 1.35);
	BOOST_CHECK_EQUAL(axis->getBinHighEdge(4), 1.35);
}

BOOST_AUTO_TEST_SUITE_END() // NonUniformSampling