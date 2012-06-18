// Created 7-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// UniformSampling class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

struct UniformSamplingFixture
{
    UniformSamplingFixture() {
    	axis.reset(new lk::UniformSampling(0.,10.,7));
    }
    ~UniformSamplingFixture() {
    }   
    lk::AbsBinningCPtr axis;
};

BOOST_FIXTURE_TEST_SUITE( UniformSampling, UniformSamplingFixture )

BOOST_AUTO_TEST_CASE( shouldThrowBinningErrorInConstructorWithInvalidParameters ) {
	// maxValue <= minValue and nSamples > 1, should throw BinningError
	BOOST_CHECK_THROW(lk::UniformSampling(105.7,105.7,2), lk::BinningError);
	BOOST_CHECK_THROW(lk::UniformSampling(105.7,.511,2), lk::BinningError);
	// maxValue != minValue and nSamples == 1, should throw BinningError
	BOOST_CHECK_THROW(lk::UniformSampling(.511,105.7,1), lk::BinningError);
	BOOST_CHECK_THROW(lk::UniformSampling(105.7,.511,1), lk::BinningError);
}

BOOST_AUTO_TEST_CASE( shouldReturnCorrectBinIndexOrThrowError ) {
	// check correct bin
	BOOST_CHECK_EQUAL(axis->getBinIndex(8.3333), 5);
	// check precision
	BOOST_CHECK_EQUAL(axis->getBinIndex(8.33349), 5);
	BOOST_CHECK_THROW(axis->getBinIndex(8.333), lk::BinningError);
	BOOST_CHECK_THROW(axis->getBinIndex(8.3335), lk::BinningError);
	// value out of range, should throw BinningError
	BOOST_CHECK_THROW(axis->getBinIndex(10.1), lk::BinningError);
	BOOST_CHECK_THROW(axis->getBinIndex(-3), lk::BinningError);
	// value out of range but within precision, should be OK
	BOOST_CHECK_EQUAL(axis->getBinIndex(10.0001), 6);
	BOOST_CHECK_EQUAL(axis->getBinIndex(-0.0001), 0);
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
	BOOST_CHECK_EQUAL(axis->getBinCenter(4), 20./3);
	BOOST_CHECK_EQUAL(axis->getBinLowEdge(4), 20./3);
	BOOST_CHECK_EQUAL(axis->getBinHighEdge(4), 20./3);
}

BOOST_AUTO_TEST_SUITE_END() // UniformSampling