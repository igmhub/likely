// Created 27-Aug-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// ExactQuantileAccumulator class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/assign.hpp>

#include "likely/likely.h"

namespace lk = likely;

struct ExactQuantileAccumulatorFixture
{
    ExactQuantileAccumulatorFixture() {
    }
    ~ExactQuantileAccumulatorFixture() { }
    lk::ExactQuantileAccumulator q;
};

BOOST_FIXTURE_TEST_SUITE( ExactQuantileAccumulator, ExactQuantileAccumulatorFixture )

BOOST_AUTO_TEST_CASE( calculateCorrectQuantiles ) {
	q.accumulate(3);
	q.accumulate(6);
	q.accumulate(7);
	q.accumulate(8);
	q.accumulate(8);
	q.accumulate(9);
	q.accumulate(10);
	q.accumulate(13);
	q.accumulate(15);
	q.accumulate(16);
	q.accumulate(20);
	BOOST_REQUIRE_EQUAL(q.getQuantile(0.25),7);
	BOOST_REQUIRE_EQUAL(q.getQuantile(0.5),9);
	BOOST_REQUIRE_EQUAL(q.getQuantile(0.75),15);
}

BOOST_AUTO_TEST_CASE( shouldThrowErrorWhenAccumulatingNegativeWeight ) {
	BOOST_CHECK_THROW(q.accumulate(3,-2),lk::RuntimeError);
}

BOOST_AUTO_TEST_CASE( shouldThrowErrorWhenGetQuantileIsCalledBeforeAnyAccumulation ) {
	BOOST_CHECK_THROW(q.getQuantile(.5),lk::RuntimeError);
}

BOOST_AUTO_TEST_CASE( shouldThrowErrorWhenIfQuantileProbabilityIsOutOfRange ) {
	q.accumulate(3);
	BOOST_CHECK_THROW(q.getQuantile(-.1),lk::RuntimeError);
	BOOST_CHECK_THROW(q.getQuantile(1.01),lk::RuntimeError);
}

BOOST_AUTO_TEST_SUITE_END() // ExactQuantileAccumulator