// Created 6-Jun-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// BinnedData class unit tests.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "likely/likely.h"
namespace lk = likely;

struct BinnedDataFixture
{
    BinnedDataFixture() {
    	std::vector<double> bins(4);
    	bins[0] = 0; bins[1] = 0.25; bins[2] = 0.35; bins[3] = 1;    
    	lk::AbsBinningCPtr
        	axis1(new lk::UniformBinning(0.,1.,3)),
        	axis2(new lk::UniformSampling(0.,1.,3)),
        	axis3(new lk::NonUniformBinning(bins));
        binnedData.reset(new lk::BinnedData(axis1, axis2, axis3));
    }
    ~BinnedDataFixture() { }   
    boost::shared_ptr<const lk::BinnedData> binnedData;
};

BOOST_FIXTURE_TEST_SUITE( BinnedData, BinnedDataFixture )

BOOST_AUTO_TEST_CASE( shouldThrowBinningErrorInConstructorWithInvalidParameters ) {
	BOOST_CHECK_EQUAL(1, 1);
}

// clone, =, swap
// +=, add
// isCongruent
// getIndex (binIndices, values)
// getIndexAtOffset, getOffsetForIndex
// getBinIndices
// getBinWidths
// hasData, getData, setData, addData


BOOST_AUTO_TEST_SUITE_END() // BinnedData