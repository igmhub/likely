// Created 27-Aug-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// A quantile accumulator test program.

#include <iostream>

#include "likely/likely.h"
 
#include "boost/program_options.hpp"

namespace lk = likely;

int main(int argc, char **argv) {
    try {
		// Initialize quantile accumulators
		lk::QuantileAccumulator median(.5);
		lk::QuantileAccumulator s1(0.841345);
		lk::QuantileAccumulator s2(0.97725);
		lk::QuantileAccumulator s3(0.99865);
		lk::ExactQuantileAccumulator exact;
		// Initialize random number generator
	    lk::RandomPtr random(new lk::Random());
        random->setSeed(100);
	    // Accumulate samples
	    double ndecades = 2.5;
	    double spacing = .2;
		int ntrials = 100;
	    for(int j = 0; j < ndecades/spacing; ++j) {
	    	long nsamples = pow(10,1+j*spacing);
	    	// Prepare quantile result accumulators
	    	lk::WeightedAccumulator approx0,approx1,approx2,approx3,exact0,exact1,exact2,exact3;
	    	for(int k = 0; k < ntrials; ++k) {
			    for(int i = 0; i < nsamples; ++i){
			    	double sample(random->getNormal());
			    	median.accumulate(sample);
			    	s1.accumulate(sample);
			    	s2.accumulate(sample);
			    	s3.accumulate(sample);
			    	exact.accumulate(sample);
			    }
			    approx0.accumulate(median.getQuantile());
			    approx1.accumulate(s1.getQuantile());
			    approx2.accumulate(s2.getQuantile());
			    approx3.accumulate(s3.getQuantile());
			    exact0.accumulate(exact.getQuantile(0.5));
			    exact1.accumulate(exact.getQuantile(0.841345));
			    exact2.accumulate(exact.getQuantile(0.97725));
			    exact3.accumulate(exact.getQuantile(0.99865));
			}
			std::cout << nsamples
				<< " " << approx0.mean() << " " << approx0.error()
				<< " " << 1 - approx1.mean() << " " << approx1.error()
				<< " " << 2 - approx2.mean() << " " << approx2.error()
				<< " " << 3 - approx3.mean() << " " << approx3.error()
				<< " " << exact0.mean() << " " << exact0.error()
				<< " " << 1 - exact1.mean() << " " << exact1.error()
				<< " " << 2 - exact2.mean() << " " << exact2.error()
				<< " " << 3 - exact3.mean() << " " << exact3.error()
				<< std::endl;
		}
	}
	catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

	return 0;
}