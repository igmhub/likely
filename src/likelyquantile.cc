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
        random->setSeed(42);
	    // Accumulate samples
	    double sample;
	    for(int j = 0; j < 20; ++j) {
	    	long maxSamples = pow(10,1+j*.2);
		    for(int i = 0; i < maxSamples; ++i){
		    	sample = random->getNormal();
		    	median.accumulate(sample);
		    	s1.accumulate(sample);
		    	s2.accumulate(sample);
		    	s3.accumulate(sample);
		    	exact.accumulate(sample);
		    }
			std::cout << maxSamples
				<< " " << 0-median.getQuantile() 
				<< " " << 1-s1.getQuantile() 
				<< " " << 2-s2.getQuantile() 
				<< " " << 3-s3.getQuantile() 
				<< std::endl;
			std::cout << maxSamples
				<< " " << 0-exact.getQuantile(.5) 
				<< " " << 1-exact.getQuantile(0.841345) 
				<< " " << 2-exact.getQuantile(0.97725) 
				<< " " << 3-exact.getQuantile(0.99865) 
				<< std::endl;
		}
	}
	catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

	return 0;
}