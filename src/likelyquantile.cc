// Created 27-Aug-2011 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// A weighted sum accumulator test program.

#include <iostream>

#include "likely/likely.h"
 
#include "boost/program_options.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"

namespace lk = likely;
namespace po = boost::program_options;

class ExactQuantileAccumulator {

public:
	ExactQuantileAccumulator(double quantileProbability) 
	: _quantileProbability(quantileProbability) {
		_weightedCount = 0;
	};
	~ExactQuantileAccumulator() {};
	// Accumulates one (possibly weighted) sample value;
	void accumulate(double value, double weight = 1) {
		valueWeightPairs.insert(std::pair<double, double>(value, weight));
		_weightedCount += weight;
	};
	// Returns the number of weighted samples accumulated.
	int count() const {
		return valueWeightPairs.size();
	};
	// Returns the quantile value based on the samples accumulated so far.
	double getQuantile() const {
		return getQuantile(_quantileProbability);
	};
	// Returns the quantile value based on the samples accumulated so far.
	double getQuantile(double quantileProbability) const {
		double _weightedSoFar = 0;
		std::multiset<std::pair<double,double> >::iterator it;
		for(it = valueWeightPairs.begin(); it != valueWeightPairs.end(); it++) {
			_weightedSoFar += (*it).second;
			if( _weightedSoFar/_weightedCount >= quantileProbability) break;
		}
		return (*it).first;
	};
private:
	std::multiset<std::pair<double,double> > valueWeightPairs;
	double _quantileProbability;
	double _weightedCount;

}; // ExactQuantileAccumulator

int main(int argc, char **argv) {
	uint32_t seed;
	int n;
	po::options_description cli("Quantile accumulator test program.");
	cli.add_options()
		("help,h", "Print this info and exits.")
		("seed", po::value<uint32_t>(&seed)->default_value(0), "Random number generator seed.")
		("nsamples,n", po::value<int>(&n)->default_value(1000), "Number of samples to accumulate.")
		;
	// Parse command line options
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
	catch(std::exception const &e) {
		std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
		return -1;
	}
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }

    try {
		// Initialize quantile accumulator
		lk::QuantileAccumulator median(.5);
		lk::QuantileAccumulator oneSig(0.841345);
		lk::QuantileAccumulator twoSig(0.97725);
		lk::QuantileAccumulator threeSig(0.99865);
		ExactQuantileAccumulator exact(.5);

		// Initialize random seed
	    if(0 == seed) seed = static_cast<uint32_t>(std::time(0));
	    srand(seed);
	    boost::mt19937 uniform(seed);
	    // Normal distribution
	    double sigma(1);
	    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
	        normal(uniform, boost::normal_distribution<>(0,sigma));
	    // Accumulate trials
	    for(int i = 0; i < n; ++i){
	    	double sample(normal());
	    	median.accumulate(sample);
	    	oneSig.accumulate(sample);
	    	twoSig.accumulate(sample);
	    	threeSig.accumulate(sample);
	    	exact.accumulate(sample);
	    }
		// 1,2,3-sigma quantiles
		std::cout << "Accumulated " << median.count() 
			<< " samples with median " << median.getQuantile() 
			<< ", one sigma at " << oneSig.getQuantile() 
			<< ", two sigma at " << twoSig.getQuantile() 
			<< ", three sigma at " << threeSig.getQuantile() 
			<< std::endl;
		// 1,2,3-sigma quantiles
		std::cout << "Exact Accumulated " << exact.count() 
			<< " samples with median " << exact.getQuantile(.5) 
			<< ", one sigma at " << exact.getQuantile(0.841345) 
			<< ", two sigma at " << exact.getQuantile(0.97725) 
			<< ", three sigma at " << exact.getQuantile(0.99865) 
			<< std::endl;
	}
	catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }

	return 0;
}