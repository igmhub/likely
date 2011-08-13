// Created 08-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_RANDOM
#define LIKELY_RANDOM

#include "boost/random/mersenne_twister.hpp"
#include "boost/function.hpp"

namespace likely {
	class Random {
	public:
		Random();
        void setSeed(int seedValue);
        // Returns a double-precision value uniformly sampled from [0,1)
        double getUniform();
        // Returns a double-precision value with mean 0 and RMS 1.
        double getNormal();
        // Returns a single-precision value uniformly sampled from [0,1) using
        // an inline coding of SFMT (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/)
        float getFastUniform();
        // Returns a reference to this object's internal generator, so that it
        // can be used for other distributions. This should only be used on the
        // global shared instance.
        boost::mt19937 &getGenerator();
        // Returns the global shared Random instance.
        static Random &instance();
	private:
        boost::mt19937 _generator;
        boost::function<double ()> _uniform, _gauss;
	}; // Random
	
	inline void Random::setSeed(int seedValue) { _generator.seed(seedValue); }
    inline double Random::getUniform() { return _uniform(); }
    inline double Random::getNormal() { return _gauss(); }
    inline boost::mt19937 &Random::getGenerator() { return _generator; }
	
} // likely

#endif // LIKELY_RANDOM
