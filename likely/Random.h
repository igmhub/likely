// Created 08-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_RANDOM
#define LIKELY_RANDOM

#include "boost/random/mersenne_twister.hpp"
#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

#include <cstddef>

namespace likely {
	class Random {
	public:
		Random();
        void setSeed(int seedValue);
        // Returns a double-precision value uniformly sampled from [0,1).
        double getUniform();
        // Returns a double-precision value with mean 0 and RMS 1.
        double getNormal();
        // Returns a single-precision value uniformly sampled from [0,1) using
        // an inline coding of SFMT (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/)
        float getFastUniform();
        // Fills the specified array with double-precision values uniformly sampled from [0,1)
        // using the specified seed (that is independent of the seed used by getUniform and
        // getNormal).
        static void fillArrayUniform(double *array, std::size_t size, int seed);
        // Returns a shared array filled with single-precision values normally distributed with
        // mean 0 and RMS 1 using the specified seed (that is independent of the seed used
        // by getUniform and getNormal). The number nrandom of numbers generated must be a
        // multiple of 4 that is >= 624.
        static boost::shared_array<float> fillArrayNormal(std::size_t nrandom, int seed);
        // Returns a reference to this object's internal generator, so that it
        // can be used for other distributions. This should only be used on the
        // global shared instance.
        boost::mt19937 &getGenerator();
        // Returns the global shared Random instance.
        static Random &instance();
	private:
        boost::mt19937 _generator;
        boost::function<double ()> _uniform, _gauss;
        static const double _ziggurat_ytab[128], _ziggurat_wtab[128];
        static const uint32_t _ziggurat_ktab[128];
	}; // Random
	
    inline double Random::getUniform() { return _uniform(); }
    inline double Random::getNormal() { return _gauss(); }
    inline boost::mt19937 &Random::getGenerator() { return _generator; }
	
    // Allocates an array with the 128-bit alignment required by the Random::fillArrayX methods
    // where size is in bytes.
    void *allocateAlignedArray(std::size_t byteSize);
    
    // Wraps the result of calling allocateAlignedArray into a smart array pointer with a
    // custom deleter that calls free.
    boost::shared_array<float> allocateAlignedFloatArray(std::size_t size);
    boost::shared_array<double> allocateAlignedDoubleArray(std::size_t size);
    
} // likely

#endif // LIKELY_RANDOM
