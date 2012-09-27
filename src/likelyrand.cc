// Created 13-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates and benchmarks the Random class.

#include "likely/Random.h"
#include "likely/WeightedAccumulator.h"
#include "likely/QuantileAccumulator.h"

#include "boost/format.hpp"

#include <sys/resource.h>
#include <cmath>
#include <iostream>
#include <cassert>

#define BENCHMARK_ASSIGN(RESULT,METHOD,ARGS) {\
getrusage(RUSAGE_SELF,&before);\
RESULT = random->METHOD ARGS;\
getrusage(RUSAGE_SELF,&after); \
std::cout << results % #METHOD % (1e3*elapsed(before,after)/repeat); \
}

#define BENCHMARK_LOOP(METHOD) {\
getrusage(RUSAGE_SELF,&before);\
for(int i = 0; i < repeat; ++i) random->METHOD();\
getrusage(RUSAGE_SELF,&after); \
std::cout << results % #METHOD % (1e3*elapsed(before,after)/repeat); \
}

namespace lk = likely;

// Returns the number of elapsed microseconds from before to after.
double elapsed(struct timeval const &before, struct timeval const &after) {
    return (after.tv_sec - before.tv_sec)*1e6 + (after.tv_usec - before.tv_usec);
}
double elapsed(struct rusage const &before, struct rusage const &after) {
    return elapsed(before.ru_utime,after.ru_utime) + elapsed(before.ru_stime,after.ru_stime);
}

int main(int argc, char **argv) {
    double rval;
    struct rusage before,after;
    boost::format results("%22s: %6.3f nanosecs/call\n");
    int repeat(1<<24); //  = 16M
    std::cout << "Testing generation of " << repeat << " random numbers." << std::endl;
    
    lk::RandomPtr random = lk::Random::instance();
    random->setSeed(1234);
    
    // Random integer tests
    {
        // Test getInteger
        int size(6),ntrials(100000);
        std::vector<int> counts(size,0);
        for(int trial = 0; trial < ntrials; ++trial) {
            counts[random->getInteger(1,size)-1]++;
        }
        std::cout << "integer counts:";
        for(int k = 0; k < size; ++k) std::cout << ' ' << counts[k];
        std::cout << std::endl;
    }
    {
        // Test partialShuffle
        int ssize(10),size(3),ntrials(10000);
        std::vector<int> sample;
        std::vector<std::vector<int> > counts;
        for(int i = 0; i < size; ++i) counts.push_back(std::vector<int>(ssize,0));
        for(int i = 0; i < ssize; ++i) sample.push_back(i);
        for(int trial = 0; trial < ntrials; ++trial) {
            random->partialShuffle(sample,size);
            for(int i = 0; i < size; ++i) counts[i][sample[i]]++;
        }
        for(int i = 0; i < size; ++i) {
            std::cout << "shuffle[" << i << "]";
            for(int j = 0; j < ssize; ++j) {
                std::cout << ' ' << counts[i][j];
            }
            std::cout << std::endl;
        }
    }
    {
        // Test sampleWithReplacement
        int ssize(10), size(15), ntrials(1000);
        std::vector<int> sample(ssize,0);
        std::vector<int> counts(ssize,0);
        for(int trial = 0; trial < ntrials; ++trial) {
            random->sampleWithReplacement(sample,size);
            for(int i = 0; i < ssize; ++i) counts[i] += sample[i];
        }
        std::cout << "bs counts";
        for(int i = 0; i < ssize; ++i) std::cout << ' ' << counts[i];
        std::cout << std::endl;
    }
    
    // Random real benchmarks
    BENCHMARK_LOOP(getUniform);
    BENCHMARK_LOOP(getNormal);
    BENCHMARK_LOOP(getFastUniform);

    int seed = 123;
    {
        std::size_t nrandom(repeat);
        boost::shared_array<double> dbuffer;
        BENCHMARK_ASSIGN(dbuffer,fillDoubleArrayUniform,(nrandom));
        lk::WeightedAccumulator stats;
        lk::QuantileAccumulator median, q90(0.9);
        for(int i = 0; i < repeat; ++i) {
            double value(dbuffer[i]);
            if(i < 10) std::cout << i << ' ' << value << std::endl;
            stats.accumulate(value);
            median.accumulate(value);
            q90.accumulate(value);
        }
        std::cout << "mean = " << stats.mean() << ", variance = " << stats.variance() << std::endl;
        std::cout << "median = " << median.getQuantile() << ", 90% quantile = "
            << q90.getQuantile() << std::endl;
    }    
    {
        std::size_t nrandom(repeat);
        boost::shared_array<float> fbuffer;
        BENCHMARK_ASSIGN(fbuffer,fillFloatArrayNormal,(nrandom));
        lk::WeightedAccumulator stats;
        lk::QuantileAccumulator median, sigma(1-0.5*0.317310508);
        for(int i = 0; i < repeat; ++i) {
            double value(fbuffer[i]);
            if(i < 10) std::cout << i << ' ' << value << std::endl;
            stats.accumulate(value);
            median.accumulate(value);
            sigma.accumulate(value);
        }
        std::cout << "mean = " << stats.mean() << ", variance = " << stats.variance() << std::endl;
        std::cout << "median = " << median.getQuantile() << ", 1-sigma quantile = "
            << sigma.getQuantile() << std::endl;
    }
    {
        std::size_t nrandom(repeat);
        boost::shared_array<double> dbuffer;
        BENCHMARK_ASSIGN(dbuffer,fillDoubleArrayNormal,(nrandom));
        lk::WeightedAccumulator stats;
        lk::QuantileAccumulator median, sigma(1-0.5*0.317310508);
        for(int i = 0; i < repeat; ++i) {
            double value(dbuffer[i]);
            if(i < 10) std::cout << i << ' ' << value << std::endl;
            stats.accumulate(value);
            median.accumulate(value);
            sigma.accumulate(value);
        }
        std::cout << "mean = " << stats.mean() << ", variance = " << stats.variance() << std::endl;
        std::cout << "median = " << median.getQuantile() << ", 1-sigma quantile = "
            << sigma.getQuantile() << std::endl;
    }
    
    // Compare the requested and returned nrandom for each array-filling method.
    for(int size = 1; size < 700; size += 25) {
        double sum;
        std::size_t nrand1 = size,nrand2 = size,nrand3 = size;
        {
            boost::shared_array<double> buf = random->fillDoubleArrayUniform(nrand1);
            for(int k = 0; k < nrand1; ++k) sum += buf[k];
        }
        {
            boost::shared_array<double> buf = random->fillDoubleArrayNormal(nrand2);
            for(int k = 0; k < nrand2; ++k) sum += buf[k];
        }
        {
            boost::shared_array<float> buf = random->fillFloatArrayNormal(nrand3);
            for(int k = 0; k < nrand3; ++k) sum += buf[k];
        }
        std::cout << size << ' ' << nrand1 << ' ' << nrand2 << ' ' << nrand3 << std::endl;
    }
}
