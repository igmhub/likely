// Created 13-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates and benchmarks the Random class.

#include "likely/likely.h"

#include "boost/format.hpp"

#include <sys/resource.h>
#include <cmath>
#include <iostream>
#include <cassert>

#define BENCHMARK_ASSIGN(RESULT,METHOD,ARGS) {\
getrusage(RUSAGE_SELF,&before);\
RESULT = random.METHOD ARGS;\
getrusage(RUSAGE_SELF,&after); \
std::cout << results % #METHOD % (1e3*elapsed(before,after)/repeat); \
}

#define BENCHMARK_LOOP(METHOD) {\
getrusage(RUSAGE_SELF,&before);\
for(int i = 0; i < repeat; ++i) random.METHOD();\
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
    
    lk::Random &random(lk::Random::instance());
    random.setSeed(1234);

    BENCHMARK_LOOP(getUniform);
    BENCHMARK_LOOP(getNormal);
    BENCHMARK_LOOP(getFastUniform);

    {
        boost::shared_array<double> dbuffer;
        BENCHMARK_ASSIGN(dbuffer,fillDoubleArrayUniform,(repeat,123));
    }    
    {
        boost::shared_array<float> fbuffer;
        BENCHMARK_ASSIGN(fbuffer,fillFloatArrayNormal,(repeat,123));
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
        boost::shared_array<double> dbuffer;
        BENCHMARK_ASSIGN(dbuffer,fillDoubleArrayNormal,(repeat,123));
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
}
