// Created 13-Aug-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates and benchmarks the Random class.

#include "likely/likely.h"

#include "boost/format.hpp"

#include <sys/resource.h>
#include <cmath>
#include <iostream>

#define BENCHMARK(METHOD) {\
getrusage(RUSAGE_SELF,&before); \
for(int i = 0; i < repeat; ++i) { rval = random.METHOD(); } \
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
    lk::Random &random(lk::Random::instance());
    int repeat(1000000);
    double rval;
    struct rusage before,after;
    boost::format results("%20s: %.3f nanosecs/call\n");
    
    BENCHMARK(getUniform);
    BENCHMARK(getNormal);
}
