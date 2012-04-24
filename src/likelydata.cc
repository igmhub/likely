// Created 24-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// Demonstates the CovarianceMatrix class.

#include "likely/likely.h"

#include <iostream>
#include <sys/resource.h>

namespace lk = likely;

// Returns the number of elapsed microseconds from before to after.
double elapsed(struct timeval const &before, struct timeval const &after) {
    return (after.tv_sec - before.tv_sec)*1e6 + (after.tv_usec - before.tv_usec);
}
double elapsed(struct rusage const &before, struct rusage const &after) {
    return elapsed(before.ru_utime,after.ru_utime) + elapsed(before.ru_stime,after.ru_stime);
}

int main(int argc, char **argv) {
    return 0;
}
