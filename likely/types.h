// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TYPES
#define LIKELY_TYPES

#include "boost/function.hpp"

namespace likely {

    typedef double const *Parameters;
    typedef boost::function<double (Parameters)> Function;

} // likely

#endif // LIKELY_TYPES
