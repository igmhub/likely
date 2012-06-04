// Created 11-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/function.h"

#include "boost/bind.hpp"

template <class P> likely::GenericFunctionPtr likely::createFunctionPtr(boost::shared_ptr<P> pimpl) {
    GenericFunctionPtr fptr(new GenericFunction(boost::bind(&P::operator(),pimpl,_1)));
    return fptr;
}
