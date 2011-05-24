// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslErrorHandler.h"

#include "boost/format.hpp"

#include <iostream>

namespace local = likely;

local::GslErrorHandler::GslErrorHandler(std::string const &context)
: _original(gsl_set_error_handler(_handle))
{
    _contextStack.push(context);
}

local::GslErrorHandler::~GslErrorHandler() {
    _contextStack.pop();
    gsl_set_error_handler(_original);
}

void local::GslErrorHandler::_handle(
const char *reason,const char *file,int line,int gsl_errno) {
    static boost::format messageFormat("%s <GSL error at line %d of %s> %s\n");
    std::cerr << boost::str(messageFormat % _contextStack.top() % line % file % reason);
}

std::stack<std::string> local::GslErrorHandler::_contextStack = std::stack<std::string>();