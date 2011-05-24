// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/GslErrorHandler.h"

namespace local = likely;

local::GslErrorHandler::GslErrorHandler() {
    _original = gsl_set_error_handler(_handle);
}

local::GslErrorHandler::~GslErrorHandler() {
    gsl_set_error_handler(_original);
}

void local::GslErrorHandler::_handle(
const char *reason,const char *file,int line,int gsl_errno) {
}
