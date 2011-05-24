// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_GSL_ERROR_HANDLER
#define LIKELY_GSL_ERROR_HANDLER

#include "gsl/gsl_errno.h"

namespace likely {
	class GslErrorHandler {
	public:
	    // Installs a new error handler, remembering the original handler.
		GslErrorHandler();
		// Restores the error handler that was in place when this object was created.
		virtual ~GslErrorHandler();
	private:
		// Handles an error by throwing a RuntimeError with a descriptive message.
		static void _handle(const char *reason,const char *file,int line,int gsl_errno);
        gsl_error_handler_t const *_original;
	}; // GslErrorHandler	
} // likely

#endif // LIKELY_GSL_ERROR_HANDLER
