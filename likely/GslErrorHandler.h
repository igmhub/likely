// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_GSL_ERROR_HANDLER
#define LIKELY_GSL_ERROR_HANDLER

#include "gsl/gsl_errno.h"

#include <string>
#include <stack>

namespace likely {
	class GslErrorHandler {
	public:
	    // Installs a new context-specific error handler, remembering the original handler.
		GslErrorHandler(std::string const &context);
		// Restores the error handler that was in place when this object was created.
		virtual ~GslErrorHandler();
	private:
		// Remembers the original handler
        gsl_error_handler_t *_original;
		// Handles an error by throwing a RuntimeError with a descriptive message.
		static void _handle(const char *reason,const char *file,int line,int gsl_errno);
        // Keeps track of the current context in case of nested handlers.
        static std::stack<std::string> &getContextStack();
	}; // GslErrorHandler	
} // likely

#endif // LIKELY_GSL_ERROR_HANDLER
