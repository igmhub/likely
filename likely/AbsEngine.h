// Created 28-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_ENGINE
#define LIKELY_ABS_ENGINE

#include "likely/types.h"

#include "boost/utility.hpp"
#include "boost/function.hpp"

#include <string>
#include <map>
#include <utility>

namespace likely {
    class FunctionMinimum;
	class AbsEngine : public boost::noncopyable {
	public:
		AbsEngine();
		virtual ~AbsEngine();

		// Declare our dynamic entry point for findMinimum.
		typedef boost::function<FunctionMinimumPtr
		    (Parameters const &pInitial, Parameters const &pErrors)> MinimumFinder;
        MinimumFinder minimumFinder;

	    // Declare a global registry for creating engines by name.
        typedef boost::function<AbsEngine* (Function, int, std::string const&)> Factory;
        typedef std::map<std::string, Factory> Registry;
        static Registry &getRegistry();
        
	}; // AbsEngine
	
	// Parses a method name of the form <engine>::<algorithm> or throws a RuntimeError.
    typedef std::pair<std::string,std::string> ParsedMethodName;
    static ParsedMethodName parseMethodName(std::string const &methodName);

    // Finds a minimum of the specified function starting from the initial parameters
    // provided, with steps sizes scaled to the error estimates provided. Returns
    // a smart pointer to a FunctionMinimum object or else throws a RuntimeError.
    // Uses the the method specified by name.
	FunctionMinimumPtr findMinimum(Function f, Parameters const &initial,
        Parameters const &errors, std::string const &methodName);
	
} // likely

#endif // LIKELY_ABS_ENGINE
