// Created 28-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_ENGINE
#define LIKELY_ABS_ENGINE

#include "likely/types.h"

#include "boost/utility.hpp"

#include <string>
#include <map>
#include <utility>

namespace likely {
    class FunctionMinimum;
	class AbsEngine : public boost::noncopyable {
	public:
		AbsEngine();
		virtual ~AbsEngine();

	    // Declares a global registry for creating engines by name.
        typedef boost::function<AbsEngine* (FunctionPtr, int, std::string const&)> Factory;
        typedef std::map<std::string, Factory> Registry;
        static Registry &getRegistry();
        
		// Declares our dynamic entry point for findMinimum.
		typedef boost::function<FunctionMinimumPtr
		    (Parameters const &pInitial, Parameters const &pErrors, double, long)>
		    MinimumFinder;
        MinimumFinder minimumFinder;

        /*
        typedef boost::function<double (Parameters const &pInitial,
            Parameters const &pErrors,Parameters &pFinal, Covariance &covariance)>
            MinimumAndCovarianceFinder;
        typedef boost::function<void (Parameters const &pValues, Gradients &gValues)>
            FunctionGradientCalculator;
        */
        
	}; // AbsEngine
	
	// Parses a method name of the form <engine>::<algorithm> or throws a RuntimeError.
    typedef std::pair<std::string,std::string> ParsedMethodName;
    static ParsedMethodName parseMethodName(std::string const &methodName);

    // Finds a minimum of the specified function starting from the initial parameters
    // provided, with steps sizes scaled to the error estimates provided. Returns
    // a smart pointer to a FunctionMinimum object or else throws a RuntimeError.
    // Uses the the method specified by name. The precision parameter determines how
    // precisely the algorithm will attempt to locate the minimum. Its exact definition
    // is algorithm dependent but a smaller value will generally require more evaluations
    // and provide a more precise minimum. Use a positive value for maxIterations to
    // request a maximum number of times that the function is called.
	FunctionMinimumPtr findMinimum(FunctionPtr f, Parameters const &initial,
        Parameters const &errors, std::string const &methodName,
        double precision = 1e-3, long maxIterations = 0);
	
} // likely

#endif // LIKELY_ABS_ENGINE
