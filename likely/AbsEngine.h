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
		// Returns the number of function evaluations by this engine.
        long getEvalCount() const;
        // Returns the number of gradient evaluations by this engine.
        long getGradCount() const;
        
        // The findMinimum functions can access our internals.
        friend FunctionMinimumPtr findMinimum(FunctionPtr,
    	    Parameters const &, Parameters const &, std::string const &, double, long);
        friend FunctionMinimumPtr findMinimum(FunctionPtr, GradientCalculatorPtr,
    	    Parameters const &, Parameters const &, std::string const &, double, long);

    protected:
        // Subclass API for managing evaluation counts.
        void incrementEvalCount() const;
        void incrementGradCount() const;

		// Declares our dynamic entry point for findMinimum.
		typedef boost::function<FunctionMinimumPtr
		    (Parameters const &pInitial, Parameters const &pErrors, double, long)>
		    MinimumFinder;
        MinimumFinder minimumFinder;

	    // Declares global registries for creating engines by name.
        typedef boost::function<AbsEngine* (FunctionPtr, GradientCalculatorPtr,
            int, std::string const&)> EngineFactory;
        typedef std::map<std::string, EngineFactory> EngineRegistry;
        static EngineRegistry &getEngineRegistry();
        
        /*
        typedef boost::function<double (Parameters const &pInitial,
            Parameters const &pErrors,Parameters &pFinal, Covariance &covariance)>
            MinimumAndCovarianceFinder;
        typedef boost::function<void (Parameters const &pValues, Gradients &gValues)>
            FunctionGradientCalculator;
        */
        
    private:        
        mutable long _evalCount, _gradCount;

	}; // AbsEngine
	
    inline long AbsEngine::getEvalCount() const { return _evalCount; }
    inline long AbsEngine::getGradCount() const { return _gradCount; }
    inline void AbsEngine::incrementEvalCount() const { _evalCount++; }
    inline void AbsEngine::incrementGradCount() const { _gradCount++; }
	
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
	FunctionMinimumPtr findMinimum(FunctionPtr f,
	    Parameters const &initial, Parameters const &errors, std::string const &methodName,
        double precision = 1e-3, long maxIterations = 0);
    // Same as above, but using a gradient calculator.
	FunctionMinimumPtr findMinimum(FunctionPtr f, GradientCalculatorPtr gc,
	    Parameters const &initial, Parameters const &errors, std::string const &methodName,
        double precision = 1e-3, long maxIterations = 0);
        
    // Remembers the number of function and gradient evaluations from the last call
    // to findMinimum.
    extern long lastMinEvalCount, lastMinGradCount;
	
} // likely

#endif // LIKELY_ABS_ENGINE
