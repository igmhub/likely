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
		
        typedef boost::function<AbsEngine* (Function, int, std::string const&)> Factory;
        typedef std::map<std::string, Factory> Registry;
        static Registry &getRegistry();

        /**
		typedef boost::function<double (Parameters const &pInitial,
		    Parameters const &pErrors, Parameters &pFinal)> MinimumFinder;
		//typedef std::pair<EngineFactory,MinimumFinder> 
        
        MinimumFinder getMinimumFinder(std::string const &methodName) const;
        **/
        
	private:	    
	    
	}; // AbsEngine
	
    typedef std::pair<std::string,std::string> ParsedMethodName;
    static ParsedMethodName parseMethodName(std::string const &methodName);

	FunctionMinimum findMinimum(Function f, Parameters const &initial,
        Parameters const &errors, std::string const &methodName);
	
} // likely

#endif // LIKELY_ABS_ENGINE
