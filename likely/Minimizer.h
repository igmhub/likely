// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_MINIMIZER
#define LIKELY_MINIMIZER

#include "likely/types.h"

#include "boost/function.hpp"
#include "boost/smart_ptr.hpp"

#include <string>
#include <map>

namespace likely {
    class AbsEngine;
	class Minimizer {
	public:
		Minimizer(Function f, int nPar, std::string const &method);
		virtual ~Minimizer();
        Parameters minimize(Parameters const& initial, Parameters const &errors);        

        typedef boost::function<AbsEngine* (Function, int)> MethodFactory;
        void registerMethod(std::string const& methodName, MethodFactory methodFactory);

        typedef std::map<std::string,MethodFactory> Registry;
        static Registry &getRegistry();

	private:
        boost::scoped_ptr<AbsEngine> _engine;
	}; // Minimizer
} // likely

#endif // LIKELY_MINIMIZER
