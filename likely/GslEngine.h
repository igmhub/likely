// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_GSL_ENGINE
#define LIKELY_GSL_ENGINE

#include "likely/types.h"
#include "likely/AbsEngine.h"

#include "gsl/gsl_multimin.h"

#include <string>
#include <stack>
#include <utility>

namespace likely {
	class GslEngine : public AbsEngine {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters.
		GslEngine(Function f, int nPar, std::string const &algorithm);
		virtual ~GslEngine();
		// Evaluates the engine's function for the specified input parameter values.
        //double operator()(Parameters const& pValues) const;
        // Runs a simplex minimization using the specified initial parameter values
        // and error estimates.
        typedef const gsl_multimin_fminimizer_type *Method;
        void minimize(Method method, Parameters const &initial, Parameters const &errors,
            double minSize = 1e-3, int maxIterations = 1000);
	private:
        int _nPar;
        Function _f;
        gsl_multimin_function _func;
        // Global C-style callback that evaluates the top function on its stack.
        static double _evaluate(const gsl_vector *v, void *params);
        // Maintain a function stack.
        typedef std::pair<Function,Parameters> Binding;
        static std::stack<Binding> &getFunctionStack();
        // Registers our named methods.
        static bool registerGslEngineMethods();
        static bool _registered;
	}; // GslEngine
} // likely

#endif // LIKELY_GSL_ENGINE
