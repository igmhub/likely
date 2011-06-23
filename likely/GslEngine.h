// Created 24-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_GSL_ENGINE
#define LIKELY_GSL_ENGINE

#include "likely/types.h"
#include "likely/AbsEngine.h"

#include "gsl/gsl_multimin.h"

#include <string>
#include <stack>

namespace likely {
    // Implements GSL multidimensional minimization algorithms. For details, see:
    // http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html
	class GslEngine : public AbsEngine {
	public:
	    // Creates a new engine for the specified function of the specified number
	    // of parameters.
		GslEngine(FunctionPtr f, int nPar, std::string const &algorithm);
        GslEngine(FunctionPtr f, GradientCalculatorPtr g, int nPar,
            std::string const &algorithm);
		virtual ~GslEngine();
        // Performs a minimization without derivatives, using the specified initial
        // parameter values and error estimates.
        typedef const gsl_multimin_fminimizer_type *fMethod;
        FunctionMinimumPtr minimize(fMethod method,
            Parameters const &initial, Parameters const &errors,
            double prec, long maxIterations);
        // Performs a minimization with derivatives, using the specified initial
        // parameter values and error estimates.
        typedef const gsl_multimin_fdfminimizer_type *fdfMethod;
        FunctionMinimumPtr minimizeWithGradient(fdfMethod method,
            Parameters const &initial, Parameters const &errors,
            double prec, long maxIterations, double lineMinTol);        
	private:
        int _nPar;
        FunctionPtr _f;
        Parameters _params;
        GradientCalculatorPtr _gc;
        Gradient _grad;
        gsl_multimin_function _func;
        gsl_multimin_function_fdf _funcWithGradient;
        // Global C-style callbacks that evaluate the top engine on the stack.
        static double _evaluate(const gsl_vector *v, void *p);
        static void _evaluateGradient(const gsl_vector *v, void *p, gsl_vector *g);
        static void _evaluateBoth(const gsl_vector *v, void *p, double *f, gsl_vector *g);
        // Maintains a function stack that the global callbacks use to determine with
        // function to invoke when they are called.
        static std::stack<GslEngine*> &_getEngineStack();
        static GslEngine* _useTopEngine(const gsl_vector *v);
	}; // GslEngine

    // Registers our named methods.
    void registerGslEngineMethods();

} // likely

#endif // LIKELY_GSL_ENGINE
