// Created 10-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_FIT_PARAMETER_STATISTICS
#define LIKELY_FIT_PARAMETER_STATISTICS

#include "likely/FitParameter.h"
#include "likely/types.h"

#include "boost/smart_ptr.hpp"

#include <iosfwd>
#include <vector>
#include <string>

namespace likely {
    class WeightedAccumulator;
    class CovarianceAccumulator;
    class ExactQuantileAccumulator;
    // Accumulates fit parameter value statistics.
	class FitParameterStatistics {
	public:
	    // Creates a new statistics accumulator for values of the specified fit parameters.
		FitParameterStatistics(FitParameters const &params);
		virtual ~FitParameterStatistics();
		// Returns the number of free parameters we are keeping statistics for.
        int getNFreeParameters() const;
        // Returns the number of times our statistics have been successfully updated.
        int getNUpdates() const;
        // Updates our statistics using the specified parameter values and function value.
        void update(Parameters pvalues, double fval);
        // Prints our statistics to the specified output stream.
        void printToStream(std::ostream &out, std::string const &formatSpec = "%12.6f") const;
	private:
        int _nfree, _nupdates;
        Parameters _baseline;
        boost::scoped_array<WeightedAccumulator> _stats;
        boost::scoped_array<ExactQuantileAccumulator> _quantiles;
        boost::scoped_ptr<CovarianceAccumulator> _accumulator;
        std::vector<std::string> _labels;
	}; // FitParameterStatistics
	
    inline int FitParameterStatistics::getNFreeParameters() const { return _nfree; }
    inline int FitParameterStatistics::getNUpdates() const { return _nupdates; }
	
} // likely

#endif // LIKELY_FIT_PARAMETER_STATISTICS
