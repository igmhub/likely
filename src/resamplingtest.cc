// Created 16-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying resampling methods.

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <cassert>
#include <iostream>

namespace lk = likely;
namespace po = boost::program_options;

class Model : public lk::FitModel {
public:
    Model(int ndim) : lk::FitModel("Multidimensional Gaussian"), _ndim(ndim) {
        defineParameter("sigma",1,0.1);
    }
    double evaluate(std::vector<double> point, lk::Parameters const &values) {
        // Calculate the norm of point
        assert(point.size() == _ndim);
        double norm(0);
        for(int i = 0; i < _ndim; ++i) norm += point[i]*point[i];
        // Lookup the parameter values we should use.
        updateParameterValues(values);
        resetParameterValuesChanged();
        double sigma = getParameterValue("sigma");
        // Return the model evaluated with these parameters at the specified point
        return std::exp(-norm/(2*sigma));
    }
private:
    int _ndim;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int ndim,nbin,ndata,seed;
    double range,sigma0;
    po::options_description cli("Resampling methods test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("ndim", po::value<int>(&ndim)->default_value(2),"Number of dimensions for binning data.")
        ("nbin", po::value<int>(&nbin)->default_value(10),"Number of data bins in each dimension.")
        ("range", po::value<double>(&range)->default_value(2),
            "Data spans [-range,+range] in each dimension.")
        ("ndata", po::value<int>(&ndata)->default_value(10),"Number of random datasets to generate.")
        ("sigma0", po::value<double>(&sigma0)->default_value(1),
            "True value of sigma used to generate random datasets.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed for generating initial parameter values.")
        ;

    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose"));
    if(ndim <= 0) {
        std::cerr << "Expected ndim > 0." << std::endl;
        return -1;
    }
    if(nbin <= 0) {
        std::cerr << "Expected nbin > 0." << std::endl;
        return -1;
    }
    if(range <= 0) {
        std::cerr << "Expected range > 0." << std::endl;
        return -1;
    }
    if(sigma0 <= 0) {
        std::cerr << "Expected sigma0 > 0." << std::endl;
        return -1;
    }

    try {
        // Initialize our random engine.
        lk::Random &random = lk::Random::instance();
        random.setSeed(seed);
        
        // Create the model.
        Model model(ndim);
        
        // Create a prototype dataset.
        std::vector<lk::AbsBinningCPtr> axes;
        for(int i = 0; i < ndim; ++i) {
            lk::AbsBinningCPtr axis(new lk::UniformBinning(-range,+range,nbin));
            axes.push_back(axis);
        }
        lk::BinnedDataPtr prototype(new lk::BinnedData(axes));
        
        // Fill each bin of the prototype dataset with the model evaluated with sigma=sigma0
        std::vector<double> point(ndim), params;
        params.push_back(sigma0);
        for(int index = 0; index < prototype->getNBinsTotal(); ++index) {
            prototype->getBinCenters(index,point);
            prototype->setData(index,model.evaluate(point,params));
        }
        
        // Create our resampler.
        lk::BinnedDataResampler resampler(seed);
        
        // Generate random datasets.
        for(int i = 0; i < ndata; ++i) {
            // Clone our prototype.
            lk::BinnedDataPtr data(prototype->clone());
            // Generate a random covariance matrix with determinant 1 for this dataset.
            lk::CovarianceMatrixPtr covariance =
                lk::generateRandomCovariance(data->getNBinsWithData(),seed);
            data->setCovarianceMatrix(covariance);
            // Sample the covariance to generate random offset for each bin.
            boost::shared_array<double> offsets = covariance->sample(1,seed);
            int nextOffset(0);
            for(lk::BinnedData::IndexIterator iter = data->begin(); iter != data->end(); ++iter) {
                data->addData(*iter,offsets[nextOffset++]);
            }
            // Add this dataset to our resampler.
            resampler.addObservation(data);
        }
        std::cout << "Resampling " << resampler.getNObservations() << " observations." << std::endl;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
