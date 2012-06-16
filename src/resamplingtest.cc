// Created 16-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying resampling methods.

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/foreach.hpp"
#include "boost/bind.hpp"
#include "boost/ref.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

namespace lk = likely;
namespace po = boost::program_options;

void dumpParameterValues(lk::Parameters const &params, std::ostream &out) {
    BOOST_FOREACH(double pvalue, params) {
        out << ' ' << pvalue;
    }
    out << std::endl;
}

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

class Fitter {
public:
    Fitter(Model &model) : _model(model) { }
    // Fits the specified data to our model.
    lk::FunctionMinimumPtr fit(lk::BinnedDataCPtr data) {
        _data = data;
        lk::FunctionPtr fptr(new lk::Function(*this));
        return _model.findMinimum(fptr,"mn2::vmetric");
    }
    // Generates MCMC samples (and updates fmin)
    void mcmc(lk::FunctionMinimumPtr fmin, int nmc, int nskip, std::ostream &out) {
        lk::FunctionPtr fptr(new lk::Function(*this));
        lk::FitParameters params(fmin->getFitParameters());
        lk::MarkovChainEngine engine(fptr,lk::GradientCalculatorPtr(),params,"saunter");
        lk::Parameters pvalues;
        int ntrial(nmc*nskip);
        lk::MarkovChainEngine::Callback callback = boost::bind(dumpParameterValues,_1,boost::ref(out));
        engine.generate(fmin,ntrial,ntrial,callback,nskip);
    }
    // Returns chiSquare/2 for the specified model parameter values.
    double operator()(lk::Parameters const &params) const {
        assert(params.size() == _model.getNParameters());
        std::vector<double> point, prediction;
        for(lk::BinnedData::IndexIterator iter = _data->begin(); iter != _data->end(); ++iter) {
            int index(*iter);
            _data->getBinCenters(index,point);
            prediction.push_back(_model.evaluate(point,params));
        }
        return 0.5*_data->chiSquare(prediction);
    }
private:
    Model &_model;
    lk::BinnedDataCPtr _data;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int ndim,nbin,ndata,nmc,nskip,ntrial,seed;
    double range,sigma0,covScale;
    po::options_description cli("Resampling methods test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("ndim", po::value<int>(&ndim)->default_value(2),"Number of dimensions for binning data.")
        ("nbin", po::value<int>(&nbin)->default_value(10),"Number of data bins in each dimension.")
        ("range", po::value<double>(&range)->default_value(2),
            "Data spans [-range,+range] in each dimension.")
        ("ndata", po::value<int>(&ndata)->default_value(10),"Number of random datasets to generate.")
        ("nmc", po::value<int>(&nmc)->default_value(100), "Number of MCMC samples to generate.")
        ("nskip", po::value<int>(&nskip)->default_value(10), "Number of MCMC trials to skip per sample.")
        ("ntrial", po::value<int>(&ntrial)->default_value(100), "Number of bootstrap trials to fit.")
        ("sigma0", po::value<double>(&sigma0)->default_value(1),
            "True value of sigma used to generate random datasets.")
        ("cov-scale", po::value<double>(&covScale)->default_value(100),
            "Scale of random dataset covariances to generate.")
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
    assert(ndim > 0);
    assert(nbin > 0);
    assert(ndata > 0);
    assert(nmc > 0);
    assert(nskip > 0);
    assert(ntrial > 0);

    try {
        // Create the model and fitter we will use.
        Model model(ndim);
        Fitter fitter(model);
        
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
        int size(prototype->getNBinsWithData());
        
        // Create our resampler.
        lk::BinnedDataResampler resampler(seed);
        
        // Generate random datasets.
        std::cout << "Generating " << ndata << " datasets with " << size << " bins..." << std::endl;
        lk::CovarianceMatrixPtr covariance;
        for(int i = 0; i < ndata; ++i) {
            // Clone our prototype.
            lk::BinnedDataPtr data(prototype->clone());
            // Generate a random covariance matrix with determinant 1 for this dataset.
            if(0 == i) covariance = lk::generateRandomCovariance(size,seed,covScale);
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
        
        // Fit the combined data.
        lk::BinnedDataCPtr combined = resampler.combined();
        lk::FunctionMinimumPtr combinedFit = fitter.fit(combined);
        combinedFit->printToStream(std::cout);
        
        // Dump the likelihood curve for the combined data.
        {
            std::ofstream out("likely.dat");
            std::vector<double> params(1);
            for(int i = 0; i < 1000; ++i) {
                params[0] = 2e-3*(i+1);
                double nll= fitter(params);
                double prob = std::exp(-nll/size);
                out << params[0] << ' ' << prob << std::endl;
            }
            out.close();
        }
        
        // Generate MCMC samples from the combined likelihood.
        {
            std::cout << "Generating " << nmc << " MCMC samples with skip " << nskip << "..." << std::endl;
            std::ofstream out("mcmc.dat");
            out << "sigma" << std::endl;
            fitter.mcmc(combinedFit,nmc,nskip,out);
            out.close();
        }
        // Loop over bootstrap trials.
        {
            std::cout << "Generating " << ntrial << " bootstrap trials..." << std::endl;
            lk::Parameters fitted;
            std::ofstream out("bstrials.dat");
            out << "sigma" << std::endl;
            for(int i = 0; i < ntrial; ++i) {
                lk::BinnedDataCPtr sample = resampler.bootstrap(ndata,false);
                lk::FunctionMinimumPtr sampleFit = fitter.fit(sample);
                dumpParameterValues(sampleFit->getParameters(true),out);
                if((i+1)%1000==0) std::cout << "...completed " << (i+1) << " trials" << std::endl;
            }
            out.close();
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
