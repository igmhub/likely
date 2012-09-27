// Created 16-Jun-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A test program for studying resampling methods.

#include "likely/BinnedData.h"
#include "likely/FitParameter.h"
#include "likely/MarkovChainEngine.h"
#include "likely/Random.h"
#include "likely/UniformBinning.h"
#include "likely/BinnedDataResampler.h"
#include "likely/CovarianceMatrix.h"
#include "likely/FunctionMinimum.h"

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

double model(std::vector<double> point, double sigma) {
    double norm(0);
    BOOST_FOREACH(double coord, point) {
        norm += coord*coord;
    }
    return std::exp(-norm/(2*sigma));
}

void updateCount(int &count, double sigma0, lk::Parameters const &params, std::ostream *out) {
    double sigma = params[0];
    if(sigma < sigma0) count++;
    if(out) *out << sigma << std::endl;
}

class Fitter {
public:
    Fitter() { }
    // Fits the specified data to our model.
    lk::FunctionMinimumPtr fit(lk::BinnedDataCPtr data) {
        _data = data;
        lk::FunctionPtr fptr(new lk::Function(*this));
        lk::FitParameters params;
        params.push_back(lk::FitParameter("sigma",1,0.1));
        return lk::findMinimum(fptr,params,"mn2::vmetric");
    }
    // Generates MCMC samples (and updates fmin)
    double mcmc(lk::FunctionMinimumPtr fmin, int nmc, int nskip, double sigma0, std::ostream *out) {
        lk::FunctionPtr fptr(new lk::Function(*this));
        lk::FitParameters params(fmin->getFitParameters());
        lk::MarkovChainEngine engine(fptr,lk::GradientCalculatorPtr(),params,"saunter");
        lk::Parameters pvalues;
        int ntrial(nmc*nskip);
        int count = 0;
        lk::MarkovChainEngine::Callback callback = boost::bind(updateCount,boost::ref(count),sigma0,_1,out);
        engine.generate(fmin,ntrial,ntrial,callback,nskip);
        return (double)count/nmc;
    }
    // Returns chiSquare/2 for the specified model parameter values.
    double operator()(lk::Parameters const &params) const {
        double sigma = params[0];
        std::vector<double> point, prediction;
        for(lk::BinnedData::IndexIterator iter = _data->begin(); iter != _data->end(); ++iter) {
            _data->getBinCenters(*iter,point);
            prediction.push_back(model(point,sigma));
        }
        return 0.5*_data->chiSquare(prediction);
    }
private:
    lk::BinnedDataCPtr _data;
};

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    int ndim,nbin,ndata,nmc,nskip,ntrial,nexpt,seed;
    double range,sigma0,covScale,covRescale;
    po::options_description cli("Resampling methods test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("ndim", po::value<int>(&ndim)->default_value(2),"Number of dimensions for binning data.")
        ("nbin", po::value<int>(&nbin)->default_value(5),"Number of data bins in each dimension.")
        ("range", po::value<double>(&range)->default_value(2),
            "Data spans [-range,+range] in each dimension.")
        ("ndata", po::value<int>(&ndata)->default_value(15),"Number of random datasets to generate.")
        ("nmc", po::value<int>(&nmc)->default_value(100), "Number of MCMC samples to generate.")
        ("nskip", po::value<int>(&nskip)->default_value(10), "Number of MCMC trials to skip per sample.")
        ("ntrial", po::value<int>(&ntrial)->default_value(100), "Number of bootstrap trials to fit.")
        ("nexpt", po::value<int>(&nexpt)->default_value(1), "Number of experiments to perform.")
        ("sigma0", po::value<double>(&sigma0)->default_value(1),
            "True value of sigma used to generate random datasets.")
        ("cov-scale", po::value<double>(&covScale)->default_value(1),
            "Scale of random dataset covariances to generate.")
        ("cov-rescale", po::value<double>(&covRescale)->default_value(1),
            "Amount to rescale covariance used in likelihood fits.")
        ("same-cov", "All datasets generated with same covariance matrix.")
        ("bs-combined", "Bootstrap samples are assigned the combined covariance.")
        ("unweighted", "Bootstrap samples are unweighted.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed for generating initial parameter values.")
        ("dump", "Dumps analysis results to text files.")
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
    bool verbose(vm.count("verbose")), dump(vm.count("dump")), sameCov(vm.count("same-cov")),
        bsCombined(vm.count("bs-combined")), unweighted(vm.count("unweighted"));
    assert(ndim > 0);
    assert(nbin > 0);
    assert(ndata > 0);
    assert(nmc > 0);
    assert(nskip > 0);
    assert(ntrial > 0);
    assert(nexpt > 0);
    if(nexpt > 1) {
        verbose = false;
        dump = false;
        std::cout << "mc bs" << std::endl;
    }

    try {
        // Initialize our random generator.
        lk::Random::instance()->setSeed(seed);
        
        // Create the fitter we will use.
        Fitter fitter;
        
        // Create a prototype dataset.
        std::vector<lk::AbsBinningCPtr> axes;
        for(int i = 0; i < ndim; ++i) {
            lk::AbsBinningCPtr axis(new lk::UniformBinning(-range,+range,nbin));
            axes.push_back(axis);
        }
        lk::BinnedDataPtr prototype(new lk::BinnedData(axes));
        
        // Fill each bin of the prototype dataset with the model evaluated with sigma=sigma0
        std::vector<double> point(ndim);
        for(int index = 0; index < prototype->getNBinsTotal(); ++index) {
            prototype->getBinCenters(index,point);
            prototype->setData(index,model(point,sigma0));
        }
        int size(prototype->getNBinsWithData());
        
        // Loop over experiments
        for(int expt = 0; expt < nexpt; ++expt) {
        
            // Create our resampler.
            lk::BinnedDataResampler resampler, unweightedResampler;
        
            // Generate random datasets.
            if(verbose) {
                std::cout << "Generating " << ndata << " datasets with " << size << " bins..." << std::endl;
            }
            lk::CovarianceMatrixPtr covariance;
            for(int i = 0; i < ndata; ++i) {
                // Clone our prototype.
                lk::BinnedDataPtr data(prototype->clone());
                // Generate a random covariance matrix with determinant 1 for this dataset.
                if(0 == i || !sameCov) {
                    covariance = lk::generateRandomCovariance(size,covScale);
                    //covariance = lk::createDiagonalCovariance(size,0.1);
                }
                // Sample the covariance to generate random offset for each bin.
                boost::shared_array<double> offsets = covariance->sample(1);
                int nextOffset(0);
                for(lk::BinnedData::IndexIterator iter = data->begin(); iter != data->end(); ++iter) {
                    data->addData(*iter,offsets[nextOffset++]);
                }
                // The data is associated with a possibly rescaled covariance.
                if(covRescale != 1) covariance->applyScaleFactor(covRescale);
                data->setCovarianceMatrix(covariance);                
                // Add this dataset to our resampler.
                resampler.addObservation(data);
                if(unweighted) {
                    lk::BinnedDataPtr copy(data->clone());
                    copy->dropCovariance();
                    unweightedResampler.addObservation(copy);
                }
            }
        
            // Fit the combined data.
            lk::BinnedDataCPtr combined = resampler.combined();
            lk::FunctionMinimumPtr combinedFit = fitter.fit(combined);
            if(verbose) combinedFit->printToStream(std::cout);
            if(combinedFit->getStatus() != lk::FunctionMinimum::OK) continue;
            // Make a copy of the combined dataset's covariance.
            lk::CovarianceMatrixPtr combinedCovariance(new lk::CovarianceMatrix(*combined->getCovarianceMatrix()));

            // Dump the likelihood curve for the combined data.
            if(dump) {
                std::ofstream out("likely.dat");
                std::vector<double> params(1);
                for(int i = 0; i < 1000; ++i) {
                    params[0] = 2e-3*(i+1);
                    double nll= fitter(params);
                    double prob = std::exp(-nll);
                    out << params[0] << ' ' << prob << std::endl;
                }
                out.close();
            }
        
            // Generate MCMC samples from the combined likelihood.
            double mcfrac;
            {
                if(verbose) {
                    std::cout << "Generating " << nmc << " MCMC samples with skip "
                        << nskip << "..." << std::endl;
                }
                std::ofstream *out(0);
                if(dump) {
                    out = new std::ofstream("mcmc.dat");
                    *out << "sigma" << std::endl;
                }
                mcfrac = fitter.mcmc(combinedFit,nmc,nskip,sigma0,out);
                if(dump) out->close();
            }
            // Loop over bootstrap trials.
            double bsfrac(-1);
            int errors(0),warnings(0);
            {
                if(verbose) std::cout << "Generating " << ntrial << " bootstrap trials..." << std::endl;
                lk::Parameters fitted;
                std::ofstream *out(0);
                if(dump) {
                    out = new std::ofstream("bstrials.dat");
                    *out << "sigma" << std::endl;
                }
                int count = 0;
                for(int i = 0; i < ntrial; ++i) {
                    lk::BinnedDataPtr sample;
                    if(unweighted) {
                        sample = unweightedResampler.bootstrap(ndata,false);
                    }
                    else {
                        sample = resampler.bootstrap(ndata,!sameCov);
                    }
                    if(bsCombined) {
                        // Make sure our data vector is unweighted.
                        sample->getData(*sample->begin(),false);
                        // (Re)set the covariance to use with this data.
                        sample->setCovarianceMatrix(combinedCovariance);
                    }
                    lk::FunctionMinimumPtr sampleFit = fitter.fit(sample);
                    if(sampleFit->getStatus() == lk::FunctionMinimum::ERROR) {
                        errors++;
                        continue;
                    }
                    else if(sampleFit->getStatus() == lk::FunctionMinimum::WARNING) {
                        warnings++;
                    }
                    double sigma = sampleFit->getParameters()[0];
                    if(sigma < sigma0) count++;
                    if(dump) *out << sigma << std::endl;
                    if(verbose && (i+1)%1000==0) {
                        std::cout << "...completed " << (i+1) << " trials" << std::endl;
                    }
                }
                if(errors < ntrial) bsfrac = (double)count/(ntrial-errors);
                if(verbose) {
                    std::cout << errors << " errors, " << warnings << " warnings" << std::endl;
                }
                if(dump) out->close();
            }
            std::cout << mcfrac << ' ' << bsfrac << std::endl;
        }
    }
    catch(int e) { }
    /*
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        throw e;
    }
    */
    return 0;
}
