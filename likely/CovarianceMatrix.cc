// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"
#include "likely/Random.h"

#include "boost/format.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

// Declare bindings to BLAS,LAPACK routines we need
extern "C" {
    // http://www.netlib.org/lapack/double/dpptrf.f
    void dpptrf_(char const *uplo, int const *n, double *ap, int *info);
    // http://www.netlib.org/lapack/double/dpptri.f
    void dpptri_(char const *uplo, int const *n, double *ap, int *info);
    // http://netlib.org/blas/dspmv.f
    void dspmv_(char const *uplo, int const *n, double const *alpha, double const *ap,
        double const *x, int const *incx, double const *beta, double *y, int const *incy);
    // http://www.netlib.org/blas/dtrmm.f
    void dtrmm_(char const *side, char const *uplo, char const *transa, const char *diag,
        int const *m, int const *n, double const *alpha, double const *a, int const *lda,
        double *b, int const *ldb);        
}

namespace local = likely;

local::CovarianceMatrix::CovarianceMatrix(int size)
: _size(size), _compressed(false), _nextSeed(0)
{
    if(size <= 0) {
        throw RuntimeError("CovarianceMatrix: expected size > 0.");
    }
    _ncov = (_size*(_size+1))/2;
    // We don't actually allocate any memory at this point. Wait until this is actually
    // necessary, and we know wether to allocate _cov or _icov.
}

local::CovarianceMatrix::CovarianceMatrix(std::vector<double> packed)
: _ncov(packed.size()), _compressed(false), _nextSeed(0)
{
    if(_ncov == 0) {
        throw RuntimeError("CovarianceMatrix: expected packed size > 0.");
    }
    _size = symmetricMatrixSize(_ncov);
    // Copy elements from the input array to our internal storage now. The first
    // call to setCovariance triggers the allocation of our internal storage.
    int index(0);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            setCovariance(row,col,packed[index++]);
        }
    }
}

local::CovarianceMatrix::~CovarianceMatrix() { }

size_t local::CovarianceMatrix::getMemoryUsage() const {
    return sizeof(*this) + sizeof(double)*(
        _cov.capacity() + _icov.capacity() + _cholesky.capacity() +
        _diag.capacity() + _offdiagIndex.capacity() + _offdiagValue.capacity());
}

std::string local::CovarianceMatrix::getMemoryState() const {
    return boost::str(boost::format("[%c%c%c%c%c%c] %d") %
        _tag('M',_cov) % _tag('I',_icov) % _tag('C',_cholesky) % _tag('D',_diag) %
        _tag('Z',_offdiagIndex) % _tag('V',_offdiagValue) % getMemoryUsage());
}

char local::CovarianceMatrix::_tag(char symbol, std::vector<double> const &vector) const {
    if(0 == vector.capacity()) return '-';
    if(0 == vector.size()) return std::tolower(symbol);
    return symbol;
}

bool local::CovarianceMatrix::compress() const {
    // Are we already compressed?
    if(_compressed) return false;
    // Do we still have valid compressed data?
    if(_diag.empty()) {
        // Reserve space for the diagonal elements, which cannot be compressed.
        _diag.reserve(_size);
        // Prepare to read the inverse covariance and check if anything been allocated yet.
        if(!_readsICov()) return false;
        // Loop over the upper-diagonal (row <= col) inverse matrix elements.
        int index(0);
        double value;
        for(int col = 0; col < _size; ++col) {
            for(int row = 0; row < col; ++row) {
                if(value = _icov[index]) {
                    _offdiagIndex.push_back(index);
                    _offdiagValue.push_back(value);
                }
                index++;
            }
            _diag.push_back(_icov[index++]);
        }
    }
    // Delete anything we don't need now.
    if(!_cov.empty()) std::vector<double>().swap(_cov);
    if(!_icov.empty()) std::vector<double>().swap(_icov);
    if(!_cholesky.empty()) std::vector<double>().swap(_cholesky);
    _compressed = true;
    return true;
}

void local::CovarianceMatrix::_uncompress() const {
    // Are we already decompressed?
    if(!_compressed) return;
    assert(0 == _cov.capacity());
    assert(0 == _icov.capacity());
    assert(0 == _cholesky.capacity());
    // Decompress the inverse covariance matrix.
    std::vector<double>(_ncov,0).swap(_icov);
    for(int k = 0; k < _offdiagIndex.size(); ++k) {
        _icov[_offdiagIndex[k]] = _offdiagValue[k];
    }
    for(int k = 0; k < _size; ++k) {
        _icov[(k*(k+3))/2] = _diag[k];
    }
    // Don't delete the uncompressed matrix data in case we can re-use it
    // because no changes are made before the next call to compress().
    _compressed = false;
}

int local::symmetricMatrixIndex(int row, int col, int size) {
    if(row < 0 || col < 0 || row >= size || col >= size) {
        throw RuntimeError("symmetricMatrixIndex: row or col out of range.");
    }
    // Ensure that row <= col
    if(row > col) std::swap(row,col);
    return row+(col*(col+1))/2;
}

int local::symmetricMatrixSize(int nelem) {
    int size = std::floor(std::sqrt(8*nelem+1)/2);
    if(nelem != (size*(size+1))/2) {
        throw RuntimeError("symmetricMatrixSize: invalid number of elements.");
    }
    return size;
}

void local::choleskyDecompose(std::vector<double> &matrix, int size) {
    static char uplo('U');
    static int info(0);
    if(0 == size) size = symmetricMatrixSize(matrix.size());
    dpptrf_(&uplo,&size,&matrix[0],&info);
    if(0 != info) {
        info = 0;
        throw RuntimeError("choleskyDecomposition: matrix is not positive definite.");
    }
}

void local::invertCholesky(std::vector<double> &matrix, int size) {
    static char uplo('U');
    static int info(0);
    if(0 == size) size = symmetricMatrixSize(matrix.size());
    dpptri_(&uplo,&size,&matrix[0],&info);
    if(0 != info) {
        info = 0;
        throw RuntimeError("invertCholesky: symmetric matrix inversion failed.");
    }
} 

void local::symmetricMatrixMultiply(std::vector<double> const &matrix,
std::vector<double> const &vector, std::vector<double> &result) {
    static char uplo('U');
    static int incr(1);
    static double alpha(1),beta(0);
    int size(vector.size());
    if(matrix.size() != (size*(size+1))/2) {
        throw RuntimeError("symmetricMatrixMultiply: incompatible matrix and vector sizes.");
    }
    // size result correctly (but do not need to zero elements since beta=0)
    std::vector<double>(size).swap(result);
    // See http://netlib.org/blas/dspmv.f
    dspmv_(&uplo,&size,&alpha,&matrix[0],&vector[0],&incr,&beta,&result[0],&incr);
}

void local::CovarianceMatrix::_changesCov() {
    _uncompress();
    // Any cached compressed matrix data is now invalid so delete it.
    if(!_diag.empty()) {
        // TODO: use resize(0) instead?
        std::vector<double>().swap(_diag);
        std::vector<double>().swap(_offdiagIndex);
        std::vector<double>().swap(_offdiagValue);
    }
    // Do we have a matrix to change?
    if(_cov.empty()) {
        // Have we allocated anything yet?
        if(_icov.empty()) {
            // Allocate a covariance matrix initialized to zero.
            std::vector<double>(_ncov,0).swap(_cov);
        }
        else {
            // Try to invert the existing inverse covariance in place. This will throw a
            // RuntimeError in case the existing inverse covariance is only partially filled in.
            choleskyDecompose(_icov,_size);
            invertCholesky(_icov,_size);
            // Remove the existing inverse covariance (by swapping with _cov), since it will
            // become invalid after we update the the covariance.
            _cov.swap(_icov);
            // Remove any existing Cholesky decomposition since it will become invalid.
            // TODO: use resize(0) instead here?
            if(!_cholesky.empty()) std::vector<double>().swap(_cholesky);
        }
    }
    else {
        // Delete any inverse covariance and Cholesky decomposition.
        if(!_icov.empty()) std::vector<double>().swap(_icov);
        if(!_cholesky.empty()) std::vector<double>().swap(_cholesky);
    }
    assert(!_cov.empty());
    assert(0 == _icov.capacity());
    assert(0 == _cholesky.capacity());
}

void local::CovarianceMatrix::_changesICov() {
    _uncompress();
    // Any cached compressed matrix data is now invalid so delete it.
    if(!_diag.empty()) {
        // TODO: use resize(0) instead?
        std::vector<double>().swap(_diag);
        std::vector<double>().swap(_offdiagIndex);
        std::vector<double>().swap(_offdiagValue);
    }
    // Do we have a matrix to change?
    if(_icov.empty()) {
        // Have we allocated anything yet?
        if(_cov.empty()) {
            // Allocate an inverse covariance matrix initialized to zero.
            std::vector<double>(_ncov,0).swap(_icov);
        }
        else {
            // Try to invert the existing covariance in place. This will throw a
            // RuntimeError in case the existing covariance is only partially filled in.
            if(_cholesky.empty()) {
                choleskyDecompose(_cov,_size);
                invertCholesky(_cov,_size);
                // Remove the existing covariance (by swapping with _icov), since it will
                // become invalid after we update the the covariance.
                _icov.swap(_cov);
            }
            else {
                // Use the cached Cholesky decomposition to invert the covariance. This
                // will clobber _cholesky but it is now invalid anyway.
                invertCholesky(_cholesky,_size);
                // Remove the existing _cholesky by swapping with _icov.
                _icov.swap(_cholesky);
                // Remove the existing _cov.
                std::vector<double>().swap(_cov);
            }
        }
    }
    else {
        // Delete and covariance and Cholesky decomposition.
        if(!_cov.empty()) std::vector<double>().swap(_cov);
        if(!_cholesky.empty()) std::vector<double>().swap(_cholesky);
    }
    assert(!_icov.empty());
    assert(0 == _cov.capacity());
    assert(0 == _cholesky.capacity());
}

bool local::CovarianceMatrix::_readsCov() const {
    _uncompress();
    // Do we have a covariance matrix allocated yet?
    if(_cov.empty()) {
        if(_icov.empty()) {
            // Nothing has been allocated yet.
            return false;
        }
        else {
            // Try to invert the existing inverse covariance into _cov. This will throw a
            // RuntimeError in case the existing inverse covariance is only partially filled in.
            _cov = _icov;
            choleskyDecompose(_cov,_size);
            // (we don't bother keeping the Cholesky decomposition of the inverse covariance)
            invertCholesky(_cov,_size);
        }
    }
    return true;
}

bool local::CovarianceMatrix::_readsICov() const {
    _uncompress();
    // Do we have an inverse covariance matrix allocated yet?
    if(_icov.empty()) {
        if(_cov.empty()) {
            // Nothing has been allocated yet.
            return false;
        }
        else {
            // Try to invert the existing covariance into _icov. This will throw a
            // RuntimeError in case the existing inverse covariance is only partially filled in.
            if(_cholesky.empty()) {
                // Calculate and save the covariance Cholesky decomposition.
                _icov = _cov;
                choleskyDecompose(_icov,_size);
                _cholesky = _icov;
            }
            else {
                // Use the cached Cholesky decomposition.
                _icov = _cholesky;
            }
            invertCholesky(_icov,_size);
        }
    }
    return true;
}

void local::CovarianceMatrix::_readsCholesky() const {
    // Make sure we have a packed Cholesky decomposition available.
    if(_cholesky.empty()) {
        if(!_readsCov()) {
            throw RuntimeError(
                "CovarianceMatrix: invalid Cholesky decomposition (no elements set yet).");
        }
        _cholesky = _cov;
        choleskyDecompose(_cholesky,_size);
    }    
}

double local::CovarianceMatrix::getCovariance(int row, int col) const {
    // Calculate the index corresponding to (row,col). This will throw a RuntimeError
    // in case of an invalid address, before we go any further.
    int index(symmetricMatrixIndex(row,col,_size));
    // Prepare to read from the covariance matrix, and return zero if nothing has
    // been allocated yet.
    if(!_readsCov()) return 0;
    return _cov[index];
}

double local::CovarianceMatrix::getInverseCovariance(int row, int col) const {
    // Calculate the index corresponding to (row,col). This will throw a RuntimeError
    // in case of an invalid address, before we go any further.
    int index(symmetricMatrixIndex(row,col,_size));
    // Prepare to read from the inverse covariance matrix, and return zero if nothing has
    // been allocated yet.
    if(!_readsICov()) return 0;
    return _icov[index];
}

void local::CovarianceMatrix::setCovariance(int row, int col, double value) {
    if(row == col && value <= 0) {
        throw RuntimeError("CovarianceMatrix: diagonal elements must be > 0.");
    }
    // Calculate the index corresponding to (row,col). This will throw a RuntimeError
    // in case of an invalid address, before we actually change anything.
    int index(symmetricMatrixIndex(row,col,_size));
    // Prepare to change the covariance matrix, which might throw a RuntimeError
    // if icov elements have already been set, but icov is not invertible.
    _changesCov();
    // Finally, set the new value here.
    _cov[index] = value;
}

void local::CovarianceMatrix::setInverseCovariance(int row, int col, double value) {
    if(row == col && value <= 0) {
        throw RuntimeError("CovarianceMatrix: diagonal elements must be > 0.");
    }
    // Calculate the index corresponding to (row,col). This will throw a RuntimeError
    // in case of an invalid address, before we actually change anything.
    int index(symmetricMatrixIndex(row,col,_size));
    // Prepare to change the inverse covariance matrix, which might throw a RuntimeError
    // if cov elements have already been set, but cov is not invertible.
    _changesICov();
    // Finally, set the new value here.
    _icov[index] = value;        
}

void local::CovarianceMatrix::multiplyByCovariance(std::vector<double> &vector) const {
    _readsCov();
    std::vector<double> result;
    symmetricMatrixMultiply(_cov,vector,result);
    vector.swap(result);
}

void local::CovarianceMatrix::multiplyByInverseCovariance(std::vector<double> &vector) const {
    _readsICov();
    std::vector<double> result;
    symmetricMatrixMultiply(_icov,vector,result);
    vector.swap(result);
}

double local::CovarianceMatrix::chiSquare(std::vector<double> const &delta) const {
    std::vector<double> icovDelta = delta;
    multiplyByInverseCovariance(icovDelta);
    double result(0);
    for(int k = 0; k < delta.size(); ++k) {
        result += delta[k]*icovDelta[k];
    }
    return result;
}

double local::CovarianceMatrix::sample(std::vector<double> &delta, Random *random) const {
    // Use the default generator if none was specified.
    if(0 == random) random = &Random::instance();
    // Clear the vector and prepare for filling.
    delta.resize(0);
    delta.reserve(_size);
    // Fill a vector delta' with uncorrelated random numbers and calculate the
    // -log(L) = delta.Cinv.delta/2 = (Linv.delta).(Linv.delta)/2 where delta'=Linv.delta
    // is the vector of uncorrelated values generated below and L is the lower-diagonal
    // Cholesky decomposition of the covariance matrix.
    double nll(0);
    std::vector<double> deltap;
    deltap.reserve(_size);
    for(int k = 0; k < _size; ++k) {
        double r(random->getNormal());
        deltap.push_back(r);
        nll += r*r;
    }
    _readsCholesky();
    // Add correlations via L.delta
    int index(0);
    for(int i = 0; i < _size; ++i) {
        double result(0);
        for(int j = 0; j <= i; ++j) {
            result += _cholesky[index++]*deltap[j];
        }
        delta.push_back(result);
    }
    return nll/2;
}

boost::shared_array<double> local::CovarianceMatrix::sample(int nsample, int seed) const {
    if(nsample <= 0) {
        throw RuntimeError("CovarianceMatrix: expected nsample > 0.");
    }
    _readsCholesky();
    // Temporarily transpose and expand the packed Cholesky matrix.
    boost::shared_array<double> expanded(new double[_size*_size]);
    double *ptr(&_cholesky[0]);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            // BLAS expects column-major ordering, i.e., with row increasing fastest.
            // Since we are transposing, col increases fastest here. We do not need to
            // initialize the upper diagonal elements since they will never be used.
            expanded[row*_size + col] = *ptr++;
        }
    }
    // Generate double-precision normally distributed (but uncorrelated) random numbers.
    std::size_t nrandom(nsample*_size), ngen(nrandom);
    if(0 == seed) seed = ++_nextSeed;
    boost::shared_array<double> array = Random::fillDoubleArrayNormal(ngen,seed);
    // Consider this array to be a rectangular matrix M of dimensions _size x nsample and
    // calculate (expanded).(M) to obtain a new matrix of dimensions _size x nsample
    // containing correlated residual vectors of length _size in each of its nsample columns.
    // The column-major ordering of BLAS means that the transformed residuals vectors
    // will be consecutive in memory.
    double alpha(1);
    char side = 'L', uplo = 'L', transa = 'N', diag = 'N';
    dtrmm_(&side,&uplo,&transa,&diag,&_size,&nsample,&alpha,
        expanded.get(),&_size,array.get(),&_size);
    return array;
}

void local::CovarianceMatrix::printToStream(std::ostream &os, std::string format) const {
    boost::format indexFormat("%5d"),valueFormat(format);
    for(int row = 0; row < _size; ++row) {
        os << (indexFormat % row);
        for(int col = 0; col <= row; ++col) {
            os << ' ' << (valueFormat % getCovariance(row,col));
        }
        os << std::endl;
    }
}

void local::CovarianceMatrix::addInverse(CovarianceMatrix const &other, double weight) {
    if(weight <= 0) {
        throw RuntimeError("CovarianceMatrix::addInverse: expected weight > 0.");
    }
    if(other.getSize() != _size) {
        throw RuntimeError("CovarianceMatrix::addInverse: incompatible sizes.");
    }
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            double otherValue = other.getInverseCovariance(row,col);
            double myValue = getInverseCovariance(row,col);
            setInverseCovariance(row,col,myValue+weight*otherValue);
        }
    }
}
