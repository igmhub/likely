// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"

#include <cassert>

// Declare bindings to BLAS,LAPACK routines we need
extern "C" {
    // http://www.netlib.org/lapack/double/dpptrf.f
    void dpptrf_(char const *uplo, int const *n, double *ap, int *info);
    // http://www.netlib.org/lapack/double/dpptri.f
    void dpptri_(char const *uplo, int const *n, double *ap, int *info);
    // http://netlib.org/blas/dspmv.f
    void dspmv_(char const *uplo, int const *n, double const *alpha, double const *ap,
        double const *x, int const *incx, double const *beta, double *y, int const *incy);
    // http://www.netlib.org/blas/dsymm.f
    void dsymm_(char const *side, char const *uplo, int const *m, int const *n,
        double const *alpha, double const *a, int const *lda, double const *b,
        int const *ldb, double const *beta, double *c, int const *ldc);
}

namespace local = likely;

local::CovarianceMatrix::CovarianceMatrix(int size)
: _size(size), _compressed(false)
{
    if(size <= 0) {
        throw RuntimeError("CovarianceMatrix: expected size > 0.");
    }
    _ncov = (_size*(_size+1))/2;
    // We don't actually allocate any memory at this point. Wait until this is actually
    // necessary, and we know wether to allocate _cov or _icov.
}

local::CovarianceMatrix::~CovarianceMatrix() { }

void local::CovarianceMatrix::compress() const {
    if(_compressed) return;
    _compressed = true;
}

void local::CovarianceMatrix::uncompress() const {
    if(!_compressed) return;
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

void local::choleskyDecompose(std::vector<double> &matrix) {
    static char uplo('U');
    static int info(0);
    int n(matrix.size());
    dpptrf_(&uplo,&n,&matrix[0],&info);
    if(0 != info) {
        info = 0;
        throw RuntimeError("choleskyDecomposition: matrix is not positive definite.");
    }    
}

void local::invertCholesky(std::vector<double> &matrix) {
    static char uplo('U');
    static int info(0);
    int n(matrix.size());
    dpptri_(&uplo,&n,&matrix[0],&info);
    if(0 != info) {
        info = 0;
        throw RuntimeError("invertCholesky: symmetric matrix inversion failed.");
    }    
} 

void local::CovarianceMatrix::_changesCov() {
    // Do we have a valid matrix to change?
    if(_cov.empty()) {
        // Have we allocated anything yet?
        if(_icov.empty()) {
            // Allocate a covariance matrix initialized to zero.
            std::vector<double>(_ncov,0).swap(_cov);
        }
        else {
            // Try to invert the existing inverse covariance in place. This will throw a
            // RuntimeError in case the existing inverse covariance is only partially filled in.
            choleskyDecompose(_icov);
            invertCholesky(_icov);
            // Remove the existing inverse covariance (by swapping with _cov), since it will
            // become invalid after we update the the covariance.
            _cov.swap(_icov);
            // Remove any existing Cholesky decomposition since it will become invalid.
            if(!_cholesky.empty()) std::vector<double>().swap(_cholesky);
        }
    }
    assert(!_cov.empty());
    assert(_icov.empty());
    assert(_cholesky.empty());
}

double local::CovarianceMatrix::getCovariance(int row, int col) const {
    // Calculate the index corresponding to (row,col). This will throw a RuntimeError
    // in case of an invalid address, before we go any further.
    int index(symmetricMatrixIndex(row,col,_size));
    // Do we have a covariance matrix allocated yet?
    if(_cov.empty()) {
        if(_icov.empty()) {
            // Nothing has been allocated yet, so return zero since that is our
            // declared initial state.
            return 0;
        }
        else {
            // Try to invert the existing inverse covariance into _cov. This will throw a
            // RuntimeError in case the existing inverse covariance is only partially filled in.
            _cov = _icov;
            choleskyDecompose(_cov);
            _cholesky = _cov;
            invertCholesky(_cov);            
        }
    }
    return _cov[index];
}

double local::CovarianceMatrix::getInverseCovariance(int row, int col) const {
    return 0;
}

void local::CovarianceMatrix::setCovariance(int row, int col, double value) {
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
    
}
