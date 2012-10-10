// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/CovarianceMatrix.h"
#include "likely/RuntimeError.h"
#include "likely/Random.h"

#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/smart_ptr.hpp"

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
    // http://www.netlib.org/blas/dsyrk.f
    void dsyrk_(char const *uplo, char const *trans, int const *n, int const *k,
        double const *alpha, double const *a, int const *lda, double const *beta,
        double *c, int const *ldc);
    // http://www.netlib.org/lapack/double/dspevd.f
    void dspevd_(char const *jobz, char const *uplo, int const *n, double *ap, double *w,
        double *z, int const *ldz, double *work, int const *lwork, int *iwork,
        int const *liwork, int *info);
}

namespace local = likely;

local::CovarianceMatrix::CovarianceMatrix(int size)
: _size(size), _compressed(false), _logDeterminant(0)
{
    if(size <= 0) {
        throw RuntimeError("CovarianceMatrix: expected size > 0.");
    }
    _ncov = (_size*(_size+1))/2;
    // We don't actually allocate any memory at this point. Wait until this is actually
    // necessary, and we know wether to allocate _cov or _icov.
}

local::CovarianceMatrix::CovarianceMatrix(std::vector<double> packed)
: _ncov(packed.size()), _compressed(false), _logDeterminant(0)
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

local::CovarianceMatrix& local::CovarianceMatrix::operator=(CovarianceMatrix other) {
    swap(*this,other);
    return *this;
}

void local::swap(CovarianceMatrix& a, CovarianceMatrix& b) {
    // Enable argument-dependent lookup (ADL)
    using std::swap;
    swap(a._size,b._size);
    swap(a._ncov,b._ncov);
    swap(a._logDeterminant,b._logDeterminant);
    swap(a._compressed,b._compressed);
    swap(a._cov,b._cov);
    swap(a._icov,b._icov);
    swap(a._cholesky,b._cholesky);
    swap(a._diag,b._diag);
    swap(a._offdiagIndex,b._offdiagIndex);
    swap(a._offdiagValue,b._offdiagValue);
}

size_t local::CovarianceMatrix::getMemoryUsage() const {
    return sizeof(*this) + sizeof(double)*(
        _cov.capacity() + _icov.capacity() + _cholesky.capacity() +
        _diag.capacity() + _offdiagIndex.capacity() + _offdiagValue.capacity());
}

std::string local::CovarianceMatrix::getMemoryState() const {
    return boost::str(boost::format("[%c%c%c%c%c%c%c] %d") %
        _tag('M',_cov) % _tag('I',_icov) % _tag('C',_cholesky) % (_logDeterminant == 0 ? '-':'L') %
        _tag('D',_diag) % _tag('Z',_offdiagIndex) % _tag('V',_offdiagValue) % getMemoryUsage());
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
    // Don't delete the compressed matrix data in case we can re-use it
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

double local::choleskyDecompose(std::vector<double> &matrix, int size) {
    static char uplo('U');
    static int info(0);
    if(0 == size) size = symmetricMatrixSize(matrix.size());
    dpptrf_(&uplo,&size,&matrix[0],&info);
    if(0 != info) {
        info = 0;
        throw RuntimeError("choleskyDecomposition: matrix is not positive definite.");
    }
    // Calculate and the product of diagonal Cholesky matrix elements squared.
    double logdet(0);
    for(int index = 0; index < size; ++index) {
        logdet += 2*std::log(matrix[symmetricMatrixIndex(index,index,size)]);
    }
    return logdet;
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

void local::matrixSquare(std::vector<double> const &matrix, std::vector<double> &result,
bool transposeLeft, int size) {
    static char uplo('U');
    static int info(0);
    static double alpha(1),beta(0);
    // Calculate the matrix size, if necessary.
    if(0 == size) size = symmetricMatrixSize(matrix.size());
    // Calculate Mt.M or M.Mt ?
    char trans = transposeLeft ? 'T' : 'N';
    boost::shared_array<double> unpackedResult(new double [size*size]);
    dsyrk_(&uplo,&trans,&size,&size,&alpha,&matrix[0],&size,&beta,unpackedResult.get(),&size);
    // Pack the result back into 'U' format.
    result.resize(0);
    result.reserve((size*(size+1))/2);
    for(int col = 0; col < size; ++col) {
        int base(col*size);
        for(int row = 0; row <= col; ++row) {
            result.push_back(unpackedResult[base++]);
        }
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

void local::symmetricMatrixEigenSolve(std::vector<double> const &matrix,
std::vector<double> &eigenvalues, std::vector<double> &eigenvectors, int size) {
    static char jobz('V'), uplo('U');
    static int info(0);
    // Calculate the matrix size if it was not provided.
    if(0 == size) size = symmetricMatrixSize(matrix.size());
    // Allocate space for the eigenvalues and vectors.
    eigenvalues.resize(size), eigenvectors.resize(size*size);
    {
        // copy the input matrix since the algorithm overwrites it
        std::vector<double> matrixCopy(matrix);
        // allocate temporory workspaces
        int workSize(1+6*size+size*size), iworkSize(3+5*size);
        boost::scoped_array<double> work(new double[workSize]);
        boost::scoped_array<int> iwork(new int[iworkSize]);
        dspevd_(&jobz,&uplo,&size,&matrixCopy[0],&eigenvalues[0],&eigenvectors[0],&size,
            &work[0],&workSize,&iwork[0],&iworkSize,&info);
        if(0 != info) {
            throw RuntimeError("symmetricMatrixEigenSolve: failed with info = " +
                boost::lexical_cast<std::string>(info));
            info = 0;
        }
        // cleanup temporary storage by closing this scope
    }   
}

void local::CovarianceMatrix::prune(std::set<int> const &keep) {
    int newSize(keep.size());
    if(newSize == getSize()) return;
    
    if(!_readsCov()) {
        throw RuntimeError("CovarianceMatrix::prune: no elements have been set.");
    }
    _changesCov();
    
    std::set<int>::const_iterator nextOldCol(keep.begin());
    for(int newCol = 0; newCol < newSize; ++newCol) {
        int oldCol = *nextOldCol++;        
        std::set<int>::const_iterator nextOldRow(keep.begin());
        for(int newRow = 0; newRow <= newCol; ++newRow) {
            int oldRow = *nextOldRow++;
            int newIndex = symmetricMatrixIndex(newRow, newCol, newSize);            
            int oldIndex = symmetricMatrixIndex(oldRow, oldCol, getSize());

            //std::cout << "prune: " << newCol << ',' << newRow << ' ' << oldCol << ',' << oldRow
            //    << ' ' << newIndex << ',' << oldIndex << std::endl;
            assert(oldIndex >= newIndex);

            _cov[newIndex] = _cov[oldIndex];
        }
    }
    _size = newSize;
    _ncov = (newSize*(newSize+1))/2;
    _cov.resize(_ncov);

    assert(0 == _icov.capacity());
    assert(0 == _cholesky.capacity());
    assert(0 == _diag.capacity());
    assert(0 == _offdiagIndex.capacity());
    assert(0 == _offdiagValue.capacity());
}

void local::CovarianceMatrix::_changesCov() {
    _uncompress();
    // Any cached determinant is now invalid.
    _logDeterminant = 0;
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
            _logDeterminant = -choleskyDecompose(_icov,_size);
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
    // Any cached determinant is now invalid.
    _logDeterminant = 0;
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
                _logDeterminant = +choleskyDecompose(_cov,_size);
                // No need to save this Cholesky decomposition since it will be invalid soon.
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
            _logDeterminant = -choleskyDecompose(_cov,_size);
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
                _logDeterminant = +choleskyDecompose(_icov,_size);
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
        _logDeterminant = +choleskyDecompose(_cholesky,_size);
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

local::CovarianceMatrix &local::CovarianceMatrix::setCovariance(int row, int col, double value) {
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
    // Return a self reference to allow chaining.
    return *this;
}

local::CovarianceMatrix &local::CovarianceMatrix::setInverseCovariance(int row, int col, double value) {
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
    // Return a self reference to allow chaining.
    return *this;
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

double local::CovarianceMatrix::chiSquareModes(std::vector<double> const &delta,
std::vector<double> &eigenvalues, std::vector<double> &eigenvectors,
std::vector<double> &chi2modes) const {
    // Solve our eigensystem for Cinv
    // TODO: if only C is available, solve its eigensystem instead, remembering to transform
    // lambda -> 1/lambda and that eigenvalue ordering is therefore reversed.
    _readsICov();
    symmetricMatrixEigenSolve(_icov,eigenvalues,eigenvectors,_size);
    // Loop over eigenmodes of Cinv
    double chi2(0);
    chi2modes.resize(0);
    chi2modes.reserve(_size);
    for(int i = 0; i < _size; ++i) {
        // Calculate the dot product of eigenvector i with delta
        double dotprod(0);
        for(int j = 0; j < _size; ++j) {
            dotprod += eigenvectors[i*_size + j]*delta[j];
        }
        // Calculate and save the contribution to chi2 due to this eigenmode.
        double chi2i = dotprod*dotprod*eigenvalues[i];
        // Replace lambda with 1/lambda so that we return decreasing eigenvalues of cov
        // instead of increasing eigenvalues of icov.
        eigenvalues[i] = 1/eigenvalues[i];
        chi2modes.push_back(chi2i);
        chi2 += chi2i;
    }
    return chi2;
}

double local::CovarianceMatrix::sample(std::vector<double> &delta, RandomPtr random) const {
    // Use the default generator if none was specified.
    if(!random) random = Random::instance();
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

void local::CovarianceMatrix::replaceWithTripleProduct(CovarianceMatrix const &other) {
    if(other.getSize() != _size) {
        throw RuntimeError("CovarianceMatrix::addInverse: incompatible sizes.");
    }
    // Any cached compressed matrix data is now invalid so delete it.
    if(!_diag.empty()) {
        // TODO: use resize(0) instead?
        std::vector<double>().swap(_diag);
        std::vector<double>().swap(_offdiagIndex);
        std::vector<double>().swap(_offdiagValue);
    }
    // Instead of calculating C -> A.Cinv.A we calculate Cinv -> Ainv.C.Ainv using:
    //
    //   Ainv.C.Ainv = Ainv.U*.U.Ainv = (U.Ainv)*.(U.Ainv)
    //
    // where U is the upper-diagonal Cholesky decomposition of C that we store in _cholesky.
    // First, calculate the elements of U = _cholesky, if necessary.
    _readsCholesky();
    
    // Free up any _cov or _icov storage now, before we allocate new temporary storage.
    if(!_cov.empty()) std::vector<double>().swap(_cov);
    if(!_icov.empty()) std::vector<double>().swap(_icov);

    // Next, multiply U.Ainv using the BLAS DTRMM routine which is optimized for the
    // upper triangular form of U, but not optimized for the symmetry of Ainv.
    // DTRMM needs both matrices to be unpacked first. Do this in two separate loops
    // so we can free the _cholesky memory before allocating the second temporary array.
    boost::shared_array<double> unpackedCholesky(new double [_size*_size]);
    double *choleskyPtr(&_cholesky[0]);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            unpackedCholesky[col*_size + row] = *choleskyPtr++;
        }
    }
    std::vector<double>().swap(_cholesky);
    boost::shared_array<double> unpackedOther(new double [_size*_size]);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row < col; ++row) {
            unpackedOther[row*_size + col] = unpackedOther[col*_size + row] =
                other.getInverseCovariance(row,col);
        }
        unpackedOther[col*_size + col] = other.getInverseCovariance(col,col);
    }

    double alpha(1);
    char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
    dtrmm_(&side,&uplo,&transa,&diag,&_size,&_size,&alpha,
        unpackedCholesky.get(),&_size,unpackedOther.get(),&_size);

    // Now calculate B*.B where B = U.Ainv is the contents of unpackedOther. Use the
    // BLAS DSYRK routine which knows that the result is symmetric. Save the result
    // into the same memory that we allocated above for the unpackedCholesky decomposition.
    double *unpackedResult = unpackedCholesky.get();
    char trans = 'T';
    double beta(0);
    dsyrk_(&uplo,&trans,&_size,&_size,&alpha,unpackedOther.get(),&_size,
        &beta,unpackedResult,&_size);
    
    // Finally, pack the result back into our inverse covariance.
    _icov.resize(0);
    _icov.reserve(_ncov);
    for(int col = 0; col < _size; ++col) {
        for(int row = 0; row <= col; ++row) {
            _icov.push_back(unpackedResult[col*_size + row]);
        }
    }
}

local::CovarianceMatrixPtr local::generateRandomCovariance(int size, double scale, RandomPtr random) {
    if(size <= 0) {
        throw RuntimeError("generateRandomCovariance: expected size > 0.");
    }
    if(scale <= 0) {
        throw RuntimeError("generateRandomCovariance: expected scale > 0.");
    }
    // Use the default generator if none was specified.
    if(!random) random = Random::instance();
    // Initialize the storage we will need.
    int sizeSq(size*size);
    CovarianceMatrixPtr C(new CovarianceMatrix(size));
    boost::shared_array<double> M;
    std::vector<double> MtM(sizeSq);
    // Loop over trials to generate a positive-definite random matrix.
    int ntrials(0), maxtrials(10);
    size_t nrandom(sizeSq);
    while(ntrials++ < maxtrials) {
        // Generate a random matrix M
        M = random->fillDoubleArrayUniform(nrandom);
        // Offset [0,1) to [-0.5,+0.5). The range used here is irrelevant
        // since we will be rescaling to get the desired determinant. However, the choice of a
        // uniform distribution does determine the distribution properties of the generated
        // covariance matrices. Might want to provide an option for e.g, Gaussian instead?
        for(int index = 0; index < sizeSq; ++index) M[index] -= 0.5;

        // Mt.M is positive definite iff M is invertible (i.e., has full rank and no zero singular values)
        // At this point, we can either calculate the singular values with BLAS DGESVD or go ahead and
        // calculate Mt.M and then check that it can be Cholesky decomposed. We do the latter since
        // the Cholesky decomposition will be cached and is potentially useful.
    
        // Multiply Mt.M to get a symmetric matrix that will be positive definite if M is invertible.
        char uplo = 'U', trans = 'T';
        double alpha(1),beta(0);
        dsyrk_(&uplo,&trans,&size,&size,&alpha,&M[0],&size,&beta,&MtM[0],&size);
    
        // Copy the upper triangle of MtM into a new CovarianceMatrix.
        // Loop over elements.
        for(int col = 0; col < size; ++col) {
            for(int row = 0; row <= col; ++row) {
                C->setCovariance(row,col,MtM[col*size + row]);
            }
        }
        // Calculate the re-scaling factor required to get the requested determinant. This
        // will throw a RuntimeError in case our original M was not invertible.
        try {
            double logdet = C->getLogDeterminant();
            double rescale = scale*std::exp(-logdet/size);
            C->applyScaleFactor(rescale);
            break;
        }
        catch(RuntimeError const &e) {
            // Oh well, try again.
        }
    }
    if(ntrials == maxtrials) {
        // This is really unlikely, so something is probably wrong.
        throw RuntimeError("generateRandomCovariance: failed after maxtrials.");
    }
    return C;
}

boost::shared_array<double> local::CovarianceMatrix::sample(int nsample, RandomPtr random) const {
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
    // Use the default generator if none was specified.
    if(!random) random = Random::instance();
    // Generate double-precision normally distributed (but uncorrelated) random numbers.
    std::size_t nrandom(nsample*_size), ngen(nrandom);
    boost::shared_array<double> array = random->fillDoubleArrayNormal(ngen);
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

void local::CovarianceMatrix::printToStream(std::ostream &os, bool normalized, std::string format,
std::vector<std::string> const &labels) const {
    if(labels.size() > 0 && labels.size() != _size) {
        throw RuntimeError("CovarianceMatrix::printToStream: unexpected number of labels.");
    }
    boost::format indexFormat("%5d :"),labelFormat("%20s :"),valueFormat(format);
    for(int row = 0; row < _size; ++row) {
        os << ((labels.size() > 0) ? (labelFormat % labels[row]) : (indexFormat % row));
        for(int col = 0; col <= row; ++col) {
            double value = getCovariance(row,col);
            if(normalized) {
                if(row == col) {
                    value = std::sqrt(value);
                }
                else {
                    value /= std::sqrt(getCovariance(row,row)*getCovariance(col,col));
                }
            }
            os << ' ' << (valueFormat % value);
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
    if(other.isCompressed()) {
        _changesICov();
        for(int k = 0; k < _size; ++k) {
            _icov[(k*(k+3))/2] += weight*other._diag[k];
        }
        for(int k = 0; k < other._offdiagIndex.size(); ++k) {
            _icov[other._offdiagIndex[k]] += weight*other._offdiagValue[k];
        }
    }
    else {
        for(int col = 0; col < _size; ++col) {
            for(int row = 0; row <= col; ++row) {
                double otherValue = other.getInverseCovariance(row,col);
                double myValue = getInverseCovariance(row,col);
                setInverseCovariance(row,col,myValue+weight*otherValue);
            }
        }
    }
}

int local::CovarianceMatrix::getNElements() const {
    // Prepare to read from the covariance matrix, and return zero if nothing has
    // been allocated yet.
    if(!_readsCov()) return 0;
    // Loop over all elements.
    int nelem(0);
    for(int index = 0; index < _ncov; ++index) {
        if(_cov[index] != 0) nelem++;
    }
    return nelem;
}

bool local::CovarianceMatrix::isPositiveDefinite() const {
    // Positive definiteness is equivalent to having a valid Cholesky decomposition for
    // either C or Cinv. Our implementation of getLogDeterminant() already keeps track
    // of whether we have a valid decomposition and efficiently triggers a new decomposition
    // when necessary, so use that.
    try {
        getLogDeterminant();
        return true;
    }
    catch(RuntimeError const &e) {
        return false;
    }
}

double local::CovarianceMatrix::getLogDeterminant() const {
    // Only do the minimum work necessary...
    if(0 == _logDeterminant) {
        _uncompress();
        // If we don't have a cached value then we have at most one of _icov or _cov,
        // but not both. Do a Cholesky decomposition of whatever we have.
        assert(_cholesky.empty());
        assert(_cov.empty() || _icov.empty());
        if(!_cov.empty()) {
            // Calculate and save the covariance Cholesky decomposition now.
            _cholesky = _cov;
            _logDeterminant = +choleskyDecompose(_cholesky,_size);
        }
        else if(!_icov.empty()) {
            // Calculate the inverse covariance Cholesky decomposition now.
            _cholesky = _icov;
            _logDeterminant = -choleskyDecompose(_cholesky,_size);
            // Don't keep this decomposition, since this was the inverse.
            std::vector<double>().swap(_cholesky);
        }
        else {
            throw RuntimeError("CovarianceMatrix::getLogDeterminant: no elements have been set.");
        }
    }
    assert(_logDeterminant != 0);
    return _logDeterminant;
}

void local::CovarianceMatrix::applyScaleFactor(double scaleFactor) {
    if(scaleFactor <= 0) {
        throw RuntimeError("CovarianceMatrix::applyScaleFactor: expected scaleFactor > 0.");
    }
    // We could actually do this on a compressed object - maybe later...
    _uncompress();
    // Transform whatever vectors we have using the appropriate scale.
    if(!_cov.empty()) {
        double scale(scaleFactor);
        for(int index = 0; index < _ncov; ++index) _cov[index] *= scale;
    }
    if(!_icov.empty()) {
        double scale(1/scaleFactor);
        for(int index = 0; index < _ncov; ++index) _icov[index] *= scale;
    }
    if(!_cholesky.empty()) {
        double scale(std::sqrt(scaleFactor));
        for(int index = 0; index < _ncov; ++index) _cholesky[index] *= scale;
    }
    if(_logDeterminant != 0) _logDeterminant += _size*std::log(scaleFactor);
}

local::CovarianceMatrixPtr local::createDiagonalCovariance(int size, double diagonalValue) {
    if(size <= 0) {
        throw RuntimeError("createDiagonalCovariance: expected size > 0.");
    }
    if(diagonalValue <= 0) {
        throw RuntimeError("createDiagonalCovariance: expected diagonalValue > 0.");
    }
    CovarianceMatrixPtr C(new CovarianceMatrix(size));
    for(int k = 0; k < size; ++k) C->setCovariance(k,k,diagonalValue);
    return C;
}

local::CovarianceMatrixPtr local::createDiagonalCovariance(std::vector<double> diagonalValues) {
    int size(diagonalValues.size());
    CovarianceMatrixPtr C(new CovarianceMatrix(size));
    for(int k = 0; k < size; ++k) C->setCovariance(k,k,diagonalValues[k]);
    return C;
}
