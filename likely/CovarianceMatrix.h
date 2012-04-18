// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_COVARIANCE_MATRIX
#define LIKELY_COVARIANCE_MATRIX

#include <vector>

namespace likely {
    // Represents an abstract interface to a covariance matrix.
	class CovarianceMatrix {
	public:
	    // Creates a new size-by-size covariance matrix with all elements initialized to zero.
	    // Throws a RuntimeError if size <= 0. Note that the matrix created by this constructor
	    // is not valid until sufficient elements have been set to make it positive definite.
		explicit CovarianceMatrix(int size);
		virtual ~CovarianceMatrix();
		// Returns the fixed size of this covariance matrix.
        int getSize() const;
        // Returns the specified (inverse) covariance matrix element or throws a RuntimeError.
        // (row,col) and (col,row) give the same result by construction.
        double getCovariance(int row, int col) const;
        double getInverseCovariance(int row, int col) const;
        // Sets the specified (inverse) covariance matrix element or throws a RuntimeError.
        // Row and column indices should be in the range [0,size-1]. Setting any element with
        // row != col will also set the symmetric element in the matrix.
        void setCovariance(int row, int col, double value);
        void setInverseCovariance(int row, int col, double value);
        // Requests that this covariance matrix be compressed to reduce its memory usage,
        // if possible. Returns immediately if we are already compressed.
        void compress();
        // Undoes any compression. Returns immediately if we are already uncompressed.
        void uncompress();
        // Returns true if this covariance matrix is currently compressed.
        bool isCompressed() const;
	private:
        // Prepares to change at least one element of the _cov array.
        void _changesCov();
        int _size, _ncov;
        bool _compressed;
        mutable std::vector<double> _cov, _icov, _cholesky;
	}; // CovarianceMatrix
	
    inline int CovarianceMatrix::getSize() const { return _size; }
    
    inline bool CovarianceMatrix::isCompressed() const { return _compressed; }

    // Returns the array offset index for the BLAS packed symmetric matrix format
    // described at http://www.netlib.org/lapack/lug/node123.html or throws a
    // RuntimeError for invalid row or col inputs.
    int symmetricMatrixIndex(int row, int col, int size);
    // Performs a Cholesky decomposition in place of a symmetric positive definite matrix
    // or throws a RuntimeError if the matrix is not positive definite. The input matrix
    // is assumed to be in the BLAS packed format implied by packedMatrixIndex(row,col).
    static void choleskyDecompose(std::vector<double> &matrix);
    // Inverts a symmetric positive definite matrix in place, or throws a RuntimeError.
    // The input matrix should already be Cholesky decomposed and in the BLAS packed format
    // implied by packedMatrixIndex(row,col), e.g. by first calling _choleskyDecompose(matrix).
    static void invertCholesky(std::vector<double> &matrix);

} // likely

#endif // LIKELY_COVARIANCE_MATRIX
