// Created 17-Apr-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_ABS_COVARIANCE_MATRIX
#define LIKELY_ABS_COVARIANCE_MATRIX

namespace likely {
    // Represents an abstract interface to a covariance matrix.
	class AbsCovarianceMatrix {
	public:
		virtual ~AbsCovarianceMatrix();
		// Returns the fixed size of this covariance matrix.
        int getSize() const;
        // Requests that this covariance matrix be compressed to reduce its memory usage,
        // if possible. Returns immediately if we are already compressed.
        void compress();
        // Undoes any compression. Returns immediately if we are already uncompressed.
        void uncompress();
        // Returns true if this covariance matrix is currently compressed.
        bool isCompressed() const;
	protected:
	    // Creates a new size-by-size covariance matrix. Throws a RuntimeError if size <= 0.
		explicit AbsCovarianceMatrix(int size);
		// Does the actual work of compressing a matrix, if possible.
        virtual void doCompress() = 0;
        // Does the actual work of undoing any compression performed by doCompress().
        virtual void doUncompress() = 0;
	private:
        int _size;
        bool _compressed;
	}; // AbsCovarianceMatrix
	
    inline int AbsCovarianceMatrix::getSize() const { return _size; }
    
    inline void AbsCovarianceMatrix::compress() {
        if(!_compressed) { doCompress(); _compressed = true; }
    }
    inline void AbsCovarianceMatrix::uncompress() {
        if(_compressed) { doUncompress(); _compressed = false; }
    }
    
    inline bool AbsCovarianceMatrix::isCompressed() const { return _compressed; }

} // likely

#endif // LIKELY_ABS_COVARIANCE_MATRIX
