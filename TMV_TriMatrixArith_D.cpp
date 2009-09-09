
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  //
  // MultMM
  //

  template <class T, class Ta> void RowMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    const size_t N = B.colsize();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	B.row(i) += A.row(i,i+1,N) * B.Rows(i+1,N);
	B.row(i) *= alpha;
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin();
      for(size_t i=0; i<N; ++i,++Aii) {
	B.row(i) *= *Aii;
	B.row(i) += A.row(i,i+1,N) * B.Rows(i+1,N);
	B.row(i) *= alpha;
      }
    }
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    const size_t N = B.colsize();

    if (A.isunit()) 
      for(size_t j=0; j<N; ++j) B.Rows(0,j) += A.col(j,0,j) ^ B.row(j);
    else {
      CVIter<Ta> Ajj = A.diag().begin();
      for(size_t j=0; j<N; ++j,++Ajj) {
	B.Rows(0,j) += A.col(j,0,j) ^ B.row(j);
	B.row(j) *= *Ajj;
      }
    }
    B *= alpha;
  }

  template <class T, class Ta> void RowMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    const int N = B.colsize();

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i) {
	B.row(i) += A.row(i,0,i) * B.Rows(0,i);
	B.row(i) *= alpha;
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin()+N-1;
      for(int i=N-1; i>=0; --i,--Aii) {
	B.row(i) *= *Aii;
	B.row(i) += A.row(i,0,i) * B.Rows(0,i);
	B.row(i) *= alpha;
      }
    }
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    const int N = B.colsize();

    if (A.isunit())
      for(int j=N-1; j>=0; --j)
	B.Rows(j+1,N) += (A.col(j,j+1,N) ^ B.row(j));
    else {
      CVIter<Ta> Ajj = A.diag().begin()+N-1;
      for(int j=N-1; j>=0; --j,--Ajj) {
	B.Rows(j+1,N) += (A.col(j,j+1,N) ^ B.row(j));
	B.row(j) *= *Ajj;
      }
    }
    B *= alpha;
  }

  // In TMV_TriMatrixArith_A.cpp
  template <class T, class Ta> extern void MultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x);
  template <class T, class Ta> extern void MultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x);

  template <class T, class Ta> inline void NonBlasMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    if (B.isrm()) {
      if (A.iscm()) ColMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    } else {
      for(size_t j=0;j<B.rowsize();++j) MultEqMV(A,B.col(j));
      if (alpha != T(1)) B *= alpha;
    }
  }

  template <class T, class Ta> inline void NonBlasMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    if (B.isrm()) {
      if (A.iscm()) ColMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    } else {
      for(size_t j=0;j<B.rowsize();++j) MultEqMV(A,B.col(j));
      if (alpha != T(1)) B *= alpha;
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasMultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasMultEqMM(alpha,A,B); }
  template <class T, class Ta> inline void BlasMultEqMM(T alpha,
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasMultEqMM(alpha,A,B); }
  inline void BlasMultEqMM(double alpha, const GenUpperTriMatrix<double>& A, 
      const MatrixView<double>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_dtrmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasLower : CblasUpper,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(double alpha, const GenLowerTriMatrix<double>& A, 
      const MatrixView<double>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_dtrmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasUpper : CblasLower,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(complex<double> alpha,
      const GenUpperTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasLower : CblasUpper;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ztrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
  inline void BlasMultEqMM(complex<double> alpha,
      const GenLowerTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasUpper : CblasLower;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ztrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
#ifndef NOFLOAT
  inline void BlasMultEqMM(float alpha, const GenUpperTriMatrix<float>& A, 
      const MatrixView<float>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_strmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasLower : CblasUpper,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(float alpha, const GenLowerTriMatrix<float>& A, 
      const MatrixView<float>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_strmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasUpper : CblasLower,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(complex<float> alpha,
      const GenUpperTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasLower : CblasUpper;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ctrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
  inline void BlasMultEqMM(complex<float> alpha,
      const GenLowerTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasUpper : CblasLower;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ctrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta> inline void MultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    if (B.isconj()) 
      MultEqMM(CONJ(alpha),A.QuickConjugate(),B.QuickConjugate());
    else 
#ifdef BLAS
      if ( (A.isrm() || A.iscm()) && (B.isrm() || B.iscm()) )
	BlasMultEqMM(alpha,A,B);
      else
#endif
	NonBlasMultEqMM(alpha,A,B);
  }

  template <class T, class Ta> inline void MultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    if (B.isconj()) 
      MultEqMM(CONJ(alpha),A.QuickConjugate(),B.QuickConjugate());
    else 
#ifdef BLAS
      if ( (A.isrm() || A.iscm()) && (B.isrm() || B.iscm()) )
	BlasMultEqMM(alpha,A,B);
      else
#endif
	NonBlasMultEqMM(alpha,A,B);
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha==T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	Matrix<T,ColMajor> tempB = B;
	MultEqMM(alpha,A,tempB.QuickView());
	if (beta != T(1)) C *= beta;
	C += tempB;
      } else {
	if (beta == T(0)) {
	  C = B;
	  MultEqMM(alpha,A,C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  if (beta != T(1)) C *= beta;
	  C += tempB;
	}
      } 
    }
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha==T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	Matrix<T,ColMajor> tempB = B;
	MultEqMM(alpha,A,tempB.QuickView());
	if (beta != T(1)) C *= beta;
	C += tempB;
      } else {
	if (beta == T(0)) {
	  C = B;
	  MultEqMM(alpha,A,C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  if (beta != T(1)) C *= beta;
	  C += tempB;
	}
      } 
    }
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha==T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	Matrix<T,ColMajor> tempB = B;
	MultEqMM(alpha,A,tempB.QuickView());
	if (beta != T(1)) C *= beta;
	C += tempB;
      } else {
	if (beta == T(0)) {
	  C = B;
	  MultEqMM(alpha,A,C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  if (beta != T(1)) C *= beta;
	  C += tempB;
	}
      } 
    }
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha==T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	Matrix<T,ColMajor> tempB = B;
	MultEqMM(alpha,A,tempB.QuickView());
	if (beta != T(1)) C *= beta;
	C += tempB;
      } else {
	if (beta == T(0)) {
	  C = B;
	  MultEqMM(alpha,A,C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  if (beta != T(1)) C *= beta;
	  C += tempB;
	}
      } 
    }
  }

  template <class T, class Ta> void RowMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	B.row(i,i+1,N) += A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	B.row(i,i+1,N) *= alpha;
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin();
      TMVAssert(!B.isunit());
      for(size_t i=0; i<N; ++i,++Aii) {
	B.row(i,i,N) *= *Aii;
	B.row(i,i+1,N) += A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	B.row(i,i,N) *= alpha;
      }
    }
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit()) {
	for(size_t j=1; j<N; ++j) {
	  B.SubMatrix(0,j,j+1,N) += A.col(j,0,j) ^ B.row(j,j+1,N);
	  B.col(j,0,j) += A.col(j,0,j);
	}
      }
      else {
	for(size_t j=1; j<N; ++j) 
	  B.SubMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
      }
    } else {
      CVIter<Ta> Ajj = A.diag().begin();
      TMVAssert(!B.isunit());
      for(size_t j=0; j<N; ++j,++Ajj) {
	B.SubMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
	B.row(j,j,N) *= (*Ajj);
      }
    } 
    B *= alpha;
  }

  template <class T, class Ta> inline void MultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    if (B.isconj()) MultEqMM(CONJ(alpha),A.QuickConjugate(),
	B.QuickConjugate());
    else if (B.isrm()) {
      if (A.iscm()) ColMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    } else {
      if (B.isunit()) {
	// Then alpha = 1 and A.isunit
	for(size_t j=1;j<B.size();++j) {
	  MultEqMV(A.SubTriMatrix(0,j),B.col(j,0,j));
	  B.col(j,0,j) += A.col(j,0,j);
	}
      }
      else {
	for(size_t j=0;j<B.size();++j)
	  MultEqMV(A.SubTriMatrix(0,j+1),B.col(j,0,j+1));
	if (alpha != T(1)) B *= alpha;
      }
    }
  }

  template <class T, class Ta> void RowMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i) {
	B.row(i,0,i) += A.row(i,0,i) * B.SubTriMatrix(0,i);
	B.row(i,0,i) *= alpha;
      }
      if (!B.isunit()) B.diag() *= alpha;
    }
    else {
      CVIter<Ta> Aii = A.diag().begin()+N-1;
      TMVAssert(!B.isunit());
      for(size_t i=N-1; i>0; --i,--Aii) {
	B.row(i,0,i+1) *= *Aii;
	B.row(i,0,i) += A.row(i,0,i) * B.SubTriMatrix(0,i);
	B.row(i,0,i+1) *= alpha;
      }
    }
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit()) {
	for(int j=N-1; j>=0; --j) {
	  B.SubMatrix(j+1,N,0,j) += A.col(j,j+1,N) ^ B.row(j,0,j);
	  B.col(j,j+1,N) += A.col(j,j+1,N);
	}
      }
      else {
	for(int j=N-1; j>=0; --j) 
	  B.SubMatrix(j+1,N,0,j+1) += A.col(j,j+1,N) ^ B.row(j,0,j+1);
      }
    } else {
      CVIter<Ta> Ajj = A.diag().begin()+N-1;
      TMVAssert(!B.isunit());
      for(int j=N-1; j>=0; --j,--Ajj) {
	B.SubMatrix(j+1,N,0,j+1) += A.col(j,j+1,N) ^ B.row(j,0,j+1);
	B.row(j,0,j+1) *= (*Ajj);
      }
    } 
    B *= alpha;
  }

  template <class T, class Ta> inline void MultEqMM(T alpha,
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    if (B.isconj()) MultEqMM(CONJ(alpha),A.QuickConjugate(),
	B.QuickConjugate());
    else if (B.isrm()) {
      if (A.iscm()) ColMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    } else {
      if (B.isunit()) {
	// Then alpha = 1 and A.isunit
	for(size_t j=0;j<B.size()-1;++j) {
	  MultEqMV(A.SubTriMatrix(j+1,B.size()),B.col(j,j+1,B.size()));
	  B.col(j,j+1,B.size()) += A.col(j,j+1,B.size());
	}
      }
      else {
	for(size_t j=0;j<B.size();++j)
	  MultEqMV(A.SubTriMatrix(j,B.size()),B.col(j,j,B.size()));
	if (alpha != T(1)) B *= alpha;
      }
    }
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const UpperTriMatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    TMVAssert(A.size() == C.size());
    TMVAssert(A.size() == B.size());
    TMVAssert(!C.isunit() || (A.isunit() && B.isunit() && alpha == T(1)) );

    if (C.size() > 0) {
      if (alpha==T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (beta != T(0) || (C.SameStorageAs(A) && C.SameStorageAs(B))) {
	if (B.isrm()) {
	  if (B.dt() == UnitDiag) {
	    UpperTriMatrix<T,UnitDiag,RowMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    if (beta != T(1)) C *= beta;
	    C += tempB;
	  } else {
	    UpperTriMatrix<T,NonUnitDiag,RowMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    if (beta != T(1)) C *= beta;
	    C += tempB;
	  }
	} else {
	  if (B.dt() == UnitDiag) {
	    UpperTriMatrix<T,UnitDiag,ColMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    if (beta != T(1)) C *= beta;
	    C += tempB;
	  } else {
	    UpperTriMatrix<T,NonUnitDiag,ColMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    if (beta != T(1)) C *= beta;
	    C += tempB;
	  }
	}
      } else if (C.SameStorageAs(A)) {
	C = A;
	MultEqMM(alpha,B.QuickTranspose(),C.QuickTranspose());
      } else {
	C = B;
	MultEqMM(alpha,A,C);
      }
    }
  }

#define InstFile "TMV_TriMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


