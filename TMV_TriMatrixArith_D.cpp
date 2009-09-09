
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t MM_BLOCKSIZE = 64;
#endif

  // MJ: Compare with ATLAS code - try to speed up.

  //
  // MultMM
  //

  template <class T, class Ta> void RowMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RowMultEq Upper alpha = "<<alpha<<endl;
#endif
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
	B.row(i) = *Aii * B.row(i) + A.row(i,i+1,N) * B.Rows(i+1,N);
	B.row(i) *= alpha;
      }
    }
  }

  template <class T, class Ta> void OPMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"OPMultEq Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    const size_t N = B.colsize();

    if (A.isunit()) 
      for(size_t j=0; j<N; ++j) 
	B.Rows(0,j) += A.col(j,0,j) ^ B.row(j);
    else {
      CVIter<Ta> Ajj = A.diag().begin();
      for(size_t j=0; j<N; ++j,++Ajj) {
	B.Rows(0,j) += A.col(j,0,j) ^ B.row(j);
	B.row(j) *= *Ajj;
      }
    }
    B *= alpha;
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"ColMultEq Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    for(size_t j=0;j<B.rowsize();++j) 
      B.col(j) = alpha * A * B.col(j);
  }

  template <class T, class Ta> void RowMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RowMultEq Lower alpha = "<<alpha<<endl;
#endif
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
	// B.row(i) = alpha * (*Aii * B.row(i) + A.row(i,0,i) * B.Rows(0,i));
	T aa = *Aii;
	if (alpha != T(1)) aa *= alpha;
	MultMV(alpha,B.Rows(0,i).Transpose(),A.row(i,0,i),aa,B.row(i));
      }
    }
  }

  template <class T, class Ta> void OPMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"OPMultEq Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);
    const int N = B.colsize();

    if (A.isunit())
      for(int j=N-1; j>=0; --j)
	B.Rows(j+1,N) += A.col(j,j+1,N) ^ B.row(j);
    else {
      CVIter<Ta> Ajj = A.diag().begin()+N-1;
      for(int j=N-1; j>=0; --j,--Ajj) {
	B.Rows(j+1,N) += A.col(j,j+1,N) ^ B.row(j);
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

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"ColMultEq Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    for(size_t j=0;j<B.rowsize();++j) 
      B.col(j) = alpha * A * B.col(j);
  }

  template <class T, class Ta> inline void NonBlasMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    if (B.isrm()) 
      if (A.iscm()) OPMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
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

    if (B.isrm()) 
      if (A.iscm()) OPMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
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
    TMVAssert(B.ct() == NonConj);
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
    TMVAssert(B.ct() == NonConj);
#ifdef BLAS
    if ( (A.isrm() || A.iscm()) && (B.isrm() || B.iscm()) )
      BlasMultEqMM(alpha,A,B);
    else
#endif
      NonBlasMultEqMM(alpha,A,B);
  }

  template <class T, class Ta, class Tb> void RowAddMultMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowAddMult Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    const size_t N = C.colsize();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	C.row(i) += alpha * A.row(i,i+1,N) * B.Rows(i+1,N);
	C.row(i) += alpha * B.row(i);
      }
    }
    else {
      for(size_t i=0; i<N; ++i) 
	C.row(i) += alpha * A.row(i,i,N) * B.Rows(i,N);
    }
  }

  template <class T, class Ta, class Tb> void OPAddMultMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPAddMult Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    const size_t N = C.colsize();

    if (A.isunit()) 
      for(size_t j=0; j<N; ++j) {
	C.Rows(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
	C.row(j) += alpha * B.row(j);
      }
    else {
      for(size_t j=0; j<N; ++j) 
	C.Rows(0,j+1) += alpha * A.col(j,0,j+1) ^ B.row(j);
    }
  }

  template <class T, class Ta, class Tb> void ColAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColAddMult Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    for(size_t j=0;j<C.rowsize();++j) 
      C.col(j) += alpha * A * B.col(j);
  }

  template <class T, class Ta, class Tb> void RowAddMultMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowAddMult Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    const size_t N = C.colsize();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	C.row(i) += alpha * A.row(i,0,i) * B.Rows(0,i);
	C.row(i) += alpha * B.row(i);
      }
    }
    else {
      for(size_t i=0; i<N; ++i) 
	C.row(i,0,i+1) += alpha * A.row(i,0,i+1) * B.Rows(0,i+1);
    }
  }

  template <class T, class Ta, class Tb> void OPAddMultMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPAddMult Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    const size_t N = C.colsize();

    if (A.isunit())
      for(size_t j=0; j<N; ++j) {
	C.Rows(j+1,N) += alpha * A.col(j,j+1,N) ^ B.row(j);
	C.row(j) += alpha * B.row(j);
      }
    else {
      for(size_t j=0; j<N; ++j) 
	C.Rows(j,N) += alpha * A.col(j,j,N) ^ B.row(j);
    }
  }

  // In TMV_TriMatrixArith_A.cpp
  template <class T, class Ta, class Tx> extern void AddMultMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y);
  template <class T, class Ta, class Tx> extern void AddMultMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y);

  template <class T, class Ta, class Tb> void ColAddMultMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColAddMult Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    for(size_t j=0;j<C.rowsize();++j) 
      C.col(j) += alpha * A * B.col(j);
  }

  template <class T, class Ta, class Tb> inline void AddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct() == NonConj);

    if (A.isrm() && C.isrm()) RowAddMultMM(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPAddMultMM(alpha,A,B,C);
    else ColAddMultMM(alpha,A,B,C);
  }

  template <class T, class Ta, class Tb> inline void AddMultMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct() == NonConj);

    if (A.isrm() && C.isrm()) RowAddMultMM(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPAddMultMM(alpha,A,B,C);
    else ColAddMultMM(alpha,A,B,C);
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (B.isrm()) {
	Matrix<T,RowMajor> tempB = B.Cols(j,j2);
	MultEqMM(alpha,A,tempB.View());
	C.Cols(j,j2) *= beta;
	C.Cols(j,j2) += tempB;
      } else {
	Matrix<T,ColMajor> tempB = B.Cols(j,j2);
	MultEqMM(alpha,A,tempB.View());
	C.Cols(j,j2) *= beta;
	C.Cols(j,j2) += tempB;
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    if (B.isrm()) {
      Matrix<T,RowMajor> tempB = B;
      MultEqMM(alpha,A,tempB.View());
      C *= beta;
      C += tempB;
    } else {
      Matrix<T,ColMajor> tempB = B;
      MultEqMM(alpha,A,tempB.View());
      C *= beta;
      C += tempB;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha==T(0)) 
	C *= beta;
      else if (A.SameStorageAs(C)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } 
      else if (C.SameStorageAs(B)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      } 
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    if (B.isrm()) {
      Matrix<T,RowMajor> tempB = B;
      MultEqMM(alpha,A,tempB.View());
      C *= beta;
      C += tempB;
    } else {
      Matrix<T,ColMajor> tempB = B;
      MultEqMM(alpha,A,tempB.View());
      C *= beta;
      C += tempB;
    }
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (B.isrm()) {
	Matrix<T,RowMajor> tempB = B.Cols(j,j2);
	MultEqMM(alpha,A,tempB.View());
	C.Cols(j,j2) *= beta;
	C.Cols(j,j2) += tempB;
      } else {
	Matrix<T,ColMajor> tempB = B.Cols(j,j2);
	MultEqMM(alpha,A,tempB.View());
	C.Cols(j,j2) *= beta;
	C.Cols(j,j2) += tempB;
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha==T(0)) 
	C *= beta;
      else if (A.SameStorageAs(C)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } 
      else if (C.SameStorageAs(B)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void DoMultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct() == NonConj);

    const size_t N = A.size();

    if (beta == T(0)) {
      C = B;
      for(size_t j=0;j<N;++j) {
	C.col(j,0,j) = A.SubMatrix(0,j,j,N) * C.col(j,j,N);
	C.col(j,j,N) = A.SubTriMatrix(j,N) * C.col(j,j,N);
      }
      C *= alpha;
    } else {
      C *= beta;
      if (B.isunit()) {
	VIter<T> Cjj = C.diag().begin();
	if (A.isunit()) {
	  for(size_t j=0;j<N;++j,++Cjj) {
	    Vector<T> temp = A.col(j,0,j);
	    temp += A.SubMatrix(0,j,j+1,N) * B.col(j,j+1,N);
	    C.col(j,0,j) += alpha * temp;
	    *Cjj += alpha*(T(1) + A.row(j,j+1,N) * B.col(j,j+1,N));
	    C.col(j,j+1,N) += alpha * A.SubTriMatrix(j+1,N) * B.col(j,j+1,N);
	  }
	} else {
	  CVIter<Ta> Ajj = A.diag().begin();
	  for(size_t j=0;j<N;++j,++Ajj,++Cjj) {
	    Vector<T> temp = A.col(j,0,j);
	    temp += A.SubMatrix(0,j,j+1,N) * B.col(j,j+1,N);
	    C.col(j,0,j) += alpha * temp;
	    *Cjj += alpha*(*Ajj + A.row(j,j+1,N) * B.col(j,j+1,N));
	    C.col(j,j+1,N) += alpha * A.SubTriMatrix(j+1,N) * B.col(j,j+1,N);
	  }
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  C.col(j,0,j) += alpha * A.SubMatrix(0,j,j,N) * B.col(j,j,N);
	  C.col(j,j,N) += alpha * A.SubTriMatrix(j,N) * B.col(j,j,N);
	}
      }
    } 
  }

  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    const size_t N = A.size();

    if (N==0) return;
    else if (alpha == T(0)) 
      C *= beta;
    else if (A.SameStorageAs(C) || B.SameStorageAs(C)) {
      if (C.isrm()) {
	Matrix<T,RowMajor> tempC(C);
	DoMultMM(alpha,A,B,beta,tempC.View());
	C = tempC;
      } else {
	Matrix<T,ColMajor> tempC(C);
	DoMultMM(alpha,A,B,beta,tempC.View());
	C = tempC;
      }
    } 
    else if (C.isconj()) 
      DoMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),CONJ(beta),
	  C.Conjugate());
    else
      DoMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct() == NonConj);

    const size_t N = A.size();

    if (beta == T(0)) {
      C = B;
      for(size_t j=0;j<N;++j) {
	C.col(j,j+1,N) = A.SubMatrix(j+1,N,0,j+1) * C.col(j,0,j+1);
        C.col(j,0,j+1) = A.SubTriMatrix(0,j+1) * C.col(j,0,j+1);
      }
      C *= alpha;
    } else {
      C *= beta;
      if (B.isunit()) {
	VIter<T> Cjj = C.diag().begin();
	if (A.isunit()) {
	  for(size_t j=0;j<N;++j,++Cjj) {
	    Vector<T> temp = A.col(j,j+1,N);
	    temp += A.SubMatrix(j+1,N,0,j) * B.col(j,0,j);
	    C.col(j,j+1,N) += alpha * temp;
	    *Cjj += alpha*(T(1) + A.row(j,0,j) * B.col(j,0,j));
	    C.col(j,0,j) += alpha * A.SubTriMatrix(0,j) * B.col(j,0,j);
	  }
	} else {
	  CVIter<Ta> Ajj = A.diag().begin();
	  for(size_t j=0;j<N;++j,++Ajj,++Cjj) {
	    Vector<T> temp = A.col(j,j+1,N);
	    temp += A.SubMatrix(j+1,N,0,j) * B.col(j,0,j);
	    C.col(j,j+1,N) += alpha * temp;
	    *Cjj += alpha*(*Ajj + A.row(j,0,j) * B.col(j,0,j));
	    C.col(j,0,j) += alpha * A.SubTriMatrix(0,j) * B.col(j,0,j);
	  }
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  C.col(j,j+1,N) += alpha * A.SubMatrix(j+1,N,0,j+1)*B.col(j,0,j+1);
	  C.col(j,0,j+1) += alpha * A.SubTriMatrix(0,j+1) * B.col(j,0,j+1);
	}
      }
    } 
  }

  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    const size_t N = A.size();

    if (N==0) return;
    else if (alpha == T(0)) 
      C *= beta;
    else if (A.SameStorageAs(C) || B.SameStorageAs(C)) {
      if (C.isrm()) {
	Matrix<T,RowMajor> tempC(C);
	DoMultMM(alpha,A,B,beta,tempC.View());
	C = tempC;
      } else {
	Matrix<T,ColMajor> tempC(C);
	DoMultMM(alpha,A,B,beta,tempC.View());
	C = tempC;
      }
    }
    else if (C.isconj()) 
      DoMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),CONJ(beta),
	  C.Conjugate());
    else 
      DoMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta> void RowMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RowMultEqMM Upper/Upper alpha = "<<alpha<<endl;
    Matrix<T> B2 = alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> B0 = B;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	B.row(i,i+1,N) += A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	if (!B.isunit()) 
	  B.row(i,i,N) *= alpha;
	else { TMVAssert(alpha == T(1)); }
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin();
      TMVAssert(!B.isunit());
      VIter<T> Bii = B.diag().begin();
      for(size_t i=0; i<N; ++i,++Aii,++Bii) {
	// B.row(i,i+1,N) = alpha * (*Aii * B.row(i,i+1,N) +
	//     A.row(i,i+1,N) * B.SubTriMatrix(i+1,N));
	T aa = *Aii;
	if (alpha != T(1)) aa *= alpha;
	MultMV(alpha,B.SubTriMatrix(i+1,N).Transpose(),A.row(i,i+1,N),
	    aa,B.row(i,i+1,N));
	*Bii *= aa;
      }
    }
#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*(Norm(A)+Norm(B))) {
      cerr<<"RowMultEqMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"--> B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void OPMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"OPMultEqMM Upper/Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit()) {
	TMVAssert(alpha == T(1));
	for(size_t j=1; j<N; ++j) {
	  B.SubMatrix(0,j,j+1,N) += A.col(j,0,j) ^ B.row(j,j+1,N);
	  B.col(j,0,j) += A.col(j,0,j);
	}
      }
      else {
	for(size_t j=0; j<N; ++j) 
	  B.SubMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
	B *= alpha;
      }
    } else {
      CVIter<Ta> Ajj = A.diag().begin();
      TMVAssert(!B.isunit());
      for(size_t j=0; j<N; ++j,++Ajj) {
	B.SubMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
	B.row(j,j,N) *= *Ajj;
      }
      B *= alpha;
    } 
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"ColMultEqMM Upper/Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    if (B.isunit()) {
      // Then alpha = 1 and A.isunit
      for(size_t j=1;j<B.size();++j) {
	B.col(j,0,j) = A.SubTriMatrix(0,j) * B.col(j,0,j);
	B.col(j,0,j) += A.col(j,0,j);
      }
    }
    else {
      for(size_t j=0;j<B.size();++j) 
	B.col(j,0,j+1) = alpha * A.SubTriMatrix(0,j+1) * B.col(j,0,j+1);
    }
  }

  template <class T, class Ta> inline void MultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);
    TMVAssert(B.size() > 0);

    if (B.isrm()) 
      if (A.iscm()) OPMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
  }

  template <class T, class Ta> void RowMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RowMultEqMM Lower/Lower: alpha = "<<alpha<<endl;
#endif
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
      if (!B.isunit()) 
	B.diag() *= alpha;
    }
    else {
      CVIter<Ta> Aii = A.diag().begin()+N-1;
      TMVAssert(!B.isunit());
      VIter<T> Bii = B.diag().begin()+N-1;
      for(size_t i=N-1; i>0; --i,--Aii,--Bii) {
	// B.row(i,0,i) = alpha * (*Aii * B.row(i,0,i) +
	//     A.row(i,0,i) * B.SubTriMatrix(0,i));
	T aa = *Aii;
	if (alpha != T(1)) aa *= alpha;
	MultMV(alpha,B.SubTriMatrix(0,i).Transpose(),A.row(i,0,i),
	    aa,B.row(i,0,i));
	*Bii *= aa;
      }
    }
  }

  template <class T, class Ta> void OPMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"OPMultEqMM Lower/Lower: alpha = "<<alpha<<endl;
#endif
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
	B.row(j,0,j+1) *= *Ajj;
      }
    } 
    B *= alpha;
  }

  template <class T, class Ta> void ColMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"ColMultEqMM Lower/Lower: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();
    if (B.isunit()) {
      // Then alpha = 1 and A.isunit
      for(size_t j=0;j<B.size()-1;++j) {
	B.col(j,j+1,N) = A.SubTriMatrix(j+1,N) * B.col(j,j+1,N);
	B.col(j,j+1,N) += A.col(j,j+1,N);
      }
    }
    else {
      for(size_t j=0;j<B.size();++j) 
	B.col(j,j,N) = alpha * A.SubTriMatrix(j,N) * B.col(j,j,N);
    }
  }

  template <class T, class Ta> inline void MultEqMM(T alpha,
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);
    TMVAssert(B.size() > 0);

    if (B.isrm()) 
      if (A.iscm()) OPMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
  }


  template <class T, class Ta, class Tb> void RowAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowAddMultMM Upper/Upper\n";
    Matrix<T> C2 = Matrix<T>(C)+alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = C.size();

    if (A.isunit()) {
      if (B.isunit()) {
	VIter<T> Cii = C.diag().begin();
	for(size_t i=0; i<N; ++i) {
	  C.row(i,i+1,N) += alpha * A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	  C.row(i,i+1,N) += alpha * B.row(i,i+1,N);
	  *Cii += alpha;
	}
      } else {
	for(size_t i=0; i<N; ++i) {
	  C.row(i,i+1,N) += alpha * A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	  C.row(i,i,N) += alpha * B.row(i,i,N);
	}
      }
    } else {
      for(size_t i=0; i<N; ++i) 
	C.row(i,i,N) += alpha * A.row(i,i,N) * B.SubTriMatrix(i,N);
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"RowAddMultMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void OPAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPAddMult Upper/Upper: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = C.size();

    if (A.isunit()) {
      if (B.isunit()) {
	VIter<T> Cjj = C.diag().begin();
	for(size_t j=0; j<N; ++j,++Cjj) {
	  C.SubMatrix(0,j,j+1,N) += alpha * A.col(j,0,j) ^ B.row(j,j+1,N);
	  C.col(j,0,j) += alpha * A.col(j,0,j);
	  C.row(j,j+1,N) += alpha * B.row(j,j+1,N);
	  *Cjj += alpha;
	}
      }
      else {
	for(size_t j=0; j<N; ++j) {
	  C.SubMatrix(0,j,j,N) += alpha * A.col(j,0,j) ^ B.row(j,j,N);
	  C.row(j,j,N) += alpha * B.row(j,j,N);
	}
      }
    } else {
      if (B.isunit()) {
	for(size_t j=0; j<N; ++j) {
	  C.SubMatrix(0,j+1,j+1,N) += alpha * A.col(j,0,j+1)^B.row(j,j+1,N);
	  C.row(j,0,j+1) += alpha * A.col(j,0,j+1);
	}
      } else {
	for(size_t j=0; j<N; ++j) 
	  C.SubMatrix(0,j+1,j,N) += alpha * A.col(j,0,j+1) ^ B.row(j,j,N);
      }
    } 
  }

  template <class T, class Ta, class Tb> void ColAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColAddMult Upper/Upper: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = C.size();

    if (B.isunit()) {
      if (A.isunit()) {
	VIter<T> Cjj = C.diag().begin();
	for(size_t j=0;j<N;++j,++Cjj) {
	  C.col(j,0,j) += alpha * A.SubTriMatrix(0,j) * B.col(j,0,j);
	  C.col(j,0,j) += alpha * A.col(j,0,j);
	  *Cjj += alpha;
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  C.col(j,0,j) += alpha * A.SubTriMatrix(0,j) * B.col(j,0,j);
	  C.col(j,0,j+1) += alpha * A.col(j,0,j+1);
	}
      }
    } else {
      for(size_t j=0;j<N;++j)
	C.col(j,0,j+1) += alpha * A.SubTriMatrix(0,j+1) * B.col(j,0,j+1);
    }
  }

  template <class T, class Ta, class Tb> inline void AddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    if (A.isrm() && C.isrm()) RowAddMultMM(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPAddMultMM(alpha,A,B,C);
    else ColAddMultMM(alpha,A,B,C);
  }
  
  template <class T, class Ta, class Tb> void TempMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const UpperTriMatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    if (B.isrm()) {
      if (B.dt() == UnitDiag) {
	UpperTriMatrix<T,UnitDiag,RowMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      } else {
	UpperTriMatrix<T,NonUnitDiag,RowMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      }
    } else {
      if (B.dt() == UnitDiag) {
	UpperTriMatrix<T,UnitDiag,ColMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      } else {
	UpperTriMatrix<T,NonUnitDiag,ColMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      }
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const UpperTriMatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == C.size());
    TMVAssert(A.size() == B.size());
    TMVAssert(!C.isunit() || (A.isunit() && B.isunit() && alpha == T(1)) );

    if (C.size() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha==T(0)) 
	C *= beta;
      else if (C.SameStorageAs(A)) 
	if (beta == T(0) && !C.SameStorageAs(B)) {
	  C = A;
	  MultEqMM(alpha,B.Transpose(),C.Transpose());
	}
	else 
	  TempMultMM(alpha,A,B,beta,C);
      else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      }
      else if (C.SameStorageAs(B)) 
	TempMultMM(alpha,A,B,beta,C);
      else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


