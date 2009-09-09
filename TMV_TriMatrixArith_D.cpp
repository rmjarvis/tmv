
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

  // MJ: Write blocked versions?

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
	MultXV(alpha,B.row(i));
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin();
      for(size_t i=0; i<N; ++i,++Aii) {
	MultXV(*Aii,B.row(i));
	B.row(i) += A.row(i,i+1,N) * B.Rows(i+1,N);
	MultXV(alpha,B.row(i));
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
	Rank1Update(T(1),A.col(j,0,j),B.row(j),B.Rows(0,j));
    else {
      CVIter<Ta> Ajj = A.diag().begin();
      for(size_t j=0; j<N; ++j,++Ajj) {
	Rank1Update(T(1),A.col(j,0,j),B.row(j),B.Rows(0,j));
	MultXV(*Ajj,B.row(j));
      }
    }
    MultXM(alpha,B);
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

    for(size_t j=0;j<B.rowsize();++j) MultEqMV(A,B.col(j));
    MultXM(alpha,B);
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
	MultXV(alpha,B.row(i));
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin()+N-1;
      for(int i=N-1; i>=0; --i,--Aii) {
	MultXV(*Aii,B.row(i));
	B.row(i) += A.row(i,0,i) * B.Rows(0,i);
	MultXV(alpha,B.row(i));
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
	Rank1Update(T(1),A.col(j,j+1,N),B.row(j),B.Rows(j+1,N));
    else {
      CVIter<Ta> Ajj = A.diag().begin()+N-1;
      for(int j=N-1; j>=0; --j,--Ajj) {
	Rank1Update(T(1),A.col(j,j+1,N),B.row(j),B.Rows(j+1,N));
	MultXV(*Ajj,B.row(j));
      }
    }
    MultXM(alpha,B);
  }

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

    for(size_t j=0;j<B.rowsize();++j) MultEqMV(A,B.col(j));
    MultXM(alpha,B);
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

  template <class T, class Ta, class Tb> void RowAddMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowAddMultEq Upper alpha = "<<alpha<<endl;
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
	AddVV(alpha,B.row(i),C.row(i));
      }
    }
    else {
      for(size_t i=0; i<N; ++i) 
	C.row(i) += alpha * A.row(i,i,N) * B.Rows(i,N);
    }
  }

  template <class T, class Ta, class Tb> void OPAddMultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPAddMultEq Upper alpha = "<<alpha<<endl;
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
	Rank1Update(alpha,A.col(j,0,j),B.row(j),C.Rows(0,j));
	AddVV(alpha,B.row(j),C.row(j));
      }
    else {
      for(size_t j=0; j<N; ++j) 
	Rank1Update(alpha,A.col(j,0,j+1),B.row(j),C.Rows(0,j+1));
    }
  }

  template <class T, class Ta, class Tb> void ColAddMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColAddMultEq Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    for(size_t j=0;j<C.rowsize();++j) 
      AddMultEqMV(alpha,A,B.col(j),C.col(j));
  }

  template <class T, class Ta, class Tb> void RowAddMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowAddMultEq Lower alpha = "<<alpha<<endl;
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
	C.row(i,0,i) += alpha * A.row(i,0,i) * B.Rows(0,i);
	AddVV(alpha,B.row(i,0,i),C.row(i,0,i));
      }
    }
    else {
      for(size_t i=0; i<N; ++i) 
	C.row(i,0,i+1) += alpha * A.row(i,0,i+1) * B.Rows(0,i+1);
    }
  }

  template <class T, class Ta, class Tb> void OPAddMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPAddMultEq Lower alpha = "<<alpha<<endl;
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
	Rank1Update(alpha,A.col(j,j+1,N),B.row(j),C.Rows(j+1,N));
	AddVV(alpha,B.row(j),C.row(j));
      }
    else {
      for(size_t j=0; j<N; ++j) 
	Rank1Update(alpha,A.col(j,j,N),B.row(j),C.Rows(j,N));
    }
  }

  template <class T, class Ta, class Tb> void ColAddMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColAddMultEq Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    for(size_t j=0;j<C.rowsize();++j) 
      AddMultEqMV(alpha,A,B.col(j),C.col(j));
  }

  // In TMV_TriMatrixArith_A.cpp
  template <class T, class Ta, class Tx> extern void AddMultEqMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y);
  template <class T, class Ta, class Tx> extern void AddMultEqMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y);

  template <class T, class Ta, class Tb> inline void AddMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
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

    if (A.isrm() && C.isrm()) RowAddMultEqMM(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPAddMultEqMM(alpha,A,B,C);
    else ColAddMultEqMM(alpha,A,B,C);
  }

  template <class T, class Ta, class Tb> inline void AddMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
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

    if (A.isrm() && C.isrm()) RowAddMultEqMM(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPAddMultEqMM(alpha,A,B,C);
    else ColAddMultEqMM(alpha,A,B,C);
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
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
      if (alpha==T(0)) {
	MultXM(beta,C);
      } else if (C.SameStorageAs(A)) {
	if (B.isrm()) {
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	}
      } else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } else if (C.SameStorageAs(B)) { 
	if (B.isrm()) {
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	}
      } else {
	MultXM(beta,C);
	AddMultEqMM(alpha,A,B,C);
      } 
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.0001) {
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

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
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
      if (alpha==T(0)) {
	MultXM(beta,C);
      } else if (C.SameStorageAs(A)) {
	if (B.isrm()) {
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	}
      } else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } else if (C.SameStorageAs(B)) { 
	if (B.isrm()) {
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	}
      } else {
	MultXM(beta,C);
	AddMultEqMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.0001) {
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

  template <class T, class Ta, class Tb> inline void DoMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = A.size();

    if (beta == T(0)) {
      C = B;
      for(size_t j=0;j<N;++j) {
	C.col(j,0,j) = A.SubMatrix(0,j,j,N)*C.col(j,j,N);
	MultEqMV(A.SubTriMatrix(j,N),C.col(j,j,N));
      }
      MultXM(alpha,C);
    } else {
      MultXM(beta,C);
      if (B.isunit()) {
	for(size_t j=0;j<N;++j) {
	  Vector<T> temp = A.col(j,0,j);
	  temp += A.SubMatrix(0,j,j+1,N)*B.col(j,j+1,N);
	  AddVV(alpha,temp,C.col(j,0,j));
	  C(j,j) += alpha*(A(j,j) + A.row(j,j+1,N)*B.col(j,j+1,N));
	  AddMultEqMV(alpha,A.SubTriMatrix(j+1,N),B.col(j,j+1,N),
	      C.col(j,j+1,N));
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  C.col(j,0,j) += alpha*A.SubMatrix(0,j,j,N)*B.col(j,j,N);
	  AddMultEqMV(alpha,A.SubTriMatrix(j,N),B.col(j,j,N),C.col(j,j,N));
	}
      }
    } 
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    if (A.size() > 0) {
      if (alpha==T(0)) {
	MultXM(beta,C);
      } else if (C.SameStorageAs(A)) {
	if (beta == T(0) && !C.SameStorageAs(B)) {
	  DoMultMM(alpha,B.QuickTranspose(),A.QuickTranspose(),beta,
	      C.QuickTranspose());
	} else if (C.isrm()) {
	  Matrix<T,RowMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	} else {
	  Matrix<T,ColMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	}
      } else if (C.SameStorageAs(B)) { 
	if (beta == T(0)) {
	  DoMultMM(alpha,A,B,beta,C);
	} else if (C.isrm()) {
	  Matrix<T,RowMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	} else {
	  Matrix<T,ColMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	}
      } else {
	if (C.isrm())
	  DoMultMM(alpha,B.QuickTranspose(),A.QuickTranspose(),beta,
	      C.QuickTranspose());
	else 
	  DoMultMM(alpha,A,B,beta,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.0001) {
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

  template <class T, class Ta, class Tb> inline void DoMultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = A.size();

    if (beta == T(0)) {
      C = B;
      for(size_t j=0;j<N;++j) {
	C.col(j,j+1,N) = A.SubMatrix(j+1,N,0,j+1)*C.col(j,0,j+1);
	MultEqMV(A.SubTriMatrix(0,j+1),C.col(j,0,j+1));
      }
      MultXM(alpha,C);
    } else {
      MultXM(beta,C);
      if (B.isunit()) {
	for(size_t j=0;j<N;++j) {
	  Vector<T> temp = A.col(j,j+1,N);
	  temp += A.SubMatrix(j+1,N,0,j)*B.col(j,0,j);
	  AddVV(alpha,temp,C.col(j,j+1,N));
	  C(j,j) += alpha*(A(j,j) + A.row(j,0,j)*B.col(j,0,j));
	  AddMultEqMV(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  C.col(j,j+1,N) += alpha*A.SubMatrix(j+1,N,0,j+1)*B.col(j,0,j+1);
	  AddMultEqMV(alpha,A.SubTriMatrix(0,j+1),B.col(j,0,j+1),
	      C.col(j,0,j+1));
	}
      }
    } 
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
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

    const size_t N = A.size();

    if (N > 0) {
      if (alpha==T(0)) {
	MultXM(beta,C);
      } else if (C.SameStorageAs(A)) {
	if (beta == T(0) && !C.SameStorageAs(B)) {
	  DoMultMM(alpha,B.QuickTranspose(),A.QuickTranspose(),beta,
	      C.QuickTranspose());
	} else if (C.isrm()) {
	  Matrix<T,RowMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	} else {
	  Matrix<T,ColMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	}
      } else if (C.SameStorageAs(B)) { 
	if (beta == T(0)) {
	  DoMultMM(alpha,A,B,beta,C);
	} else if (C.isrm()) {
	  Matrix<T,RowMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	} else {
	  Matrix<T,ColMajor> tempC = C;
	  MultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	}
      } else {
	if (C.isrm())
	  DoMultMM(alpha,B.QuickTranspose(),A.QuickTranspose(),beta,
	      C.QuickTranspose());
	else 
	  DoMultMM(alpha,A,B,beta,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.0001) {
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
	if (!B.isunit()) MultXV(alpha,B.row(i,i,N));
	else { TMVAssert(alpha == T(1)); }
      }
    }
    else {
      CVIter<Ta> Aii = A.diag().begin();
      TMVAssert(!B.isunit());
      for(size_t i=0; i<N; ++i,++Aii) {
	MultXV(*Aii,B.row(i,i,N));
	B.row(i,i+1,N) += A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	MultXV(alpha,B.row(i,i,N));
      }
    }
#ifdef XDEBUG
    if (Norm(B2-B) > 1.e-4) {
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
	  Rank1Update(T(1),A.col(j,0,j),B.row(j,j+1,N),B.SubMatrix(0,j,j+1,N));
	  AddVV(T(1),A.col(j,0,j),B.col(j,0,j));
	}
      }
      else {
	for(size_t j=0; j<N; ++j) 
	  Rank1Update(T(1),A.col(j,0,j),B.row(j,j,N),B.SubMatrix(0,j,j,N));
	MultXM(alpha,B);
      }
    } else {
      CVIter<Ta> Ajj = A.diag().begin();
      TMVAssert(!B.isunit());
      for(size_t j=0; j<N; ++j,++Ajj) {
	Rank1Update(T(1),A.col(j,0,j),B.row(j,j,N),B.SubMatrix(0,j,j,N));
	MultXV(*Ajj,B.row(j,j,N));
      }
      MultXM(alpha,B);
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
	MultEqMV(A.SubTriMatrix(0,j),B.col(j,0,j));
	AddVV(T(1),A.col(j,0,j),B.col(j,0,j));
      }
    }
    else {
      for(size_t j=0;j<B.size();++j)
	MultEqMV(A.SubTriMatrix(0,j+1),B.col(j,0,j+1));
      MultXM(alpha,B);
    }
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
    else if (B.isrm()) 
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
	MultXV(alpha,B.row(i,0,i));
      }
      if (!B.isunit()) MultXV(alpha,B.diag());
    }
    else {
      CVIter<Ta> Aii = A.diag().begin()+N-1;
      TMVAssert(!B.isunit());
      for(size_t i=N-1; i>0; --i,--Aii) {
	MultXV(*Aii,B.row(i,0,i+1));
	B.row(i,0,i) += A.row(i,0,i) * B.SubTriMatrix(0,i);
	MultXV(alpha,B.row(i,0,i+1));
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
	  Rank1Update(T(1),A.col(j,j+1,N),B.row(j,0,j),B.SubMatrix(j+1,N,0,j));
	  AddVV(T(1),A.col(j,j+1,N),B.col(j,j+1,N));
	}
      }
      else {
	for(int j=N-1; j>=0; --j) 
	  Rank1Update(T(1),A.col(j,j+1,N),B.row(j,0,j+1),
	      B.SubMatrix(j+1,N,0,j+1));
      }
    } else {
      CVIter<Ta> Ajj = A.diag().begin()+N-1;
      TMVAssert(!B.isunit());
      for(int j=N-1; j>=0; --j,--Ajj) {
	Rank1Update(T(1),A.col(j,j+1,N),B.row(j,0,j+1),
	    B.SubMatrix(j+1,N,0,j+1));
	MultXV(*Ajj,B.row(j,0,j+1));
      }
    } 
    MultXM(alpha,B);
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

    if (B.isunit()) {
      // Then alpha = 1 and A.isunit
      for(size_t j=0;j<B.size()-1;++j) {
	MultEqMV(A.SubTriMatrix(j+1,B.size()),B.col(j,j+1,B.size()));
	AddVV(T(1),A.col(j,j+1,B.size()),B.col(j,j+1,B.size()));
      }
    }
    else {
      for(size_t j=0;j<B.size();++j)
	MultEqMV(A.SubTriMatrix(j,B.size()),B.col(j,j,B.size()));
      MultXM(alpha,B);
    }
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
    else if (B.isrm()) 
      if (A.iscm()) OPMultEqMM(alpha,A,B);
      else RowMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
  }

  template <class T, class Ta, class Tb> void RowAddMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowAddMultEqMM Upper/Upper\n";
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
	for(size_t i=0; i<N; ++i) {
	  AddMultEqMV(alpha,B.SubTriMatrix(i+1,N).QuickTranspose(),
	      A.row(i,i+1,N),C.row(i,i+1,N));
	  AddVV(alpha,B.row(i,i+1,N),C.row(i,i+1,N));
	  C(i,i) += alpha;
	}
      } else {
	for(size_t i=0; i<N; ++i) {
	  AddMultEqMV(alpha,B.SubTriMatrix(i+1,N).QuickTranspose(),
	      A.row(i,i+1,N),C.row(i,i+1,N));
	  AddVV(alpha,B.row(i,i,N),C.row(i,i,N));
	}
      }
    } else {
      for(size_t i=0; i<N; ++i) {
	AddMultEqMV(alpha,B.SubTriMatrix(i,N).QuickTranspose(),
	    A.row(i,i,N),C.row(i,i,N));
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 1.e-4) {
      cerr<<"RowAddMultEqMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void OPAddMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPAddMultEq Upper/Upper: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = C.size();

    if (A.isunit()) {
      if (B.isunit()) {
	for(size_t j=0; j<N; ++j) {
	  Rank1Update(alpha,A.col(j,0,j),B.row(j,j+1,N),C.SubMatrix(0,j,j+1,N));
	  AddVV(alpha,A.col(j,0,j),C.col(j,0,j));
	  AddVV(alpha,B.row(j,j+1,N),C.row(j,j+1,N));
	  C(j,j) += alpha;
	}
      }
      else {
	for(size_t j=0; j<N; ++j) {
	  Rank1Update(alpha,A.col(j,0,j),B.row(j,j,N),C.SubMatrix(0,j,j,N));
	  AddVV(alpha,B.row(j,j,N),C.row(j,j,N));
	}
      }
    } else {
      if (B.isunit()) {
	for(size_t j=0; j<N; ++j) {
	  Rank1Update(alpha,A.col(j,0,j+1),B.row(j,j+1,N),
	      C.SubMatrix(0,j+1,j+1,N));
	  AddVV(alpha,A.col(j,0,j+1),C.row(j,0,j+1));
	}
      } else {
	for(size_t j=0; j<N; ++j) {
	  Rank1Update(alpha,A.col(j,0,j+1),B.row(j,j,N),
	      C.SubMatrix(0,j+1,j,N));
	}
      }
    } 
  }

  template <class T, class Ta, class Tb> void ColAddMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColAddMultEq Upper/Upper: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = C.size();

    if (B.isunit()) {
      if (A.isunit()) {
	for(size_t j=0;j<N;++j) {
	  AddMultEqMV(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));
	  AddVV(alpha,A.col(j,0,j),C.col(j,0,j));
	  C(j,j) += alpha;
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  AddMultEqMV(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));
	  AddVV(alpha,A.col(j,0,j+1),C.col(j,0,j+1));
	}
      }
    } else {
      for(size_t j=0;j<N;++j)
	AddMultEqMV(alpha,A.SubTriMatrix(0,j+1),B.col(j,0,j+1),
	    C.col(j,0,j+1));
    }
  }

  template <class T, class Ta, class Tb> inline void AddMultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));

    if (B.isconj()) AddMultEqMM(CONJ(alpha),A.QuickConjugate(),
	B.QuickConjugate(),C.QuickConjugate());
    else if (A.isrm() && C.isrm()) RowAddMultEqMM(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPAddMultEqMM(alpha,A,B,C);
    else ColAddMultEqMM(alpha,A,B,C);
  }

  template <class T, class Ta, class Tb> inline void MultMM(const T alpha,
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
    if (!(!C.isunit() || (A.isunit() && B.isunit() && alpha == T(1)) ))
    { 
      cerr<<"alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<endl;
      cerr<<"B = "<<Type(B)<<endl;
      cerr<<"C = "<<Type(C)<<endl;
    }
    TMVAssert(!C.isunit() || (A.isunit() && B.isunit() && alpha == T(1)) );

    if (C.size() > 0) {
      if (alpha==T(0)) {
	MultXM(beta,C);
      } else if (C.SameStorageAs(A)) {
	if (beta == T(0) && !C.SameStorageAs(B)) {
	  C = A;
	  MultEqMM(alpha,B.QuickTranspose(),C.QuickTranspose());
	} else if (B.isrm()) {
	  if (B.dt() == UnitDiag) {
	    UpperTriMatrix<T,UnitDiag,RowMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  } else {
	    UpperTriMatrix<T,NonUnitDiag,RowMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  }
	} else {
	  if (B.dt() == UnitDiag) {
	    UpperTriMatrix<T,UnitDiag,ColMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  } else {
	    UpperTriMatrix<T,NonUnitDiag,ColMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  }
	}
      } else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } else if (C.SameStorageAs(B)) { 
	if (B.isrm()) {
	  if (B.dt() == UnitDiag) {
	    UpperTriMatrix<T,UnitDiag,RowMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  } else {
	    UpperTriMatrix<T,NonUnitDiag,RowMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  }
	} else {
	  if (B.dt() == UnitDiag) {
	    UpperTriMatrix<T,UnitDiag,ColMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  } else {
	    UpperTriMatrix<T,NonUnitDiag,ColMajor> tempB(B);
	    MultEqMM(alpha,A,tempB.QuickView());
	    MultXM(beta,C);
	    AddMM(T(1),tempB,T(1),C);
	  }
	}
      } else {
	MultXM(beta,C);
	AddMultEqMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.0001) {
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


