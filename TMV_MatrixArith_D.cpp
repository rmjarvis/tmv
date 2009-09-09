
#include "TMV.h"

namespace tmv {

  //
  // MultMM
  //

  template <class T, class Ta, class Tb> inline void RowMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    for(size_t i=0;i<C.colsize();++i) 
      MultMV(alpha,B.Transpose(),A.row(i),beta,C.row(i));
  }

  template <class T, class Ta, class Tb> inline void ColMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    for(size_t j=0;j<C.rowsize();++j) MultMV(alpha,A,B.col(j),beta,C.col(j));
  }

  template <class T, class Ta, class Tb> inline void OPMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    if (beta == T(0)) C.Zero();
    else if (beta != T(1)) C *= beta;
    for(size_t k=0;k<A.rowsize();++k) Rank1Update(alpha,A.col(k),B.row(k),C);
  }

  template <class T, class Ta, class Tb> inline void NonBlasMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    // We want to make the inner loops as efficient as possible.
    // As mentioned above, we want the inner loops to be unit stride
    // or with long vectors.

    if (A.isrm() && C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (B.iscm() && C.iscm()) ColMultMM(alpha,A,B,beta,C);
    else if (A.iscm() && B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else {
      const size_t M = C.colsize();
      const size_t N = C.rowsize();
      const size_t K = A.rowsize();
      if (M < N && M < K) RowMultMM(alpha,A,B,beta,C);
      else if (N < M && N < K) ColMultMM(alpha,A,B,beta,C);
      else if (K < M && K < N) OPMultMM(alpha,A,B,beta,C);
      else if (M < N) RowMultMM(alpha,A,B,beta,C);
      else if (A.isrm() || C.isrm()) RowMultMM(alpha,A,B,beta,C);
      else if (B.iscm() || C.iscm()) ColMultMM(alpha,A,B,beta,C);
      else if (A.iscm() || B.isrm()) OPMultMM(alpha,A,B,beta,C);
      else ColMultMM(alpha,A,B,beta,C);
    }
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> inline void BlasMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  { NonBlasMultMM(alpha,A,B,beta,C); }
  template <> inline void BlasMultMM(
      double alpha, const GenMatrix<double>& A, const GenMatrix<double>& B,
      double beta, const MatrixView<double>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());

    cblas_dgemm(C.isrm() ? CblasRowMajor : CblasColMajor,
	(A.isrm()==C.isrm())?CblasNoTrans:CblasTrans,
	(B.isrm()==C.isrm())?CblasNoTrans:CblasTrans,
	C.colsize(),C.rowsize(),A.rowsize(),
	alpha,A.cptr(),A.isrm()?A.stepi():A.stepj(),
	B.cptr(),B.isrm()?B.stepi():B.stepj(),
	beta,C.ptr(),C.isrm()?C.stepi():C.stepj());
  }
  template <> inline void BlasMultMM(
      complex<double> alpha, const GenMatrix<complex<double> >& A,
      const GenMatrix<complex<double> >& B,
      complex<double> beta, const MatrixView<complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());

    if ((A.stor() == C.stor() && A.isconj()) ||
        (B.stor() == C.stor() && B.isconj()) )
      NonBlasMultMM(alpha,A,B,beta,C);
    else {
      CBLAS_ORDER order = C.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_TRANSPOSE Atran = 
	A.isconj() ?  CblasConjTrans :
	(A.stor()==C.stor() ? CblasNoTrans : CblasTrans);
      CBLAS_TRANSPOSE Btran = 
	B.isconj() ?  CblasConjTrans :
	(B.stor()==C.stor() ? CblasNoTrans : CblasTrans);
      cblas_zgemm(order,Atran,Btran,
	  C.colsize(),C.rowsize(),A.rowsize(),
	  &alpha,A.cptr(),A.isrm()?A.stepi():A.stepj(),
	  B.cptr(),B.isrm()?B.stepi():B.stepj(),
	  &beta,C.ptr(),C.isrm()?C.stepi():C.stepj());
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMM(
      float alpha, const GenMatrix<float>& A, const GenMatrix<float>& B,
      float beta, const MatrixView<float>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());

    cblas_sgemm(C.isrm() ? CblasRowMajor : CblasColMajor,
	(A.isrm()==C.isrm())?CblasNoTrans:CblasTrans,
	(B.isrm()==C.isrm())?CblasNoTrans:CblasTrans,
	C.colsize(),C.rowsize(),A.rowsize(),
	alpha,A.cptr(),A.isrm()?A.stepi():A.stepj(),
	B.cptr(),B.isrm()?B.stepi():B.stepj(),
	beta,C.ptr(),C.isrm()?C.stepi():C.stepj());
  }
  template <> inline void BlasMultMM(
      complex<float> alpha, const GenMatrix<complex<float> >& A,
      const GenMatrix<complex<float> >& B,
      complex<float> beta, const MatrixView<complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());

    if ((A.stor() == C.stor() && A.isconj()) ||
        (B.stor() == C.stor() && B.isconj()) )
      NonBlasMultMM(alpha,A,B,beta,C);
    else {
      CBLAS_ORDER order = C.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_TRANSPOSE Atran = 
	A.isconj() ? CblasConjTrans :
	(A.stor()==C.stor() ? CblasNoTrans : CblasTrans);
      CBLAS_TRANSPOSE Btran = 
	B.isconj() ? CblasConjTrans :
	(B.stor()==C.stor() ? CblasNoTrans : CblasTrans);
      cblas_cgemm(order,Atran,Btran,
	  C.colsize(),C.rowsize(),A.rowsize(),
	  &alpha,A.cptr(),A.isrm()?A.stepi():A.stepj(),
	  B.cptr(),B.isrm()?B.stepi():B.stepj(),
	  &beta,C.ptr(),C.isrm()?C.stepi():C.stepj());
    }
  }
#endif //NOFLOAT
#ifdef LAP
  template <> inline void BlasMultMM(
      const complex<double> alpha, const GenMatrix<complex<double> >& A,
      const GenMatrix<double>& B, const complex<double> beta, 
      const MatrixView<complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);

    if (A.iscm() && B.iscm() && C.iscm() && 
	beta == double(0) && B.IsSquare() && !A.isconj()) {
      int m = C.colsize();
      int n = C.rowsize();
      int lda = A.stepj();
      int ldb = B.stepj();
      int ldc = C.stepj();
      int lwork = 2*m*n;
      double* rwork = LAP_DWork(lwork);
      zlacrm(&m,&n,LAP_Complex(A.cptr()),&lda,
	  const_cast<double*>(B.cptr()),&ldb,
	  LAP_Complex(C.ptr()),&ldc,rwork);
      if (alpha != double(1)) MultXM(alpha,C);
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
  template <> inline void BlasMultMM(
      const complex<double> alpha, const GenMatrix<double>& A,
      const GenMatrix<complex<double> >& B, const complex<double> beta, 
      const MatrixView<complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);

    if (A.isrm() && B.isrm() && C.isrm() && 
	beta == double(0) && A.IsSquare() && !B.isconj()) {
      int m = C.rowsize();
      int n = C.colsize();
      int lda = A.stepi();
      int ldb = B.stepi();
      int ldc = C.stepi();
      int lwork = 2*m*n;
      double* rwork = LAP_DWork(lwork);
      zlacrm(&m,&n,LAP_Complex(B.cptr()),&ldb,
	  const_cast<double*>(A.cptr()),&lda,
	  LAP_Complex(C.ptr()),&ldc,rwork);
      if (alpha != double(1)) MultXM(alpha,C);
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMM(
      const complex<float> alpha, const GenMatrix<complex<float> >& A,
      const GenMatrix<float>& B, const complex<float> beta, 
      const MatrixView<complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);

    if (A.iscm() && B.iscm() && C.iscm() && 
	beta == float(0) && B.IsSquare() && !A.isconj()) {
      int m = C.colsize();
      int n = C.rowsize();
      int lda = A.stepj();
      int ldb = B.stepj();
      int ldc = C.stepj();
      int lwork = 2*m*n;
      float* rwork = LAP_SWork(lwork);
      clacrm(&m,&n,LAP_Complex(A.cptr()),&lda,
	  const_cast<float*>(B.cptr()),&ldb,
	  LAP_Complex(C.ptr()),&ldc,rwork);
      if (alpha != float(1)) MultXM(alpha,C);
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
  template <> inline void BlasMultMM(
      const complex<float> alpha, const GenMatrix<float>& A,
      const GenMatrix<complex<float> >& B, const complex<float> beta, 
      const MatrixView<complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);

    if (A.isrm() && B.isrm() && C.isrm() && 
	beta == float(0) && A.IsSquare() && !B.isconj()) {
      int m = C.rowsize();
      int n = C.colsize();
      int lda = A.stepi();
      int ldb = B.stepi();
      int ldc = C.stepi();
      int lwork = 2*m*n;
      float* rwork = LAP_SWork(lwork);
      clacrm(&m,&n,LAP_Complex(B.cptr()),&ldb,
	  const_cast<float*>(A.cptr()),&lda,
	  LAP_Complex(C.ptr()),&ldc,rwork);
      if (alpha != float(1)) MultXM(alpha,C);
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
#endif // NOFLOAT
#endif // LAP
#endif // BLAS

  template <class T, class Ta, class Tb> void DoMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

    if (C.isconj()) DoMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	CONJ(beta),C.Conjugate());
    else
#ifdef BLAS
      if ( (A.isrm() || A.iscm()) && 
	  (B.isrm() || B.iscm()) &&
	  (C.isrm() || C.iscm()) )
	BlasMultMM(alpha,A,B,beta,C);
      else
#endif
	NonBlasMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0) C.Zero(); 
      else if (alpha == T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	if (C.SameAs(A) && !C.SameStorageAs(B)) {
	  for(size_t i=0;i<C.colsize();++i) {
	    MultMV(alpha,Transpose(B),C.row(i),beta,C.row(i));
	  }
	} else {
	  if (C.isrm()) {
	    Matrix<T,RowMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else {
	    Matrix<T,ColMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  }
	}
      } else if (C.SameStorageAs(B)) {
	if (C.SameAs(B)) {
	  for(size_t j=0;j<C.rowsize();++j) {
	    MultMV(alpha,A,C.col(j),beta,C.col(j));
	  }
	} else {
	  if (C.isrm()) {
	    Matrix<T,RowMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else {
	    Matrix<T,ColMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  }
	}
      } else DoMultMM(alpha, A, B, beta, C);
    }
  }

#define InstFile "TMV_MatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


