
#include "TMV.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t MM_BLOCKSIZE = 64;
#endif

  // MJ: Look at Atlas code, and try to mimic structure to make this faster.
 
  //
  // MultMM
  //

  template <class T, class Ta, class Tb> void RowMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowMultMM\n";
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    for(size_t i=0;i<C.colsize();++i) 
      // C.row(i) = beta*C.row(i) + alpha * A.row(i) * B;
      // (Arith for this doesn't inline all the way yet.)
      MultMV(alpha,B.Transpose(),A.row(i),beta,C.row(i));
  }

  template <class T, class Ta, class Tb> void OPMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPMultMM\n";
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    C *= beta;
    for(size_t k=0;k<A.rowsize();++k) 
      C += alpha * (A.col(k) ^ B.row(k));
  }

  template <class T, class Ta, class Tb> void NonBlasMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"NonBlockMultMM\n";
#endif
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
    else if (B.iscm() && C.iscm()) 
      RowMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
    else if (A.iscm() && B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else {
      const size_t M = C.colsize();
      const size_t N = C.rowsize();
      const size_t K = A.rowsize();
      if (M < N && M < K) RowMultMM(alpha,A,B,beta,C);
      else if (N < M && N < K) 
	RowMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
      else if (K < M && K < N) OPMultMM(alpha,A,B,beta,C);
      else if (M < N) RowMultMM(alpha,A,B,beta,C);
      else if (A.isrm() || C.isrm()) RowMultMM(alpha,A,B,beta,C);
      else if (B.iscm() || C.iscm()) 
	RowMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
      else if (A.iscm() || B.isrm()) OPMultMM(alpha,A,B,beta,C);
      else RowMultMM(alpha,A,B,beta,C);
    }
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> void BlasMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  { NonBlasMultMM(alpha,A,B,beta,C); }
  template <> void BlasMultMM(
      double alpha, const GenMatrix<double>& A, const GenMatrix<double>& B,
      double beta, const MatrixView<double>& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM double\n";
#endif
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
  template <> void BlasMultMM(
      complex<double> alpha, const GenMatrix<complex<double> >& A,
      const GenMatrix<complex<double> >& B,
      complex<double> beta, const MatrixView<complex<double> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM c double\n";
#endif
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
  template <> void BlasMultMM(
      float alpha, const GenMatrix<float>& A, const GenMatrix<float>& B,
      float beta, const MatrixView<float>& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM float\n";
#endif
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
  template <> void BlasMultMM(
      complex<float> alpha, const GenMatrix<complex<float> >& A,
      const GenMatrix<complex<float> >& B,
      complex<float> beta, const MatrixView<complex<float> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM c float\n";
#endif
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
  template <> void BlasMultMM(
      const complex<double> alpha, const GenMatrix<complex<double> >& A,
      const GenMatrix<double>& B, const complex<double> beta, 
      const MatrixView<complex<double> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM lap c double\n";
#endif
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
      C *= alpha;
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
  template <> void BlasMultMM(
      const complex<double> alpha, const GenMatrix<double>& A,
      const GenMatrix<complex<double> >& B, const complex<double> beta, 
      const MatrixView<complex<double> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM lap2 c double\n";
#endif
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
      C *= alpha;
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
#ifndef NOFLOAT
  template <> void BlasMultMM(
      const complex<float> alpha, const GenMatrix<complex<float> >& A,
      const GenMatrix<float>& B, const complex<float> beta, 
      const MatrixView<complex<float> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM lap c float\n";
#endif
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
      C *= alpha;
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
  template <> void BlasMultMM(
      const complex<float> alpha, const GenMatrix<float>& A,
      const GenMatrix<complex<float> >& B, const complex<float> beta, 
      const MatrixView<complex<float> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM lap2 c float\n";
#endif
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
      C *= alpha;
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
#endif // NOFLOAT
#endif // LAP
#endif // BLAS

  template <class T, class Ta, class Tb> inline void DoMultMM(const T alpha,
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
    TMVAssert(C.ct() == NonConj);

#ifdef BLAS
    if ( ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) && 
	((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0)) && 
	((C.isrm() && C.stepi()>0) || (C.iscm() && C.stepj()>0)) )
      BlasMultMM(alpha,A,B,beta,C);
    else
#endif
      NonBlasMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C   via a temporary
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    } else {
      Matrix<T,ColMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    }
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * C + beta * C
  {
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (C.isrm()) {
	Matrix<T,RowMajor> tempB = B.Cols(j,j2);
	DoMultMM(alpha,A,tempB,beta,C.Cols(j,j2));
      } else {
	Matrix<T,ColMajor> tempB = B.Cols(j,j2);
	DoMultMM(alpha,A,tempB,beta,C.Cols(j,j2));
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    Matrix<T> B2 = B;
    B2 *= alpha;
    Matrix<T> A2 = A;
    Matrix<T> C3(C.colsize(),C.rowsize());
    if (C.colsize()>0 && C.rowsize()>0)
      if (A.rowsize()==0) C3.Zero();
      else DoMultMM(T(1),A2,B2,T(0),C3.View());
    Matrix<T> C2 = C;
    C2 *= beta;
    C2 += C3;
    Matrix<T> C0 = C;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (A.rowsize() == 0 || alpha == T(0)) 
	C *= beta;
      else if (C.SameStorageAs(A)) 
	if (C.SameStorageAs(B)) 
	  FullTempMultMM(alpha,A,B,beta,C);
	else if (C.stepi() == A.stepi() && C.stepj() == A.stepj())
	  BlockTempMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else if (C.SameStorageAs(B))
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else
	DoMultMM(alpha,A,B,beta,C);
    }
      
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*max(RealType(T)(1),Norm(C)+Norm(B)+Norm(A))) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


