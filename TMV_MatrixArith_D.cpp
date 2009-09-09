
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

  template <class T, class Ta, class Tb> inline void RowMultMM(const T alpha,
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
      MultMV(alpha,B.Transpose(),A.row(i),beta,C.row(i));
  }

  template <class T, class Ta, class Tb> inline void OPMultMM(const T alpha,
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

  template <class T, class Ta, class Tb> inline void NonBlasMultMM(
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
    TMVAssert(!C.isrm());
    // We want to make the inner loops as efficient as possible.
    // We want the inner loops to be unit stride or with long vectors.

    if (B.iscm() && C.iscm()) 
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
      else if (A.isrm()) RowMultMM(alpha,A,B,beta,C);
      else if (B.iscm() || C.iscm()) 
	RowMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
      else if (A.iscm() || B.isrm()) OPMultMM(alpha,A,B,beta,C);
      else RowMultMM(alpha,A,B,beta,C);
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
    TMVAssert(C.iscm());

    int m = C.colsize();
    int n = C.rowsize();
    int k = A.rowsize();
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    int ldc = C.stepj();
    BLASNAME(dgemm) (BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
	B.iscm()?BLASCH_NT:BLASCH_T,
	BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	BLASV(beta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> inline void BlasMultMM(
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
    TMVAssert(C.iscm());

    if (A.iscm() && A.isconj()) {
      Matrix<complex<double> > AA = alpha*A;
      return BlasMultMM(complex<double>(1),AA,B,beta,C);
    } else if (B.iscm() && B.isconj()) {
      Matrix<complex<double> > BB = alpha*B;
      return BlasMultMM(complex<double>(1),A,BB,beta,C);
    } else {
      int m = C.colsize();
      int n = C.rowsize();
      int k = A.rowsize();
      int lda = A.isrm()?A.stepi():A.stepj();
      int ldb = B.isrm()?B.stepi():B.stepj();
      int ldc = C.stepj();
      BLASNAME(zgemm) (BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
	  BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	  BLASP(&beta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMM(
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
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    int m = C.colsize();
    int n = C.rowsize();
    int k = A.rowsize();
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    int ldc = C.stepj();
    BLASNAME(sgemm) (BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
	B.iscm()?BLASCH_NT:BLASCH_T,
	BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	BLASV(beta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> inline void BlasMultMM(
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
    TMVAssert(C.iscm());

    if (A.iscm() && A.isconj()) {
      Matrix<complex<float> > AA = alpha*A;
      return BlasMultMM(complex<float>(1),AA,B,beta,C);
    } else if (B.iscm() && B.isconj()) {
      Matrix<complex<float> > BB = alpha*B;
      return BlasMultMM(complex<float>(1),A,BB,beta,C);
    } else {
      int m = C.colsize();
      int n = C.rowsize();
      int k = A.rowsize();
      int lda = A.isrm()?A.stepi():A.stepj();
      int ldb = B.isrm()?B.stepi():B.stepj();
      int ldc = C.stepj();
      BLASNAME(cgemm) (BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
	  BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	  BLASP(&beta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
  }
#endif //NOFLOAT
#ifdef ELAP
  template <> inline void BlasMultMM(
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
    TMVAssert(C.iscm());

    if (A.iscm() && B.iscm() && 
	beta == double(0) && B.IsSquare() && !A.isconj()) {
      int m = C.colsize();
      int n = C.rowsize();
      int lda = A.stepj();
      int ldb = B.stepj();
      int ldc = C.stepj();
#ifndef LAPNOWORK
      int lwork = 2*m*n;
      double* rwork = LAP_DWork(lwork);
#endif
      LAPNAMEX(zlacrm) (LAPCM LAPV(m),LAPV(n),LAPP(A.cptr()),LAPV(lda),
	  LAPP(B.cptr()),LAPV(ldb),LAPP(C.ptr()),LAPV(ldc) LAPWK(rwork));
      C *= alpha;
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMM(
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
    TMVAssert(C.iscm());

    if (A.iscm() && B.iscm() &&
	beta == float(0) && B.IsSquare() && !A.isconj()) {
      int m = C.colsize();
      int n = C.rowsize();
      int lda = A.stepj();
      int ldb = B.stepj();
      int ldc = C.stepj();
#ifndef LAPNOWORK
      int lwork = 2*m*n;
      float* rwork = LAP_SWork(lwork);
#endif
      LAPNAMEX(clacrm) (LAPCM LAPV(m),LAPV(n),LAPP(A.cptr()),LAPV(lda),
	  LAPP(B.cptr()),LAPV(ldb),LAPP(C.ptr()),LAPV(ldc) LAPWK(rwork));
      C *= alpha;
    } else NonBlasMultMM(alpha,A,B,beta,C);
  }
#endif // NOFLOAT
#endif // ELAP
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
    TMVAssert(!C.isrm());

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tb())))
      BlasMultMM(alpha,A,B,beta,C);
    else if (!(C.iscm() && C.stepj()>0)) {
      Matrix<T,ColMajor> C2 = beta*C;
      DoMultMM(alpha,A,B,T(1),C2.View());
      C = C2;
    } else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
      if (IMAG(alpha) == RealType(T)(0)) {
	Matrix<Ta,ColMajor> A2 = REAL(alpha)*A;
	DoMultMM(T(1),A2,B,beta,C);
      } else {
	Matrix<T,ColMajor> A2 = alpha*A;
	DoMultMM(T(1),A2,B,beta,C);
      }
    } else if (!((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
      if (IMAG(alpha) == RealType(T)(0)) {
	Matrix<Tb,ColMajor> B2 = REAL(alpha)*B;
	DoMultMM(T(1),A,B2,beta,C);
      } else {
	Matrix<T,ColMajor> B2 = alpha*B;
	DoMultMM(T(1),A,B2,beta,C);
      }
    } else {
      BlasMultMM(alpha,A,B,beta,C);
    }
#else
    NonBlasMultMM(alpha,A,B,beta,C);
#endif
  }

  template <class T, class Ta, class Tb> inline void FullTempMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C   via a temporary
  {
    Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
    DoMultMM(T(1),A,B,T(0),C2.View());
    C = alpha * C2 + beta * C;
  }

  template <class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * C + beta * C
  {
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (C.isrm()) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM(T(1),B2.Transpose(),A.Transpose(),beta,
	      C.Cols(j,j2).Transpose());
	} else {
	  Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM(T(1),B2.Transpose(),A.Transpose(),beta,
	      C.Cols(j,j2).Transpose());
	}
      } else {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	} else {
	  Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	}
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> B2 = B;
    B2 *= alpha;
    Matrix<T> A2 = A;
    Matrix<T,ColMajor> C3(C.colsize(),C.rowsize());
    if (C.colsize()>0 && C.rowsize()>0)
      if (A.rowsize()==0) C3.Zero();
      else DoMultMM(T(1),A2,B2,T(0),C3.View());
    Matrix<T> C2 = C;
    C2 *= beta;
    C2 += C3;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0)) 
	C *= beta;
      else if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (C.isrm()) MultMM(alpha,B.Transpose(),A.Transpose(),
	  beta,C.Transpose());
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
    if (Norm(C2-C) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(C0)))) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
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


