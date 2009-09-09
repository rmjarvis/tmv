
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_MM_BLOCKSIZE 64
#define SYM_MM_BLOCKSIZE2 32
#endif

  //
  // MultMM
  //

  template <class T, class Ta, class Tb> inline void RRowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RRowMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();
    for(size_t j=0;j<N;++j) {
      C.row(j) = beta * C.row(j) + alpha * A.row(j,0,j+1) * B.Rows(0,j+1);
      C.Rows(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
    }
  }

  template <class T, class Ta, class Tb> inline void CRowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CRowMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();
    for(int j=N-1;j>=0;--j) {
      C.row(j) = beta * C.row(j) + alpha * A.row(j,j,N) * B.Rows(j,N);
      C.Rows(j+1,N) += alpha * A.col(j,j+1,N) ^ B.row(j);
    }
  }

  template <class T, class Ta, class Tb> inline void RowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (A.iscm()) CRowMultMM(alpha,A,B,beta,C);
    else RRowMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> inline void ColMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    for(size_t j=0;j<C.rowsize();++j) 
      C.col(j) = beta * C.col(j) + alpha * A * B.col(j);
  }

  template <class T, class Ta, class Tb> inline void RecursiveMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"NonBlasMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();
    if (N <= SYM_MM_BLOCKSIZE2) {
      if (B.isrm() && C.isrm()) RowMultMM(alpha,A,B,beta,C);
      else if (B.iscm() && C.iscm()) ColMultMM(alpha,A,B,beta,C);
      else if (C.colsize() < C.rowsize()) RowMultMM(alpha,A,B,beta,C);
      else ColMultMM(alpha,A,B,beta,C);
    } else {
      size_t k = N/2;
      const size_t nb = SYM_MM_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;

      // [ A00 A10t ] [ B0 ] = [ A00 B0 + A10t B1 ]
      // [ A10 A11  ] [ B1 ]   [ A10 B0 + A11 B1  ]

      ConstSymMatrixView<Ta> A00 = A.SubSymMatrix(0,k);
      ConstSymMatrixView<Ta> A11 = A.SubSymMatrix(k,N);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      ConstMatrixView<Tb> B0 = B.Rows(0,k);
      ConstMatrixView<Tb> B1 = B.Rows(k,N);
      MatrixView<T> C0 = C.Rows(0,k);
      MatrixView<T> C1 = C.Rows(k,N);

      RecursiveMultMM(alpha,A00,B0,beta,C0);
      RecursiveMultMM(alpha,A11,B1,beta,C1);
      C1 += alpha * A10 * B0;
      if (A.issym())
	C0 += alpha * A10.Transpose() * B1;
      else
	C0 += alpha * A10.Adjoint() * B1;
    }
  }

  template <class T, class Ta, class Tb> inline void NonBlasMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"NonBlasMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

    if (A.uplo() == Upper)
      if (A.isherm()) NonBlasMultMM(alpha,A.Adjoint(),B,beta,C);
      else NonBlasMultMM(alpha,A.Transpose(),B,beta,C);
    else if (C.isconj())
      NonBlasMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
    else RecursiveMultMM(alpha,A,B,beta,C);
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> inline void BlasMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  { NonBlasMultMM(alpha,A,B,beta,C); }
  template <> inline void BlasMultMM(
      double alpha, const GenSymMatrix<double>& A, const GenMatrix<double>& B,
      double beta, const MatrixView<double>& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM double\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    BLASNAME(dsymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
	A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(B.cptr()),BLASV(ldb),BLASV(beta),
	BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> inline void BlasMultMM(
      complex<double> alpha, const GenSymMatrix<complex<double> >& A,
      const GenMatrix<complex<double> >& B,
      complex<double> beta, const MatrixView<complex<double> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM c double\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==C.ct());
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    if (A.issym())
      BLASNAME(zsymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
	  A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.cptr()),BLASV(ldb),BLASP(&beta),
	  BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    else
      BLASNAME(zhemm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
	  A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.cptr()),BLASV(ldb),BLASP(&beta),
	  BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMM(
      float alpha, const GenSymMatrix<float>& A, const GenMatrix<float>& B,
      float beta, const MatrixView<float>& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM float\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    BLASNAME(ssymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
	A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(B.cptr()),BLASV(ldb),BLASV(beta),
	BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> inline void BlasMultMM(
      complex<float> alpha, const GenSymMatrix<complex<float> >& A,
      const GenMatrix<complex<float> >& B,
      complex<float> beta, const MatrixView<complex<float> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM c float\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==C.ct());
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    if (A.issym())
      BLASNAME(csymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
	  A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.cptr()),BLASV(ldb),BLASP(&beta),
	  BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    else
      BLASNAME(chemm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
	  A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.cptr()),BLASV(ldb),BLASP(&beta),
	  BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
#endif //NOFLOAT
#endif // BLAS

  template <class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    //cerr<<"Start DoMultMM: A = "<<Type(A)<<endl;
    //cerr<<"B = "<<Type(B)<<endl;
    //cerr<<"C = "<<Type(C)<<endl;
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tb())))
      BlasMultMM(alpha,A,B,beta,C);
    else if (A.isrm())
      DoMultMM(alpha,A.issym()?A.Transpose():A.Adjoint(),B,beta,C);
    else if (A.isconj())
      DoMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),CONJ(beta),
	  C.Conjugate());
    else if (!((C.isrm() && C.stepi()>0) || (C.iscm() && C.stepj()>0)) ||
	(C.iscm() && C.isconj()) || (C.isrm() && C.isconj()==A.issym())) {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM(T(1),A,B,T(0),C2.View());
      AddMM(alpha,C2,beta,C);
    } else if (!(A.iscm() && A.stepj()>0)) {
      if (IMAG(alpha) == RealType(T)(0)) {
	if (A.isherm()) {
	  if (A.uplo() == Upper) {
	    HermMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  } else {
	    HermMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  }
	} else {
	  if (A.uplo() == Upper) {
	    SymMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  } else {
	    SymMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  }
	}
      } else {
	if (!A.issym()) {
	  if (A.uplo() == Upper) {
	    HermMatrix<T,Upper,ColMajor> A2 = alpha*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  } else {
	    HermMatrix<T,Lower,ColMajor> A2 = alpha*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  }
	} else {
	  if (A.uplo() == Upper) {
	    SymMatrix<T,Upper,ColMajor> A2 = alpha*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  } else {
	    SymMatrix<T,Lower,ColMajor> A2 = alpha*A;
	    DoMultMM(T(1),A2,B,beta,C);
	  }
	}
      }
    } else if ((B.stor() != C.stor()) || B.isconj() != C.isconj() || 
	!((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
      if (IMAG(alpha) == RealType(T)(0)) {
	if (C.isconj()) {
	  if (C.iscm()) {
	    Matrix<Tb,ColMajor> B2 = REAL(alpha)*B.Conjugate();
	    DoMultMM(T(1),A,B2.Conjugate(),beta,C);
	  } else {
	    Matrix<Tb,RowMajor> B2 = REAL(alpha)*B.Conjugate();
	    DoMultMM(T(1),A,B2.Conjugate(),beta,C);
	  }
	} else {
	  if (C.iscm()) {
	    Matrix<Tb,ColMajor> B2 = REAL(alpha)*B;
	    DoMultMM(T(1),A,B2,beta,C);
	  } else {
	    Matrix<Tb,RowMajor> B2 = REAL(alpha)*B;
	    DoMultMM(T(1),A,B2,beta,C);
	  }
	}
      } else {
	if (C.isconj()) {
	  if (C.iscm()) {
	    Matrix<T,ColMajor> B2 = CONJ(alpha)*B.Conjugate();
	    DoMultMM(T(1),A,B2.Conjugate(),beta,C);
	  } else {
	    Matrix<T,RowMajor> B2 = CONJ(alpha)*B.Conjugate();
	    DoMultMM(T(1),A,B2.Conjugate(),beta,C);
	  }
	} else {
	  if (C.iscm()) {
	    Matrix<T,ColMajor> B2 = alpha*B;
	    DoMultMM(T(1),A,B2,beta,C);
	  } else {
	    Matrix<T,RowMajor> B2 = alpha*B;
	    DoMultMM(T(1),A,B2,beta,C);
	  }
	}
      }
    } else {
      BlasMultMM(alpha,A,B,beta,C);
    }
#else
    NonBlasMultMM(alpha,A,B,beta,C);
#endif
  }

  template <class T, class Ta, class Tb> inline void FullTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = beta*C+alpha*C2;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = beta*C+alpha*C2;
    }
  }

  template <class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    for(size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+SYM_MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
	if (C.isrm()) {
	  Matrix<Tb,RowMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	} else {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	}
      } else {
	if (C.isrm()) {
	  Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
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
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    //cerr<<"Start MultMM: alpha = "<<alpha<<", beta = "<<beta<<endl;
    //cerr<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<B.cptr()<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<C.cptr()<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0+beta*C0;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha == T(0)) 
	C *= beta;
      else if (A.SameStorageAs(C)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (B.SameStorageAs(C)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else DoMultMM(alpha, A, B, beta, C);
    }

#ifdef XDEBUG
    //cerr<<"Done: C = "<<C<<endl;
    if (Norm(C-C2) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
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

  template <class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = A.size();

    for(size_t j=0;j<N;) {
      size_t j2 = min(N,j+SYM_MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
	if (C.isrm()) {
	  Matrix<Tb,RowMajor> B2(N,j2-j);
	  B2.Rows(0,j) = REAL(alpha) * B.SubMatrix(0,j,j,j2);
	  B2.Rows(j,j2) = REAL(alpha) * B.SubSymMatrix(j,j2);
	  B2.Rows(j2,N) = REAL(alpha) * B.SubMatrix(j2,N,j,j2);
	  DoMultMM(T(1),A,B2.View(),beta,C.Cols(j,j2));
	} else {
	  Matrix<Tb,ColMajor> B2(N,j2-j);
	  B2.Rows(0,j) = REAL(alpha) * B.SubMatrix(0,j,j,j2);
	  B2.Rows(j,j2) = REAL(alpha) * B.SubSymMatrix(j,j2);
	  B2.Rows(j2,N) = REAL(alpha) * B.SubMatrix(j2,N,j,j2);
	  DoMultMM(T(1),A,B2.View(),beta,C.Cols(j,j2));
	}
      } else {
	if (C.isrm()) {
	  Matrix<T,RowMajor> B2(N,j2-j);
	  B2.Rows(0,j) = alpha * B.SubMatrix(0,j,j,j2);
	  B2.Rows(j,j2) = alpha * B.SubSymMatrix(j,j2);
	  B2.Rows(j2,N) = alpha * B.SubMatrix(j2,N,j,j2);
	  DoMultMM(T(1),A,B2.View(),beta,C.Cols(j,j2));
	} else {
	  Matrix<T,ColMajor> B2(N,j2-j);
	  B2.Rows(0,j) = alpha * B.SubMatrix(0,j,j,j2);
	  B2.Rows(j,j2) = alpha * B.SubSymMatrix(j,j2);
	  B2.Rows(j2,N) = alpha * B.SubMatrix(j2,N,j,j2);
	  DoMultMM(T(1),A,B2.View(),beta,C.Cols(j,j2));
	}
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> inline void FullTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      BlockTempMultMM(T(1),A,B,T(0),C2.View());
      C = alpha * C2 + beta * C;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      BlockTempMultMM(T(1),A,B,T(0),C2.View());
      C = alpha * C2 + beta * C;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());
#ifdef XDEBUG
    //cerr<<"Start MultMM: alpha = "<<alpha<<", beta = "<<beta<<endl;
    //cerr<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<B.cptr()<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<C.cptr()<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0+beta*C0;
#endif

    if (A.size() > 0) {
      if (alpha == T(0)) 
	C *= beta;
      else if (A.SameStorageAs(C) || B.SameStorageAs(C))
	FullTempMultMM(alpha,A,B,beta,C);
      else BlockTempMultMM(alpha, A, B, beta, C);
    }

#ifdef XDEBUG
    //cerr<<"Done: C = "<<C<<endl;
    if (Norm(C-C2) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
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

#define InstFile "TMV_SymMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


