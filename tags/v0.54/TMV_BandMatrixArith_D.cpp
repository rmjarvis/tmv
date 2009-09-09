
#include "TMV.h"
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define MM_BLOCKSIZE TMV_BLOCKSIZE
#else
#define MM_BLOCKSIZE 64
#endif

  //
  // MultMM
  //

  template <class T, class Ta, class Tb> inline void RowMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RowMultMM: alpha,beta = "<<alpha<<','<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
#endif

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    size_t j1=0;
    size_t k=A.nlo();
    size_t j2=A.nhi()+1;

    for(size_t i=0;i<C.colsize(); ++i) {
      // C.row(i) = beta * C.row(i) + alpha * A.row(i,j1,j2) * B.Rows(j1,j2);
      MultMV(alpha,B.Rows(j1,j2).Transpose(),A.row(i,j1,j2),
	  beta,C.row(i));
      if (k>0) --k; else ++j1;
      if (j2<A.rowsize()) ++j2;
      else if (j1==A.rowsize()) {
	C.Rows(i+1,C.colsize()) *= beta;
	break;
      }
    }
#ifdef XDEBUG
    //cerr<<"Done: C = "<<C<<endl;
#endif
  }

  template <class T, class Ta, class Tb> inline void OPMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"OPMultMM: alpha,beta = "<<alpha<<','<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
#endif

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    size_t i1=0;
    size_t k=A.nhi();
    size_t i2=A.nlo()+1;

    C *= beta;
    for(size_t j=0;j<A.rowsize();++j) {
      C.Rows(i1,i2) += alpha * A.col(j,i1,i2) ^ B.row(j);
      if (k>0) --k; else ++i1;
      if (i2<A.colsize()) ++i2;
      else if (i1==A.colsize()) break;
    }
#ifdef XDEBUG
    //cerr<<"Done: C = "<<C<<endl;
#endif
  }

  template <class T, class Ta, class Tb> inline void ColMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColMultMM: alpha,beta = "<<alpha<<','<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
#endif

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    for(size_t j=0;j<B.rowsize();j++)
      // C.col(j) = beta * C.col(j) + alpha * A * B.col(j);
      MultMV(alpha,A,B.col(j),beta,C.col(j));
#ifdef XDEBUG
    //cerr<<"Done: C = "<<C<<endl;
#endif
  }

  template <bool add, class T, class Ta, class Tb> 
    inline void NonLapTriDiagMultMM(
	const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
	const MatrixView<T>& C)
    {
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      TMVAssert(A.rowsize() > 0);
      TMVAssert(C.rowsize() > 0);
      TMVAssert(C.colsize() > 0);
      TMVAssert(A.ct()==NonConj);
      TMVAssert(A.nlo() == 1);
      TMVAssert(A.nhi() == 1);
      TMVAssert(A.isdm());

      const size_t N = A.diag().size();
      const size_t M = A.rowsize()>A.colsize() ? N : N-1;

      const Ta* di = A.cptr();
      const Ta* dui = A.diag(1).cptr();
      const Ta* dli = A.diag(-1).cptr()-1;

      for(size_t i=0;i<N;++i,++di,++dui,++dli) {
	if (add) C.row(i) += *di*B.row(i);
	else C.row(i) = *di*B.row(i);
	if (i>0) C.row(i) += *dli*B.row(i-1);
	if (i<M) C.row(i) += *dui*B.row(i+1);
      }
      if (A.colsize() > A.rowsize()) {
	if (add) C.row(N) += *dli*B.row(N-1);
	else C.row(N) = *dli*B.row(N-1);
      }
    }

#ifdef ELAP
  template <class T, class Ta, class Tb> inline void LapTriDiagMultMM(
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B, int beta,
      const MatrixView<T>& C)
  {
    if (beta) NonLapTriDiagMultMM<true>(A,B,C);
    else NonLapTriDiagMultMM<false>(A,B,C);
  }
  template <> inline void LapTriDiagMultMM(
      const GenBandMatrix<double>& A, const GenMatrix<double>& B, int beta,
      const MatrixView<double>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.nlo() == 1);
    TMVAssert(A.nhi() == 1);
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());

    int n = A.colsize();
    int nrhs = B.rowsize();
    double a(1);
    int ldB = B.stepj();
    double b(beta);
    int ldC = C.stepj();
    LAPNAMEX(dlagtm) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(b),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
  template <> inline void LapTriDiagMultMM(
      const GenBandMatrix<complex<double> >& A, 
      const GenMatrix<complex<double> >& B, int beta,
      const MatrixView<complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.nlo() == 1);
    TMVAssert(A.nhi() == 1);
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());
    TMVAssert(B.isconj() == C.isconj());

    int n = A.colsize();
    int nrhs = B.rowsize();
    double a(1);
    int ldB = B.stepj();
    double b(beta);
    int ldC = C.stepj();
    LAPNAMEX(zlagtm) (LAPCM B.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(b),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
#ifndef NOFLOAT
  template <> inline void LapTriDiagMultMM(
      const GenBandMatrix<float>& A, const GenMatrix<float>& B, int beta,
      const MatrixView<float>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.nlo() == 1);
    TMVAssert(A.nhi() == 1);
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());

    int n = A.colsize();
    int nrhs = B.rowsize();
    float a(1);
    int ldB = B.stepj();
    float b(beta);
    int ldC = C.stepj();
    LAPNAMEX(slagtm) (LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(b),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
  template <> inline void LapTriDiagMultMM(
      const GenBandMatrix<complex<float> >& A, 
      const GenMatrix<complex<float> >& B, int beta,
      const MatrixView<complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.nlo() == 1);
    TMVAssert(A.nhi() == 1);
    TMVAssert(A.isdm());
    TMVAssert(A.IsSquare());
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());
    TMVAssert(B.isconj() == C.isconj());

    int n = A.colsize();
    int nrhs = B.rowsize();
    float a(1);
    int ldB = B.stepj();
    float b(beta);
    int ldC = C.stepj();
    LAPNAMEX(clagtm) (LAPCM B.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
	LAPV(a),LAPP(A.cptr()+A.stepj()),
	LAPP(A.cptr()),LAPP(A.cptr()+A.stepi()),
	LAPP(B.cptr()),LAPV(ldB),LAPV(b),LAPP(C.ptr()),LAPV(ldC) LAP1);
  }
#endif // FLOAT
#endif // ELAP

  template <bool add, class T, class Ta, class Tb> inline void DoTriDiagMultMM(
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef ELAP
    if (A.IsSquare() && B.iscm() && C.iscm() && B.isconj() == C.isconj())
      LapTriDiagMultMM(A,B,add?1:0,C);
    else
#endif
      NonLapTriDiagMultMM<add>(A,B,C);
  }

  template <class T, class Ta, class Tb> inline void TriDiagMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (beta == T(0)) {
      if (alpha == T(1) && A.isdm()) {
	if (A.isconj()) 
	  DoTriDiagMultMM<false>(A.Conjugate(),B.Conjugate(),C.Conjugate());
	else
	  DoTriDiagMultMM<false>(A,B,C);
      } else if (IMAG(alpha) == RealType(T)(0)) {
	BandMatrix<Ta,DiagMajor> A1 = REAL(alpha)*A;
	DoTriDiagMultMM<false>(A1,B,C);
      } else {
	BandMatrix<T,DiagMajor> A1 = alpha*A;
	DoTriDiagMultMM<false>(A1,B,C);
      }
    } else {
      if (beta != T(1)) C *= beta;
      if (alpha == T(1) && A.isdm()) {
	if (A.isconj()) 
	  DoTriDiagMultMM<true>(A.Conjugate(),B.Conjugate(),C.Conjugate());
	else
	  DoTriDiagMultMM<true>(A,B,C);
      } else if (IMAG(alpha) == RealType(T)(0)) {
	BandMatrix<Ta,DiagMajor> A1 = REAL(alpha)*A;
	DoTriDiagMultMM<true>(A1,B,C);
      } else {
	BandMatrix<Ta,DiagMajor> A1 = alpha*A;
	DoTriDiagMultMM<true>(A1,B,C);
      }
    }
  }

  // MJ: Put in a recursive block calculation here.  (Also in Band*Band)
  template <class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (A.isrm() && C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (A.iscm() && B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else if (B.iscm() && C.iscm()) ColMultMM(alpha,A,B,beta,C);
    else if (A.nlo() == 1 && A.nhi() == 1) {
      if (A.isconj())
	TriDiagMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	    CONJ(beta),C.Conjugate());
      else
	TriDiagMultMM(alpha,A,B,beta,C);
    }
    else if (C.colsize() < C.rowsize()) RowMultMM(alpha,A,B,beta,C);
    else ColMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> inline void FullTempMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = alpha*C2+beta*C;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = alpha*C2+beta*C;
    }
  }

  template <class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
	if (C.isrm()) {
	  Matrix<Tb,RowMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	} else  {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	}
      } else {
	if (C.isrm()) {
	  Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	} else  {
	  Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM(T(1),A,B2,beta,C.Cols(j,j2));
	}
      }
      j=j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

#ifdef XDEBUG
    //cerr<<"Start MultMM:\n";
    //cerr<<"alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"beta = "<<beta<<endl;
    //cerr<<"C = "<<Type(C);
    //if (beta != T(0)) cerr<<"  "<<C;
    //cerr<<endl;
    Matrix<Tb> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C;
    C2 *= beta;
    C2 += alpha*A0*B0;
    //cerr<<"C2 = "<<C2<<endl;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0)) {
	C *= beta;
      } else if (A.rowsize() > A.colsize()+A.nhi()) {
	MultMM(alpha,A.Cols(0,A.colsize()+A.nhi()),
	    B.Rows(0,A.colsize()+A.nhi()),beta,C);
      } else if (A.colsize() > A.rowsize()+A.nlo()) {
	MultMM(alpha,A.Rows(0,A.rowsize()+A.nlo()),
	    B,beta,C.Rows(0,A.rowsize()+A.nlo()));
	C.Rows(A.rowsize()+A.nlo(),A.colsize()) *= beta;
      } else if (C.isconj()) {
	MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	    CONJ(beta),C.Conjugate());
      } else if (C.SameStorageAs(A)) {
	FullTempMultMM(alpha,A,B,beta,C);
      } else if (C.SameStorageAs(B)) {
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else 
	  FullTempMultMM(alpha,A,B,beta,C);
      } else {
	DoMultMM(alpha, A, B, beta, C);
      }
    }
#ifdef XDEBUG
    //cerr<<"C -> "<<C<<endl;
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

#define InstFile "TMV_BandMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv
