
#include "TMV.h"
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t MM_BLOCKSIZE = 64;
#endif

  //
  // MultMM
  //

  template <class T, class Ta, class Tb> void RowMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
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

  template <class T, class Ta, class Tb> void OPMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
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

  template <class T, class Ta, class Tb> void ColMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
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

  // MJ: Lap lagtm is Tridiagonal * Matrix
  // MJ: Put in a recursive block calculation here.  (Also in Band*Band)
  template <class T, class Ta, class Tb> inline void DoMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
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
    else if (C.colsize() < C.rowsize()) RowMultMM(alpha,A,B,beta,C);
    else ColMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
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

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
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
    Matrix<Tb> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C;
    C2 *= beta;
    C2 += alpha*A0*B0;
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
    if (Norm(C2-C) > 0.001*max(RealType(T)(1),Norm(C)+Norm(B)+Norm(A))) {
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
