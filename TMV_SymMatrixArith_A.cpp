
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  template <class T, class Ta, class Tx> inline void DoUnitAMultMV(
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    const size_t N = A.size();
    y = beta*y + A.LowerTri() * x;

    if (N > 1)
      y.SubVector(0,N-1) += A.UpperTri().OffDiag() * x.SubVector(1,N);
  }

  template <class T, class Ta, class Tx> inline void UnitAMultMV(
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    // Check for 0's in the beginning or end of x:
    //      [ A11 A12 A13 ] [ 0 ]          [ A12 ]
    // y += [ A21 A22 A23 ] [ x ] --> y += [ A22 ] x
    //      [ A31 A32 A33 ] [ 0 ]          [ A32 ]

    const size_t N = x.size(); // == A.size()
    size_t j2 = N;
    for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
    if (j2 == 0) {
      y *= beta;
      return;
    }
    size_t j1 = 0;
    for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) DoUnitAMultMV(A,x,beta,y);
    else {
      if (j1 > 0)
	MultMV(T(1),A.SubMatrix(0,j1,j1,j2),x.SubVector(j1,j2),
	    beta,y.SubVector(0,j1));
      TMVAssert(j1 != j2);
      DoUnitAMultMV(A.SubSymMatrix(j1,j2),x.SubVector(j1,j2),
	  beta,y.SubVector(j1,j2));
      if (j2 < N)
	MultMV(T(1),A.SubMatrix(j2,N,j1,j2),x.SubVector(j1,j2),
	    beta,y.SubVector(j2,N));
    }
  }

  template <class T, class Ta, class Tx> inline void NonBlasMultMV(
      const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(y.size() > 0);
    TMVAssert(!x.SameStorageAs(y));

    if (A.uplo() == Upper) 
      if (A.isherm()) NonBlasMultMV(alpha,A.Adjoint(),x,beta,y);
      else NonBlasMultMV(alpha,A.Transpose(),x,beta,y);
    else if (y.isconj())
      NonBlasMultMV(CONJ(alpha),A.Conjugate(),x.Conjugate(),
	  CONJ(beta),y.Conjugate());
    else {
      if (x.step() != 1) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Vector<Tx> xx = REAL(alpha)*x;
	  if (y.step() != 1) {
	    Vector<T> yy(y.size());
	    UnitAMultMV(A,xx,T(0),yy.View());
	    AddVV(T(1),yy,beta,y);
	  }
	  else 
	    UnitAMultMV(A,xx,beta,y);
	} else {
	  Vector<T> xx = alpha*x;
	  if (y.step()!=1) {
	    Vector<T> yy(y.size());
	    UnitAMultMV(A,xx,T(0),yy.View());
	    AddVV(T(1),yy,beta,y);
	  }
	  else
	    UnitAMultMV(A,xx,beta,y);
	}
      } else if (y.step()!=1 || alpha!=T(1)) {
	Vector<T> yy(y.size());
	UnitAMultMV(A,x,T(0),yy.View());
	AddVV(alpha,yy,beta,y);
      } else {
	TMVAssert(alpha == T(1));
	TMVAssert(y.step() == 1);
	TMVAssert(x.step() == 1);
	UnitAMultMV(A,x,beta,y);
      }
    }
  }

#ifdef BLAS
  template <class T, class Ta, class Tx> inline void BlasMultMV(
      const T alpha, const GenSymMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  { NonBlasMultMV(alpha,A,x,beta,y); }
  template <> inline void BlasMultMV(const double alpha,
      const GenSymMatrix<double>& A, const GenVector<double>& x,
      const double beta, const VectorView<double>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));

    int n = A.size();
    int lda = A.stepj();
    int xs = x.step();
    int ys = y.step();
    BLASNAME(dsymv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(beta),
	BLASP(y.ptr()),BLASV(ys) BLAS1);
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<complex<double> >& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));

    if (A.isherm()) {
      int n = A.size();
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      BLASNAME(zhemv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&beta),
	  BLASP(y.ptr()),BLASV(ys) BLAS1);
    } else {
#ifdef ELAP
      int n = A.size();
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      LAPNAMEX(zsymv) (LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
	  LAPV(n),LAPP(&alpha),LAPP(A.cptr()),LAPV(lda),
	  LAPP(x.cptr()),LAPV(xs),LAPP(&beta),
	  LAPP(y.ptr()),LAPV(ys) LAP1);
#else
      NonBlasMultMV(alpha,A,x,beta,y);
#endif
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMV(const float alpha,
      const GenSymMatrix<float>& A, const GenVector<float>& x,
      const float beta, const VectorView<float>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));
    TMVAssert(A.isherm());

    int n = A.size();
    int lda = A.stepj();
    int xs = x.step();
    int ys = y.step();
    BLASNAME(ssymv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(beta),
	BLASP(y.ptr()),BLASV(ys) BLAS1);
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<complex<float> >& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));

    if (A.isherm()) {
      int n = A.size();
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      BLASNAME(chemv) (BLASCM A.uplo() == Upper?BLASCH_UP:BLASCH_LO, 
	  BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&beta),
	  BLASP(y.ptr()),BLASV(ys) BLAS1);
    } else {
#ifdef ELAP
      int n = A.size();
      int lda = A.stepj();
      int xs = x.step();
      int ys = y.step();
      LAPNAMEX(csymv) (LAPCM A.uplo()==Upper ? LAPCH_UP : LAPCH_LO,
	  LAPV(n),LAPP(&alpha),LAPP(A.cptr()),LAPV(lda),
	  LAPP(x.cptr()),LAPV(xs),LAPP(&beta),
	  LAPP(y.ptr()),LAPV(ys) LAP1);
#else
      NonBlasMultMV(alpha,A,x,beta,y);
#endif
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta, class Tx> inline void DoMultMV(
      const T alpha, const GenSymMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!x.SameStorageAs(y));

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tx())))
      BlasMultMV(alpha,A,x,beta,y);
    else if (A.isconj()) DoMultMV(CONJ(alpha),A.Conjugate(),x.Conjugate(),
	CONJ(beta),y.Conjugate());
    else if (A.isrm()) {
      if (A.isherm()) DoMultMV(alpha,A.Adjoint(),x,beta,y);
      else DoMultMV(alpha,A.Transpose(),x,beta,y);
    } else {
      if (A.iscm()&&A.stepj()>0) {
	if (!y.isconj() && y.step() > 0) { 
	  if (!x.isconj() && x.step() > 0)
	    BlasMultMV(alpha,A,x,beta,y);
	  else {
	    Vector<T> xx = alpha*x;
	    BlasMultMV(T(1),A,xx,beta,y);
	  }
	} else {
	  Vector<T> yy(y.size());
	  if (!x.isconj() && x.step() > 0) {
	    BlasMultMV(T(1),A,x,T(0),yy.View());
	    AddVV(alpha,yy,beta,y);
	  } else {
	    Vector<T> xx = alpha*x;
	    BlasMultMV(T(1),A,xx,T(0),yy.View());
	    AddVV(T(1),yy,beta,y);
	  }
	}
      } else {
	if (IMAG(alpha) == RealType(T)(0)) {
	  if (A.isherm()) {
	    if (A.uplo() == Upper) {
	      HermMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    } else {
	      HermMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    }
	  } else {
	    if (A.uplo() == Upper) {
	      SymMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    } else {
	      SymMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    }
	  }
	} else {
	  if (A.isherm()) {
	    if (A.uplo() == Upper) {
	      HermMatrix<T,Upper,ColMajor> A2 = alpha*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    } else {
	      HermMatrix<T,Lower,ColMajor> A2 = alpha*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    }
	  } else {
	    if (A.uplo() == Upper) {
	      SymMatrix<T,Upper,ColMajor> A2 = alpha*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    } else {
	      SymMatrix<T,Lower,ColMajor> A2 = alpha*A;
	      DoMultMV(T(1),A2,x,beta,y);
	    }
	  }
	}
      }
    }
#else
    NonBlasMultMV(alpha,A,x,beta,y);
#endif
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
#ifdef XDEBUG
    //cerr<<"Start MultMV: alpha, beta = "<<alpha<<"  "<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
    Matrix<T> A0 = A;
    Vector<T> x0 = x;
    Vector<T> y0 = y;
    Vector<T> y2 = alpha*A0*x0+beta*y0;
#endif

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) {
	y *= beta;
      } else if (x.SameStorageAs(y)) {
	Vector<T> yy(y.size());
	DoMultMV(T(1),A,x,T(0),yy.View());
	AddVV(alpha,yy,beta,y);
      } else {
	DoMultMV(alpha,A,x,beta,y);
      } 
    }
#ifdef XDEBUG
    //cerr<<"--> y = "<<y<<endl;
    //cerr<<"y2 = "<<y2<<endl;
    if (Norm(y-y2) > 0.001*(abs(alpha)*Norm(A0)*Norm(x0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(y0)))) {
      cerr<<"MultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"--> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


