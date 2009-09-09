
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  template <class T, class Ta, class Tx> void DoUnitAMultMV(
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    const size_t N = A.size();
    y = beta*y + A.LowerTri() * x;

    if (N > 1)
      y.SubVector(0,N-1) += A.UpperTri().OffDiag() * x.SubVector(1,N);
  }

  template <class T, class Ta, class Tx> void UnitAMultMV(
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

  template <class T, class Ta, class Tx> void NonBlasMultMV(
      const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(!x.SameStorageAs(y));

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

#ifdef BLAS
  template <class T, class Ta, class Tx> void BlasMultMV(
      const T alpha, const GenSymMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  { NonBlasMultMV(alpha,A,x,beta,y); }
  template <> void BlasMultMV(const double alpha,
      const GenSymMatrix<double>& A, const GenVector<double>& x,
      const double beta, const VectorView<double>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(!x.SameStorageAs(y));
    TMVAssert(A.isherm());

    cblas_dsymv(A.isrm()?CblasRowMajor:CblasColMajor,CblasLower,
	A.size(),alpha,A.cptr(),A.isrm()?A.stepi():A.stepj(),
	x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<complex<double> >& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(!x.SameStorageAs(y));

    if (A.isherm()) {
      cblas_zhemv(A.isconj()==A.isrm()?CblasColMajor:CblasRowMajor,
	  A.isconj()?CblasUpper:CblasLower,
	  A.size(),&alpha,A.cptr(),
	  A.isrm()?A.stepi():A.stepj(),
	  x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    } else {
#ifdef LAP
      char uplo = A.iscm() ? 'L' : 'U';
      int n = A.size();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      int incx = x.step();
      int incy = y.step();
      zsymv(&uplo,&n,LAP_Complex(&alpha),LAP_Complex(A.cptr()),&lda,
	  LAP_Complex(x.cptr()),&incx,LAP_Complex(&beta),
	  LAP_Complex(y.ptr()),&incy);
#else
      NonBlasMultMV(alpha,A,x,beta,y);
#endif
    }
  }
  template <> void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<double>& A,
      const GenVector<double>& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      if (REAL(beta) != double(1)) y.Imag() *= REAL(beta);
    } else {
      Vector<double> y1 = A*x;
      AddVV(alpha,y1,beta,y);
    }
  }
  template <> void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<double>& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else {
      Vector<complex<double> > y1(y.size());
      y1.Real() = A*x.Real();
      y1.Imag() = A*x.Imag();
      AddVV(alpha,y1,beta,y);
    }
  }
  template <> void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<complex<double> >& A,
      const GenVector<double>& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    Vector<complex<double> > x1 = x;
    BlasMultMV(alpha,A,x1,beta,y);
  }
#ifndef NOFLOAT
  template <> void BlasMultMV(const float alpha,
      const GenSymMatrix<float>& A, const GenVector<float>& x,
      const float beta, const VectorView<float>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(!x.SameStorageAs(y));
    TMVAssert(A.isherm());
    
    cblas_ssymv(A.isrm()?CblasRowMajor:CblasColMajor,CblasLower,
	A.size(),alpha,A.cptr(),A.isrm()?A.stepi():A.stepj(),
	x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<complex<float> >& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(!x.SameStorageAs(y));

    if (A.isherm()) {
      cblas_chemv(A.isconj()==A.isrm()?CblasColMajor:CblasRowMajor,
	  A.isconj()?CblasUpper:CblasLower,
	  A.size(),&alpha,A.cptr(),
	  A.isrm()?A.stepi():A.stepj(),
	  x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    } else {
#ifdef LAP
      char uplo = A.iscm() ? 'L' : 'U';
      int n = A.size();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      int incx = x.step();
      int incy = y.step();
      csymv(&uplo,&n,LAP_Complex(&alpha),LAP_Complex(A.cptr()),&lda,
	  LAP_Complex(x.cptr()),&incx,LAP_Complex(&beta),
	  LAP_Complex(y.ptr()),&incy);
#else
      NonBlasMultMV(alpha,A,x,beta,y);
#endif
    }
  }
  template <> void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<float>& A,
      const GenVector<float>& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      if (REAL(beta) != float(1)) y.Imag() *= REAL(beta);
    } else {
      Vector<float> y1 = A*x;
      AddVV(alpha,y1,beta,y);
    }
  }
  template <> void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<float>& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else {
      Vector<complex<float> > y1(y.size());
      y1.Real() = A*x.Real();
      y1.Imag() = A*x.Imag();
      AddVV(alpha,y1,beta,y);
    }
  }
  template <> void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<complex<float> >& A,
      const GenVector<float>& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    Vector<complex<float> > x1 = x;
    BlasMultMV(alpha,A,x1,beta,y);
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
    if (y.isconj()) DoMultMV(CONJ(alpha),A.Conjugate(),x.Conjugate(),
	CONJ(beta),y.Conjugate());
    else if (A.uplo() == Upper) {
      if (A.isherm()) DoMultMV(alpha,A.Adjoint(),x,beta,y);
      else DoMultMV(alpha,A.Transpose(),x,beta,y);
    } else {
#ifdef BLAS
      if (((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0))) {
	if (!A.isconj() && y.step() > 0) { 
	  if (!x.isconj() && x.step() > 0 && !x.SameStorageAs(y))
	    BlasMultMV(alpha,A,x,beta,y);
	  else {
	    Vector<T> xx = alpha*x;
	    BlasMultMV(T(1),A,xx,beta,y);
	  }
	} else {
	  Vector<T> yy(y.size());
	  if (A.isconj()) {
	    if (x.isconj() && x.step() > 0) {
	      BlasMultMV(T(1),A.Conjugate(),x.Conjugate(),T(0),yy.View());
	      AddVV(alpha,yy.Conjugate(),beta,y);
	    } else {
	      Vector<T> xx = CONJ(alpha)*x.Conjugate();
	      BlasMultMV(T(1),A.Conjugate(),xx,T(0),yy.View());
	      AddVV(T(1),yy.Conjugate(),beta,y);
	    }
	  } else {
	    if (!x.isconj() && x.step() > 0) {
	      BlasMultMV(T(1),A,x,T(0),yy.View());
	      AddVV(alpha,yy,beta,y);
	    } else {
	      Vector<T> xx = alpha*x;
	      BlasMultMV(T(1),A,xx,T(0),yy.View());
	      AddVV(T(1),yy,beta,y);
	    }
	  }
	}
      }
      else {
	if (A.isherm()) {
	  HermMatrix<Ta,Upper,RowMajor> A2 = alpha*A;
	  DoMultMV(T(1),A2,x,beta,y);
	} else {
	  SymMatrix<Ta,Upper,RowMajor> A2 = alpha*A;
	  DoMultMV(T(1),A2,x,beta,y);
	}
      }
#else
      NonBlasMultMV(alpha,A,x,beta,y);
#endif
    }
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
    Vector<T> y2 = alpha*Matrix<T>(A)*x+beta*y;
    Vector<T> y0 = y;
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
    if (Norm(y-y2) > 0.001*max(RealType(T)(1),Norm(y))) {
      cerr<<"MultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
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


