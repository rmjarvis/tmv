
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  // In TMV_TriMatrixArith_A.cpp:
  template <bool rm, bool ca, bool ua, class T, class Ta> 
    extern void DoRowMultEqMV(
	const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x);
  template <bool cm, bool ca, bool ua, class T, class Ta> 
    extern void DoColMultEqMV(
	const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x);
  template <bool rm, bool cx, bool ca, class T, class Ta, class Tx> 
    extern void DoRowAddMultMV(
	const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y);
  template <bool cm, bool cx, bool ca, class T, class Ta, class Tx> 
    extern void DoColAddMultMV(
	const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y);
  template <bool rm, bool cx, bool ca, class T, class Ta, class Tx> 
    extern void DoRowAddMultMV(
	const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y);
  template <bool cm, bool cx, bool ca, class T, class Ta, class Tx> 
    extern void DoColAddMultMV(
	const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y);

  template <bool b0, bool rm, bool cx, class T, class Ta, class Tx>
    inline void RowMultMV(const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
#ifdef XDEBUG
      //cerr<<"RowMultMV\n";
      //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
#endif
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());
      TMVAssert(rm == A.isrm());

      const size_t N = A.size();

      if (b0) {
	y = x;
	if (IsComplex(T()) && A.isconj())
	  DoRowMultEqMV<rm,true,false>(A.LowerTri(),y);
	else
	  DoRowMultEqMV<rm,false,false>(A.LowerTri(),y);
      } else {
	if (IsComplex(T()) && A.isconj())
	  DoRowAddMultMV<rm,cx,true>(A.LowerTri(),x,y);
	else
	  DoRowAddMultMV<rm,cx,false>(A.LowerTri(),x,y);
      }

      if (N > 1)
	if (IsComplex(T()) && A.isconj() == A.issym())
	  DoColAddMultMV<rm,cx,true>(A.UpperTri().OffDiag(),
	      x.SubVector(1,N),y.SubVector(0,N-1));
	else
	  DoColAddMultMV<rm,cx,false>(A.UpperTri().OffDiag(),
	      x.SubVector(1,N),y.SubVector(0,N-1));
    }

  template <bool b0, bool cm, bool cx, class T, class Ta, class Tx>
    inline void ColMultMV(const GenSymMatrix<Ta>& A,
	const GenVector<Tx>& x, const VectorView<T>& y)
    {
#ifdef XDEBUG
      //cerr<<"ColMultMV\n";
      //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
#endif
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());
      TMVAssert(cm == A.iscm());

      const size_t N = A.size();

      if (b0) {
	y = x;
	if (IsComplex(T()) && A.isconj())
	  DoColMultEqMV<cm,true,false>(A.LowerTri(),y);
	else
	  DoColMultEqMV<cm,false,false>(A.LowerTri(),y);
      } else {
	if (IsComplex(T()) && A.isconj())
	  DoColAddMultMV<cm,cx,true>(A.LowerTri(),x,y);
	else
	  DoColAddMultMV<cm,cx,false>(A.LowerTri(),x,y);
      }

      if (N > 1) 
	if (IsComplex(T()) && A.isconj() == A.issym())
	  DoRowAddMultMV<cm,cx,true>(A.UpperTri().OffDiag(),
	      x.SubVector(1,N),y.SubVector(0,N-1));
	else
	  DoRowAddMultMV<cm,cx,false>(A.UpperTri().OffDiag(),
	      x.SubVector(1,N),y.SubVector(0,N-1));
    }

  template <bool b0, bool cx, class T, class Ta, class Tx>
    void UnitAMultMV1(
	const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step() == 1);
      TMVAssert(y.step() == 1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());

      if (A.isrm())
	RowMultMV<b0,true,cx>(A,x,y);
      else if (A.iscm())
	ColMultMV<b0,true,cx>(A,x,y);
      else 
	RowMultMV<b0,false,cx>(A,x,y);
    }

  template <bool b0, bool cx, class T, class Ta, class Tx>
    extern void UnitAMultMV1(const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y);

  template <bool b0, bool cx, class T, class Ta, class Tx> void UnitAMultMV(
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    // Check for 0's in the beginning or end of x:
    //      [ A11 A12 A13 ] [ 0 ]          [ A12 ]
    // y += [ A21 A22 A23 ] [ x ] --> y += [ A22 ] x
    //      [ A31 A32 A33 ] [ 0 ]          [ A32 ]

    const size_t N = x.size(); // == A.size()
    size_t j2 = N;
    for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
    if (j2 == 0) {
      if (b0) y.Zero();
      return;
    }
    size_t j1 = 0;
    for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) UnitAMultMV1<b0,cx>(A,x,y);
    else {
      if (j1 > 0)
	UnitAMultMV1<b0,cx>(A.SubMatrix(0,j1,j1,j2),x.SubVector(j1,j2),
	    y.SubVector(0,j1));
      TMVAssert(j1 != j2);
      UnitAMultMV1<b0,cx>(A.SubSymMatrix(j1,j2),x.SubVector(j1,j2),
	  y.SubVector(j1,j2));
      if (j2 < N)
	UnitAMultMV1<b0,cx>(A.SubMatrix(j2,N,j1,j2),x.SubVector(j1,j2),
	    y.SubVector(j2,N));
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

    const RealType(T) one(1);
    const RealType(T) zero(0);

    if (x.step() != 1) {
      if (IMAG(alpha) == zero) {
	Vector<Tx> xx = REAL(alpha)*x;
	if (y.step() != 1) {
	  Vector<T> yy(y.size());
	  UnitAMultMV<true,false>(A,xx,yy.View());
	  AddVV(T(1),yy,beta,y);
	}
	else if (beta == zero)
	  UnitAMultMV<true,false>(A,xx,y);
	else  {
	  if (beta != one) y *= beta;
	  UnitAMultMV<false,false>(A,xx,y);
	}
      } else {
	Vector<T> xx = alpha*x;
	if (y.step()!=1) {
	  Vector<T> yy(y.size());
	  UnitAMultMV<true,false>(A,xx,yy.View());
	  AddVV(T(1),yy,beta,y);
	}
	else if (beta == zero)
	    UnitAMultMV<true,false>(A,xx,y);
	else {
	  if (beta != one) y *= beta;
	  UnitAMultMV<false,false>(A,xx,y);
	}
      }
    } else if (y.step()!=1 || alpha!=one) {
      Vector<T> yy(y.size());
      if (x.isconj())
	UnitAMultMV<true,true>(A,x,yy.View());
      else
	UnitAMultMV<true,false>(A,x,yy.View());
      AddVV(alpha,yy,beta,y);
    } else {
      TMVAssert(alpha == T(1));
      TMVAssert(y.step() == 1);
      TMVAssert(x.step() == 1);
      if (beta == zero)
	if (x.isconj())
	  UnitAMultMV<true,true>(A,x,y);
	else
	  UnitAMultMV<true,false>(A,x,y);
      else {
	if (beta != one) y *= beta;
	if (x.isconj())
	  UnitAMultMV<false,true>(A,x,y);
	else
	  UnitAMultMV<false,false>(A,x,y);
      }
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
    TMVAssert(A.uplo() == Lower);
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
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isherm());

    cblas_zhemv(A.isconj()==A.isrm()?CblasColMajor:CblasRowMajor,
	A.isconj()?CblasUpper:CblasLower,
	A.size(),&alpha,A.cptr(),
	A.isrm()?A.stepi():A.stepj(),
	x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
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
    TMVAssert(A.uplo() == Lower);
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
    TMVAssert(A.isherm());
    TMVAssert(!x.SameStorageAs(y));

    cblas_chemv(A.isconj()==A.isrm()?CblasColMajor:CblasRowMajor,
	A.isconj()?CblasUpper:CblasLower,
	A.size(),&alpha,A.cptr(),
	A.isrm()?A.stepi():A.stepj(),
	x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
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
    } else 
#ifdef BLAS
      if (((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0)) && A.isherm()) {
	if (y.step() > 0) { 
	  if (!x.isconj() && x.step() > 0 && !x.SameStorageAs(y))
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
      }
      else
#endif
	NonBlasMultMV(alpha,A,x,beta,y);
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


