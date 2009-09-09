
#include "TMV.h"

//#define XDEBUG

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool cx, bool cy, bool y1, bool cm, class T, class Tx, class Ty> 
    void ColRank1Update(
	const GenVector<Tx>& x, const GenVector<Ty>& y,
	const MatrixView<T>& A)
    {
#ifdef XDEBUG
      //cerr<<"Row Rank1\n";
#endif
      TMVAssert(A.colsize() == x.size());
      TMVAssert(A.rowsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(!A.isrm());
      TMVAssert(x.step() == 1);
      TMVAssert(cm == A.iscm());
      TMVAssert(y1 == (y.step() == 1));
      TMVAssert(cx == x.isconj());
      TMVAssert(cy == y.isconj());

      const Ty* yj = y.cptr();
      const Tx*const xptr = x.cptr();
      T* Acolj = A.ptr();
      const int sj = A.stepj();
      const int si = (cm ? 1 : A.stepi());
      const int ys = y.step();
      const size_t M = A.colsize();
      const size_t N = A.rowsize();

      for (size_t j=N; j>0; --j,(y1?++yj:yj+=ys),Acolj+=sj) if (*yj!=Ty(0)) {
	T* Aij = Acolj;
	const Tx* xi = xptr;
	for (size_t i=M; i>0; --i,++xi,(cm?++Aij:Aij+=si)) 
	  *Aij += (cx ? CONJ(*xi) : *xi) * (cy ? CONJ(*yj) : *yj);
      }
    }

  template <bool cx, bool y1, class T, class Tx, class Ty> 
    inline void UnitARank1Update(
	const GenVector<Tx>& x,
	const GenVector<Ty>& y, const MatrixView<T>& A)
    {
      TMVAssert(A.colsize() == x.size());
      TMVAssert(A.rowsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(!A.isrm());
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(y1 == (y.step() == 1));

      if (A.iscm()) 
	if (y.isconj())
	  ColRank1Update<cx,true,y1,true>(x,y,A);
	else
	  ColRank1Update<cx,false,y1,true>(x,y,A);
      else
	if (y.isconj())
	  ColRank1Update<cx,true,y1,false>(x,y,A);
	else
	  ColRank1Update<cx,false,y1,false>(x,y,A);
    }

  template <class T, class Tx, class Ty> void NonBlasRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(!A.isrm());

    if (x.step() != 1 || alpha != T(1)) {
      if (x.step() == 1 && y.size() < x.size()) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Vector<Ty> yy = REAL(alpha)*y;
	  if (x.isconj())
	    UnitARank1Update<true,true>(x,yy,A);
	  else
	    UnitARank1Update<false,true>(x,yy,A);
	} else {
	  Vector<T> yy = alpha*y;
	  if (x.isconj())
	    UnitARank1Update<true,true>(x,yy,A);
	  else
	    UnitARank1Update<false,true>(x,yy,A);
	}
      } else {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Vector<Tx> xx = REAL(alpha)*x;
	  if (y.step() == 1)
	    UnitARank1Update<false,true>(xx,y,A);
	  else
	    UnitARank1Update<false,false>(xx,y,A);
	} else {
	  Vector<T> xx = alpha*x;
	  if (y.step() == 1)
	    UnitARank1Update<false,true>(xx,y,A);
	  else
	    UnitARank1Update<false,false>(xx,y,A);
	}
      }
    } else {
      if (x.isconj())
	if (y.step() == 1)
	  UnitARank1Update<true,true>(x,y,A);
	else
	  UnitARank1Update<true,false>(x,y,A);
      else
	if (y.step() == 1)
	  UnitARank1Update<false,true>(x,y,A);
	else
	  UnitARank1Update<false,false>(x,y,A);
    }
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  { NonBlasRank1Update(alpha,x,y,A); }
  template <> void BlasRank1Update(
      const double alpha, const GenVector<double>& x,
      const GenVector<double>& y, const MatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 double\n";
#endif
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    cblas_dger(CblasColMajor, A.colsize(),A.rowsize(),alpha,
	x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank1Update(
      const complex<double> alpha, const GenVector<complex<double> >& x, 
      const GenVector<complex<double> >& y,
      const MatrixView<complex<double> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 c double\n";
#endif
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
    TMVAssert(A.iscm());

    if (x.isconj())
      cblas_zgerc(CblasRowMajor,A.rowsize(),A.colsize(),&alpha,
	  y.cptr(),y.step(),x.cptr(),x.step(),A.ptr(),A.stepj()); 
    else if (y.isconj())
      cblas_zgerc(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
    else
      cblas_zgeru(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
#ifndef NOFLOAT
  template <> void BlasRank1Update(
      const float alpha, const GenVector<float>& x,
      const GenVector<float>& y, const MatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 float\n";
#endif
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    cblas_sger(CblasColMajor, A.colsize(),A.rowsize(),alpha,
	x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank1Update(
      const complex<float> alpha, const GenVector<complex<float> >& x, 
      const GenVector<complex<float> >& y,
      const MatrixView<complex<float> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 c float\n";
#endif
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
    TMVAssert(A.iscm());

    if (x.isconj())
      cblas_cgerc(CblasRowMajor,A.rowsize(),A.colsize(),&alpha,
	  y.cptr(),y.step(),x.cptr(),x.step(),A.ptr(),A.stepj()); 
    else if (y.isconj())
      cblas_cgerc(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
    else
      cblas_cgeru(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx, class Ty> void Rank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
    // A = A + alpha * x * yT
  {
#ifdef XDEBUG
    Matrix<T> A2 = A;
    Vector<T> x2 = x;
    Vector<T> y2 = y;
    for(size_t i=0;i<x.size();i++) for(size_t j=0;j<y.size();j++) 
      A2(i,j) += alpha*x2(i)*y2(j);
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    if (alpha != T(0) && A.colsize() > 0 && A.rowsize() > 0) {
      if (A.isconj()) 
	Rank1Update(CONJ(alpha),x.Conjugate(),y.Conjugate(),A.Conjugate());
      else if (A.isrm())
	Rank1Update(alpha,y,x,A.Transpose());
#ifdef BLAS
      else if ((A.iscm() && A.stepj()>0)) {
	if (x.step() < 0) {
	  if (y.step() < 0) {
	    if (x.size() <= y.size()) {
	      Vector<Tx> xx = alpha*x;
	      Vector<Ty> yy = y;
	      BlasRank1Update(T(1),xx,yy,A);
	    } else {
	      Vector<Tx> xx = x;
	      Vector<Ty> yy = alpha*y;
	      BlasRank1Update(T(1),xx,yy,A);
	    }
	  } else {
	    Vector<Tx> xx = alpha*x;
	    BlasRank1Update(T(1),xx,y,A);
	  }
	} else {
	  if (y.step() < 0) {
	    Vector<Ty> yy = alpha*y;
	    BlasRank1Update(T(1),x,yy,A);
	  } else if (x.isconj() && y.isconj()) {
	    if (x.size() <= y.size()) {
	      Vector<Tx> xx = alpha*x;
	      BlasRank1Update(T(1),xx,y,A);
	    } else {
	      Vector<Ty> yy = alpha*y;
	      BlasRank1Update(T(1),x,yy,A);
	    }
	  } else {
	    BlasRank1Update(alpha,x,y,A);
	  }
	}
      }
#endif
      else NonBlasRank1Update(alpha,x,y,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*max(RealType(T)(1),Norm(A))) {
      cerr<<"Rank1Update: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<"  step = "<<y.step()<<"  "<<y<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


