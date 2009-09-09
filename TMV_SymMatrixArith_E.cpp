
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool ua, bool cx, bool ha, bool rm, class T, class Tx, class Ta> 
    void RowRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
#ifdef XDEBUG
      //cerr<<"Row Rank1\n";
#endif
      TMVAssert(A.size() == x.size());
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(rm == A.isrm());
      TMVAssert(ua == (alpha == Ta(1)));
      TMVAssert(x.step()==1);
      TMVAssert(rm == A.isrm());
      TMVAssert(cx == x.isconj());
      TMVAssert(ha == A.isherm());

      const Tx* xi = x.cptr();
      const Tx*const x0 = x.cptr();
      T* Arowi = A.ptr();
      const size_t si = A.stepi();
      const size_t sj = A.stepj();
      const size_t N = A.size();

      for (size_t i=0;i<N;++i,++xi,Arowi+=si) if (*xi != Tx(0)) {
	T* Aij = Arowi;
	const Tx* xj = x0;

	if (ua || IsReal(alpha)) {
	  Tx axi = cx ? CONJ(*xi) : *xi;
	  if (!ua) axi *= REAL(alpha);
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	    *Aij += axi * (ha ? (cx?*xj:CONJ(*xj)) : (cx?CONJ(*xj):*xj));
	} else {
	  T axi = alpha * (cx ? CONJ(*xi) : *xi);
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	    *Aij += axi * (ha ? (cx?*xj:CONJ(*xj)) : (cx?CONJ(*xj):*xj));
	}
      }
    }

  template <bool ua, bool cx, bool ha, bool cm, class T, class Tx, class Ta> 
    void ColRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
#ifdef XDEBUG
      //cerr<<"Col Rank1\n";
#endif
      TMVAssert(x.step()==1);
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(ua == (alpha == Ta(1)));
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(cm == A.iscm());
      TMVAssert(ha == A.isherm());

      const Tx* xj = x.cptr();
      T* Ajj = A.ptr();
      const size_t si = cm ? 1 : A.stepi();
      const size_t ds = A.stepj()+si;
      const size_t N = A.size();

      for (size_t j=N;j>0;--j,++xj,Ajj+=ds) if (*xj!=Tx(0)) {
	// A.col(j,j,N) += *xj * x.SubVector(j,N);
	T* Aij = Ajj;
	const Tx* xi = xj;
	if (ua || IsReal(alpha)) {
	  Tx axj = (ha!=cx) ? CONJ(*xj) : *xj;
	  if (!ua) axj *= REAL(alpha);
	  for(size_t i=j;i>0;--i,++xi,(cm?++Aij:Aij+=si))
	    *Aij += axj * (cx ? CONJ(*xi) : *xi);
	} else {
	  T axj = alpha * ((ha!=cx) ? CONJ(*xj) : *xj);
	  for(size_t i=j;i>0;--i,++xi,(cm?++Aij:Aij+=si))
	    *Aij += axj * (cx ? CONJ(*xi) : *xi);
	}
      }
    }

  template <bool ua, bool cx, class T, class Ta, class Tx>
    inline void DoRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(ua == (alpha == Ta(1)));
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());

      if (A.isherm()) 
	if (A.isrm()) RowRank1Update<ua,cx,true,true>(alpha,x,A);
	else if (A.iscm()) ColRank1Update<ua,cx,true,true>(alpha,x,A);
	else RowRank1Update<ua,cx,true,false>(alpha,x,A);
      else
	if (A.isrm()) RowRank1Update<ua,cx,false,true>(alpha,x,A);
	else if (A.iscm()) ColRank1Update<ua,cx,false,true>(alpha,x,A);
	else RowRank1Update<ua,cx,false,false>(alpha,x,A);
    }

  template <class T, class Tx> void NonBlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.step() != 1) {
      Vector<T> xx = x;
      if (IMAG(alpha) == RealType(T)(0)) 
	if (REAL(alpha) == RealType(T)(1)) 
	  DoRank1Update<true,false>(REAL(alpha),xx,A);
	else
	  DoRank1Update<false,false>(REAL(alpha),xx,A);
      else 
	DoRank1Update<false,false>(alpha,xx,A);
    } 
    else 
      if (x.isconj())
	if (IMAG(alpha) == RealType(T)(0)) 
	  if (REAL(alpha) == RealType(T)(1)) 
	    DoRank1Update<true,true>(REAL(alpha),x,A);
	  else
	    DoRank1Update<false,true>(REAL(alpha),x,A);
	else 
	  DoRank1Update<false,true>(alpha,x,A);
      else
	if (IMAG(alpha) == RealType(T)(0)) 
	  if (REAL(alpha) == RealType(T)(1)) 
	    DoRank1Update<true,false>(REAL(alpha),x,A);
	  else
	    DoRank1Update<false,false>(REAL(alpha),x,A);
	else 
	  DoRank1Update<false,false>(alpha,x,A);
  }

#ifdef BLAS
  template <class T, class Tx> void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  { NonBlasRank1Update(alpha,x,A); }
  template <> void BlasRank1Update(
      const double alpha, const GenVector<double>& x,
      const SymMatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 double\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(x.step() > 0);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.isherm());

    if (A.isrm())
      cblas_dsyr(CblasRowMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepi()); 
    else
      cblas_dsyr(CblasColMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank1Update(
      const complex<double> alpha, const GenVector<complex<double> >& x, 
      const SymMatrixView<complex<double> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 c double\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(IMAG(alpha)==double(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(x.step() > 0);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.isherm());

    if (A.isrm()) 
      cblas_zher(CblasRowMajor,CblasLower,A.size(),REAL(alpha),
	  x.cptr(),x.step(),A.ptr(),A.stepi()); 
    else 
      cblas_zher(CblasColMajor,CblasLower,A.size(),REAL(alpha),
	  x.cptr(),x.step(),A.ptr(),A.stepj()); 
  }
#ifndef NOFLOAT
  template <> void BlasRank1Update(
      const float alpha, const GenVector<float>& x,
      const SymMatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 float\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(x.step() > 0);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.isherm());

    if (A.isrm())
      cblas_ssyr(CblasRowMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepi()); 
    else
      cblas_ssyr(CblasColMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank1Update(
      const complex<float> alpha, const GenVector<complex<float> >& x, 
      const SymMatrixView<complex<float> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 c float\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != float(0));
    TMVAssert(IMAG(alpha)==float(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(x.step() > 0);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.isherm());

    TMVAssert(IMAG(alpha)==float(0));
    if (A.isrm()) 
      cblas_cher(CblasRowMajor,CblasLower,A.size(),REAL(alpha),
	  x.cptr(),x.step(),A.ptr(),A.stepi()); 
    else 
      cblas_cher(CblasColMajor,CblasLower,A.size(),REAL(alpha),
	  x.cptr(),x.step(),A.ptr(),A.stepj()); 
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx> void Rank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
    // A = A + alpha * x * xT
  {
#ifdef XDEBUG
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += (alpha*x^x.Conjugate());
    else 
      A2 += (alpha*x^x);
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank1Update(CONJ(alpha),x.Conjugate(),A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.issym()) Rank1Update(alpha,x,A.Transpose());
	else Rank1Update(alpha,x,A.Adjoint());
      }
#ifdef BLAS
      else if ( ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) &&
	  A.isherm()) {
	if (!x.isconj() && x.step() > 0)
	  BlasRank1Update(alpha,x,A);
	else {
	  Vector<Tx> xx = x;
	  BlasRank1Update(alpha,xx,A);
	}
      }
#endif
      else NonBlasRank1Update(alpha,x,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*max(RealType(T)(1),Norm(A))) {
      cerr<<"Rank1Update: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


