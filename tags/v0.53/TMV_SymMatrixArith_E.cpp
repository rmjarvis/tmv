
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool a1, bool cx, bool ha, bool rm, bool add, class T, class Tx, class Ta> 
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
      TMVAssert(a1 == (alpha == Ta(1)));
      TMVAssert(x.step()==1);
      TMVAssert(rm == A.isrm());
      TMVAssert(cx == x.isconj());
      TMVAssert(ha == A.isherm());

      const size_t si = A.stepi();
      const size_t sj = A.stepj();
      const size_t N = A.size();
      const Tx* xi = x.cptr();
      const Tx*const x0 = x.cptr();

      T A00;
      if (*xi == Tx(0)) A00 = T(0);
      else if (ha) A00 = alpha * NORM(*xi);
      else if (cx) A00 = alpha * CONJ(*xi * *xi);
      else A00 = alpha * (*xi * *xi);
      ++xi;
      T* Arowi = A.ptr()+si;

      for (size_t i=1;i<N;++i,++xi,Arowi+=si) if (*xi != Tx(0)) {
	// A.row(i,0,i+1) += ax * x.SubVector(0,i+1);
	T* Aij = Arowi;
	const Tx* xj = x0;

	if (a1 || IsReal(alpha)) {
	  Tx axi = cx ? CONJ(*xi) : *xi;
	  if (!a1) axi *= REAL(alpha);
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
	    const T temp = axi * (ha==cx ? *xj : CONJ(*xj));
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	} else {
	  T axi = alpha * (cx ? CONJ(*xi) : *xi);
	  for(size_t j=i+1;j>0;--j,++xj,(rm?++Aij:Aij+=sj)) {
	    const T temp = axi * (ha==cx ? *xj : CONJ(*xj));
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	}
      }
      if (add) *A.ptr() += A00; 
      else *A.ptr() = A00;
      // Do this at end in case A.ptr = x.ptr, so don't mess up x
    }

  template <bool a1, bool cx, bool ha, bool cm, bool add, class T, class Tx, class Ta> 
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
      TMVAssert(a1 == (alpha == Ta(1)));
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());
      TMVAssert(cm == A.iscm());
      TMVAssert(ha == A.isherm());

      const size_t si = cm ? 1 : A.stepi();
      const size_t ds = A.stepj()+si;
      const size_t N = A.size();
      const Tx* xj = x.cptr()+N-1;
      T* Ajj = A.ptr()+(N-1)*ds;

      for (size_t jj=N,Nmj=1;jj>0;--jj,++Nmj,--xj,Ajj-=ds) if (*xj!=Tx(0)) {
	// Nmj = N-j
	// A.col(j,j,N) += *xj * x.SubVector(j,N);
	T* Aij = Ajj;
	const Tx* xi = xj;
	if (a1 || IsReal(alpha)) {
	  Tx axj = (ha!=cx) ? CONJ(*xj) : *xj;
	  if (!a1) axj *= REAL(alpha);
	  for(size_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
	    const T temp = axj * (cx ? CONJ(*xi) : *xi);
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	} else {
	  T axj = alpha * ((ha!=cx) ? CONJ(*xj) : *xj);
	  for(size_t i=Nmj;i>0;--i,++xi,(cm?++Aij:Aij+=si)) {
	    const T temp = axj * (cx ? CONJ(*xi) : *xi);
	    if (add) *Aij += temp;
	    else *Aij = temp;
	  }
	}
      }
    }

  template <bool a1, bool cx, bool add, class T, class Ta, class Tx>
    inline void DoRank1Update(const Ta alpha,
	const GenVector<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.size());
      TMVAssert(x.size() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(a1 == (alpha == Ta(1)));
      TMVAssert(x.step() == 1);
      TMVAssert(cx == x.isconj());

      if (A.isherm()) 
	if (A.isrm()) RowRank1Update<a1,cx,true,true,add>(alpha,x,A);
	else if (A.iscm()) ColRank1Update<a1,cx,true,true,add>(alpha,x,A);
	else RowRank1Update<a1,cx,true,false,add>(alpha,x,A);
      else
	if (A.isrm()) RowRank1Update<a1,cx,false,true,add>(alpha,x,A);
	else if (A.iscm()) ColRank1Update<a1,cx,false,true,add>(alpha,x,A);
	else RowRank1Update<a1,cx,false,false,add>(alpha,x,A);
    }

  template <bool add, class T, class Tx> void NonBlasRank1Update(
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
	  DoRank1Update<true,false,add>(REAL(alpha),xx,A);
	else
	  DoRank1Update<false,false,add>(REAL(alpha),xx,A);
      else 
	DoRank1Update<false,false,add>(alpha,xx,A);
    } 
    else 
      if (x.isconj())
	if (IMAG(alpha) == RealType(T)(0)) 
	  if (REAL(alpha) == RealType(T)(1)) 
	    DoRank1Update<true,true,add>(REAL(alpha),x,A);
	  else
	    DoRank1Update<false,true,add>(REAL(alpha),x,A);
	else 
	  DoRank1Update<false,true,add>(alpha,x,A);
      else
	if (IMAG(alpha) == RealType(T)(0)) 
	  if (REAL(alpha) == RealType(T)(1)) 
	    DoRank1Update<true,false,add>(REAL(alpha),x,A);
	  else
	    DoRank1Update<false,false,add>(REAL(alpha),x,A);
	else 
	  DoRank1Update<false,false,add>(alpha,x,A);
  }

#ifdef BLAS
  template <class T, class Tx> void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  { NonBlasRank1Update<true>(alpha,x,A); }
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

    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==double(0));
      if (A.isrm())
	cblas_zher(CblasRowMajor,CblasLower,A.size(),REAL(alpha),
	    x.cptr(),x.step(),A.ptr(),A.stepi());
      else
	cblas_zher(CblasColMajor,CblasLower,A.size(),REAL(alpha),
	    x.cptr(),x.step(),A.ptr(),A.stepj());
    } else {
#ifdef LAP
      char uplo = A.iscm() ? 'L' : 'U';
      int n = A.size();
      int incx = x.step();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      zsyr(&uplo,&n,LAP_Complex(&alpha),LAP_Complex(x.cptr()),&incx,
	  LAP_Complex(A.ptr()),&lda);
#else
      NonBlasRank1Update<true>(alpha,x,A);
#endif
    }
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

    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==double(0));
      if (A.isrm())
	cblas_cher(CblasRowMajor,CblasLower,A.size(),REAL(alpha),
	    x.cptr(),x.step(),A.ptr(),A.stepi());
      else
	cblas_cher(CblasColMajor,CblasLower,A.size(),REAL(alpha),
	    x.cptr(),x.step(),A.ptr(),A.stepj());
    } else {
#ifdef LAP
      char uplo = A.iscm() ? 'L' : 'U';
      int n = A.size();
      int incx = x.step();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      csyr(&uplo,&n,LAP_Complex(&alpha),LAP_Complex(x.cptr()),&incx,
	  LAP_Complex(A.ptr()),&lda);
#else
      NonBlasRank1Update<true>(alpha,x,A);
#endif
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx> void Rank1Update(
      const T alpha, const GenVector<Tx>& x, const int beta,
      const SymMatrixView<T>& A)
    // A = beta*A + alpha * x * xT
  {
#ifdef XDEBUG
    //cerr<<"Start Rank1Update:\n";
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
    //cerr<<"alpha,beta = "<<alpha<<"  "<<beta<<endl;
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += (alpha*x^x.Conjugate());
    else 
      A2 += (alpha*x^x);
    //cerr<<"A2 = "<<A2<<endl;
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(beta == 0 || beta == 1);
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank1Update(CONJ(alpha),x.Conjugate(),beta,A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.issym()) Rank1Update(alpha,x,beta,A.Transpose());
	else Rank1Update(alpha,x,beta,A.Adjoint());
      }
#ifdef BLAS
      else if ( ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
	if (x.isconj() || x.step() < 0 || x.cptr() ==  (Tx*)(A.ptr())) {
	  Vector<Tx> xx = x;
	  if (beta == 0) A.Zero();
	  BlasRank1Update(alpha,xx,A);
	}
	else {
	  if (beta == 0) A.Zero();
	  BlasRank1Update(alpha,x,A);
	}
      }
#endif
      else if (beta == 0)
	NonBlasRank1Update<false>(alpha,x,A);
      else
	NonBlasRank1Update<true>(alpha,x,A);
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


