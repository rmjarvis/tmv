
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  template <class T1, class Ta, class Tx, class T2, class Ty>
    void RowMajorMultMV(const T1 alpha, 
	const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const T2 beta, const VectorView<Ty>& y)
    {
#ifdef XDEBUG
      //cerr<<"RowMajorMultMV\n";
#endif
      TMVAssert(A.isrm());
      TMVAssert(x.step()==1);
      TMVAssert(A.rowsize()==x.size());
      TMVAssert(A.colsize()==y.size());
      TMVAssert(alpha != T1(0));
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);

      VIt<Ty,Step,NonConj> yit = y.begin();
      const Ta* Aptr = A.cptr();
      for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	Ty temp;
	if (A.isconj())
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(x.begin()),A.rowsize());
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(x.begin()),A.rowsize());
	else
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(x.begin()),A.rowsize());
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(x.begin()),A.rowsize());
	if (beta == T2(0)) {
	  if (alpha == T1(1)) *yit = temp;
	  else if (alpha == T1(-1)) *yit = -temp;
	  else *yit = alpha * temp;
	} else {
	  if (beta != T2(1)) *yit *= beta;
	  if (alpha == T1(1)) *yit += temp;
	  else if (alpha == T1(-1)) *yit -= temp;
	  else *yit += alpha * temp;
	}
      }
    }

  template <class T, class Ta, class Tx> void DoRowMajorMultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(y.ct()==NonConj);
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	RowMajorMultMV(REAL(alpha),A,x,REAL(beta),y);
      else
	RowMajorMultMV(REAL(alpha),A,x,beta,y);
    else 
      if (IMAG(beta) == RealType(T)(0))
	RowMajorMultMV(alpha,A,x,REAL(beta),y);
      else
	RowMajorMultMV(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void RowMultMV(const T alpha, 
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
#ifdef XDEBUG
    //cerr<<"RowMultMV\n";
#endif
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    VIt<T,Step,NonConj> yit = y.begin();
    for(size_t i=0; i<A.colsize(); ++i,++yit) {
      T temp = MultVV(A.row(i),x);
      if (beta == T()) {
	if (alpha == T(1)) *yit = temp;
	else *yit = alpha * temp;
      } else {
	if (beta != T(1)) *yit *= beta;
	if (alpha != T(1)) temp *= alpha;
	*yit += temp;
      }
    }
  }

  template <class T1, class Ta, class Tx, class T2, class Ty>
    void ColMajorMultMV(const T1 alpha, const GenMatrix<Ta>& A,
	const GenVector<Tx>& x, const T2 beta, const VectorView<Ty>& y)
    {
#ifdef XDEBUG
      //cerr<<"ColMajorMultMV\n";
#endif
      TMVAssert(A.iscm());
      TMVAssert(y.step()==1);
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(alpha != T1(0));
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);

      const Ta* Aptr = A.cptr();
      CVIter<Tx> xj = x.begin();
      VIt<Ty,Unit,NonConj> y0 = y.begin();

      if (beta == T2(0)) y.Zero();
      else if (beta != T2(1))
	DoMultXV(beta,y0,y.size());
      for(size_t j=0; j<A.rowsize(); ++j,++xj,Aptr+=A.stepj()) {
	if (*xj != Tx(0)) {
	  Ty ax = *xj;
	  if (alpha != T1(1)) ax *= alpha;
	  if (IMAG(ax) == RealType(Ty)(0))
	    if (A.isconj())
	      DoAddVV(REAL(ax),CVIt<Ta,Unit,Conj>(Aptr,1),y0,y.size());
	    else
	      DoAddVV(REAL(ax),CVIt<Ta,Unit,NonConj>(Aptr,1),y0,y.size());
	  else
	    if (A.isconj())
	      DoAddVV(ax,CVIt<Ta,Unit,Conj>(Aptr,1),y0,y.size());
	    else
	      DoAddVV(ax,CVIt<Ta,Unit,NonConj>(Aptr,1),y0,y.size());
	}
      }
    }

  template <class T, class Ta, class Tx> void DoColMajorMultMV(const T alpha, 
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	ColMajorMultMV(REAL(alpha),A,x,REAL(beta),y);
      else
	ColMajorMultMV(REAL(alpha),A,x,beta,y);
    else
      if (IMAG(beta) == RealType(T)(0))
	ColMajorMultMV(alpha,A,x,REAL(beta),y);
      else
	ColMajorMultMV(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void ColMultMV(const T alpha, 
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
#ifdef XDEBUG
    //cerr<<"ColMultMV\n";
#endif
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct()==NonConj);

    CVIter<Tx> xj = x.begin();
    MultXV(beta,y);
    for(size_t j=0; j<A.rowsize(); ++j,++xj) {
      if (*xj != Tx(0)) {
	T ax = *xj;
	if (alpha != T(1)) ax *= alpha;
	AddVV(ax,A.col(j),y);
      }
    }
  }

  template <class T, class Ta, class Tx> void NonBlasMultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    // The main performance issues are to maximize the length
    // of the vectors in the innermost loop and to use
    // unit stride vectors (which are faster to pipeline).
    // It seems that the step size is the more important concern
    // generally, so I check that one first.
    //
    // If col(A) and y are unit step, the column major algorithm is faster.
    // If row(A) and x are unit step, The row major algorithm is faster.
    // Otherwise I look at the lengths of the vectors involved.
    //
    // If the length of x is longer than y, the row version is faster.
    // If the length of y is longer than x, the column version is faster.
    //
    // If they are the same size, I go back to the step sizes, requiring
    // just 1 of the two vectors to be unit step.
    //

    if (A.isrm() && x.step()==1) 
      DoRowMajorMultMV(alpha,A,x,beta,y);
    else if (A.iscm() && y.step()==1)
      DoColMajorMultMV(alpha,A,x,beta,y);
    else if ( (A.isrm() || x.step()==1) && (x.size() > y.size()) )
      RowMultMV(alpha,A,x,beta,y);
    else if ( (A.iscm() || y.step()==1) && (x.size() < y.size()) )
      ColMultMV(alpha,A,x,beta,y);
    else if (A.isrm()) RowMultMV(alpha,A,x,beta,y);
    else if (A.iscm()) ColMultMV(alpha,A,x,beta,y);
    else if (x.step()==1) RowMultMV(alpha,A,x,beta,y);
    else if (y.step()==1) ColMultMV(alpha,A,x,beta,y);
    else RowMultMV(alpha,A,x,beta,y);
  }

#ifdef BLAS
  template <class T, class Ta, class Tx> inline void BlasMultMV(
      const T alpha, const GenMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  { NonBlasMultMV(alpha,A,x,beta,y); }
  template <> inline void BlasMultMV(const double alpha,
      const GenMatrix<double>& A, const GenVector<double>& x,
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
    if (A.isrm())
      cblas_dgemv(CblasRowMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),alpha,A.cptr(), A.stepi(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_dgemv(CblasColMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),alpha,A.cptr(), A.stepj(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenMatrix<complex<double> >& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    if (x.isconj()) NonBlasMultMV(alpha,A,x,beta,y);
    else if (A.isconj())
      if (A.isrm()) 
	cblas_zgemv(CblasColMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),&alpha,A.cptr(), A.stepi(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_zgemv(CblasRowMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),&alpha,A.cptr(), A.stepj(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else
      if (A.isrm())
	cblas_zgemv(CblasRowMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),&alpha,A.cptr(), A.stepi(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_zgemv(CblasColMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),&alpha,A.cptr(), A.stepj(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenMatrix<double>& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      if (x.isconj())
	BlasMultMV(-REAL(alpha),A,x.Conjugate().Imag(),REAL(beta),y.Imag());
      else
	BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenMatrix<double>& A,
      const GenVector<double>& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      MultXV(REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMV( 
      const float alpha,
      const GenMatrix<float>& A, const GenVector<float>& x,
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
    if (A.isrm()) 
      cblas_sgemv(CblasRowMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),alpha,A.cptr(), A.stepi(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_sgemv(CblasColMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),alpha,A.cptr(), A.stepj(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenMatrix<complex<float> >& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    if (x.isconj()) NonBlasMultMV(alpha,A,x,beta,y);
    else if (A.isconj())
      if (A.isrm()) 
	cblas_cgemv(CblasColMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),&alpha,A.cptr(), A.stepi(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_cgemv(CblasRowMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),&alpha,A.cptr(), A.stepj(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else
      if (A.isrm())
	cblas_cgemv(CblasRowMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),&alpha,A.cptr(), A.stepi(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_cgemv(CblasColMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),&alpha,A.cptr(), A.stepj(),
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenMatrix<float>& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      if (x.isconj())
	BlasMultMV(-REAL(alpha),A,x.Conjugate().Imag(),REAL(beta),y.Imag());
      else
	BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenMatrix<float>& A,
      const GenVector<float>& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      MultXV(REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta, class Tx> inline void DoMultMV(
      const T alpha, const GenMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    if (y.isconj()) DoMultMV(CONJ(alpha),A.Conjugate(),x.Conjugate(),
	CONJ(beta),y.Conjugate());
    else 
#ifdef BLAS
      if (((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0)) && 
	  x.step()>0 && y.step()>0 )
	BlasMultMV(alpha,A,x,beta,y);
      else
#endif
	NonBlasMultMV(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
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
    Vector<T> x2 = x;
    x2 *= alpha;
    Matrix<T> A2 = A;
    Vector<T> y3(y.size());
    if (x.size()>0 && y.size()>0) 
      DoMultMV(T(1),A2,x2,T(0),y3.View());
    else y3.Zero();
    Vector<T> y2 = y;
    y2 *= beta;
    y2 += y3;
    Vector<T> y0 = y;
#endif

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) {
	MultXV(beta,y);
      } else if (x.SameStorageAs(y)) {
	Vector<Tx> temp = x;
	DoMultMV(alpha,A,temp,beta,y);
      } else {
	DoMultMV(alpha,A,x,beta,y);
      } 
    }
#ifdef XDEBUG
    if (Norm(y2-y) > 0.0001) {
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

#define InstFile "TMV_MatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


