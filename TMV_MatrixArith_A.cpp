
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"

namespace tmv {

  //
  // MultMV
  //

  template <ConjItType Ca, class T1, class Ta, class Tx, class T2, ConjItType Cx, class Ty>
    void RowMajorMultMV(const T1 alpha, 
	const GenMatrix<Ta>& A, CVIt<Tx,Unit,Cx> xit,
	const T2 beta, VIt<Ty,Step,NonConj> yit)
    {
      TMVAssert(A.isrm());
      TMVAssert(A.ct()==Ca);
      TMVAssert(alpha != T1(0));

      const Ta* Aptr = A.cptr();
      if (beta == T2(0)) {
	if (alpha == T1(1)) {
	    for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) 
	      DoMultVV(*yit,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	}
	else if (alpha == T1(-1)) {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	      DoMultVV(*yit,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit = -*yit;
	  }
	}
	else {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	      DoMultVV(*yit,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit *= alpha;
	  }
	}
      }
      else if (beta == T2(1)) {
	if (alpha == T1(1)) {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	    Ty temp;
	      DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit += temp;
	  }
	}
	else if (alpha == T1(-1)) {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	    Ty temp;
	      DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit -= temp;
	  }
	}
	else {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	    Ty temp;
	      DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit += alpha * temp;
	  }
	}
      }
      else {
	if (alpha == T1(1)) {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	    Ty temp;
	      DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit = temp + beta * (*yit);
	  }
	}
	else if (alpha == T1(-1)) {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	    Ty temp;
	      DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit = beta * (*yit) - temp;
	  }
	}
	else {
	  for(size_t i=0; i<A.colsize(); ++i,++yit,Aptr+=A.stepi()) {
	    Ty temp;
	      DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,A.rowsize());
	    *yit = alpha * temp + beta * (*yit);
	  }
	}
      }
    }

  template <class T1, class Ta, class Tx, class T2, class Ty> 
    void DoRowMajorMultMV2(
	const T1 alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const T2 beta, const VectorView<Ty>& y)
    {
      if (A.isconj())
	if (x.isconj()) 
	  RowMajorMultMV<Conj>(alpha,A,CVIt<Tx,Unit,Conj>(x.begin()),beta,
	      VIt<Ty,Step,NonConj>(y.begin()));
	else 
	  RowMajorMultMV<Conj>(alpha,A,CVIt<Tx,Unit,NonConj>(x.begin()),beta,
	      VIt<Ty,Step,NonConj>(y.begin()));
      else
	if (x.isconj()) 
	  RowMajorMultMV<NonConj>(alpha,A,CVIt<Tx,Unit,Conj>(x.begin()),beta,
	      VIt<Ty,Step,NonConj>(y.begin()));
	else 
	  RowMajorMultMV<NonConj>(alpha,A,CVIt<Tx,Unit,NonConj>(x.begin()),beta,
	      VIt<Ty,Step,NonConj>(y.begin()));
    }

  template <class T, class Ta, class Tx> void DoRowMajorMultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(y.ct()==NonConj);
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	DoRowMajorMultMV2(REAL(alpha),A,x,REAL(beta),y);
      else
	DoRowMajorMultMV2(REAL(alpha),A,x,beta,y);
    else 
      if (IMAG(beta) == RealType(T)(0))
	DoRowMajorMultMV2(alpha,A,x,REAL(beta),y);
      else
	DoRowMajorMultMV2(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void RowMultMV(const T alpha, 
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    VIt<T,Step,NonConj> yit = y.begin();
    if (beta == T(0)) {
      if (alpha == T(1)) 
	for(size_t i=0; i<A.colsize(); ++i,++yit) *yit = A.row(i) * x;
      else 
	for(size_t i=0; i<A.colsize(); ++i,++yit) 
	  *yit = alpha * (A.row(i) * x);
    } else if (beta == T(1)) {
      if (alpha == T(1)) 
	for(size_t i=0; i<A.colsize(); ++i,++yit) *yit += A.row(i) * x;
      else
	for(size_t i=0; i<A.colsize(); ++i,++yit) 
	  *yit += alpha * (A.row(i) * x);
    } else {
      if (alpha == T(1)) 
	for(size_t i=0; i<A.colsize(); ++i,++yit) 
	  *yit = A.row(i) * x + beta * (*yit);
      else 
	for(size_t i=0; i<A.colsize(); ++i,++yit) 
	  *yit = alpha * (A.row(i) * x) + beta * (*yit);
    }
  }

  template <ConjItType Ca, class T1, class Ta, class Tx, class T2, class Ty>
    void ColMajorMultMV(const T1 alpha, const GenMatrix<Ta>& A,
	const GenVector<Tx>& x, const T2 beta, const VectorView<Ty>& y)
    {
      TMVAssert(A.iscm());
      TMVAssert(y.step()==1);
      TMVAssert(A.colsize() == y.size());
      TMVAssert(alpha != T1(0));
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct()==Ca);
      TMVAssert(y.ct()==NonConj);

      const Ta* Aptr = A.cptr();
      CVIter<Tx> xit = x.begin();
      VIt<Ty,Unit,NonConj> yit = y.begin();

      if (beta == T2(0)) {
	y.Zero();
	for(size_t j=0; j<A.rowsize(); ++j,++xit,Aptr+=A.stepj()) {
	  if (*xit != Tx(0))
	    DoAddVV(*xit,CVIt<Ta,Unit,Ca>(Aptr,1),yit,y.size());
	}
	if (alpha != T1(1)) 
	  DoMultXV(alpha,VIt<Ty,Unit,NonConj>(y.begin()),y.size());
      } else {
	if (alpha == T1(1)) {
	  if (beta != T2(1)) y *= beta;
	  for(size_t j=0; j<A.rowsize(); ++j,++xit,Aptr+=A.stepj()) {
	    if (*xit != Tx(0))
	      DoAddVV(*xit,CVIt<Ta,Unit,Ca>(Aptr,1),yit,y.size());
	  }
	} else if (alpha == T1(-1)) {
	  if (beta != T2(1)) y *= beta;
	  for(size_t j=0; j<A.rowsize(); ++j,++xit,Aptr+=A.stepj()) {
	    if (*xit != Tx(0))
	      DoAddVV(-*xit,CVIt<Ta,Unit,Ca>(Aptr,1),yit,y.size());
	  }
	} else {
	  // Requires temporary
	  Vector<Ty> betay = beta*y;
	  ColMajorMultMV<Ca>(alpha,A,x,RealType(Ty)(0),y);
	  DoAddVV(RealType(Ty)(1),CVIt<Ty,Unit,NonConj>(betay.begin()),
	      VIt<Ty,Unit,NonConj>(y.begin()),y.size());
	}
      }
    }

  template <class T1, class Ta, class Tx, class T2, class Ty> 
    void DoColMajorMultMV2(const T1 alpha, const GenMatrix<Ta>& A,
	const GenVector<Tx>& x, const T2 beta, const VectorView<Ty>& y)
    {
      if (A.isconj()) ColMajorMultMV<Conj>(alpha,A,x,beta,y);
      else ColMajorMultMV<NonConj>(alpha,A,x,beta,y);
    }

  template <class T, class Ta, class Tx> void DoColMajorMultMV(const T alpha, 
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	DoColMajorMultMV2(REAL(alpha),A,x,REAL(beta),y);
      else
	DoColMajorMultMV2(REAL(alpha),A,x,beta,y);
    else
      if (IMAG(beta) == RealType(T)(0))
	DoColMajorMultMV2(alpha,A,x,REAL(beta),y);
      else
	DoColMajorMultMV2(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void ColMultMV(const T alpha, 
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    CVIter<Tx> xit = x.begin();
    if (beta == T(0)) {
      y.Zero();
      for(size_t j=0; j<A.rowsize(); ++j,++xit) y += (*xit) * A.col(j);
      if (alpha != T(1)) y *= alpha;
    } else {
      if (alpha == T(1)) {
	if (beta != T(1)) y *= beta;
	for(size_t j=0; j<A.rowsize(); ++j,++xit) y += (*xit) * A.col(j);
      } else if (alpha == T(-1)) {
	if (beta != T(1)) y *= beta;
	for(size_t j=0; j<A.rowsize(); ++j,++xit) y -= (*xit) * A.col(j);
      } else {
	// Requires temporary
	Vector<T> betay = beta*y;
	ColMultMV(alpha,A,x,T(0),y);
	y += betay;
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
      if (REAL(beta) != double(1)) y.Imag() *= REAL(beta);
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
      if (REAL(beta) != float(1)) y.Imag() *= REAL(beta);
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
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<"  "<<y<<endl;

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) MultXV(beta,y);
      else if (x.SameStorageAs(y)) {
	Vector<Tx> temp = x;
	DoMultMV(alpha,A,temp,beta,y);
      } else {
	DoMultMV(alpha,A,x,beta,y);
      } 
    }
  }

#define InstFile "TMV_MatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


