
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"

//#define XDEBUG

namespace tmv {

  // 
  // Rank1Update
  //

  template <class T, class T1, class Tx, class Ty> void RowMajorRank1Update(
      const T1 alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"RowMajor Rank1\n";
#endif
    TMVAssert(A.isrm());
    TMVAssert(y.step()==1);
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);

    CVIter<Tx> xit = x.begin();
    T* Aptr = A.ptr();

    for (size_t i=0;i<A.colsize();++i,++xit,Aptr+=A.stepi()) {
      if (*xit != Tx(0)) {
	T ax = *xit;
	if (alpha != T1(1)) ax *= alpha;
	// A.row(i) += ax * y;
	if (IMAG(ax) == RealType(T)(0))
	  if (y.isconj())
	    DoAddVV(REAL(ax),CVIt<Ty,Unit,Conj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		y.size());
	  else
	    DoAddVV(REAL(ax),CVIt<Ty,Unit,NonConj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		y.size());
	else
	  if (y.isconj())
	    DoAddVV(ax,CVIt<Ty,Unit,Conj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		y.size());
	  else
	    DoAddVV(ax,CVIt<Ty,Unit,NonConj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		y.size());
      }
    }
  }

  template <class T, class Tx, class Ty> inline void DoRowMajorRank1Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const MatrixView<T>& A)
  {
    if (IMAG(alpha) == RealType(T)(0))
      RowMajorRank1Update(REAL(alpha),x,y,A);
    else
      RowMajorRank1Update(alpha,x,y,A);
  }

  template <class T, class Tx, class Ty> void RowRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Row Rank1\n";
#endif
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);

    CVIter<Tx> xit = x.begin();
    for (size_t i=0;i<A.colsize();++i,++xit) {
      T ax = *xit;
      if (alpha != T(1)) ax *= alpha;
      A.row(i) += ax * y;
    }
  }

  template <class T, class T1, class Tx, class Ty> void ColMajorRank1Update(
      const T1 alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"ColMajor Rank1\n";
#endif
    TMVAssert(A.iscm());
    TMVAssert(x.step()==1);
    TMVAssert(A.colsize() == x.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);

    CVIter<Ty> yit = y.begin();
    T* Aptr = A.ptr();
    for (size_t j=0;j<A.rowsize();++j,++yit,Aptr+=A.stepj()) {
      if (*yit != Ty(0)) {
	T ay = *yit;
	if (alpha != T1(1)) ay *= alpha;
	// A.col(j) += ay * x;
	if (IMAG(ay) == RealType(T)(0)) 
	  if (x.isconj())
	    DoAddVV(REAL(ay),CVIt<Tx,Unit,Conj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		x.size());
	  else
	    DoAddVV(REAL(ay),CVIt<Tx,Unit,NonConj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		x.size());
	else
	  if (x.isconj())
	    DoAddVV(ay,CVIt<Tx,Unit,Conj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		x.size());
	  else
	    DoAddVV(ay,CVIt<Tx,Unit,NonConj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
		x.size());
      }
    }
  }

  template <class T, class Tx, class Ty> inline void DoColMajorRank1Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const MatrixView<T>& A)
  {
    if (IMAG(alpha) == RealType(T)(0))
      ColMajorRank1Update(REAL(alpha),x,y,A);
    else
      ColMajorRank1Update(alpha,x,y,A);
  }

  template <class T, class Tx, class Ty> void ColRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Col Rank1\n";
#endif
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);

    CVIter<Ty> yit = y.begin();
    for (size_t j=0;j<A.rowsize();++j,++yit) {
      T ay = *yit;
      if (alpha != T(1)) ay *= alpha;
      A.col(j) += ay * x;
    }
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
    // The performance issues here are the essentially the same as
    // for MultMV, so see those comments for more about these choices.
    //
    if (A.isrm() && y.step()==1) 
      DoRowMajorRank1Update(alpha,x,y,A);
    else if (A.iscm() && x.step()==1) 
      DoColMajorRank1Update(alpha,x,y,A);
    else if ( (A.isrm() || y.step()==1) && y.size() < x.size() )
      RowRank1Update(alpha,x,y,A);
    else if ( (A.iscm() || x.step()==1) && x.size() < y.size() )
      ColRank1Update(alpha,x,y,A);
    else if (A.isrm()) RowRank1Update(alpha,x,y,A);
    else if (A.iscm()) ColRank1Update(alpha,x,y,A);
    else if (x.step()==1) ColRank1Update(alpha,x,y,A);
    else RowRank1Update(alpha,x,y,A);
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> inline void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  { NonBlasRank1Update(alpha,x,y,A); }
  template <> inline void BlasRank1Update(
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

    if (A.isrm())
      cblas_dger(CblasRowMajor, A.colsize(),A.rowsize(),alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
    else
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
    if (x.isconj())
      if (y.isconj())
	NonBlasRank1Update(alpha,x,y,A);
      else
	if (A.isrm())
	  cblas_zgerc(CblasColMajor,A.rowsize(),A.colsize(),&alpha,
	      y.cptr(),y.step(),x.cptr(),x.step(),A.ptr(),A.stepi()); 
	else
	  cblas_zgerc(CblasRowMajor,A.rowsize(),A.colsize(),&alpha,
	      y.cptr(),y.step(),x.cptr(),x.step(),A.ptr(),A.stepj()); 
    else 
      if (y.isconj())
	if (A.isrm())
	  cblas_zgerc(CblasRowMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
	else
	  cblas_zgerc(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
      else
	if (A.isrm())
	  cblas_zgeru(CblasRowMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
	else
	  cblas_zgeru(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
#ifndef NOFLOAT
  template <> inline void BlasRank1Update(
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

    if (A.isrm())
      cblas_sger(CblasRowMajor, A.colsize(),A.rowsize(),alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
    else
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
    if (x.isconj())
      if (y.isconj())
	NonBlasRank1Update(alpha,x,y,A);
      else
	if (A.isrm())
	  cblas_cgerc(CblasColMajor,A.rowsize(),A.colsize(),&alpha,
	      y.cptr(),y.step(),x.cptr(),x.step(),A.ptr(),A.stepi()); 
	else
	  cblas_cgerc(CblasRowMajor,A.rowsize(),A.colsize(),&alpha,
	      y.cptr(),y.step(),x.cptr(),x.step(),A.ptr(),A.stepj()); 
    else 
      if (y.isconj())
	if (A.isrm())
	  cblas_cgerc(CblasRowMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
	else
	  cblas_cgerc(CblasColMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
      else
	if (A.isrm())
	  cblas_cgeru(CblasRowMajor,A.colsize(),A.rowsize(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
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
#ifdef BLAS
      else if ((A.isrm() || A.iscm()) && x.step()>0 && y.step()>0)
	BlasRank1Update(alpha,x,y,A);
#endif
      else NonBlasRank1Update(alpha,x,y,A);
    }
  
#ifdef XDEBUG
    if (Norm(A2-A) > 0.001*Norm(A)) {
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


