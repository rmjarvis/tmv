
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  template <class T1, class Ta, class Tx, class T2, class Ty>
    void RowMajorMultMV(const T1 alpha, 
	const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
	const T2 beta, const VectorView<Ty>& y)
    {
#ifdef XDEBUG
      //cerr<<"RowMajorMultMV\n";
#endif
      TMVAssert(A.isrm());
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(A.size()==x.size());
      TMVAssert(A.size()==y.size());
      TMVAssert(alpha != T1(0));
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);

      const size_t N = A.size();
      const Tx* xj = x.cptr();
      VIt<Ty,Unit,NonConj> yj = y.begin();
      const Ta* Aptr = A.cptr();
      for(size_t j=0; j<N; ++j,++xj,++yj,Aptr+=A.stepi()) {
	Ty temp;
	// temp = A.row(j,0,j+1) * x.SubVector(0,j+1);
	if (A.isconj())
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(x.begin()),j+1);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(x.begin()),j+1);
	else
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(x.begin()),j+1);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(x.begin()),j+1);
	if (alpha != T1(1)) temp *= alpha;
	if (beta == T2(0)) *yj = temp;
	else {
	  if (beta != T2(1)) *yj *= beta;
	  *yj += temp;
	}
	if (j>0 && *xj != Tx(0)) {
	  Ty ax = x.isconj() ? CONJ(*xj) : *xj;
	  if (alpha != T1(1)) ax *= alpha;
	  // y.SubVector(0,j) += ax * A.col(j,0,j);
	  if (A.isconj() != A.isherm())
	    DoAddVV(ax,CVIt<Ta,Unit,Conj>(Aptr,1),
		VIt<Ty,Unit,NonConj>(y.begin()),j);
	  else
	    DoAddVV(ax,CVIt<Ta,Unit,NonConj>(Aptr,1),
		VIt<Ty,Unit,NonConj>(y.begin()),j);
	}
      }
    }

  template <class T, class Ta, class Tx> void DoRowMajorMultMV(
      const T alpha, const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
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
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
#ifdef XDEBUG
    //cerr<<"RowMultMV\n";
#endif
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct()==NonConj);

    const size_t N = A.size();
    CVIter<Tx> xj = x.begin();
    VIt<T,Step,NonConj> yj = y.begin();
    for(size_t j=0; j<N; ++j,++xj,++yj) {
      T temp = A.row(j,0,j+1) * x.SubVector(0,j+1);
      if (alpha != T(1)) temp *= alpha;
      if (beta == T(0)) *yj = temp;
      else {
	if (beta != T(1)) *yj *= beta;
	*yj += temp;
      }
      if (*xj != Tx(0)) {
	T ax = *xj;
	if (alpha != T(1)) ax *= alpha;
	y.SubVector(0,j) += ax * A.col(j,0,j);
      }
    }
  }

  template <class T1, class Ta, class Tx, class T2, class Ty>
    void ColMajorMultMV(const T1 alpha, const GenSymMatrix<Ta>& A,
	const GenVector<Tx>& x, const T2 beta, const VectorView<Ty>& y)
    {
#ifdef XDEBUG
      //cerr<<"ColMajorMultMV\n";
#endif
      TMVAssert(A.iscm());
      TMVAssert(A.uplo() == Lower);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(alpha != T1(0));
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);

      const size_t N = A.size();
      const size_t ds = A.stepj()+1;
      const Tx* xj = x.cptr()+N-1;
      VIt<Ty,Unit,NonConj> yj = y.begin()+N-1;
      const Ta* Aptr = A.cptr()+(N-1)*ds;
      size_t len = 1;
      for(int j=N-1; j>=0; --j,--xj,--yj,Aptr-=ds,++len) {
	Ty temp;
	// temp = A.row(j,j,N) * x.SubVector(j,N);
	if (A.isconj() != A.isherm())
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(xj,1),len);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(xj,1),len);
	else
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(xj,1),len);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(xj,1),len);
	if (alpha != T1(1)) temp *= alpha;
	if (beta == T2(0)) *yj = temp;
	else {
	  if (beta != T2(1)) *yj *= beta;
	  *yj += temp;
	}
	if (len>1 && *xj != Tx(0)) {
	  Ty ax = x.isconj() ? CONJ(*xj) : *xj;
	  if (alpha != T1(1)) ax *= alpha;
	  // y.SubVector(j+1,N) += ax * A.col(j,j+1,N);
	  if (A.isconj())
	    DoAddVV(ax,CVIt<Ta,Unit,Conj>(Aptr+1,1),
		VIt<Ty,Unit,NonConj>(yj+1),len-1);
	  else
	    DoAddVV(ax,CVIt<Ta,Unit,NonConj>(Aptr+1,1),
		VIt<Ty,Unit,NonConj>(yj+1),len-1);
	}
      }
    }

  template <class T, class Ta, class Tx> void DoColMajorMultMV(const T alpha, 
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
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
      const GenSymMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
#ifdef XDEBUG
    //cerr<<"ColMultMV\n";
#endif
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct()==NonConj);

    const size_t N = A.size();
    CVIter<Tx> xj = x.begin()+N-1;
    VIt<T,Step,NonConj> yj = y.begin()+N-1;
    for(int j=N-1; j>=0; --j,--xj,--yj) {
      T temp = A.row(j,j,N) * x.SubVector(j,N);
      if (alpha != T(1)) temp *= alpha;
      if (beta == T(0)) *yj = temp;
      else {
	if (beta != T(1)) *yj *= beta;
	*yj += temp;
      }
      if (*xj != Tx(0)) {
	T ax = *xj;
	if (alpha != T(1)) ax *= alpha;
	y.SubVector(j+1,N) += ax * A.col(j,j+1,N);
      }
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
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.step() == 1 && y.step() == 1) {
      if (A.isrm()) DoRowMajorMultMV(alpha,A,x,beta,y);
      else if (A.iscm()) DoColMajorMultMV(alpha,A,x,beta,y);
      else RowMultMV(alpha,A,x,beta,y);
    } else {
      if (A.isrm()) RowMultMV(alpha,A,x,beta,y);
      else if (A.iscm()) ColMultMV(alpha,A,x,beta,y);
      else RowMultMV(alpha,A,x,beta,y);
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
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    if (A.isrm())
      cblas_dsymv(CblasRowMajor,CblasLower,
	  A.size(),alpha,A.cptr(),A.stepi(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_dsymv(CblasColMajor,CblasLower,
	  A.size(),alpha,A.cptr(),A.stepj(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<complex<double> >& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.isconj()) NonBlasMultMV(alpha,A,x,beta,y);
    else if (A.isherm())
      if (A.isconj())
	if (A.isrm()) 
	  cblas_zhemv(CblasColMajor,CblasUpper,
	      A.size(),&alpha,A.cptr(),A.stepi(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
	else
	  cblas_zhemv(CblasRowMajor,CblasUpper,
	      A.size(),&alpha,A.cptr(),A.stepj(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	if (A.isrm())
	  cblas_zhemv(CblasRowMajor,CblasLower,
	      A.size(),&alpha,A.cptr(),A.stepi(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
	else
	  cblas_zhemv(CblasColMajor,CblasLower,
	      A.size(),&alpha,A.cptr(),A.stepj(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<double>& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      if (x.isconj())
	BlasMultMV(-REAL(alpha),A,x.Conjugate().Imag(),REAL(beta),y.Imag());
      else
	BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenSymMatrix<double>& A,
      const GenVector<double>& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      y.Imag() *= REAL(beta);
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMV(const float alpha,
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
    
    if (A.isrm())
      cblas_ssymv(CblasRowMajor,CblasLower,
	  A.size(),alpha,A.cptr(),A.stepi(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_ssymv(CblasColMajor,CblasLower,
	  A.size(),alpha,A.cptr(),A.stepj(),
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
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
    TMVAssert(A.uplo() == Lower);

    if (x.isconj()) NonBlasMultMV(alpha,A,x,beta,y);
    else if (A.isherm())
      if (A.isconj())
	if (A.isrm()) 
	  cblas_chemv(CblasColMajor,CblasUpper,
	      A.size(),&alpha,A.cptr(),A.stepi(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
	else
	  cblas_chemv(CblasRowMajor,CblasUpper,
	      A.size(),&alpha,A.cptr(),A.stepj(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	if (A.isrm())
	  cblas_chemv(CblasRowMajor,CblasLower,
	      A.size(),&alpha,A.cptr(),A.stepi(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
	else
	  cblas_chemv(CblasColMajor,CblasLower,
	      A.size(),&alpha,A.cptr(),A.stepj(),
	      x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<float>& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      if (x.isconj())
	BlasMultMV(-REAL(alpha),A,x.Conjugate().Imag(),REAL(beta),y.Imag());
      else
	BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenSymMatrix<float>& A,
      const GenVector<float>& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      y.Imag() *= REAL(beta);
    } else NonBlasMultMV(alpha,A,x,beta,y);
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
      if (((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0)) && 
	  x.step()>0 && y.step()>0 )
	BlasMultMV(alpha,A,x,beta,y);
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
	Vector<Tx> temp = x;
	DoMultMV(alpha,A,temp,beta,y);
      } else {
	DoMultMV(alpha,A,x,beta,y);
      } 
    }
#ifdef XDEBUG
    if (Norm(y2-y) > 0.001*Norm(y)) {
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


