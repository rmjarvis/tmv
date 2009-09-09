
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  template <bool b0, bool cx, bool ca, bool rm,  class T, class Ta, class Tx>
    void RowMultMV(const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      TMVAssert(A.rowsize()==x.size());
      TMVAssert(A.colsize()==y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());
      TMVAssert(ca == A.isconj());
      TMVAssert(rm == A.isrm());

      const int si = A.stepi();
      const int sj = (rm ? 1 : A.stepj());
      const int ds = A.diagstep();
      const size_t M = A.colsize();
      const size_t N = A.rowsize();

      const Ta* Aij1 = A.cptr();
      const Tx* xj1 = x.cptr();
      T* yi = y.ptr();

      size_t k=A.nlo();
      size_t i=0;
      size_t j1=0;
      size_t j2=A.nhi()+1;
      size_t len = j2; // len = j2-j1
      for(; i<M; ++i, ++yi) {
	if (b0) *yi = T(0);

	// *yi += A.row(i,j1,j2) * x.SubVector(j1,j2);
	const Ta* Aij = Aij1;
	const Tx* xj = xj1;
	for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	  *yi += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);

	if (k>0) { --k; ++len; Aij1+=si; }
	else { ++j1; ++xj1; Aij1+=ds; }
	if (j2<N) ++j2;
	else { --len; if (j1==N) { ++i, ++yi; break; } }
      }
      if (b0) for(;i<M; ++i, ++yi) *yi = T(0);
    }

  template <bool b0, bool cx, bool ca, bool cm, class T, class Ta, class Tx>
    void ColMultMV(const GenBandMatrix<Ta>& A,
	const GenVector<Tx>& x, const VectorView<T> y)
    {
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());
      TMVAssert(ca == A.isconj());
      TMVAssert(cm == A.iscm());

      const size_t N = A.rowsize();
      const size_t M = A.colsize();
      const int si = (cm ? 1 : A.stepi());
      const int sj = A.stepj();
      const int ds = A.diagstep();

      const Ta* Ai1j = A.cptr();
      const Tx* xj = x.cptr();
      T* yi1 = y.ptr();

      size_t k=A.nhi();
      size_t i1=0;
      size_t i2=A.nlo()+1;
      size_t len = i2; // = i2-i1

      if (b0) y.Zero();

      for(size_t j=N; j>0; --j,++xj) {
	if (*xj != Tx(0)) {
	  // y.SubVector(i1,i2) += *xj * A.col(j,i1,i2);
	  const Ta* Aij = Ai1j;
	  T* yi = yi1;
	  for(size_t i=len;i>0;--i,++yi,(cm?++Aij:Aij+=si))
	    *yi += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);
	}
	if (k>0) { --k; Ai1j+=sj; ++len; }
	else { ++i1; ++yi1; Ai1j+=ds; }
	if (i2<M) ++i2; 
	else { --len; if (i1==M) break; }
      }
    }

  template <bool b0, bool cx, bool ca, bool dm, class T, class Ta, class Tx>
    void DiagMultMV(
	const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());
      TMVAssert(ca == A.isconj());
      TMVAssert(dm == A.isdm());

      const int si = A.stepi();
      const int sj = A.stepj();
      const int ds = A.diagstep();
      const int lo = A.nlo();
      const int hi = A.nhi();
      const size_t M = A.colsize();
      const size_t N = A.rowsize();

      const Ta* Ai1j1 = A.cptr() + lo*si;
      const Tx* xj1 = x.cptr();
      T* yi1 = y.ptr() + lo;

      size_t j2=min(int(M)-lo,int(N));
      size_t len=j2; // == j2-j1

      if (b0) y.Zero();

      for(int k=-A.nlo(); k<=hi; ++k) {
	// y.SubVector(i1,i2) += DiagMatrixViewOf(A.diag(k)) * 
	//     x.SubVector(j1,j2);
	const Ta* Aij = Ai1j1;
	const Tx* xj = xj1;
	T* yi = yi1;
	for(size_t i=len;i>0;--i,++yi,++xj,(dm?++Aij:Aij+=ds))
	  *yi += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);
	if (k<0) { --yi1; ++len; Ai1j1-=si; } 
	else { ++xj1; Ai1j1+=sj; }
	if (j2 < N) ++j2; else --len; 
      }
    }

  template <bool b0, bool cx, class T, class Ta, class Tx> 
    inline void UnitAMultMV1(
	const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      //cerr<<"UnitAMultMV: A = "<<A<<endl;
      //cerr<<"x = "<<x<<endl;
      //cerr<<"y = "<<y<<endl;
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(!x.SameStorageAs(y));
      TMVAssert(cx == x.isconj());

      if (A.isrm())
	if (A.isconj())
	  RowMultMV<b0,cx,true,true>(A,x,y);
	else
	  RowMultMV<b0,cx,false,true>(A,x,y);
      else if (A.iscm())
	if (A.isconj())
	  ColMultMV<b0,cx,true,true>(A,x,y);
	else
	  ColMultMV<b0,cx,false,true>(A,x,y);
      else if (A.isdm())
	if (A.isconj())
	  DiagMultMV<b0,cx,true,true>(A,x,y);
	else
	  DiagMultMV<b0,cx,false,true>(A,x,y);
      else
	if (A.isconj())
	  DiagMultMV<b0,cx,true,false>(A,x,y);
	else
	  DiagMultMV<b0,cx,false,false>(A,x,y);
    }

  template <bool b0, bool cx, class T, class Ta, class Tx> void UnitAMultMV(
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    // Check for 0's in beginning or end of x:
    // y += [ A1 A2 A3 ] [ 0 ]  -->  y += A2 x
    //                   [ x ]
    //                   [ 0 ]

    const size_t N = x.size(); // = A.rowsize()
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
      const int hi = A.nhi();
      const int lo = A.nlo();
      const size_t M = y.size(); // = A.colsize()
      // This next bit is copied from the BandMatrix Cols function
      size_t i1 = int(j1) > hi ? j1-hi : 0;
      size_t i2 = min(j2+lo,M);
      int newhi = int(j1) < hi ? hi-int(j1) : 0;
      int newlo = lo+hi-newhi;
      int newM = i2-i1;
      int newN = j2-j1;
      TMVAssert(newM > 0);
      TMVAssert(newN > 0);
      if (newhi >= int(newN)) newhi = newN-1;
      if (newlo >= int(newM)) newlo = newM-1;
      TMVAssert(A.OKSubBandMatrix(i1,i2,j1,j2,newlo,newhi,1,1));
      const Ta* p = A.cptr()+i1*A.stepi()+j1*A.stepj();
      ConstBandMatrixView<Ta> Acols(p,newM,newN,newlo,newhi,
	  A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
      UnitAMultMV1<b0,cx>(Acols,x.SubVector(j1,j2),y.SubVector(i1,i2));
      if (b0) {
	y.SubVector(0,i1).Zero();
        y.SubVector(i2,M).Zero();
      }
    }
  }

  template <class T, class Ta, class Tx> void NonBlasMultMV(
      const T alpha, const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);

#ifdef XDEBUG
    Vector<T> y2 = beta*y + alpha*Matrix<T>(A)*x;
    Vector<T> y0 = y;
#endif

    const RealType(T) one(1);
    const RealType(T) zero(0);

    if (x.step() != 1 || x.SameStorageAs(y)) {
      if (IMAG(alpha) == zero) {
	Vector<Tx> xx = REAL(alpha) * x;
	if (y.step() != 1) {
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
      } else {
	Vector<T> xx = alpha * x;
	if (y.step() != 1) {
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
    } else if (y.step() != 1 || alpha != one) {
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
      TMVAssert(!x.SameStorageAs(y));
      if (x.isconj())
	if (beta == zero)
	  UnitAMultMV<true,true>(A,x,y);
	else {
	  if (beta != one) y *= beta;
	  UnitAMultMV<false,true>(A,x,y);
	}
      else
	if (beta == zero)
	  UnitAMultMV<true,false>(A,x,y);
	else {
	  if (beta != one) y *= beta;
	  UnitAMultMV<false,false>(A,x,y);
	}
    }

#ifdef XDEBUG
    if (Norm(y2-y) > 0.001*Norm(y)) {
      cerr<<"NonBlas Band MultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"--> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

#ifdef BLAS
  template <class T, class Ta, class Tx> void DoBlasMultMV(
      const T alpha, const GenBandMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  { NonBlasMultMV(alpha,A,x,beta,y); }
  template <> void DoBlasMultMV(const double alpha,
      const GenBandMatrix<double>& A, const GenVector<double>& x,
      const double beta, const VectorView<double>& y)
  {
    TMVAssert(alpha != 0.);
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));
    if (A.isrm())
      cblas_dgbmv(CblasRowMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nlo(), A.stepi()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_dgbmv(CblasColMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nhi(), A.stepj()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> void DoBlasMultMV(const complex<double> alpha,
      const GenBandMatrix<complex<double> >& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(alpha != complex<double>(0));
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));
    if (A.isconj())
      if (A.isrm()) 
	cblas_zgbmv(CblasColMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_zgbmv(CblasRowMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else
      if (A.isrm())
	cblas_zgbmv(CblasRowMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_zgbmv(CblasColMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
  }
#ifndef NOFLOAT
  template <> void DoBlasMultMV(const float alpha,
      const GenBandMatrix<float>& A, const GenVector<float>& x,
      const float beta, const VectorView<float>& y)
  {
    TMVAssert(alpha != 0.);
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));
    if (A.isrm())
      cblas_sgbmv(CblasRowMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nlo(), A.stepi()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_sgbmv(CblasColMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nhi(), A.stepj()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> void DoBlasMultMV(const complex<float> alpha,
      const GenBandMatrix<complex<float> >& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(alpha != complex<float>(0));
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!x.SameStorageAs(y));
    if (A.isconj())
      if (A.isrm()) 
	cblas_cgbmv(CblasColMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_cgbmv(CblasRowMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else
      if (A.isrm())
	cblas_cgbmv(CblasRowMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_cgbmv(CblasColMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta, class Tx> void BlasMultMV(const T alpha,
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    if (A.rowsize() > A.colsize()) {
      size_t N1 = A.colsize() + A.nhi();
      if (A.rowsize() > N1) {
	DoBlasMultMV(alpha,A.SubBandMatrix(0,A.colsize(),0,N1,
	      A.nlo(),A.nhi()),x.SubVector(0,N1),beta,y);
      }
      else DoBlasMultMV(alpha,A,x,beta,y);
    }
    else if (A.colsize() > A.rowsize()) {
      size_t N2 = A.rowsize() + A.nlo();
      if (A.colsize() > N2) {
	DoBlasMultMV(alpha,A.SubBandMatrix(0,N2,0,A.rowsize(),
	      A.nlo(),A.nhi()),x,beta,y.SubVector(0,N2));
	y.SubVector(N2,y.size()) *= beta;
      }
      else DoBlasMultMV(alpha,A,x,beta,y);
    }
    else DoBlasMultMV(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void DoMultMV(const T alpha,
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    if (y.isconj()) DoMultMV(CONJ(alpha),A.Conjugate(),
	x.Conjugate(),CONJ(beta),y.Conjugate());
    else {
#ifdef BLAS
      if (((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0))) {
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
  }

  //
  // MultEqMV
  //

  template <bool rm, bool ca, class T, class Ta> void DoRowUpperMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMultEqMV Upper\n";
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step()==1);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(rm == A.isrm());
    TMVAssert(ca == A.isconj());

    const size_t N = x.size();
    const int sj = (rm ? 1 : A.stepj());
    const int ds = A.diagstep();

    T* xi = x.ptr();
    const Ta* Aii = A.cptr();
    size_t j2=A.nhi()+1;
    size_t len = j2-1;

    for(; len>0; ++xi,Aii+=ds) {
      // i = 0..N-2
      // x(i) = A.row(i,i,j2) * x.SubVector(i,j2);
      *xi *= (ca ? CONJ(*Aii) : *Aii);
      const T* xj = xi+1;
      const Ta* Aij = Aii+sj;
      for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	*xi += (*xj) * (ca ? CONJ(*Aij) : *Aij);

      if (j2<N) ++j2;
      else --len;
    }
    *xi *= (ca ? CONJ(*Aii) : *Aii);
  }

  template <bool rm, class T, class Ta> inline void RowUpperMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      DoRowUpperMultEqMV<rm,true>(A,x);
    else
      DoRowUpperMultEqMV<rm,false>(A,x);
  }

  template <bool cm, bool ca, class T, class Ta> void DoColUpperMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"ColMultEqMV Upper\n";
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step()==1);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(cm == A.iscm());
    TMVAssert(ca == A.isconj());

    const size_t N = x.size();
    const int si = cm ? 1 : A.stepi();
    const int sj = A.stepj();
    const int ds = A.diagstep();

    const Ta* Ai1j = A.cptr();
    T* xi1 = x.ptr();

    *xi1 *= (ca ? CONJ(*Ai1j) : *Ai1j);
    Ai1j += sj;
    const T* xj = x.cptr()+1;

    size_t k=A.nhi()-1;
    size_t len = 1;
    for(size_t j=1; j<N; ++j,++xj) {
      if (*xj != T(0)) {
	// j = 1..N-1
	// x.SubVector(i1,j) += x(j) * A.col(j,i1,j);
	const Ta* Aij = Ai1j;
	T* xi = xi1;
	for(size_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si))
	  *xi += *xj * (ca ? CONJ(*Aij) : *Aij);
	// Now Aij == Ajj, xi == xj
	// so this next statement is really *xj *= *Ajj
	*xi *= (ca ? CONJ(*Aij) : *Aij);
      }

      if (k>0) { --k; Ai1j+=sj; ++len; }
      else { ++xi1; Ai1j+=ds; }
    }
  }

  template <bool cm, class T, class Ta> inline void ColUpperMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      DoColUpperMultEqMV<cm,true>(A,x);
    else
      DoColUpperMultEqMV<cm,false>(A,x);
  }

  template <bool rm, bool ca, class T, class Ta> void DoRowLowerMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMultEqMV Lower\n";
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step()==1);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(rm == A.isrm());
    TMVAssert(ca == A.isconj());

    const size_t N = x.size();
    const int si = A.stepi();
    const int sj = (rm ? 1 : A.stepj());
    const int ds = A.diagstep();

    size_t j1 = N-1-A.nlo();
    const T* xj1 = x.cptr() + j1;
    T* xi = x.ptr() + N-1;
    const Ta* Aii = A.cptr() + (N-1)*ds;
    const Ta* Aij1 = Aii - A.nlo()*sj;
    size_t len = A.nlo();

    for(; len>0; --xi,Aii-=ds) {
      // i = N-1..1
      // x(i) = A.row(i,j1,i+1) * x.SubVector(j1,i+1);
      *xi *= (ca ? CONJ(*Aii) : *Aii);
      const Ta* Aij = Aij1;
      const T* xj = xj1;
      for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	*xi += *xj * (ca ? CONJ(*Aij) : *Aij);

      if (j1>0) { --j1; Aij1-=ds; --xj1; }
      else { --len; Aij1-=si; }
    }
    *xi *= (ca ? CONJ(*Aii) : *Aii);
  }

  template <bool rm, class T, class Ta> inline void RowLowerMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      DoRowLowerMultEqMV<rm,true>(A,x);
    else
      DoRowLowerMultEqMV<rm,false>(A,x);
  }

  template <bool cm, bool ca, class T, class Ta> void DoColLowerMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"ColMultEqMV Lower\n";
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(cm == A.iscm());
    TMVAssert(ca == A.isconj());

    const size_t N = x.size();
    const int si = cm ? 1 : A.stepi();
    const int ds = A.diagstep();

    T* xj = x.ptr() + N-1;
    const Ta* Ajj = A.cptr()+(N-1)*ds;

    *xj *= (ca ? CONJ(*Ajj) : *Ajj);
    --xj;
    Ajj -= ds;

    int k=A.nlo()-1;
    for(size_t j=N-1,len=1;j>0;--j,--xj,Ajj-=ds) {
      if (*xj!=T(0)) {
	// Actual j = N-2..0
	// x.SubVector(j+1,N) += *xj * A.col(j,j+1,N);
	T* xi = xj+1;
	const Ta* Aij = Ajj+si;
	for (size_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si))
	  *xi += *xj * (ca ? CONJ(*Aij) : *Aij);
	*xj *= (ca ? CONJ(*Ajj) : *Ajj);

      }
      if (k>0) { --k; ++len; }
    }
  }

  template <bool cm, class T, class Ta> inline void ColLowerMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      DoColLowerMultEqMV<cm,true>(A,x);
    else
      DoColLowerMultEqMV<cm,false>(A,x);
  }

  template <class T, class Ta> inline void DoUpperMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
    // x = A * x
  {
    if (A.isrm()) RowUpperMultEqMV<true>(A,x);
    else if (A.iscm()) ColUpperMultEqMV<true>(A,x);
    else RowUpperMultEqMV<false>(A,x);
  }

  template <class T, class Ta> inline void DoLowerMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
    if (A.isrm()) RowLowerMultEqMV<true>(A,x);
    else if (A.iscm()) ColLowerMultEqMV<true>(A,x);
    else RowLowerMultEqMV<false>(A,x);
  }

  template <class T, class Ta> void NonBlasUpperMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(x.ct() == NonConj);

    //     [ A11 A12  0  ] [ 0  ]   [ A12 x2 ]
    // x = [  0  A22 A23 ] [ x2 ] = [ A22 x2 ]
    //     [  0   0  A33 ] [ 0  ]   [   0    ]

    const size_t N = x.size(); // = A.size()
    size_t j2 = N;
    for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
    if (j2 == 0) return;
    size_t j1 = 0;
    for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) DoUpperMultEqMV(A,x);
    else {
      TMVAssert(j1 < j2);
      const Ta* p22 = A.cptr() + j1*A.diagstep();
      const size_t N22 = j2-j1;
      VectorView<T> x2 = x.SubVector(j1,j2);
      if (N22 > size_t(A.nhi())) {
	ConstBandMatrixView<Ta> A22(p22,N22,N22,0,A.nhi(),
	    A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	if (j1 > 0) {
	  const size_t jx = j1+A.nhi();
	  if (j1 < size_t(A.nhi())) {
	    const Ta* p12 = A.cptr() + j1*A.stepj();
	    ConstBandMatrixView<Ta> A12(p12,j1,A.nhi(),j1-1,A.nhi()-j1,
		A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	    UnitAMultMV1<true,false>(A12,x.SubVector(j1,jx),x.SubVector(0,j1));
	  } else {
	    const Ta* p12 = p22 - A.nhi()*A.stepi();
	    ConstBandMatrixView<Ta> A12(p12,A.nhi(),A.nhi(),A.nhi()-1,0,
		A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	    VectorView<T> x1x = x.SubVector(j1-A.nhi(),j1);
	    x1x = x.SubVector(j1,jx);
	    DoLowerMultEqMV(A12,x1x);
	  }
	}
	DoUpperMultEqMV(A22,x2);
      } else {
	ConstBandMatrixView<Ta> A22(p22,N22,N22,0,N22-1,
	    A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	if (j1 > 0) {
	  const size_t M12 = (j1 < size_t(A.nhi())) ? j1 : A.nhi();
	  const Ta* p12 = p22 - M12*A.stepi();
	  size_t newhi = A.nhi()-M12;
	  if (newhi >= N22) newhi = N22-1;
	  ConstBandMatrixView<Ta> A12(p12,M12,N22,M12-1,newhi,
	      A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	  UnitAMultMV1<true,false>(A12,x2,x.SubVector(j1-M12,j1));
	} 
	DoUpperMultEqMV(A22,x2);
      }
    }
  }

  template <class T, class Ta> void NonBlasLowerMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
    // x = A * x
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(x.ct() == NonConj);

    //     [ A11  0   0  ] [ 0  ]   [   0    ]
    // x = [ A21 A22  0  ] [ x2 ] = [ A22 x2 ]
    //     [  0  A32 A33 ] [ 0  ]   [ A32 x2 ]

    const size_t N = x.size(); // = A.size()
    size_t j2 = N;
    for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
    if (j2 == 0) return;
    size_t j1 = 0;
    for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) DoLowerMultEqMV(A,x);
    else {
      TMVAssert(j1 < j2);
      const Ta* p22 = A.cptr() + j1*A.diagstep();
      const size_t N22 = j2-j1;
      VectorView<T> x2 = x.SubVector(j1,j2);
      if (N22 > size_t(A.nlo())) {
	ConstBandMatrixView<Ta> A22(p22,N22,N22,A.nlo(),0,
	    A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	if (j2 < N) {
	  const size_t jx = j2-A.nlo();
	  const Ta* p32 = A.cptr() + j2*A.diagstep() - A.nlo()*A.stepj();
	  if (j2+A.nlo() > N) {
	    ConstBandMatrixView<Ta> A32(p32,N-j2,A.nlo(),0,A.nlo()-1,
		A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	    UnitAMultMV1<true,false>(A32,x.SubVector(jx,j2),x.SubVector(j2,N));
	  } else {
	    ConstBandMatrixView<Ta> A32(p32,A.nlo(),A.nlo(),0,A.nlo()-1,
		A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	    VectorView<T> x3x = x.SubVector(j2,j2+A.nlo());
	    x3x = x.SubVector(jx,j2);
	    DoUpperMultEqMV(A32,x3x);
	  }
	}
	DoLowerMultEqMV(A22,x2);
      } else {
	ConstBandMatrixView<Ta> A22(p22,N22,N22,N22-1,0,
	    A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	if (j2 < N) {
	  const Ta* p32 = p22 + N22*A.stepi();
	  const size_t M32 = (j2+A.nlo() > N) ? N-j2 : A.nlo();
	  size_t newlo = A.nlo()-N22;
	  if (newlo >= M32) newlo = M32-1;
	  ConstBandMatrixView<Ta> A32(p32,M32,N22,newlo,N22-1,
	      A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct());
	  UnitAMultMV1<true,false>(A32,x2,x.SubVector(j2,j2+M32));
	} 
	DoLowerMultEqMV(A22,x2);
      }
    }
  }

#ifdef BLAS
  template <class T, class Ta> void BlasMultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.nlo() == 0) NonBlasUpperMultEqMV(A,x);
    else NonBlasLowerMultEqMV(A,x);
  }
  template <> void BlasMultEqMV( 
      const GenBandMatrix<double>& A, const VectorView<double>& x)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() == 1);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    cblas_dtbmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	A.nlo() == 0 ? CblasUpper : CblasLower, CblasNoTrans, CblasNonUnit,
	A.colsize(), A.nlo() == 0 ? A.nhi() : A.nlo(), 
	A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(), 
	A.diagstep(), x.ptr(), x.step()); 
  }
  template <> void BlasMultEqMV(const GenBandMatrix<complex<double> >& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() == 1);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    if (A.isconj()) 
      cblas_ztbmv(A.isrm() ? CblasColMajor : CblasRowMajor,
	  A.nlo() == 0 ? CblasLower : CblasUpper, CblasConjTrans, CblasNonUnit,
	  A.colsize(), A.nlo() == 0 ? A.nhi() : A.nlo(), 
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(), 
	  A.diagstep(), x.ptr(), x.step()); 
    else
      cblas_ztbmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	  A.nlo() == 0 ? CblasUpper : CblasLower, CblasNoTrans, CblasNonUnit,
	  A.colsize(), A.nlo() == 0 ? A.nhi() : A.nlo(), 
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(), 
	  A.diagstep(), x.ptr(), x.step()); 
  }
#ifndef NOFLOAT
  template <> void BlasMultEqMV( 
      const GenBandMatrix<float>& A, const VectorView<float>& x)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() == 1);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    cblas_stbmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	A.nlo() == 0 ? CblasUpper : CblasLower, CblasNoTrans, CblasNonUnit,
	A.colsize(), A.nlo() == 0 ? A.nhi() : A.nlo(), 
	A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(), 
	A.diagstep(), x.ptr(), x.step()); 
  }
  template <> void BlasMultEqMV(const GenBandMatrix<complex<float> >& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() == 1);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);
    if (A.isconj()) 
      cblas_ctbmv(A.isrm() ? CblasColMajor : CblasRowMajor,
	  A.nlo() == 0 ? CblasLower : CblasUpper, CblasConjTrans, CblasNonUnit,
	  A.colsize(), A.nlo() == 0 ? A.nhi() : A.nlo(), 
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(), 
	  A.diagstep(), x.ptr(), x.step()); 
    else
      cblas_ctbmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	  A.nlo() == 0 ? CblasUpper : CblasLower, CblasNoTrans, CblasNonUnit,
	  A.colsize(), A.nlo() == 0 ? A.nhi() : A.nlo(), 
	  A.isrm() ? A.cptr()-A.nlo() : A.cptr()-A.nhi(), 
	  A.diagstep(), x.ptr(), x.step()); 
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta> void MultEqMV(
      const GenBandMatrix<Ta>& A, const VectorView<T>& x)
  {
#ifdef XDEBUG
    Vector<T> x2 = Matrix<T>(A) * x;
    Vector<T> x0 = x;
#endif
    TMVAssert(A.IsSquare());
    TMVAssert(A.colsize() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(A.nlo() == 0 || A.nhi() == 0);

    if (x.isconj()) MultEqMV(A.Conjugate(),x.Conjugate());
    else 
#ifdef BLAS
      if ( (A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0) )
	BlasMultEqMV(A,x);
      else
#endif
	if (A.nlo() == 0) NonBlasUpperMultEqMV(A,x);
	else NonBlasLowerMultEqMV(A,x);
#ifdef XDEBUG
    if (Norm(x-x2) > 0.001*Norm(x)) {
      cerr<<"MultEqMV: \n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tx> void MultMV(const T alpha,
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
#ifdef XDEBUG
    //cerr<<"Start Band: MultMV\n";
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<"  step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<"  step "<<y.step()<<"  "<<y<<endl;
    Vector<T> y2 = beta*y+alpha*Matrix<T>(A)*x;
    Vector<T> y0 = y;
#endif

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) {
	y *= beta;
      } else if (A.rowsize() > A.colsize()+A.nhi()) {
	MultMV(alpha,A.Cols(0,A.colsize()+A.nhi()),
	    x.SubVector(0,A.colsize()+A.nhi()),beta,y);
      } else if (A.colsize() > A.rowsize()+A.nlo()) {
	MultMV(alpha,A.Rows(0,A.rowsize()+A.nlo()),
	    x,beta,y.SubVector(0,A.rowsize()+A.nlo()));
	y.SubVector(A.rowsize()+A.nlo(),A.colsize()) *= beta;
      } else if (A.IsSquare() && (A.nlo() == 0 || A.nhi() == 0)) {
	if (A.nlo() == 0 && A.nhi() == 0)
	  MultMV(alpha,DiagMatrixViewOf(A.diag()),x,beta,y);
	else if (beta == T(0) && y.step() == 1) {
	  y = alpha * x;
	  MultEqMV(A,y);
	} else {
	  Vector<T> xx = alpha*x;
	  MultEqMV(A,xx.View());
	  AddVV(T(1),xx,beta,y);
	}
      } else DoMultMV(alpha,A,x,beta,y);
    }
#ifdef XDEBUG
    //cerr<<"->y = "<<y<<endl;
    if (Norm(y2-y) > 0.001*Norm(y)) {
      cerr<<"Band MultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"--> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


