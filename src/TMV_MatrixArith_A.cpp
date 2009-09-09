///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_Matrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_MatrixArith_A.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_VIt.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> const T* MatrixComposite<T>::cptr() const
  {
    if (!itsm.get()) {
      size_t len = this->colsize()*this->rowsize();
      itsm.reset(new T[len]);
      MatrixView<T>(itsm.get(),this->colsize(),this->rowsize(),
	  stepi(),stepj(),this->stor(),
	  NonConj,len FIRSTLAST1(itsm.get(),itsm.get()+len) ) = *this;
    }
    return itsm.get();
  }

  template <class T> int MatrixComposite<T>::stepi() const 
  { return this->isrm() ? this->rowsize() : 1; }

  template <class T> int MatrixComposite<T>::stepj() const 
  { return this->isrm() ? 1 : this->colsize(); }

  template <class T> size_t MatrixComposite<T>::ls() const 
  { return this->rowsize() * this->colsize(); }

  // 
  //
  // MultMV
  //

  // These routines are designed to work even if y has the same storage
  // as either x or the first row/column of A.
  
  // Most of the optimizations here are copied from ATLAS.
  // While the particular values of 4,8,32 used herein may not be 
  // optimal for a particular machine/data type, they are likely to be 
  // better than not blocking at all for almost all machines.
  // MJ: Haven't done the optimizing here yet.
  
  template <bool add, bool cx, bool ca, bool rm, class T, class Ta, class Tx>
    inline void RowMultMV(
	const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);
      TMVAssert(x.step() == 1);
      TMVAssert(y.step() == 1);
      TMVAssert(!SameStorage(x,y));
      TMVAssert(cx == x.isconj());
      TMVAssert(ca == A.isconj());
      TMVAssert(rm == A.isrm());

      const size_t M = A.colsize();
      const size_t N = A.rowsize();
      const int si = A.stepi();
      const int sj = (rm ? 1 : A.stepj());

      const Ta* Ai0 = A.cptr();
      const Tx*const x0 = x.cptr();
      T* yi = y.ptr();

      for(size_t i=M; i>0; --i,++yi,Ai0+=si) {
	// *yi += A.row(i) * x

	const Ta* Aij = Ai0;
	const Tx* xj = x0;
	register T temp(0);
	for(size_t j=N; j>0; --j,++xj,(rm?++Aij:Aij+=sj))
	  temp += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);

	if (add) *yi += temp;
	else *yi = temp;
      }
    }

  template <bool add, bool cx, bool ca, bool cm, class T, class Ta, class Tx> 
    inline void ColMultMV(
	const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct()==NonConj);
      TMVAssert(x.step() == 1);
      TMVAssert(y.step() == 1);
      TMVAssert(!SameStorage(x,y));
      TMVAssert(cx == x.isconj());
      TMVAssert(ca == A.isconj());
      TMVAssert(cm == A.iscm());

      //Vector<T> y2 = y;
      //Vector<T> y3 = Matrix<T,RowMajor>(A)*x;
      //if (add) y3 += y;

      const size_t M = A.colsize();
      size_t N = A.rowsize();
      const int si = (cm ? 1 : A.stepi());
      const int sj = A.stepj();

      const Ta* A0j = A.cptr();
      const Tx* xj = x.cptr();
      T*const y0 = y.ptr();

      if (!add) {
	//y2 = x(0) * A.col(0);
	if (*xj == Tx(0)) {
	  y.Zero();
	} else {
	  const Ta* Aij = A0j;
	  T* yi = y0;
	  const Tx xjval = (cx ? CONJ(*xj) : *xj);
	  for(size_t i=M; i>0; --i,++yi,(cm?++Aij:Aij+=si))
	    *yi = xjval * (ca ? CONJ(*Aij) : *Aij);
	}
	++xj; A0j+=sj; --N;
      }

      for(; N>0; --N,++xj,A0j+=sj) {
	// y += *xj * A.col(j)
	if (*xj != Tx(0)) {
	  const Ta* Aij = A0j;
	  T* yi = y0;
	  const Tx xjval = (cx ? CONJ(*xj) : *xj);
	  for(size_t i=M; i>0; --i,++yi,(cm?++Aij:Aij+=si))
	    *yi += xjval * (ca ? CONJ(*Aij) : *Aij);
	}
      }
    }

  template <bool add, bool cx, class T, class Ta, class Tx> 
    void UnitAMultMV1(
	const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      TMVAssert(A.rowsize() == x.size());
      TMVAssert(A.colsize() == y.size());
      TMVAssert(x.size() > 0);
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(x.step() == 1);
      TMVAssert(y.step() == 1);
      TMVAssert(!SameStorage(x,y));
      TMVAssert(cx == x.isconj());

      if (A.isrm()) 
	if (A.isconj())
	  RowMultMV<add,cx,true,true>(A,x,y);
	else
	  RowMultMV<add,cx,false,true>(A,x,y);
      else if (A.iscm())
	if (A.isconj())
	  ColMultMV<add,cx,true,true>(A,x,y);
	else
	  ColMultMV<add,cx,false,true>(A,x,y);
      else if ( A.rowsize() >= A.colsize() )
	if (A.isconj())
	  RowMultMV<add,cx,true,false>(A,x,y);
	else
	  RowMultMV<add,cx,false,false>(A,x,y);
      else 
	if (A.isconj())
	  ColMultMV<add,cx,true,false>(A,x,y);
	else
	  ColMultMV<add,cx,false,false>(A,x,y);
    }

  template <bool add, bool cx, class T, class Ta, class Tx> 
    inline void UnitAMultMV(
	const GenMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
#ifdef XDEBUG
      //cerr<<"Start UnitAMultMV: \n";
      //cerr<<"add = "<<add<<endl;
      //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
      Vector<Tx> x0 = x;
      Vector<T> y0 = y;
      Matrix<Ta> A0 = A;
      Vector<T> y2 = y;
      for(size_t i=0;i<y.size();i++) {
	if (add)
	  y2(i) += (A.row(i) * x0);
	else
	  y2(i) = (A.row(i) * x0);
      }
      //cerr<<"y2 = "<<y2<<endl;
#endif
      // Check for 0's in beginning or end of x:
      // y += [ A1 A2 A3 ] [ 0 ]  -->  y += A2 x
      //                   [ x ]
      //                   [ 0 ]

      const size_t N = x.size(); // = A.rowsize()
      size_t j2 = N;
      for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
      if (j2 == 0) {
	if (!add) y.Zero();
	return;
      }
      size_t j1 = 0;
      for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
      TMVAssert(j1 !=j2);
      if (j1 == 0 && j2 == N) UnitAMultMV1<add,cx>(A,x,y);
      else UnitAMultMV1<add,cx>(A.Cols(j1,j2),x.SubVector(j1,j2),y);

#ifdef XDEBUG
      //cerr<<"y => "<<y<<endl;
      if (Norm(y-y2) > 0.001*(Norm(A0)*Norm(x0)+
	    (add?Norm(y0):RealType(T)(0)))) {
	cerr<<"MultMV: \n";
	cerr<<"add = "<<add<<endl;
	cerr<<"A = "<<Type(A);
	if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
	else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
	cerr<<endl<<"x = "<<Type(x)<<" step "<<x.step();
	if (x.size() < 30) cerr<<"  "<<x0;
	cerr<<endl<<"y = "<<Type(y)<<" step "<<y.step();
	if (y.size() < 30) cerr<<"  "<<y0;
	cerr<<endl<<"Aptr = "<<A.cptr();
	cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
	if (y.size() < 200) {
	  cerr<<"--> y = "<<y<<endl;
	  cerr<<"y2 = "<<y2<<endl;
	} else {
	  size_t imax;
	  MaxAbsElement(y-y2,&imax);
	  cerr<<"y("<<imax<<") = "<<y(imax)<<endl;
	  cerr<<"y2("<<imax<<") = "<<y2(imax)<<endl;
	}
	cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
	cerr<<"Norm(x0) = "<<Norm(x0)<<endl;
	cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
	cerr<<"|A0|*|x0|+?|y0| = "<<
	  Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0))<<endl;
	cerr<<"Norm(y-y2) = "<<Norm(y-y2)<<endl;
	cerr<<"NormInf(y-y2) = "<<NormInf(y-y2)<<endl;
	cerr<<"Norm1(y-y2) = "<<Norm1(y-y2)<<endl;
	abort();
      }
#endif
    }

  template <bool add, class T, class Ta, class Tx> inline void NonBlasMultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y (+)= alpha * A * x
  {
#ifdef XDEBUG
    //cerr<<"Start MultMV: alpha = "<<alpha<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
    Vector<Tx> x0 = x;
    Vector<T> y0 = y;
    Matrix<Ta> A0 = A;
    Vector<T> y2 = y;
    for(size_t i=0;i<y.size();i++) {
      if (add)
	y2(i) += alpha * (A.row(i) * x0);
      else
	y2(i) = alpha * (A.row(i) * x0);
    }
    //cerr<<"y2 = "<<y2<<endl;
#endif
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);

    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    if (x.step() != 1 || SameStorage(x,y) ||
	(alpha != RealType(T)(1) && y.step() == 1 && M/4 >= N)) {
      // This last check is taken from the ATLAS version of this code.
      // Apparently M = 4N is the dividing line between applying alpha
      // here versus at the end when adding Ax to y
      if (IMAG(alpha) == RealType(T)(0)) {
	Vector<Tx> xx = REAL(alpha)*x;
	if (y.step()!=1) {
	  Vector<T> yy(y.size());
	  UnitAMultMV<false,false>(A,xx,yy.View());
	  if (add) y += yy;
	  else y = yy;
	} 
	else 
	  UnitAMultMV<add,false>(A,xx,y);
      } else {
	Vector<T> xx = alpha*x;
	if (y.step() != 1) {
	  Vector<T> yy(y.size());
	  UnitAMultMV<false,false>(A,xx,yy.View());
	  if (add) y += yy;
	  else y = yy;
	} 
	else 
	  UnitAMultMV<add,false>(A,xx,y);
      }
    } else if (y.step() != 1 || alpha != RealType(T)(1)) {
      Vector<T> yy(y.size());
      if (x.isconj())
	UnitAMultMV<false,true>(A,x,yy.View());
      else
	UnitAMultMV<false,false>(A,x,yy.View());
      if (add) y += alpha*yy;
      else y = alpha*yy;
    } else {
      TMVAssert(alpha == T(1));
      TMVAssert(y.step() == 1);
      TMVAssert(x.step() == 1);
      TMVAssert(!SameStorage(x,y));
      if (x.isconj())
	UnitAMultMV<add,true>(A,x,y);
      else
	UnitAMultMV<add,false>(A,x,y);
    } 
#ifdef XDEBUG
    //cerr<<"y => "<<y<<endl;
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A);
      if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
      else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
      cerr<<endl<<"x = "<<Type(x)<<" step "<<x.step();
      if (x.size() < 30) cerr<<"  "<<x0;
      cerr<<endl<<"y = "<<Type(y)<<" step "<<y.step();
      if (y.size() < 30) cerr<<"  "<<y0;
      cerr<<endl<<"Aptr = "<<A.cptr();
      cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
      if (y.size() < 200) {
	cerr<<"--> y = "<<y<<endl;
	cerr<<"y2 = "<<y2<<endl;
      } else {
	size_t imax;
	MaxAbsElement(y-y2,&imax);
	cerr<<"y("<<imax<<") = "<<y(imax)<<endl;
	cerr<<"y2("<<imax<<") = "<<y2(imax)<<endl;
      }
      cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
      cerr<<"Norm(x0) = "<<Norm(x0)<<endl;
      cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
      cerr<<"|alpha|*|A0|*|x0|+?|y0| = "<<
	ABS(alpha)*Norm(A0)*Norm(x0)+
	(add?Norm(y0):RealType(T)(0))<<endl;
      cerr<<"Norm(y-y2) = "<<Norm(y-y2)<<endl;
      cerr<<"NormInf(y-y2) = "<<NormInf(y-y2)<<endl;
      cerr<<"Norm1(y-y2) = "<<Norm1(y-y2)<<endl;
      abort();
    }
#endif
  }

#ifdef BLAS
  template <class T, class Ta, class Tx> inline void BlasMultMV(
      const T alpha, const GenMatrix<Ta>& A,
      const GenVector<Tx>& x, const int beta, const VectorView<T>& y)
  { 
    if (beta == 0) NonBlasMultMV<false>(alpha,A,x,y); 
    else NonBlasMultMV<true>(alpha,A,x,y); 
  }
#ifdef INST_DOUBLE
  template <> inline void BlasMultMV(
      const double alpha, const GenMatrix<double>& A,
      const GenVector<double>& x, const int beta, const VectorView<double>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    int xs = x.step();
    int ys = y.step();
    double xbeta(beta);

    BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),BLASP(y.ptr()),BLASV(ys)
	BLAS1);
  }
  template <> inline void BlasMultMV(
      const std::complex<double> alpha,
      const GenMatrix<std::complex<double> >& A,
      const GenVector<std::complex<double> >& x,
      const int beta, const VectorView<std::complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    if (x.isconj()
#ifndef CBLAS
	&& !(A.isconj() && A.iscm()) 
#endif
	) {
      Vector<std::complex<double> > xx = alpha*x;
      return BlasMultMV(std::complex<double>(1),A,xx,beta,y);
    } 

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    int xs = x.step();
    int ys = y.step();
    std::complex<double> xbeta(beta);
    if (A.isconj() && A.iscm()) {
#ifdef CBLAS
      std::swap(m,n);
      BLASNAME(zgemv) (BLASRM BLASCH_CT,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	  BLAS1);
#else
      std::complex<double> ca = CONJ(alpha);
      if (x.isconj()) {
	y.ConjugateSelf();
	BLASNAME(zgemv) (BLASCM BLASCH_NT,
	    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
	    BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	    BLAS1);
	y.ConjugateSelf();
      } else {
	Vector<std::complex<double> > xx = ca*x.Conjugate();
	ca = std::complex<double>(1);
	xs = 1;
	y.ConjugateSelf();
	BLASNAME(zgemv) (BLASCM BLASCH_NT,
	    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
	    BLASP(xx.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	    BLAS1);
	y.ConjugateSelf();
      }
#endif
    } else {
      BLASNAME(zgemv) (BLASCM A.isrm()?A.isconj()?BLASCH_CT:BLASCH_T:BLASCH_NT,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	  BLAS1);
    }
  }
#endif
#ifdef INST_FLOAT
  template <> inline void BlasMultMV( 
      const float alpha, const GenMatrix<float>& A,
      const GenVector<float>& x, const int beta, const VectorView<float>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    int xs = x.step();
    int ys = y.step();
    float xbeta(beta);

    BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.cptr()),BLASV(xs),BLASV(xbeta),BLASP(y.ptr()),BLASV(ys)
	BLAS1);
  }
  template <> inline void BlasMultMV(
      const std::complex<float> alpha, 
      const GenMatrix<std::complex<float> >& A,
      const GenVector<std::complex<float> >& x,
      const int beta, const VectorView<std::complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));
    if (x.isconj()
#ifndef CBLAS
	&& !(A.isconj() && A.iscm()) 
#endif
	) {
      Vector<std::complex<float> > xx = alpha*x;
      return BlasMultMV(std::complex<float>(1),A,xx,beta,y);
    }

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    int xs = x.step();
    int ys = y.step();
    std::complex<float> xbeta(beta);
    if (A.isconj() && A.iscm()) {
#ifdef CBLAS
      std::swap(m,n);
      BLASNAME(cgemv) (BLASRM BLASCH_CT,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	  BLAS1);
#else
      std::complex<float> ca = CONJ(alpha);
      if (x.isconj()) {
	y.ConjugateSelf();
	BLASNAME(cgemv) (BLASCM BLASCH_NT,
	    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
	    BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	    BLAS1);
	y.ConjugateSelf();
      } else {
	Vector<std::complex<float> > xx = ca*x.Conjugate();
	ca = std::complex<float>(1);
	xs = 1;
	y.ConjugateSelf();
	BLASNAME(cgemv) (BLASCM BLASCH_NT,
	    BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
	    BLASP(xx.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	    BLAS1);
	y.ConjugateSelf();
      }
#endif
    } else {
      BLASNAME(cgemv) (BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.cptr()),BLASV(xs),BLASP(&xbeta),BLASP(y.ptr()),BLASV(ys)
	  BLAS1);
    }
  }
#endif 
#endif // BLAS

  template <bool add, class T, class Ta, class Tx> inline void DoMultMV(
      const T alpha, const GenMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tx())))
      BlasMultMV(alpha,A,x,add?1:0,y);
    else if ((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0)) {
      if (y.step() > 0 && !SameStorage(A,y)) {
	if (x.step() > 0 && !SameStorage(x,y) && !SameStorage(A,x))
	  BlasMultMV(alpha,A,x,add?1:0,y);
	else {
	  Vector<T> xx = alpha*x;
	  BlasMultMV(T(1),A,xx,add?1:0,y);
	}
      } else {
	Vector<T> yy(y.size());
	if (x.step() > 0 && !SameStorage(A,x)) {
	  BlasMultMV(T(1),A,x,0,yy.View());
	  if (add) y += alpha*yy;
	  else y = alpha*yy;
	} else {
	  Vector<T> xx = alpha*x;
	  BlasMultMV(T(1),A,xx,0,yy.View());
	  if (add) y += yy;
	  else y = yy;
	}
      }
    }
    else {
      if (IMAG(alpha) == T(0)) {
	Matrix<Ta,RowMajor> A2 = REAL(alpha)*A;
	DoMultMV<add>(T(1),A2,x,y);
      } else {
	Matrix<T,RowMajor> A2 = alpha*A;
	DoMultMV<add>(T(1),A2,x,y);
      }
    }
#else
    NonBlasMultMV<add>(alpha,A,x,y);
#endif
  }

  template <bool add, class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y (+)= alpha * A * x
  { 
#ifdef XDEBUG
    //cerr<<"Start MultMV: alpha = "<<alpha<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
    Vector<Tx> x0 = x;
    Vector<T> y0 = y;
    Matrix<Ta> A0 = A;
    Vector<T> y2 = y;
    for(size_t i=0;i<y.size();i++) {
      if (add)
	y2(i) += alpha * (A.row(i) * x0);
      else
	y2(i) = alpha * (A.row(i) * x0);
    }
    //cerr<<"y2 = "<<y2<<endl;
#endif
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) {
	if (!add) y.Zero();
      }
      else if (y.isconj())
	DoMultMV<add>(CONJ(alpha),A.Conjugate(),x.Conjugate(),y.Conjugate());
      else DoMultMV<add>(alpha,A,x,y);
    }
      
#ifdef XDEBUG
    //cerr<<"y => "<<y<<endl;
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A);
      if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
      else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
      cerr<<endl<<"x = "<<Type(x)<<" step "<<x.step();
      if (x.size() < 30) cerr<<"  "<<x0;
      cerr<<endl<<"y = "<<Type(y)<<" step "<<y.step();
      if (y.size() < 30) cerr<<"  "<<y0;
      cerr<<endl<<"Aptr = "<<A.cptr();
      cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
      if (y.size() < 200) {
	cerr<<"--> y = "<<y<<endl;
	cerr<<"y2 = "<<y2<<endl;
      } else {
	size_t imax;
	MaxAbsElement(y-y2,&imax);
	cerr<<"y("<<imax<<") = "<<y(imax)<<endl;
	cerr<<"y2("<<imax<<") = "<<y2(imax)<<endl;
      }
      cerr<<"Norm(A0) = "<<Norm(A0)<<endl;
      cerr<<"Norm(x0) = "<<Norm(x0)<<endl;
      cerr<<"Norm(y0) = "<<Norm(y0)<<endl;
      cerr<<"|alpha|*|A0|*|x0|+?|y0| = "<<
	  ABS(alpha)*Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0))<<endl;
      cerr<<"Norm(y-y2) = "<<Norm(y-y2)<<endl;
      cerr<<"NormInf(y-y2) = "<<NormInf(y-y2)<<endl;
      cerr<<"Norm1(y-y2) = "<<Norm1(y-y2)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


