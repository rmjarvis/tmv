///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG


#include "TMV_Blas.h"
#include "tmv/TMV_MatrixArithFunc.h"
#include "TMV_MultMV.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
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

  template <bool add, bool cx, bool ca, bool rm, class T, class Ta, class Tx>
  static void RowMultMV(
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

    const int M = A.colsize();
    const int N = A.rowsize();
    const int si = A.stepi();
    const int sj = (rm ? 1 : A.stepj());

    const Ta* Ai0 = A.cptr();
    const Tx*const x0 = x.cptr();
    T* yi = y.ptr();

    for(int i=M; i>0; --i,++yi,Ai0+=si) {
      // *yi += A.row(i) * x

      const Ta* Aij = Ai0;
      const Tx* xj = x0;
      register T temp(0);
      for(int j=N; j>0; --j,++xj,(rm?++Aij:Aij+=sj))
        temp += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);

#ifdef TMVFLDEBUG
      TMVAssert(yi >= y.first);
      TMVAssert(yi < y.last);
#endif
      if (add) *yi += temp;
      else *yi = temp;
    }
  }

  template <bool add, bool cx, bool ca, bool cm, class T, class Ta, class Tx> 
  static void ColMultMV(
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

    const int M = A.colsize();
    int N = A.rowsize();
    const int si = (cm ? 1 : A.stepi());
    const int sj = A.stepj();

    const Ta* A0j = A.cptr();
    const Tx* xj = x.cptr();
    T*const y0 = y.ptr();

    if (!add) {
      if (*xj == Tx(0)) {
        y.Zero();
      } else {
        const Ta* Aij = A0j;
        T* yi = y0;
        const Tx xjval = (cx ? CONJ(*xj) : *xj);
        for(int i=M; i>0; --i,++yi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
          TMVAssert(yi >= y.first);
          TMVAssert(yi < y.last);
#endif
          *yi = xjval * (ca ? CONJ(*Aij) : *Aij);
        }
      }
      ++xj; A0j+=sj; --N;
    }

    for(; N>0; --N,++xj,A0j+=sj) {
      // y += *xj * A.col(j)
      if (*xj != Tx(0)) {
        const Ta* Aij = A0j;
        T* yi = y0;
        const Tx xjval = (cx ? CONJ(*xj) : *xj);
        for(int i=M; i>0; --i,++yi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
          TMVAssert(yi >= y.first);
          TMVAssert(yi < y.last);
#endif
          *yi += xjval * (ca ? CONJ(*Aij) : *Aij);
        }
      }
    }
  }

  template <bool add, bool cx, class T, class Ta, class Tx> 
  extern void UnitAMultMV1(
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
  static void UnitAMultMV(
      const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
#ifdef XDEBUG
    //cout<<"Start UnitAMultMV: \n";
    //cout<<"add = "<<add<<endl;
    //cout<<"A = "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"x = "<<TypeText(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cout<<"y = "<<TypeText(y)<<" step "<<y.step()<<"  "<<y<<endl;
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
    //cout<<"y2 = "<<y2<<endl;
#endif
    // Check for 0's in beginning or end of x:
    // y += [ A1 A2 A3 ] [ 0 ]  -->  y += A2 x
    //                   [ x ]
    //                   [ 0 ]

    const int N = x.size(); // = A.rowsize()
    int j2 = N;
    for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
    if (j2 == 0) {
      if (!add) y.Zero();
      return;
    }
    int j1 = 0;
    for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
    TMVAssert(j1 !=j2);
    if (j1 == 0 && j2 == N) UnitAMultMV1<add,cx>(A,x,y);
    else UnitAMultMV1<add,cx>(A.Cols(j1,j2),x.SubVector(j1,j2),y);

#ifdef XDEBUG
    //cout<<"y => "<<y<<endl;
    if (Norm(y-y2) > 0.001*(Norm(A0)*Norm(x0)+
          (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: \n";
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<TypeText(A);
      if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
      else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
      cerr<<endl<<"x = "<<TypeText(x)<<" step "<<x.step();
      if (x.size() < 30) cerr<<"  "<<x0;
      cerr<<endl<<"y = "<<TypeText(y)<<" step "<<y.step();
      if (y.size() < 30) cerr<<"  "<<y0;
      cerr<<endl<<"Aptr = "<<A.cptr();
      cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
      if (y.size() < 200) {
        cerr<<"--> y = "<<y<<endl;
        cerr<<"y2 = "<<y2<<endl;
      } else {
        int imax;
        (y-y2).MaxAbsElement(&imax);
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

  template <bool add, class T, class Ta, class Tx> static void NonBlasMultMV(
      const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  // y (+)= alpha * A * x
  {
#ifdef XDEBUG
    //cout<<"Start MultMV: alpha = "<<alpha<<endl;
    //cout<<"add = "<<add<<endl;
    //cout<<"A = "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"x = "<<TypeText(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cout<<"y = "<<TypeText(y)<<" step "<<y.step()<<"  "<<y<<endl;
    Vector<Tx> x0 = x;
    Vector<T> y0 = y;
    Matrix<Ta> A0 = A;
    Vector<T> y2 = y;
    for(int i=0;i<int(y.size());i++) {
      if (add)
        y2(i) += alpha * (A.row(i) * x0);
      else
        y2(i) = alpha * (A.row(i) * x0);
    }
    //cout<<"y2 = "<<y2<<endl;
#endif
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);

    const int M = A.colsize();
    const int N = A.rowsize();

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
    //cout<<"y => "<<y<<endl;
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
          (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<TypeText(A);
      if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
      else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
      cerr<<endl<<"x = "<<TypeText(x)<<" step "<<x.step();
      if (x.size() < 30) cerr<<"  "<<x0;
      cerr<<endl<<"y = "<<TypeText(y)<<" step "<<y.step();
      if (y.size() < 30) cerr<<"  "<<y0;
      cerr<<endl<<"Aptr = "<<A.cptr();
      cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
      if (y.size() < 200) {
        cerr<<"--> y = "<<y<<endl;
        cerr<<"y2 = "<<y2<<endl;
      } else {
        int imax;
        (y-y2).MaxAbsElement(&imax);
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
  template <class T, class Ta, class Tx> static inline void BlasMultMV(
      const T alpha, const GenMatrix<Ta>& A,
      const GenVector<Tx>& x, const int beta, const VectorView<T>& y)
  { 
    if (beta == 0) NonBlasMultMV<false>(alpha,A,x,y); 
    else NonBlasMultMV<true>(alpha,A,x,y); 
  }
#ifdef INST_DOUBLE
  template <> void BlasMultMV(
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
    const double* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    double* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    double xbeta(beta);

    //std::cout<<"Before dgemv"<<std::endl;
    //std::cout<<"A = "<<A<<std::endl;
    //std::cout<<"x = "<<x<<std::endl;
    //std::cout<<"y = "<<y<<std::endl;
    //std::cout<<"m = "<<m<<std::endl;
    //std::cout<<"n = "<<n<<std::endl;
    //std::cout<<"lda = "<<lda<<std::endl;
    //std::cout<<"xs = "<<xs<<std::endl;
    //std::cout<<"ys = "<<ys<<std::endl;
    //std::cout<<"alpha = "<<alpha<<std::endl;
    //std::cout<<"beta = "<<xbeta<<std::endl;
    //std::cout<<"aptr = "<<A.cptr()<<std::endl;
    //std::cout<<"xp = "<<xp<<std::endl;
    //std::cout<<"yp = "<<yp<<std::endl;
    //std::cout<<"NT = "<<(A.isrm()?'T':'N')<<std::endl;
    //if (A.isrm()) {
    //std::cout<<"x.size = "<<x.size()<<std::endl;
    //std::cout<<"x = ";
    //for(int i=0;i<m;i++) std::cout<<*(xp+i*xs)<<" ";
    //std::cout<<std::endl;
    //std::cout<<"y.size = "<<y.size()<<std::endl;
    //std::cout<<"y = ";
    //for(int i=0;i<n;i++) std::cout<<*(yp+i*ys)<<" ";
    //std::cout<<std::endl;
    //std::cout<<"A.size = "<<A.colsize()<<','<<A.rowsize()<<std::endl;
    //std::cout<<"A = ";
    //for(int i=0;i<n*lda;i++) std::cout<<*(A.cptr()+i)<<" ";
    //std::cout<<std::endl;
    //}
    BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
        BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
        BLASP(xp),BLASV(xs),BLASV(xbeta),BLASP(yp),BLASV(ys)
        BLAS1);
    //std::cout<<"After dgemv"<<std::endl;
  }
  template <> void BlasMultMV(
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
    const std::complex<double>* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    std::complex<double>* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    std::complex<double> xbeta(beta);
    if (A.isconj() && A.iscm()) {
#ifdef CBLAS
      __TMV_SWAP(m,n);
      BLASNAME(zgemv) (BLASRM BLASCH_CT,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
          BLAS1);
#else
      std::complex<double> ca = CONJ(alpha);
      if (x.isconj()) {
        y.ConjugateSelf();
        BLASNAME(zgemv) (BLASCM BLASCH_NT,
            BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
        y.ConjugateSelf();
      } else {
        Vector<std::complex<double> > xx = ca*x.Conjugate();
        ca = std::complex<double>(1);
        xs = 1;
        xp = xx.cptr();
        y.ConjugateSelf();
        BLASNAME(zgemv) (BLASCM BLASCH_NT,
            BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
        y.ConjugateSelf();
      }
#endif
    } else {
      BLASNAME(zgemv) (BLASCM A.isrm()?A.isconj()?BLASCH_CT:BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
          BLAS1);
    }
  }
  template <> void BlasMultMV(
      const std::complex<double> alpha,
      const GenMatrix<std::complex<double> >& A,
      const GenVector<double>& x,
      const int beta, const VectorView<std::complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    if (A.iscm()) {
      if (y.step() != 1) {
        Vector<std::complex<double> > yy(y.size());
        BlasMultMV(std::complex<double>(1),A,x,0,yy.View());
        if (beta == 0) y = alpha*yy;
        else y += alpha*yy;
      } else {
        if (beta == 0) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          int xs = x.step();
          int ys = 1;
          const double* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          double* yp = (double*) y.ptr();
          double xalpha(1);
          double xbeta(0);
          BLASNAME(dgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((double*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(xbeta),
              BLASP(yp),BLASV(ys) BLAS1);
          if (A.isconj()) y.ConjugateSelf();
          y *= alpha;
        } else if (A.isconj()) {
          Vector<std::complex<double> > yy(y.size());
          BlasMultMV(std::complex<double>(1),A.Conjugate(),x,0,yy.View());
          y += alpha*yy.Conjugate();
        } else if (IMAG(alpha) == 0.) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          int xs = x.step();
          int ys = 1;
          const double* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          double* yp = (double*) y.ptr();
          if (ys < 0) yp += (y.size()-1)*ys;
          double xalpha(REAL(alpha));
          double xbeta(1);
          BLASNAME(dgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((double*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(xbeta),
              BLASP(yp),BLASV(ys) BLAS1);
        } else {
          Vector<std::complex<double> > yy(y.size());
          BlasMultMV(std::complex<double>(1),A,x,0,yy.View());
          y += alpha*yy;
        }
      }
    } else { // A.isrm
      BlasMultMV(alpha,A,Vector<std::complex<double> >(x),beta,y);
    }
  }
  template <> void BlasMultMV(
      const std::complex<double> alpha,
      const GenMatrix<double>& A,
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

    if (beta == 0) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      int xs = 2*x.step();
      int ys = 2*y.step();
      const double* xp = (const double*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      double* yp = (double*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      double xalpha(1);
      double xbeta(beta);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(xbeta),
          BLASP(yp+1),BLASV(ys) BLAS1);
      if (x.isconj()) y.ConjugateSelf();
      y *= alpha;
    } else if (IMAG(alpha) == 0. && !x.isconj()) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      int xs = 2*x.step();
      int ys = 2*y.step();
      const double* xp = (const double*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      double* yp = (double*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      double xalpha(REAL(alpha));
      double xbeta(beta);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(xbeta),
          BLASP(yp+1),BLASV(ys) BLAS1);
    } else {
      Vector<std::complex<double> > xx = alpha*x;
      BlasMultMV(std::complex<double>(1),A,xx,1,y);
    }
  }
  template <> void BlasMultMV(
      const std::complex<double> alpha,
      const GenMatrix<double>& A,
      const GenVector<double>& x,
      const int beta, const VectorView<std::complex<double> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    int xs = x.step();
    int ys = 2*y.step();
    const double* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    double* yp = (double*) y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    double ar(REAL(alpha));
    double ai(IMAG(alpha));
    double xbeta(beta);
    if (ar == 0.) {
      if (beta == 0) y.Real().Zero();
    }
    else 
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp),BLASV(ys) BLAS1);
    if (ai == 0.) {
      if (beta == 0) y.Imag().Zero();
    }
    else
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp+1),BLASV(ys) BLAS1);
  }
#endif
#ifdef INST_FLOAT
  template <> void BlasMultMV(
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
    const float* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    float* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    float xbeta(beta);

    BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
        BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
        BLASP(xp),BLASV(xs),BLASV(xbeta),BLASP(yp),BLASV(ys)
        BLAS1);
  }
  template <> void BlasMultMV(
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
    const std::complex<float>* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    std::complex<float>* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    std::complex<float> xbeta(beta);
    if (A.isconj() && A.iscm()) {
#ifdef CBLAS
      __TMV_SWAP(m,n);
      BLASNAME(cgemv) (BLASRM BLASCH_CT,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
          BLAS1);
#else
      std::complex<float> ca = CONJ(alpha);
      if (x.isconj()) {
        y.ConjugateSelf();
        BLASNAME(cgemv) (BLASCM BLASCH_NT,
            BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
        y.ConjugateSelf();
      } else {
        Vector<std::complex<float> > xx = ca*x.Conjugate();
        ca = std::complex<float>(1);
        xs = 1;
        xp = xx.cptr();
        y.ConjugateSelf();
        BLASNAME(cgemv) (BLASCM BLASCH_NT,
            BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
        y.ConjugateSelf();
      }
#endif
    } else {
      BLASNAME(cgemv) (BLASCM A.isrm()?A.isconj()?BLASCH_CT:BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
          BLAS1);
    }
  }
  template <> void BlasMultMV(
      const std::complex<float> alpha,
      const GenMatrix<std::complex<float> >& A,
      const GenVector<float>& x,
      const int beta, const VectorView<std::complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    if (A.iscm()) {
      if (y.step() != 1) {
        Vector<std::complex<float> > yy(y.size());
        BlasMultMV(std::complex<float>(1),A,x,0,yy.View());
        if (beta == 0) y = alpha*yy;
        else y += alpha*yy;
      } else {
        if (beta == 0) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          int xs = x.step();
          int ys = 1;
          const float* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          float* yp = (float*) y.ptr();
          float xalpha(1);
          float xbeta(0);
          BLASNAME(sgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((float*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(xbeta),
              BLASP(yp),BLASV(ys) BLAS1);
          if (A.isconj()) y.ConjugateSelf();
          y *= alpha;
        } else if (A.isconj()) {
          Vector<std::complex<float> > yy(y.size());
          BlasMultMV(std::complex<float>(1),A.Conjugate(),x,0,yy.View());
          y += alpha*yy.Conjugate();
        } else if (IMAG(alpha) == 0.F) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          int xs = x.step();
          int ys = 1;
          const float* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          float* yp = (float*) y.ptr();
          float xalpha(REAL(alpha));
          float xbeta(1);
          BLASNAME(sgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((float*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(xbeta),
              BLASP(yp),BLASV(ys) BLAS1);
        } else {
          Vector<std::complex<float> > yy(y.size());
          BlasMultMV(std::complex<float>(1),A,x,0,yy.View());
          y += alpha*yy;
        }
      } 
    } else { // A.isrm
      BlasMultMV(alpha,A,Vector<std::complex<float> >(x),beta,y);
    }
  }
  template <> void BlasMultMV(
      const std::complex<float> alpha,
      const GenMatrix<float>& A,
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

    if (beta == 0) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      int xs = 2*x.step();
      int ys = 2*y.step();
      const float* xp = (const float*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      float* yp = (float*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      float xalpha(1);
      float xbeta(beta);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(xbeta),
          BLASP(yp+1),BLASV(ys) BLAS1);
      if (x.isconj()) y.ConjugateSelf();
      y *= alpha;
    } else if (IMAG(alpha) == 0.F && !x.isconj()) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      int xs = 2*x.step();
      int ys = 2*y.step();
      const float* xp = (const float*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      float* yp = (float*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      float xalpha(REAL(alpha));
      float xbeta(beta);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(xbeta),
          BLASP(yp+1),BLASV(ys) BLAS1);
    } else {
      Vector<std::complex<float> > xx = alpha*x;
      BlasMultMV(std::complex<float>(1),A,xx,1,y);
    }
  }
  template <> void BlasMultMV(
      const std::complex<float> alpha,
      const GenMatrix<float>& A,
      const GenVector<float>& x,
      const int beta, const VectorView<std::complex<float> >& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(!SameStorage(x,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    int xs = x.step();
    int ys = 2*y.step();
    const float* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    float* yp = (float*) y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    float ar(REAL(alpha));
    float ai(IMAG(alpha));
    float xbeta(beta);
    if (ar == 0.F) {
      if (beta == 0) y.Real().Zero();
    }
    else 
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp),BLASV(ys) BLAS1);
    if (ai == 0.F) {
      if (beta == 0) y.Imag().Zero();
    }
    else
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(xbeta),
          BLASP(yp+1),BLASV(ys) BLAS1);
  }
#endif
#endif // BLAS

  template <bool add, class T, class Ta, class Tx> static void DoMultMV(
      const T alpha, const GenMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);

#ifdef BLAS
    if (x.step() == 0) {
      if (x.size() <= 1) 
        DoMultMV<add>(alpha,A,
            ConstVectorView<Tx>(x.cptr(),x.size(),1,x.ct()),y);
      else 
        DoMultMV<add>(alpha,A,Vector<Tx>(x),y);
    } else if (y.step() == 0) {
      TMVAssert(y.size() <= 1);
      DoMultMV<add>(alpha,A,x,VectorView<T>(y.ptr(),y.size(),1,y.ct()));
    } else if ((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0)) {
      if (/*y.step() > 0 &&*/ !SameStorage(A,y)) {
        if (/*x.step() > 0 &&*/ !SameStorage(x,y) && !SameStorage(A,x))
          BlasMultMV(alpha,A,x,add?1:0,y);
        else {
          Vector<T> xx = alpha*x;
          BlasMultMV(T(1),A,xx,add?1:0,y);
        }
      } else {
        Vector<T> yy(y.size());
        if (/*x.step() > 0 &&*/ !SameStorage(A,x)) {
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
    cout<<"Start MultMV: alpha = "<<alpha<<endl;
    cout<<"add = "<<add<<endl;
    cout<<"A = "<<TypeText(A)<<"  "<<A<<endl;
    cout<<"x = "<<TypeText(x)<<" step "<<x.step()<<"  "<<x<<endl;
    cout<<"y = "<<TypeText(y)<<" step "<<y.step()<<"  "<<y<<endl;
    Vector<Tx> x0 = x;
    Vector<T> y0 = y;
    Matrix<Ta> A0 = A;
    Vector<T> y2 = y;
    for(int i=0;i<int(y.size());i++) {
      if (add)
        y2(i) += alpha * (A.row(i) * x0);
      else
        y2(i) = alpha * (A.row(i) * x0);
    }
    cout<<"y2 = "<<y2<<endl;
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
    cout<<"y => "<<y<<endl;
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
          (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<TypeText(A);
      if (A.rowsize() < 30 && A.colsize() < 30) cerr<<"  "<<A0;
      else cerr<<"  "<<A.colsize()<<" x "<<A.rowsize();
      cerr<<endl<<"x = "<<TypeText(x)<<" step "<<x.step();
      if (x.size() < 30) cerr<<"  "<<x0;
      cerr<<endl<<"y = "<<TypeText(y)<<" step "<<y.step();
      if (y.size() < 30) cerr<<"  "<<y0;
      cerr<<endl<<"Aptr = "<<A.cptr();
      cerr<<", xptr = "<<x.cptr()<<", yptr = "<<y.cptr()<<endl;
      if (y.size() < 200) {
        cerr<<"--> y = "<<y<<endl;
        cerr<<"y2 = "<<y2<<endl;
      } else {
        int imax;
        (y-y2).MaxAbsElement(&imax);
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

#define InstFile "TMV_MultMV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


