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
#include <iostream>


#include "TMV_Blas.h"
#include "tmv/TMV_MatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool cx, bool cy, bool cm, bool add, class T, class Tx, class Ty> 
  static void ColRank1Update(
      const GenVector<Tx>& x, const GenVector<Ty>& y,
      const MatrixView<T>& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(!A.isrm());
    TMVAssert(x.step() == 1);
    TMVAssert(cm == A.iscm());
    TMVAssert(cx == x.isconj());
    TMVAssert(cy == y.isconj());

    const Ty* yj = y.cptr();
    const Tx*const xptr = x.cptr();
    T* Acolj = A.ptr();
    const int sj = A.stepj();
    const int si = (cm ? 1 : A.stepi());
    const int ys = y.step();
    const int M = A.colsize();
    const int N = A.rowsize();

    for (int j=N; j>0; --j,yj+=ys,Acolj+=sj) {
      if (*yj!=Ty(0)) {
        T* Aij = Acolj;
        const Tx* xi = xptr;
        for (int i=M; i>0; --i,++xi,(cm?++Aij:Aij+=si)) {
          const T temp = (cx ? CONJ(*xi) : *xi) * (cy ? CONJ(*yj) : *yj);
#ifdef TMVFLDEBUG
          TMVAssert(Aij >= A.first);
          TMVAssert(Aij < A.last);
#endif
          if (add) *Aij += temp;
          else *Aij = temp;
        }
      } else if (!add) {
        T* Aij = Acolj;
        if (cm) std::fill_n(Aij,M,T(0));
        else for (int i=M; i>0; --i,Aij+=si) {
#ifdef TMVFLDEBUG
          TMVAssert(Aij >= A.first);
          TMVAssert(Aij < A.last);
#endif
          *Aij = T(0);
        }
      }
    }
  }

  template <bool cx, bool add, class T, class Tx, class Ty> 
  static void UnitARank1Update(
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

    if (A.iscm()) 
      if (y.isconj())
        ColRank1Update<cx,true,true,add>(x,y,A);
      else
        ColRank1Update<cx,false,true,add>(x,y,A);
    else
      if (y.isconj())
        ColRank1Update<cx,true,false,add>(x,y,A);
      else
        ColRank1Update<cx,false,false,add>(x,y,A);
  }

  template <bool add, class T, class Tx, class Ty> 
  static void NonBlasRank1Update(
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
          if (x.isconj()) UnitARank1Update<true,add>(x,yy,A);
          else UnitARank1Update<false,add>(x,yy,A);
        } else {
          Vector<T> yy = alpha*y;
          if (x.isconj()) UnitARank1Update<true,add>(x,yy,A);
          else UnitARank1Update<false,add>(x,yy,A);
        }
      } else {
        if (IMAG(alpha) == RealType(T)(0)) {
          Vector<Tx> xx = REAL(alpha)*x;
          UnitARank1Update<false,add>(xx,y,A);
        } else {
          Vector<T> xx = alpha*x;
          UnitARank1Update<false,add>(xx,y,A);
        }
      }
    } else {
      if (x.isconj()) UnitARank1Update<true,add>(x,y,A);
      else UnitARank1Update<false,add>(x,y,A);
    }
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> static inline void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  { NonBlasRank1Update<true>(alpha,x,y,A); }
#ifdef INST_DOUBLE
  template <> void BlasRank1Update(
      const double alpha, const GenVector<double>& x,
      const GenVector<double>& y, const MatrixView<double>& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    const double* xp = x.cptr();
    if (xs < 0) xp += (m-1)*xs;
    const double* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(alpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP(A.ptr()),BLASV(lda));
  }
  template <> void BlasRank1Update(
      const std::complex<double> alpha,
      const GenVector<std::complex<double> >& x, 
      const GenVector<std::complex<double> >& y,
      const MatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    const std::complex<double>* xp = x.cptr();
    if (xs < 0) xp += (m-1)*xs;
    const std::complex<double>* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    if (x.isconj()) {
#ifdef CBLAS
      BLASNAME(zgerc) (BLASRM BLASV(n),BLASV(m),BLASP(&alpha),
          BLASP(yp),BLASV(ys),BLASP(xp),BLASV(xs),
          BLASP(A.ptr()),BLASV(lda));
#else
      Vector<std::complex<double> > xx = alpha*x;
      xs = 1;
      xp = xx.cptr();
      std::complex<double> alpha2(1);
      BLASNAME(zgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha2),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP(A.ptr()),BLASV(lda));
#endif
    }
    else if (y.isconj())
      BLASNAME(zgerc) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP(A.ptr()),BLASV(lda));
    else
      BLASNAME(zgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP(A.ptr()),BLASV(lda));
  }
  template <> void BlasRank1Update(
      const std::complex<double> alpha,
      const GenVector<double>& x, 
      const GenVector<std::complex<double> >& y,
      const MatrixView<std::complex<double> >& A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ (yr + I yi)
    // Ar += ar * x ^ yr - ai * x ^ yi
    // Ai += ai * x ^ yr + ar * x ^ yi
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = 2*y.step();
    const double* yp = (const double*) y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    Vector<double> xx(2*x.size());
    xx.SubVector(0,xx.size(),2) = REAL(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = IMAG(alpha)*x;
    const double* xp = xx.cptr();
    double xalpha(1);
    if (ys == 0) {
      std::cout<<"x = "<<x<<std::endl;
      std::cout<<"y = "<<y<<std::endl;
      std::cout<<"A = "<<A<<std::endl;
    }
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((double*)A.ptr()),BLASV(lda));
    if (y.isconj()) {
      xx.SubVector(0,xx.size(),2) = IMAG(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = -REAL(alpha)*x;
    } else {
      xx.SubVector(0,xx.size(),2) = -IMAG(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = REAL(alpha)*x;
    }
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp+1),BLASV(ys),
        BLASP((double*)A.ptr()),BLASV(lda));
  }
  template <> void BlasRank1Update(
      const std::complex<double> alpha,
      const GenVector<std::complex<double> >& x,
      const GenVector<double>& y, 
      const MatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = y.step();
    const double* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    if (x.step() == 1 && !x.isconj() && IMAG(alpha) == 0.) {
      const double* xp = (double*) x.cptr();
      double xalpha(REAL(alpha));
      BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP((double*)A.ptr()),BLASV(lda));
    } else {
      Vector<std::complex<double> > xx = alpha*x;
      const double* xp = (double*) xx.cptr();
      double xalpha(1);
      BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP((double*)A.ptr()),BLASV(lda));
    } 
  }
  template <> void BlasRank1Update(
      const std::complex<double> alpha,
      const GenVector<double>& x, const GenVector<double>& y, 
      const MatrixView<std::complex<double> >& A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ y
    // Ar += ar * x ^ y
    // Ai += ai * x ^ y
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = y.step();
    const double* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    Vector<double> xx(2*x.size());
    xx.SubVector(0,xx.size(),2) = REAL(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = IMAG(alpha)*x;
    const double* xp = xx.cptr();
    double xalpha(1);
    if (ys == 0) {
      std::cout<<"x = "<<x<<std::endl;
      std::cout<<"y = "<<y<<std::endl;
      std::cout<<"A = "<<A<<std::endl;
    }
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((double*)A.ptr()),BLASV(lda));
  }
#endif
#ifdef INST_FLOAT
  template <> void BlasRank1Update(
      const float alpha, const GenVector<float>& x,
      const GenVector<float>& y, const MatrixView<float>& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.F);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    const float* xp = x.cptr();
    if (xs < 0) xp += (m-1)*xs;
    const float* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(alpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP(A.ptr()),BLASV(lda));
  }
  template <> void BlasRank1Update(
      const std::complex<float> alpha,
      const GenVector<std::complex<float> >& x, 
      const GenVector<std::complex<float> >& y,
      const MatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.F);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj || y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = A.colsize();
    int n = A.rowsize();
    int xs = x.step();
    int ys = y.step();
    const std::complex<float>* xp = x.cptr();
    if (xs < 0) xp += (m-1)*xs;
    const std::complex<float>* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    if (x.isconj()) {
#ifdef CBLAS
      BLASNAME(cgerc) (BLASRM BLASV(n),BLASV(m),BLASP(&alpha),
          BLASP(yp),BLASV(ys),BLASP(xp),BLASV(xs),
          BLASP(A.ptr()),BLASV(lda));
#else
      Vector<std::complex<float> > xx = alpha*x;
      xs = 1;
      xp = xx.cptr();
      std::complex<float> alpha2(1);
      BLASNAME(cgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha2),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP(A.ptr()),BLASV(lda));
#endif
    }
    else if (y.isconj())
      BLASNAME(cgerc) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP(A.ptr()),BLASV(lda));
    else
      BLASNAME(cgeru) (BLASCM BLASV(m),BLASV(n),BLASP(&alpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP(A.ptr()),BLASV(lda));
  }
  template <> void BlasRank1Update(
      const std::complex<float> alpha,
      const GenVector<float>& x, 
      const GenVector<std::complex<float> >& y,
      const MatrixView<std::complex<float> >& A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ (yr + I yi)
    // Ar += ar * x ^ yr - ai * x ^ yi
    // Ai += ai * x ^ yr + ar * x ^ yi
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.F);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = 2*y.step();
    const float* yp = (float*) y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    Vector<float> xx(2*x.size());
    xx.SubVector(0,xx.size(),2) = REAL(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = IMAG(alpha)*x;
    const float* xp = xx.cptr();
    float xalpha(1);
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((float*)A.ptr()),BLASV(lda));
    if (y.isconj()) {
      xx.SubVector(0,xx.size(),2) = IMAG(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = -REAL(alpha)*x;
    } else {
      xx.SubVector(0,xx.size(),2) = -IMAG(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = REAL(alpha)*x;
    }
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp+1),BLASV(ys),
        BLASP((float*)A.ptr()),BLASV(lda));
  }
  template <> void BlasRank1Update(
      const std::complex<float> alpha,
      const GenVector<std::complex<float> >& x,
      const GenVector<float>& y, 
      const MatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.F);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = y.step();
    const float* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    if (x.step() == 1 && !x.isconj() && IMAG(alpha) == 0.F) {
      const float* xp = (float*) x.cptr();
      float xalpha(REAL(alpha));
      BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP((float*)A.ptr()),BLASV(lda));
    } else {
      Vector<std::complex<float> > xx = alpha*x;
      const float* xp = (float*) xx.cptr();
      float xalpha(1);
      BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
          BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
          BLASP((float*)A.ptr()),BLASV(lda));
    } 
  }
  template <> void BlasRank1Update(
      const std::complex<float> alpha,
      const GenVector<float>& x, const GenVector<float>& y, 
      const MatrixView<std::complex<float> >& A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ y
    // Ar += ar * x ^ y
    // Ai += ai * x ^ y
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(alpha != 0.F);
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.iscm());

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = y.step();
    const float* yp = y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    Vector<float> xx(2*x.size());
    xx.SubVector(0,xx.size(),2) = REAL(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = IMAG(alpha)*x;
    const float* xp = xx.cptr();
    float xalpha(1);
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((float*)A.ptr()),BLASV(lda));
  }
#endif
#endif // BLAS

  template <bool add, class T, class Tx, class Ty> void Rank1Update(
      const T alpha, const GenVector<Tx>& x,
      const GenVector<Ty>& y, const MatrixView<T>& A)
  // A (+)= beta + alpha * x * yT
  {
#ifdef XDEBUG
    //cout<<"Rank1Update: alpha = "<<alpha<<endl;
    //cout<<"add = "<<add<<endl;
    //cout<<"x = "<<TypeText(x)<<"  "<<x<<endl;
    //cout<<"y = "<<TypeText(y)<<"  "<<y<<endl;
    //cout<<"A = "<<TypeText(A)<<"  "<<A<<endl;
    Vector<Tx> x0 = x;
    Vector<Ty> y0 = y;
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    for(int i=0;i<int(x.size());i++) for(int j=0;j<int(y.size());j++) 
      if (add)
        A2(i,j) += alpha*x0(i)*y0(j);
      else
        A2(i,j) = alpha*x0(i)*y0(j);
#endif

    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());

    if (A.colsize() > 0 && A.rowsize() > 0) {
      if (alpha == T(0)) {
        if (!add) A.Zero();
      } else {
        if (A.isconj()) 
          Rank1Update<add>(CONJ(alpha),x.Conjugate(),y.Conjugate(),
              A.Conjugate());
        else if (A.isrm())
          Rank1Update<add>(alpha,y,x,A.Transpose());
#ifdef BLAS
        else if (!((A.iscm() && A.stepj()>0))) {
          Matrix<T,ColMajor> A2(A);
          Rank1Update<add>(alpha,x,y,A2.View());
          A = A2;
        } else {
          if (SameStorage(x,A)) {
            if (SameStorage(y,A)) {
              if (IMAG(alpha) == RealType(T)(0)) {
                if (x.size() <= y.size()) {
                  Vector<Tx> xx = REAL(alpha)*x;
                  Vector<Ty> yy = y;
                  if (!add) A.Zero();
                  BlasRank1Update(T(1),xx,yy,A);
                } else {
                  Vector<Tx> xx = x;
                  Vector<Ty> yy = REAL(alpha)*y;
                  if (!add) A.Zero();
                  BlasRank1Update(T(1),xx,yy,A);
                }
              } else {
                if (x.size() <= y.size()) {
                  Vector<T> xx = alpha*x;
                  Vector<Ty> yy = y;
                  if (!add) A.Zero();
                  BlasRank1Update(T(1),xx,yy,A);
                } else {
                  Vector<Tx> xx = x;
                  Vector<T> yy = alpha*y;
                  if (!add) A.Zero();
                  BlasRank1Update(T(1),xx,yy,A);
                }
              }
            } else {
              if (IMAG(alpha) == RealType(T)(0)) {
                Vector<Tx> xx = REAL(alpha)*x;
                if (!add) A.Zero();
                BlasRank1Update(T(1),xx,y,A);
              } else {
                Vector<T> xx = alpha*x;
                if (!add) A.Zero();
                BlasRank1Update(T(1),xx,y,A);
              }
            }
          } else {
            if (SameStorage(A,y)) {
              if (IMAG(alpha) == RealType(T)(0)) {
                Vector<Ty> yy = REAL(alpha)*y;
                if (!add) A.Zero();
                BlasRank1Update(T(1),x,yy,A);
              } else {
                Vector<T> yy = alpha*y;
                if (!add) A.Zero();
                BlasRank1Update(T(1),x,yy,A);
              }
            } else {
              if (!add) A.Zero();
              if (x.isconj() && y.isconj()) {
                if (IMAG(alpha) == RealType(T)(0)) {
                  if (x.size() <= y.size()) {
                    Vector<Tx> xx = REAL(alpha)*x;
                    BlasRank1Update(T(1),xx,y,A);
                  } else {
                    Vector<Ty> yy = REAL(alpha)*y;
                    BlasRank1Update(T(1),x,yy,A);
                  }
                } else {
                  if (x.size() <= y.size()) {
                    Vector<T> xx = alpha*x;
                    BlasRank1Update(T(1),xx,y,A);
                  } else {
                    Vector<T> yy = alpha*y;
                    BlasRank1Update(T(1),x,yy,A);
                  }
                }
              } else {
                BlasRank1Update(alpha,x,y,A);
              }
            }
          }
        }
#else
        else NonBlasRank1Update<add>(alpha,x,y,A);
#endif
      }
    }

#ifdef XDEBUG
    //cout<<"Done Rank1Update: A->"<<A<<endl;
    if (Norm(A-A2) > 0.001*(ABS(alpha)*Norm(x0)*Norm(y0)+
          (add?Norm(A0):RealType(T)(0)))) {
      cerr<<"Rank1Update: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"x = "<<TypeText(x)<<"  step = "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<TypeText(y)<<"  step = "<<y.step()<<"  "<<y0<<endl;
      cerr<<"A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_Rank1_VVM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

