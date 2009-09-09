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


#include "tmv/TMV_AddVV.h"

namespace tmv {

  // 
  //
  // MultMV
  //

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
  template <int Si, int Sj>
  static void BlasMultMV(
      double alpha,
      const ConstMatrixView<double,Si,Sj>& A,
      const ConstVectorView<double>& x,
      double beta, VectorView<double> y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));
    TMVAssert(!SameStorage(A,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    if (lda < m) { TMVAssert(n==1); lda = m; }
    int xs = x.step();
    int ys = y.step();
    if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
    const double* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    double* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
        BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
        BLASP(xp),BLASV(xs),BLASV(beta),BLASP(yp),BLASV(ys)
        BLAS1);
  }
  template <int Si, int Sj, bool C1, bool C2>
  static void BlasMultMV(
      std::complex<double> alpha,
      const ConstMatrixView<std::complex<double>,Si,Sj,C1>& A,
      const ConstVectorView<std::complex<double>,UNKNOWN,C2>& x,
      double beta, VectorView<std::complex<double> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));
    TMVAssert(!SameStorage(A,y));

    if (x.isconj()
#ifndef CBLAS
        && !(A.isconj() && A.iscm()) 
#endif
       ) {
      const Vector<std::complex<double> > xx = alpha*x;
      return BlasMultMV(std::complex<double>(1.),A,xx.XView(),beta,y);
    } 

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    if (lda < m) { TMVAssert(n==1); lda = m; }
    int xs = x.step();
    int ys = y.step();
    if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
    const std::complex<double>* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    std::complex<double>* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    std::complex<double> xbeta = beta;
    if (A.isconj() && A.iscm()) {
#ifdef CBLAS
      TMV_SWAP(m,n);
      BLASNAME(zgemv) (BLASRM BLASCH_CT,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
          BLAS1);
#else
      std::complex<double> ca = TMV_CONJ(alpha);
      if (x.isconj()) {
        y.ConjugateSelf();
        BLASNAME(zgemv) (BLASCM BLASCH_NT,
            BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
        y.ConjugateSelf();
      } else {
        Vector<std::complex<double> > xx = ca*x.Conjugate();
        ca = std::complex<double>(1.);
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
  template <int Si, int Sj, bool C1>
  static void BlasMultMV(
      std::complex<double> alpha,
      const ConstMatrixView<std::complex<double>,Si,Sj,C1>& A,
      const ConstVectorView<double>& x,
      double beta, VectorView<std::complex<double> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    if (A.iscm()) {
      if (y.step() != 1) {
        Vector<std::complex<double> > yy(y.size());
        BlasMultMV(1.,A,x,0.,yy.View());
        if (beta == 0.) y = alpha*yy;
        else y += alpha*yy;
      } else {
        if (beta == 0.) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          if (lda < m) { TMVAssert(n==1); lda = m; }
          int xs = x.step();
          int ys = 1;
          const double* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          double* yp = (double*) y.ptr();
          double xalpha(1);
          BLASNAME(dgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((double*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(beta),
              BLASP(yp),BLASV(ys) BLAS1);
          if (A.isconj()) y.ConjugateSelf();
          y *= alpha;
        } else if (A.isconj()) {
          Vector<std::complex<double> > yy(y.size());
          BlasMultMV(1.,A.Conjugate(),x,0.,yy.View());
          y += alpha*yy.Conjugate();
        } else if (TMV_IMAG(alpha) == 0.) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          if (lda < m) { TMVAssert(n==1); lda = m; }
          int xs = x.step();
          int ys = 1;
          const double* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          double* yp = (double*) y.ptr();
          if (ys < 0) yp += (y.size()-1)*ys;
          double xalpha(TMV_REAL(alpha));
          BLASNAME(dgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((double*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(beta),
              BLASP(yp),BLASV(ys) BLAS1);
        } else {
          Vector<std::complex<double> > yy(y.size());
          BlasMultMV(std::complex<double>(1.),A,x,0.,yy.View());
          y += alpha*yy;
        }
      }
    } else { // A.isrm
      const Vector<std::complex<double> > xx = x;
      BlasMultMV(alpha,A,xx.XView(),beta,y);
    }
  }
  template <int Si, int Sj, bool C2>
  static void BlasMultMV(  
      std::complex<double> alpha,
      const ConstMatrixView<double,Si,Sj>& A,
      const ConstVectorView<std::complex<double>,UNKNOWN,C2>& x,
      double beta, VectorView<std::complex<double> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    if (beta == 0.) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      if (lda < m) { TMVAssert(n==1); lda = m; }
      int xs = 2*x.step();
      int ys = 2*y.step();
      if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
      if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
      const double* xp = (const double*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      double* yp = (double*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      double xalpha(1);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(beta),
          BLASP(yp+1),BLASV(ys) BLAS1);
      if (x.isconj()) y.ConjugateSelf();
      y *= alpha;
    } else if (TMV_IMAG(alpha) == 0. && !x.isconj()) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      if (lda < m) { TMVAssert(n==1); lda = m; }
      int xs = 2*x.step();
      int ys = 2*y.step();
      if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
      if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
      const double* xp = (const double*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      double* yp = (double*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      double xalpha(TMV_REAL(alpha));
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(beta),
          BLASP(yp+1),BLASV(ys) BLAS1);
    } else {
      const Vector<std::complex<double> > xx = alpha*x;
      BlasMultMV(std::complex<double>(1.),A,xx.XView(),1.,y);
    }
  }
  template <int Si, int Sj>
  static void BlasMultMV(
      std::complex<double> alpha,
      const ConstMatrixView<double,Si,Sj>& A,
      const ConstVectorView<double>& x,
      double beta, VectorView<std::complex<double> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    if (lda < m) { TMVAssert(n==1); lda = m; }
    int xs = x.step();
    int ys = 2*y.step();
    if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
    const double* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    double* yp = (double*) y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    double ar(TMV_REAL(alpha));
    double ai(TMV_IMAG(alpha));
    if (ar == 0.) {
      if (beta == 0.) y.Real().Zero();
    }
    else 
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp),BLASV(ys) BLAS1);
    if (ai == 0.) {
      if (beta == 0.) y.Imag().Zero();
    }
    else
      BLASNAME(dgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp+1),BLASV(ys) BLAS1);
  }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
  template <int Si, int Sj>
  static void BlasMultMV(
      float alpha,
      const ConstMatrixView<float,Si,Sj>& A,
      const ConstVectorView<float>& x,
      float beta, VectorView<float> y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));
    TMVAssert(!SameStorage(A,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    if (lda < m) { TMVAssert(n==1); lda = m; }
    int xs = x.step();
    int ys = y.step();
    if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
    const float* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    float* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
        BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
        BLASP(xp),BLASV(xs),BLASV(beta),BLASP(yp),BLASV(ys)
        BLAS1);
  }
  template <int Si, int Sj, bool C1, bool C2>
  static void BlasMultMV(
      std::complex<float> alpha,
      const ConstMatrixView<std::complex<float>,Si,Sj,C1>& A,
      const ConstVectorView<std::complex<float>,UNKNOWN,C2>& x,
      float beta, VectorView<std::complex<float> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));
    TMVAssert(!SameStorage(A,y));

    if (x.isconj()
#ifndef CBLAS
        && !(A.isconj() && A.iscm()) 
#endif
       ) {
      const Vector<std::complex<float> > xx = alpha*x;
      return BlasMultMV(std::complex<float>(1.F),A,xx.XView(),beta,y);
    } 

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    if (lda < m) { TMVAssert(n==1); lda = m; }
    int xs = x.step();
    int ys = y.step();
    if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
    const std::complex<float>* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    std::complex<float>* yp = y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    std::complex<float> xbeta = beta;
    if (A.isconj() && A.iscm()) {
#ifdef CBLAS
      TMV_SWAP(m,n);
      BLASNAME(cgemv) (BLASRM BLASCH_CT,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
          BLAS1);
#else
      std::complex<float> ca = TMV_CONJ(alpha);
      if (x.isconj()) {
        y.ConjugateSelf();
        BLASNAME(cgemv) (BLASCM BLASCH_NT,
            BLASV(m),BLASV(n),BLASP(&ca),BLASP(A.cptr()),BLASV(lda),
            BLASP(xp),BLASV(xs),BLASP(&xbeta),BLASP(yp),BLASV(ys)
            BLAS1);
        y.ConjugateSelf();
      } else {
        Vector<std::complex<float> > xx = ca*x.Conjugate();
        ca = std::complex<float>(1.F);
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
  template <int Si, int Sj, bool C1>
  static void BlasMultMV(
      std::complex<float> alpha,
      const ConstMatrixView<std::complex<float>,Si,Sj,C1>& A,
      const ConstVectorView<float>& x,
      float beta, VectorView<std::complex<float> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    if (A.iscm()) {
      if (y.step() != 1) {
        Vector<std::complex<float> > yy(y.size());
        BlasMultMV(1.F,A,x,0.F,yy.View());
        if (beta == 0.F) y = alpha*yy;
        else y += alpha*yy;
      } else {
        if (beta == 0.F) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          if (lda < m) { TMVAssert(n==1); lda = m; }
          int xs = x.step();
          int ys = 1;
          const float* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          float* yp = (float*) y.ptr();
          float xalpha(1);
          BLASNAME(sgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((float*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(beta),
              BLASP(yp),BLASV(ys) BLAS1);
          if (A.isconj()) y.ConjugateSelf();
          y *= alpha;
        } else if (A.isconj()) {
          Vector<std::complex<float> > yy(y.size());
          BlasMultMV(1.F,A.Conjugate(),x,0.F,yy.View());
          y += alpha*yy.Conjugate();
        } else if (TMV_IMAG(alpha) == 0.F) {
          int m = 2*A.colsize();
          int n = A.rowsize();
          int lda = 2*A.stepj();
          if (lda < m) { TMVAssert(n==1); lda = m; }
          int xs = x.step();
          int ys = 1;
          const float* xp = x.cptr();
          if (xs < 0) xp += (x.size()-1)*xs;
          float* yp = (float*) y.ptr();
          if (ys < 0) yp += (y.size()-1)*ys;
          float xalpha(TMV_REAL(alpha));
          BLASNAME(sgemv) (BLASCM BLASCH_NT,
              BLASV(m),BLASV(n),BLASV(xalpha),
              BLASP((float*)A.cptr()),BLASV(lda),
              BLASP(xp),BLASV(xs),BLASV(beta),
              BLASP(yp),BLASV(ys) BLAS1);
        } else {
          Vector<std::complex<float> > yy(y.size());
          BlasMultMV(std::complex<float>(1.F),A,x,0.F,yy.View());
          y += alpha*yy;
        }
      }
    } else { // A.isrm
      const Vector<std::complex<float> > xx = x;
      BlasMultMV(alpha,A,xx.XView(),beta,y);
    }
  }
  template <int Si, int Sj, bool C2>
  static void BlasMultMV(  
      std::complex<float> alpha,
      const ConstMatrixView<float,Si,Sj>& A,
      const ConstVectorView<std::complex<float>,UNKNOWN,C2>& x,
      float beta, VectorView<std::complex<float> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    if (beta == 0.F) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      if (lda < m) { TMVAssert(n==1); lda = m; }
      int xs = 2*x.step();
      int ys = 2*y.step();
      if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
      if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
      const float* xp = (const float*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      float* yp = (float*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      float xalpha(1);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(beta),
          BLASP(yp+1),BLASV(ys) BLAS1);
      if (x.isconj()) y.ConjugateSelf();
      y *= alpha;
    } else if (TMV_IMAG(alpha) == 0.F && !x.isconj()) {
      int m = A.iscm() ? A.colsize() : A.rowsize();
      int n = A.iscm() ? A.rowsize() : A.colsize();
      int lda = A.iscm() ? A.stepj() : A.stepi();
      if (lda < m) { TMVAssert(n==1); lda = m; }
      int xs = 2*x.step();
      int ys = 2*y.step();
      if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
      if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
      const float* xp = (const float*) x.cptr();
      if (xs < 0) xp += (x.size()-1)*xs;
      float* yp = (float*) y.ptr();
      if (ys < 0) yp += (y.size()-1)*ys;
      float xalpha(TMV_REAL(alpha));
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp),BLASV(ys) BLAS1);
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(xalpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp+1),BLASV(xs),BLASV(beta),
          BLASP(yp+1),BLASV(ys) BLAS1);
    } else {
      const Vector<std::complex<float> > xx = alpha*x;
      BlasMultMV(std::complex<float>(1.F),A,xx.XView(),1.F,y);
    }
  }
  template <int Si, int Sj>
  static void BlasMultMV(
      std::complex<float> alpha,
      const ConstMatrixView<float,Si,Sj>& A,
      const ConstVectorView<float>& x,
      float beta, VectorView<std::complex<float> > y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(!SameStorage(x,y));

    int m = A.iscm() ? A.colsize() : A.rowsize();
    int n = A.iscm() ? A.rowsize() : A.colsize();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    if (lda < m) { TMVAssert(n==1); lda = m; }
    int xs = x.step();
    int ys = 2*y.step();
    if (xs == 0) { TMVAssert(x.size() == 1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size() == 1); ys = 1; }
    const float* xp = x.cptr();
    if (xs < 0) xp += (x.size()-1)*xs;
    float* yp = (float*) y.ptr();
    if (ys < 0) yp += (y.size()-1)*ys;
    float ar(TMV_REAL(alpha));
    float ai(TMV_IMAG(alpha));
    if (ar == 0.F) {
      if (beta == 0.F) y.Real().Zero();
    }
    else 
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ar),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp),BLASV(ys) BLAS1);
    if (ai == 0.F) {
      if (beta == 0.F) y.Imag().Zero();
    }
    else
      BLASNAME(sgemv) (BLASCM A.isrm()?BLASCH_T:BLASCH_NT,
          BLASV(m),BLASV(n),BLASV(ai),BLASP(A.cptr()),BLASV(lda),
          BLASP(xp),BLASV(xs),BLASV(beta),
          BLASP(yp+1),BLASV(ys) BLAS1);
  }
#endif // FLOAT
#endif // BLAS

} // namespace tmv


