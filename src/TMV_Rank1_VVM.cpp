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


#include "TMV_Blas.h"
#include "tmv/TMV_Rank1_VVM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMM.h"

namespace tmv {

  // 
  // Rank1Update
  //

  template <bool add, class T, class V1, class V2, class M3>
  static void DoRank1Update(const T x, const V1& v1, const V2& v2, M3& m3)
  {
    // Check for non-unit step and x != 1, and do the necessary copies here,
    // rather than in the InlineRank1Update function.  
    // This is faster to compile, since it keeps the InlineRank1Update
    // algo path to the ones that have vstep == 1.

    typedef RealType(T) RT;
    const Scaling<1,RT> one;

    if (x == RT(0))
    { Maybe<!add>::zero(m3); }
    else if (v2.step() != 1)
    {
      typedef typename Traits2<T,typename V2::value_type>::type  PT2;
      Vector<PT2> v2c = x*v2;
      if (v1.step() != 1)
      {
        typedef typename V1::value_type T1;
        Vector<T1> v1c = v1;
        InlineRank1Update<add>(one,v1c.UnitView(),v2c.UnitView(),m3);
      }
      else
      {
        InlineRank1Update<add>(one,v1.UnitView(),v2c.UnitView(),m3);
      }
    }
    else if (v1.step() != 1)
    {
      typedef typename Traits2<T,typename V1::value_type>::type  PT1;
      Vector<PT1> v1c = x*v1;
      InlineRank1Update<add>(one,v1c.UnitView(),v2.UnitView(),m3);
    }
    else
    {
      if (x == RT(1))
        InlineRank1Update<add>(one,v1.UnitView(),v2.UnitView(),m3);
      else if (x == RT(-1))
        InlineRank1Update<add>(Scaling<-1,RT>(),v1.UnitView(),v2.UnitView(),m3);
      else if (TMV_IMAG(x) == RT(0))
        InlineRank1Update<add>(Scaling<0,RT>(TMV_REAL(x)),
            v1.UnitView(),v2.UnitView(),m3);
      else
        InlineRank1Update<add>(Scaling<0,T>(x),v1.UnitView(),v2.UnitView(),m3);
    }
  }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
  static void BlasRank1Update(
      double alpha, const ConstVectorView<double>& x,
      const ConstVectorView<double>& y,
      MatrixView<double,1> A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
  template <bool C1, bool C2>
  static void BlasRank1Update(
      std::complex<double> alpha,
      const ConstVectorView<std::complex<double>,UNKNOWN,C1>& x, 
      const ConstVectorView<std::complex<double>,UNKNOWN,C2>& y,
      MatrixView<std::complex<double>,1> A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    if (C1 && C2) {
      const Vector<std::complex<double> > xx = alpha * x;
      return BlasRank1Update(std::complex<double>(1.),xx.XView(),y,A);
    }
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
  template <bool C2>
  static void BlasRank1Update(
      std::complex<double> alpha,
      const ConstVectorView<double>& x, 
      const ConstVectorView<std::complex<double>,UNKNOWN,C2>& y,
      MatrixView<std::complex<double>,1> A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ (yr + I yi)
    // Ar += ar * x ^ yr - ai * x ^ yi
    // Ai += ai * x ^ yr + ar * x ^ yi
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
    xx.SubVector(0,xx.size(),2) = real(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = imag(alpha)*x;
    const double* xp = xx.cptr();
    double xalpha(1);
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((double*)A.ptr()),BLASV(lda));
    if (y.isconj()) {
      xx.SubVector(0,xx.size(),2) = imag(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = -real(alpha)*x;
    } else {
      xx.SubVector(0,xx.size(),2) = -imag(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = real(alpha)*x;
    }
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp+1),BLASV(ys),
        BLASP((double*)A.ptr()),BLASV(lda));
  }
  template <bool C1>
  static void BlasRank1Update(
      std::complex<double> alpha,
      const ConstVectorView<std::complex<double>,UNKNOWN,C1>& x,
      const ConstVectorView<double>& y, 
      MatrixView<std::complex<double>,1> A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
    if (x.step() == 1 && !x.isconj() && imag(alpha) == 0.) {
      const double* xp = (double*) x.cptr();
      double xalpha(real(alpha));
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
  static void BlasRank1Update(
      std::complex<double> alpha,
      const ConstVectorView<double>& x, const ConstVectorView<double>& y, 
      MatrixView<std::complex<double>,1> A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ y
    // Ar += ar * x ^ y
    // Ai += ai * x ^ y
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
    xx.SubVector(0,xx.size(),2) = real(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = imag(alpha)*x;
    const double* xp = xx.cptr();
    double xalpha(1);
    BLASNAME(dger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((double*)A.ptr()),BLASV(lda));
  }
  template <bool add>
  static void DoRank1Update(const double x,
      const ConstVectorView<double>& v1,
      const ConstVectorView<double>& v2, MatrixView<double,1> m3)
  { 
    if (!add) m3.Zero();
    BlasRank1Update(x,v1,v2,m3); 
  }
  template <bool add, class T1, bool C1, class T2, bool C2>
  static void DoRank1Update(const std::complex<double> x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2,
      MatrixView<std::complex<double>,1> m3)
  {
    if (!add) m3.Zero();
    BlasRank1Update(x,v1,v2,m3); 
  }
#endif //DOUBLE
#ifdef TMV_INST_FLOAT
  static void BlasRank1Update(
      float alpha, const ConstVectorView<float>& x,
      const ConstVectorView<float>& y,
      MatrixView<float,1> A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
  template <bool C1, bool C2>
  static void BlasRank1Update(
      std::complex<float> alpha,
      const ConstVectorView<std::complex<float>,UNKNOWN,C1>& x, 
      const ConstVectorView<std::complex<float>,UNKNOWN,C2>& y,
      MatrixView<std::complex<float>,1> A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    if (C1 && C2) {
      const Vector<std::complex<float> > xx = alpha * x;
      return BlasRank1Update(std::complex<float>(1.),xx.XView(),y,A);
    }
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
  template <bool C2>
  static void BlasRank1Update(
      std::complex<float> alpha,
      const ConstVectorView<float>& x, 
      const ConstVectorView<std::complex<float>,UNKNOWN,C2>& y,
      MatrixView<std::complex<float>,1> A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ (yr + I yi)
    // Ar += ar * x ^ yr - ai * x ^ yi
    // Ai += ai * x ^ yr + ar * x ^ yi
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    int m = 2*A.colsize();
    int n = A.rowsize();
    int xs = 1;
    int ys = 2*y.step();
    const float* yp = (const float*) y.cptr();
    if (ys < 0) yp += (n-1)*ys;
    int lda = 2*A.stepj();
    if (xs == 0) { TMVAssert(x.size()==1); xs = 1; }
    if (ys == 0) { TMVAssert(y.size()==1); ys = 1; }
    if (lda < m) { TMVAssert(n == 1); lda = m; }
    Vector<float> xx(2*x.size());
    xx.SubVector(0,xx.size(),2) = real(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = imag(alpha)*x;
    const float* xp = xx.cptr();
    float xalpha(1);
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((float*)A.ptr()),BLASV(lda));
    if (y.isconj()) {
      xx.SubVector(0,xx.size(),2) = imag(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = -real(alpha)*x;
    } else {
      xx.SubVector(0,xx.size(),2) = -imag(alpha)*x;
      xx.SubVector(1,xx.size()+1,2) = real(alpha)*x;
    }
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp+1),BLASV(ys),
        BLASP((float*)A.ptr()),BLASV(lda));
  }
  template <bool C1>
  static void BlasRank1Update(
      std::complex<float> alpha,
      const ConstVectorView<std::complex<float>,UNKNOWN,C1>& x,
      const ConstVectorView<float>& y, 
      MatrixView<std::complex<float>,1> A)
  {
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
    if (x.step() == 1 && !x.isconj() && imag(alpha) == 0.) {
      const float* xp = (float*) x.cptr();
      float xalpha(real(alpha));
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
  static void BlasRank1Update(
      std::complex<float> alpha,
      const ConstVectorView<float>& x, const ConstVectorView<float>& y, 
      MatrixView<std::complex<float>,1> A)
  {
    // A += a * x ^ y
    // (Ar + I Ai) += (ar + I ai) * x ^ y
    // Ar += ar * x ^ y
    // Ai += ai * x ^ y
    TMVAssert(A.colsize() == x.size());
    TMVAssert(A.rowsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

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
    xx.SubVector(0,xx.size(),2) = real(alpha)*x;
    xx.SubVector(1,xx.size()+1,2) = imag(alpha)*x;
    const float* xp = xx.cptr();
    float xalpha(1);
    BLASNAME(sger) (BLASCM BLASV(m),BLASV(n),BLASV(xalpha),
        BLASP(xp),BLASV(xs),BLASP(yp),BLASV(ys),
        BLASP((float*)A.ptr()),BLASV(lda));
  }
  template <bool add>
  static void DoRank1Update(const float x,
      const ConstVectorView<float>& v1,
      const ConstVectorView<float>& v2, MatrixView<float,1> m3)
  { 
    if (!add) m3.Zero();
    BlasRank1Update(x,v1,v2,m3); 
  }
  template <bool add, class T1, bool C1, class T2, bool C2>
  static void DoRank1Update(const std::complex<float> x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2,
      MatrixView<std::complex<float>,1> m3)
  {
    if (!add) m3.Zero();
    BlasRank1Update(x,v1,v2,m3); 
  }
#endif // FLOAT
#endif // BLAS

  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstRank1Update(const T3 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2, MatrixView<T3> m3)
  {
    if (m3.isrm())
    {
      MatrixView<T3,1> m3t = m3.Transpose();
      DoRank1Update<false>(x,v2,v1,m3t);
    }
    else if (m3.iscm())
    {
      MatrixView<T3,1> m3cm = m3;
      DoRank1Update<false>(x,v1,v2,m3cm);
    }
    else 
    {
      Matrix<T3,ColMajor> m3x(m3.colsize(),m3.rowsize());
      MatrixView<T3,1> m3cm = m3x.CMView();
      DoRank1Update<false>(x,v1,v2,m3cm);
      m3 = m3x;
    }
  }
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddRank1Update(const T3 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2, MatrixView<T3> m3)
  {
    if (m3.isrm())
    {
      MatrixView<T3,1> m3t = m3.Transpose();
      DoRank1Update<true>(x,v2,v1,m3t);
    }
    else if (m3.iscm())
    {
      MatrixView<T3,1> m3cm = m3;
      DoRank1Update<true>(x,v1,v2,m3cm);
    }
    else 
    {
      Matrix<T3,ColMajor> m3x(m3.colsize(),m3.rowsize());
      MatrixView<T3,1> m3cm = m3x.CMView();
      DoRank1Update<false>(RealType(T3)(1),v1,v2,m3cm);
      m3 += x * m3x;
    }
  }


#define InstFile "TMV_Rank1_VVM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


