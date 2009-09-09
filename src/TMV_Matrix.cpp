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
#include <iostream>
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_SwapM.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

  const int XX = UNKNOWN;

#if 0 // Put this in TMV_DivHelper.cpp
  template <class T> 
  void GenMatrix<T>::NewDivider() const
  {
    switch (this->GetDivType()) {
      case LU : this->SetDiv(new LUDiv<T>(*this,
                      this->IsDivInPlace())); 
                break;
      case QR : this->SetDiv(new QRDiv<T>(*this,
                      this->IsDivInPlace())); 
                break;
      case QRP : this->SetDiv(new QRPDiv<T>(*this,
                       this->IsDivInPlace())); 
                 break;
      case SV : this->SetDiv(new SVDiv<T>(*this,
                      this->IsDivInPlace())); 
                break;
      default : TMVAssert(FALSE);
    }
  }
#endif

  //
  // Copy Matrices
  //

  template <class T, bool C>
  void InstCopy(const ConstMatrixView<T,XX,XX,C>& m1, MatrixView<T> m2)
  {
    if (m2.isrm() && !m2.iscm())
    {
      InstCopy(m1.Transpose(),m2.Transpose());
    }
    else if (m2.iscm())
    {
      MatrixView<T,1> m2cm = m2;
      if (m1.isrm())
        InlineCopy(m1.RMView(),m2cm);
      else if (m1.iscm())
      {
        if (m1.CanLinearize() && m2.CanLinearize()) 
        {
          VectorView<T> m2l = m2.LinearView();
          InstCopy(m1.LinearView().XView(),m2l);
        }
        else
          InlineCopy(m1.CMView(),m2cm);
      }
      else
        InlineCopy(m1,m2cm);
    }
    else 
    {
      if (m1.isrm())
        InlineCopy(m1.RMView(),m2);
      else if (m2.iscm())
        InlineCopy(m1.CMView(),m2);
      else
        InlineCopy(m1,m2);
    }
  }

#ifdef ELAP
#ifdef TMV_INST_DOUBLE
  template <> 
  void InstCopy(
      const ConstMatrixView<double>& m1, MatrixView<double> m2)
  {
    TMVAssert(m1.cptr() != m2.cptr());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    if (ld1 < 0 || ld2 < 0) DoCopy(m1,m2);
    TMVAssert(ld1 >= m);
    TMVAssert(ld2 >= m);
    LAPNAME(dlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> 
  void InstCopy(
      const ConstMatrixView<std::complex<double> >& m1,
      MatrixView<std::complex<double> > m2)
  {
    TMVAssert(m1.cptr() != m2.cptr());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    if (ld1 < 0 || ld2 < 0) DoCopy(m1,m2);
    TMVAssert(ld1 >= m);
    TMVAssert(ld2 >= m);
    LAPNAME(zlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#ifdef TMV_INST_FLOAT
  template <> 
  void InstCopy(
      const ConstMatrixView<float>& m1, MatrixView<float> m2)
  {
    TMVAssert(m1.cptr() != m2.cptr());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    if (ld1 < 0 || ld2 < 0) DoCopy(m1,m2);
    TMVAssert(ld1 >= m);
    TMVAssert(ld2 >= m);
    LAPNAME(slacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> 
  void InstCopy(
      const ConstMatrixView<std::complex<float> >& m1,
      MatrixView<std::complex<float> > m2)
  {
    TMVAssert(m1.cptr() != m2.cptr());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    if (ld1 < 0 || ld2 < 0) DoCopy(m1,m2);
    TMVAssert(ld1 >= m);
    TMVAssert(ld2 >= m);
    LAPNAME(clacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#endif

  //
  // Swap Matrices
  //

  template <class T, bool C> 
  void InstSwap(MatrixView<T,XX,XX,C> m1, MatrixView<T> m2)
  {
    if (m2.isrm() && !m2.iscm())
    {
      InstSwap(m1.Transpose(),m2.Transpose());
    }
    else if (m2.iscm())
    {
      MatrixView<T,1> m2cm = m2.CMView();
      if (m1.isrm())
      {
        MatrixView<T,XX,1,C> m1rm = m1.RMView();
        InlineSwap(m1rm,m2cm);
      }
      else if (m1.iscm())
      {
        if (m1.CanLinearize() && m2.CanLinearize()) 
        {
          VectorView<T,XX,C> m1l = m1.LinearView();
          VectorView<T> m2l = m2.LinearView();
          InstSwap(m1l,m2l);
        }
        else
        {
          MatrixView<T,1,XX,C> m1cm = m1.CMView();
          InlineSwap(m1cm,m2cm);
        }
      }
      else
        InlineSwap(m1,m2cm);
    }
    else 
    {
      if (m1.isrm())
      {
        MatrixView<T,XX,1,C> m1rm = m1.RMView();
        InlineSwap(m1rm,m2);
      }
      else if (m1.iscm())
      {
        MatrixView<T,1,XX,C> m1cm = m1.CMView();
        InlineSwap(m1cm,m2);
      }
      else
        InlineSwap(m1,m2);
    }
  }

  
  // 
  // TransposeSelf
  //

  template <class T>
  void InstTransposeSelf(MatrixView<T> m)
  {
    if (m.isrm() && !m.iscm())
    {
      InstTransposeSelf(m.Transpose());
    }
    else if (m.iscm())
    {
      MatrixView<T,1> mcm = m;
      InlineTransposeSelf(mcm);
    }
    else 
    {
      InlineTransposeSelf(m);
    }
  }


  //
  // PermuteRows
  //

  // There is a Lapack version of this for colmajor matrices called 
  // laswp.  I don't use it, though, because it requires p to be
  // 1-based, rather than 0-based.  So it would require copying 
  // p into new storage and incrementing each value by 1.
  // Note that big a deal, but since the native TMV code is basically
  // just as fast, and this function is not usually a big fraction of the 
  // time of any algorithm anyway, it doesn't seem worth the bother.
  template <class T>
  void InstPermuteRows(MatrixView<T> m,
      const int*const p, const int i1, const int i2)
  {
    if (m.isrm())
    {
      MatrixView<T,XX,1> mrm = m;
      InlinePermuteRows(mrm,p,i1,i2);
    }
    else if (m.iscm())
    {
      MatrixView<T,1> mcm = m;
      InlinePermuteRows(mcm,p,i1,i2);
    }
    else 
    {
      InlinePermuteRows(m,p,i1,i2);
    }
  }

  template <class T>
  void InstReversePermuteRows(MatrixView<T> m,
      const int*const p, const int i1, const int i2)
  {
    if (m.isrm())
    {
      MatrixView<T,XX,1> mrm = m;
      InlineReversePermuteRows(mrm,p,i1,i2);
    }
    else if (m.iscm())
    {
      MatrixView<T,1> mcm = m;
      InlineReversePermuteRows(mcm,p,i1,i2);
    }
    else 
    {
      InlineReversePermuteRows(m,p,i1,i2);
    }
  }
 



  //
  // Norms, SumElements
  //
  
  template <class T>
  T InstSumElements(const ConstMatrixView<T>& m)
  {
    if (m.CanLinearize()) return m.LinearView().SumElements();
    else if (m.isrm()) return InlineSumElements(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineSumElements(m.CMView());
    else return InlineSumElements(m);
  }

  template <class T>
  RealType(T) InstSumAbsElements(const ConstMatrixView<T>& m)
  {
    if (m.CanLinearize()) return m.LinearView().SumAbsElements();
    else if (m.isrm()) return InlineSumAbsElements(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineSumAbsElements(m.CMView());
    else return InlineSumAbsElements(m);
  }

  template <class T>
  RealType(T) InstNormSq(const ConstMatrixView<T>& m)
  {
    if (m.CanLinearize()) return m.LinearView().NormSq();
    else if (m.isrm()) return InlineNormSq(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineNormSq(m.CMView());
    else return InlineNormSq(m);
  }

  template <class T>
  RealType(T) InstNormSq(const ConstMatrixView<T>& m, const RealType(T) scale)
  {
    if (m.CanLinearize()) return m.LinearView().NormSq();
    else if (m.isrm()) return InlineNormSq(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineNormSq(m.CMView());
    else return InlineNormSq(m);
  }

  template <class T>
  RealType(T) InstNormF(const ConstMatrixView<T>& m)
  {
    if (m.CanLinearize()) return m.LinearView().Norm2();
    else if (m.isrm()) return InlineNormF(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineNormF(m.CMView());
    else return InlineNormF(m);
  }

  template <class T>
  RealType(T) InstMaxAbsElement(const ConstMatrixView<T>& m)
  {
    if (m.CanLinearize()) return m.LinearView().MaxAbsElement();
    else if (m.isrm()) return InlineMaxAbsElement(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineMaxAbsElement(m.CMView());
    else return InlineMaxAbsElement(m);
  }

  template <class T>
  RealType(T) InstNorm1(const ConstMatrixView<T>& m)
  {
    if (m.isrm()) return InlineNormInf(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineNorm1(m.CMView());
    else return InlineNorm1(m);
  }

  template <class T>
  RealType(T) InstNormInf(const ConstMatrixView<T>& m)
  {
    if (m.isrm()) return InlineNorm1(m.Transpose().CMView()); 
    else if (m.iscm()) return InlineNormInf(m.CMView());
    else return InlineNormInf(m);
  }

#ifdef XLAP
#ifdef TMV_INST_DOUBLE
  template <> 
  double InstNormF(const ConstMatrixView<double>& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormF(m.Transpose());
    if (!m.iscm()) return InlineNormF(m);
    char c = 'F';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(dlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }
  template <> 
  double InstNormF(
      const ConstMatrixView<std::complex<double> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormF(m.Transpose());
    if (!m.iscm()) return InlineNormF(m);
    char c = 'F';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(zlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }

  template <> 
  double InstMaxAbsElement(const ConstMatrixView<double>& m)
  {
    if (m.isrm() && !m.iscm()) return InstMaxAbsElement(m.Transpose());
    if (!m.iscm()) return InlineMaxAbsElement(m);
    char c = 'M';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(dlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }
  template <> 
  double InstMaxAbsElement(
      const ConstMatrixView<std::complex<double> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstMaxAbsElement(m.Transpose());
    if (!m.iscm()) return InlineMaxAbsElement(m);
    char c = 'M';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(zlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }

  template <> 
  double InstNorm1(const ConstMatrixView<double>& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormInf(m.Transpose());
    if (!m.iscm()) return InlineNorm1(m);
    char c = '1';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(dlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }
  template <> 
  double InstNorm1(
      const ConstMatrixView<std::complex<double> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormInf(m.Transpose());
    if (!m.iscm()) return InlineNorm1(m);
    char c = '1';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(zlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }

  template <> 
  double InstNormInf(const ConstMatrixView<double>& m)
  {
    if (m.isrm() && !m.iscm()) return InstNorm1(m.Transpose());
    if (!m.iscm()) return InlineNormInf(m);
    char c = 'I';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAME(dlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
    return norm;
  }
  template <> 
  double InstNormInf(
      const ConstMatrixView<std::complex<double> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstNorm1(m.Transpose());
    if (!m.iscm()) return InlineNormInf(m);
    char c = 'I';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<double> work(new double[M]);
#endif
    double norm = LAPNAME(zlange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
    return norm;
  }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
  template <> 
  float InstNormF(const ConstMatrixView<float>& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormF(m.Transpose());
    if (!m.iscm()) return InlineNormF(m);
    char c = 'F';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(slange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }
  template <> 
  float InstNormF(
      const ConstMatrixView<std::complex<float> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormF(m.Transpose());
    if (!m.iscm()) return InlineNormF(m);
    char c = 'F';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(clange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }

  template <> 
  float InstMaxAbsElement(const ConstMatrixView<float>& m)
  {
    if (m.isrm() && !m.iscm()) return InstMaxAbsElement(m.Transpose());
    if (!m.iscm()) return InlineMaxAbsElement(m);
    char c = 'M';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(slange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }
  template <> 
  float InstMaxAbsElement(
      const ConstMatrixView<std::complex<float> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstMaxAbsElement(m.Transpose());
    if (!m.iscm()) return InlineMaxAbsElement(m);
    char c = 'M';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(clange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }

  template <> 
  float InstNorm1(const ConstMatrixView<float>& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormInf(m.Transpose());
    if (!m.iscm()) return InlineNorm1(m);
    char c = '1';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(slange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }
  template <> 
  float InstNorm1(
      const ConstMatrixView<std::complex<float> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstNormInf(m.Transpose());
    if (!m.iscm()) return InlineNorm1(m);
    char c = '1';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(clange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
    return norm;
  }

  template <> 
  float InstNormInf(const ConstMatrixView<float>& m)
  {
    if (m.isrm() && !m.iscm()) return InstNorm1(m.Transpose());
    if (!m.iscm()) return InlineNormInf(m);
    char c = 'I';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAME(slange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
    return norm;
  }
  template <> 
  float InstNormInf(
      const ConstMatrixView<std::complex<float> >& m)
  {
    if (m.isrm() && !m.iscm()) return InstNorm1(m.Transpose());
    if (!m.iscm()) return InlineNormInf(m);
    char c = 'I';
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<float> work(new float[M]);
#endif
    float norm = LAPNAME(clange) (LAPCM LAPV(c),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
    return norm;
  }
#endif // FLOAT
#endif // XLAP

#if 0
  template <class T> 
  RT GenMatrix<T>::DoNorm2() const
  {
    if (this->colsize() < this->rowsize()) return Transpose().DoNorm2();
    if (this->rowsize() == 0) return RT(0);
    Matrix<T> m = *this;
    DiagMatrix<RT> S(this->rowsize());
    SV_Decompose(m.View(),S.View(),false);
    return S(0);
  }

  template <class T> 
  RT GenMatrix<T>::DoCondition() const
  {
    if (this->colsize() < this->rowsize()) return Transpose().DoNorm2();
    if (this->rowsize() == 0) return RT(1);
    Matrix<T> m = *this;
    DiagMatrix<RT> S(this->rowsize());
    SV_Decompose(m.View(),S.View(),false);
    return S(0)/S(S.size()-1);
  }

  template <class T> 
  QuotXM<T,T> GenMatrix<T>::QInverse() const
  { return QuotXM<T,T>(T(1),*this); }
#endif

  //
  // I/O
  //

  template <class T, bool C>
  void InstWrite(std::ostream& os, const ConstMatrixView<T,XX,XX,C>& m)
  {
    if (m.isrm()) InlineWrite(os,m.RMView());
    else if (m.iscm()) InlineWrite(os,m.CMView());
    else InlineWrite(os,m);
  }

  template <class T, bool C>
  void InstWrite(std::ostream& os, const ConstMatrixView<T,XX,XX,C>& m,
      RealType(T) thresh)
  {
    if (m.isrm()) InlineWrite(os,m.RMView(),thresh);
    else if (m.iscm()) InlineWrite(os,m.CMView(),thresh);
    else InlineWrite(os,m,thresh);
  }

  template <class T, bool C>
  void InstRead(std::istream& is, MatrixView<T,XX,XX,C> m)
  {
    if (m.isrm()) 
    {
      MatrixView<T,XX,1,C> mrm = m.RMView();
      InlineRead(is,mrm);
    }
    else if (m.iscm()) 
    {
      MatrixView<T,1,XX,C> mcm = m.CMView();
      InlineRead(is,mcm);
    }
    else InlineRead(is,m);
  }

#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


