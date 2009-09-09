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
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixIO.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_SwapU.h"
#include "tmv/TMV_NormU.h"

namespace tmv {

  const int XX = UNKNOWN;

  //
  // Copy Matrices
  //

  template <class T, DiagType D1, bool C, DiagType D2>
  void InstCopy(const ConstUpperTriMatrixView<T,D1,XX,XX,C>& m1,
      UpperTriMatrixView<T,D2> m2)
  {
    if (m2.isrm())
    {
      UpperTriMatrixView<T,D2,XX,1> m2rm = m2;
      if (m1.isrm())
        InlineCopy(m1.RMView(),m2rm);
      else if (m1.iscm())
        InlineCopy(m1.RMView(),m2rm);
      else
        InlineCopy(m1,m2rm);
    }
    else if (m2.iscm())
    {
      UpperTriMatrixView<T,D2,1> m2cm = m2;
      if (m1.isrm())
        InlineCopy(m1.RMView(),m2cm);
      else if (m1.iscm())
        InlineCopy(m1.CMView(),m2cm);
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


  //
  // Swap Matrices
  //

  template <class T, DiagType D, bool C> 
  void InstSwap(UpperTriMatrixView<T,D,XX,XX,C> m1,
      UpperTriMatrixView<T,D> m2)
  {
    if (m2.isrm())
    {
      UpperTriMatrixView<T,D,XX,1> m2rm = m2;
      if (m1.isrm())
      {
        UpperTriMatrixView<T,D,XX,1,C> m1rm = m1;
        InlineSwap(m1rm,m2rm);
      }
      else if (m1.iscm())
      {
        UpperTriMatrixView<T,D,1,XX,C> m1cm = m1;
        InlineSwap(m1cm,m2rm);
      }
      else
        InlineSwap(m1,m2rm);
    }
    else if (m2.iscm())
    {
      UpperTriMatrixView<T,D,1> m2cm = m2;
      if (m1.isrm())
      {
        UpperTriMatrixView<T,D,XX,1,C> m1rm = m1;
        InlineSwap(m1rm,m2cm);
      }
      else if (m1.iscm())
      {
        UpperTriMatrixView<T,D,1,XX,C> m1cm = m1;
        InlineSwap(m1cm,m2cm);
      }
      else
        InlineSwap(m1,m2cm);
    }
    else
    {
      if (m1.isrm())
      {
        UpperTriMatrixView<T,D,XX,1,C> m1rm = m1;
        InlineSwap(m1rm,m2);
      }
      else if (m1.iscm())
      {
        UpperTriMatrixView<T,D,1,XX,C> m1cm = m1;
        InlineSwap(m1cm,m2);
      }
      else
        InlineSwap(m1,m2);
    }
  }


  //
  // Norms, SumElements
  //
  
  template <class T, DiagType D>
  T DoInstSumElements(const ConstUpperTriMatrixView<T,D>& m)
  { 
    if (m.isrm()) return InlineSumElements(m.RMView());
    else if (m.iscm()) return InlineSumElements(m.CMView());
    else return InlineSumElements(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstSumAbsElements(const ConstUpperTriMatrixView<T,D>& m)
  {
    if (m.isrm()) return InlineSumAbsElements(m.RMView());
    else if (m.iscm()) return InlineSumAbsElements(m.CMView());
    else return InlineSumAbsElements(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstNormSq(const ConstUpperTriMatrixView<T,D>& m)
  {
    if (m.isrm()) return InlineNormSq(m.RMView());
    else if (m.iscm()) return InlineNormSq(m.CMView());
    else return InlineNormSq(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstNormSq(const ConstUpperTriMatrixView<T,D>& m,
      const RealType(T) scale)
  {
    if (m.isrm()) return InlineNormSq(m.RMView());
    else if (m.iscm()) return InlineNormSq(m.CMView());
    else return InlineNormSq(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstNormF(const ConstUpperTriMatrixView<T,D>& m)
  {
    if (m.isrm()) return InlineNormF(m.RMView());
    else if (m.iscm()) return InlineNormF(m.CMView());
    else return InlineNormF(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstMaxAbsElement(const ConstUpperTriMatrixView<T,D>& m)
  {
    if (m.isrm()) return InlineMaxAbsElement(m.RMView());
    else if (m.iscm()) return InlineMaxAbsElement(m.CMView());
    else return InlineMaxAbsElement(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstNorm1(const ConstUpperTriMatrixView<T,D>& m)
  {
    if (m.isrm()) return InlineNorm1(m.RMView());
    else if (m.iscm()) return InlineNorm1(m.CMView());
    else return InlineNorm1(m);
  }

  template <class T, DiagType D>
  RealType(T) DoInstNormInf(const ConstUpperTriMatrixView<T,D>& m)
  {
    if (m.isrm()) return InlineNormInf(m.RMView());
    else if (m.iscm()) return InlineNormInf(m.CMView());
    else return InlineNormInf(m);
  }

#ifdef XLAP
#ifdef TMV_INST_DOUBLE
  template <DiagType D>
  double DoInstNormF(const ConstUpperTriMatrixView<double,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
    char c = 'F';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(dlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  double DoInstNormF(
      const ConstUpperTriMatrixView<std::complex<double>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
    char c = 'F';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(zlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }

  template <DiagType D>
  double DoInstMaxAbsElement(const ConstUpperTriMatrixView<double,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
    char c = 'M';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(dlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  double DoInstMaxAbsElement(
      const ConstUpperTriMatrixView<std::complex<double>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
    char c = 'M';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    double work;
    double norm = LAPNAME(zlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }

  template <DiagType D>
  double DoInstNorm1(const ConstUpperTriMatrixView<double,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
    char c = m.isrm() ? 'I' : '1';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAME(dlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  double DoInstNorm1(
      const ConstUpperTriMatrixView<std::complex<double>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
    char c = m.isrm() ? 'I' : '1';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAME(zlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }

  template <DiagType D>
  double DoInstNormInf(const ConstUpperTriMatrixView<double,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstNormInf(m.copy());
    char c = m.isrm() ? '1' : 'I';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAME(dlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  double DoInstNormInf(
      const ConstUpperTriMatrixView<std::complex<double>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstNormInf(m.copy());
    char c = m.isrm() ? '1' : 'I';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAME(zlantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
  template <DiagType D>
  float DoInstNormF(const ConstUpperTriMatrixView<float,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
    char c = 'F';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(slantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  float DoInstNormF(
      const ConstUpperTriMatrixView<std::complex<float>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
    char c = 'F';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(clantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }

  template <DiagType D>
  float DoInstMaxAbsElement(const ConstUpperTriMatrixView<float,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
    char c = 'M';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(slantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  float DoInstMaxAbsElement(
      const ConstUpperTriMatrixView<std::complex<float>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
    char c = 'M';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
    float work;
    float norm = LAPNAME(clantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
        LAP1 LAP1 LAP1);
    return norm;
  }

  template <DiagType D>
  float DoInstNorm1(const ConstUpperTriMatrixView<float,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
    char c = m.isrm() ? 'I' : '1';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAME(slantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  float DoInstNorm1(
      const ConstUpperTriMatrixView<std::complex<float>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
    char c = m.isrm() ? 'I' : '1';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAME(clantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }

  template <DiagType D>
  float DoInstNormInf(const ConstUpperTriMatrixView<float,D>& m)
  { 
    if (!m.iscm() && !m.isrm()) return DoInstNormInf(m.copy());
    char c = m.isrm() ? '1' : 'I';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAME(slantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }
  template <DiagType D>
  float DoInstNormInf(
      const ConstUpperTriMatrixView<std::complex<float>,D>& m)
  {
    if (!m.iscm() && !m.isrm()) return DoInstNormInf(m.copy());
    char c = m.isrm() ? '1' : 'I';
    int N = m.size();
    int M = N;
    int lda = m.iscm() ? m.stepj() : m.stepi();
    TMVAssert(lda >= m);
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAME(clantr) (LAPCM LAPV(c),
        m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
        LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
        LAP1 LAP1 LAP1);
    return norm;
  }
#endif // FLOAT
#endif // XLAP

  template <class T, DiagType D>
  T InstSumElements(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstSumElements(m); }

  template <class T, DiagType D>
  RealType(T) InstSumAbsElements(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstSumAbsElements(m); }

  template <class T, DiagType D>
  RealType(T) InstNormSq(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstNormSq(m); }

  template <class T, DiagType D>
  RealType(T) InstNormSq(const ConstUpperTriMatrixView<T,D>& m,
      const RealType(T) scale)
  { return DoInstNormSq(m,scale); }

  template <class T, DiagType D>
  RealType(T) InstNormF(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstNormF(m); }

  template <class T, DiagType D>
  RealType(T) InstMaxAbsElement(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstMaxAbsElement(m); }

  template <class T, DiagType D>
  RealType(T) InstNorm1(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstNorm1(m); }

  template <class T, DiagType D>
  RealType(T) InstNormInf(const ConstUpperTriMatrixView<T,D>& m)
  { return DoInstNormInf(m); }

#if 0
  template <class T> 
  RT GenUpperTriMatrix<T>::DoNorm2() const
  { return Matrix<T>(*this).DoNorm2(); }

  template <class T> 
  RT GenUpperTriMatrix<T>::DoCondition() const
  { return Matrix<T>(*this).DoCondition(); }

  template <class T> 
  QuotXU<T,T> GenUpperTriMatrix<T>::QInverse() const
  { return QuotXU<T,T>(T(1),*this); }
#endif

  //
  // I/O
  //

  template <class T, DiagType D, bool C>
  void InstWriteCompact(std::ostream& os,
      const ConstUpperTriMatrixView<T,D,XX,XX,C>& m)
  {
    if (m.isrm()) InlineWriteCompact(os,m.RMView());
    else if (m.iscm()) InlineWriteCompact(os,m.CMView());
    else InlineWriteCompact(os,m);
  }

  template <class T, DiagType D, bool C>
  void InstWriteCompact(std::ostream& os,
      const ConstLowerTriMatrixView<T,D,XX,XX,C>& m)
  {
    if (m.isrm()) InlineWriteCompact(os,m.RMView());
    else if (m.iscm()) InlineWriteCompact(os,m.CMView());
    else InlineWriteCompact(os,m);
  }

  template <class T, DiagType D, bool C>
  void InstWriteCompact(std::ostream& os,
      const ConstUpperTriMatrixView<T,D,XX,XX,C>& m, RealType(T) thresh)
  {
    if (m.isrm()) InlineWriteCompact(os,m.RMView(),thresh);
    else if (m.iscm()) InlineWriteCompact(os,m.CMView(),thresh);
    else InlineWriteCompact(os,m,thresh);
  }
 
  template <class T, DiagType D, bool C>
  void InstWriteCompact(std::ostream& os,
      const ConstLowerTriMatrixView<T,D,XX,XX,C>& m, RealType(T) thresh)
  {
    if (m.isrm()) InlineWriteCompact(os,m.RMView(),thresh);
    else if (m.iscm()) InlineWriteCompact(os,m.CMView(),thresh);
    else InlineWriteCompact(os,m,thresh);
  }

  template <class T, DiagType D, bool C>
  void InstRead(std::istream& is, UpperTriMatrixView<T,D,XX,XX,C> m)
  {
    if (m.isrm()) 
    {
      UpperTriMatrixView<T,D,XX,1,C> mrm = m.RMView();
      InlineRead(is,mrm);
    }
    else if (m.iscm()) 
    {
      UpperTriMatrixView<T,D,1,XX,C> mcm = m.CMView();
      InlineRead(is,mcm);
    }
    else InlineRead(is,m);
  }
  template <class T, DiagType D, bool C>
  void InstRead(std::istream& is, LowerTriMatrixView<T,D,XX,XX,C> m)
  {
    if (m.isrm()) 
    {
      LowerTriMatrixView<T,D,XX,1,C> mrm = m.RMView();
      InlineRead(is,mrm);
    }
    else if (m.iscm()) 
    {
      LowerTriMatrixView<T,D,1,XX,C> mcm = m.CMView();
      InlineRead(is,mcm);
    }
    else InlineRead(is,m);
  }

#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


