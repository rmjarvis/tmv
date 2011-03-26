///////////////////////////////////////////////////////////////////////////////
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

//#define PRINTALGO_NormU

#include "TMV_Blas.h"
#include <iostream>
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixIO.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_SwapU.h"
#include "tmv/TMV_NormU.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_ConjugateV.h"

namespace tmv {

    //
    // Copy Matrices
    //

    template <class M1, class M2>
    void DoCopyU(const M1& m1, M2& m2)
    {
        if (m1.iscm() && m2.iscm()) {
            typename M2::cmview_type m2cm = m2.cmView();
            InlineCopy(m1.cmView(),m2cm);
        } else if (m1.isrm() && m2.isrm()) {
            typename M2::rmview_type m2rm = m2.rmView();
            InlineCopy(m1.rmView(),m2rm);
        } else {
            InlineCopy(m1,m2);
        }
    }

    template <class T1, int C1, class T2>
    void InstCopy(
        const ConstUpperTriMatrixView<T1,C1>& m1, UpperTriMatrixView<T2> m2)
    {
        if (m1.isunit()) {
            if (m2.size() > 1) {
                UpperTriMatrixView<T2,NonUnitDiag> m2o = m2.offDiag();
                DoCopyU(m1.nonConj().offDiag(),m2o);
                Maybe<C1>::conjself(m2o);
            }
            if (!m2.isunit()) m2.diag().setAllTo(T2(1));
        } else {
            UpperTriMatrixView<T2,NonUnitDiag> m2n = m2.viewAsNonUnitDiag();
            TMVAssert(!m2.isunit());
            DoCopyU(m1.nonConj().viewAsNonUnitDiag(),m2n);
            Maybe<C1>::conjself(m2n);
        }
    }

    template <class T1, int C1, class T2>
    void InstAliasCopy(
        const ConstUpperTriMatrixView<T1,C1>& m1, UpperTriMatrixView<T2> m2)
    { InlineAliasCopy(m1,m2); }


    //
    // Swap Matrices
    //

    template <class T, int C> 
    void InstSwap(
        UpperTriMatrixView<T,C> m1, UpperTriMatrixView<T> m2)
    {
        if (m1.iscm() && m2.iscm()) {
            UpperTriMatrixView<T,C|ColMajor> m1cm = m1;
            UpperTriMatrixView<T,ColMajor> m2cm = m2;
            InlineSwap(m1cm,m2cm);
        } else if (m1.isrm() && m2.isrm()) {
            UpperTriMatrixView<T,C|RowMajor> m1rm = m1;
            UpperTriMatrixView<T,ColMajor> m2rm = m2;
            InlineSwap(m1rm,m2rm);
        } else {
            InlineSwap(m1,m2);
        }
    }

    template <class T, int C> 
    void InstAliasSwap(
        UpperTriMatrixView<T,C> m1, UpperTriMatrixView<T> m2)
    { InlineAliasSwap(m1,m2); }


    //
    // Norms, SumElements
    //

    template <class T>
    static T DoInstSumElements(const ConstUpperTriMatrixView<T>& m)
    { 
        if (m.iscm()) return InlineSumElements(m.cmView());
        else if (m.isrm()) return InlineSumElements(m.rmView());
        else return InlineSumElements(m);
    }

    template <class T>
    static typename ConstUpperTriMatrixView<T>::float_type DoInstSumAbsElements(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumAbsElements(m.cmView());
        else if (m.isrm()) return InlineSumAbsElements(m.rmView());
        else return InlineSumAbsElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbs2Elements(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumAbs2Elements(m.cmView());
        else if (m.isrm()) return InlineSumAbs2Elements(m.rmView());
        else return InlineSumAbs2Elements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormSq(m.cmView());
        else if (m.isrm()) return InlineNormSq(m.rmView());
        else return InlineNormSq(m);
    }

    template <class T>
    static typename ConstUpperTriMatrixView<T>::float_type DoInstNormSq(
        const ConstUpperTriMatrixView<T>& m,
        const typename ConstUpperTriMatrixView<T>::float_type scale)
    {
        if (m.iscm()) return InlineNormSq(m.cmView(),scale);
        else if (m.isrm()) return InlineNormSq(m.rmView(),scale);
        else return InlineNormSq(m,scale);
    }

    template <class T>
    static typename ConstUpperTriMatrixView<T>::float_type DoInstNormF(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormF(m.cmView());
        else if (m.isrm()) return InlineNormF(m.rmView());
        else return InlineNormF(m);
    }

    template <class T>
    static typename ConstUpperTriMatrixView<T>::float_type DoInstMaxAbsElement(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineMaxAbsElement(m.cmView());
        else if (m.isrm()) return InlineMaxAbsElement(m.rmView());
        else return InlineMaxAbsElement(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstMaxAbs2Element(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineMaxAbs2Element(m.cmView());
        else if (m.isrm()) return InlineMaxAbs2Element(m.rmView());
        else return InlineMaxAbs2Element(m);
    }

    template <class T>
    static typename ConstUpperTriMatrixView<T>::float_type DoInstNorm1(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNorm1(m.cmView());
        else if (m.isrm()) return InlineNorm1(m.rmView());
        else return InlineNorm1(m);
    }

    template <class T>
    static typename ConstLowerTriMatrixView<T>::float_type DoInstNorm1(
        const ConstLowerTriMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNorm1(m.cmView());
        else if (m.isrm()) return InlineNorm1(m.rmView());
        else return InlineNorm1(m);
    }

#ifdef XLAP
#ifdef TMV_INST_DOUBLE
    static double DoInstNormF(const ConstUpperTriMatrixView<double>& m)
    { 
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        char c = 'F';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(dlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }
    static double DoInstNormF(
        const ConstUpperTriMatrixView<std::complex<double> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        char c = 'F';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(zlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }

    static double DoInstMaxAbsElement(const ConstUpperTriMatrixView<double>& m)
    { 
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        char c = 'M';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(dlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }
    static double DoInstMaxAbsElement(
        const ConstUpperTriMatrixView<std::complex<double> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        char c = 'M';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(zlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }

    static double DoInstMaxAbs2Element(const ConstUpperTriMatrixView<double>& m)
    { return DoInstMaxAbsElement(m); }
    static double DoInstMaxAbs2Element(
        const ConstUpperTriMatrixView<std::complex<double> >& m)
    { return InlineMaxAbs2Element(m); }

    static double DoInstNorm1(const ConstUpperTriMatrixView<double>& m)
    { 
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? '1' : 'I';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(dlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
    static double DoInstNorm1(
        const ConstUpperTriMatrixView<std::complex<double> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? '1' : 'I';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(zlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static float DoInstNormF(const ConstUpperTriMatrixView<float>& m)
    { 
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        char c = 'F';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(slantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }
    static float DoInstNormF(
        const ConstUpperTriMatrixView<std::complex<float> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        char c = 'F';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(clantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }

    static float DoInstMaxAbsElement(const ConstUpperTriMatrixView<float>& m)
    { 
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        char c = 'M';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(slantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }
    static float DoInstMaxAbsElement(
        const ConstUpperTriMatrixView<std::complex<float> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        char c = 'M';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(clantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(&work)
            LAP1 LAP1 LAP1);
        return norm;
    }

    static float DoInstMaxAbs2Element(const ConstUpperTriMatrixView<float>& m)
    {  return DoInstMaxAbsElement(m); }
    static float DoInstMaxAbs2Element(
        const ConstUpperTriMatrixView<std::complex<float> >& m)
    { return InlineMaxAbs2Element(m); } 

    static float DoInstNorm1(const ConstUpperTriMatrixView<float>& m)
    { 
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? 'I' : '1';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(slantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
    static float DoInstNorm1(
        const ConstUpperTriMatrixView<std::complex<float> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? 'I' : '1';
        int N = m.size();
        int M = N;
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(clantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
#endif // FLOAT
#endif // XLAP

    template <class T>
    T InstSumElements(const ConstUpperTriMatrixView<T>& m)
    { return DoInstSumElements(m); }

    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstSumAbsElements(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstSumAbsElements(m); }

    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstSumAbs2Elements(m); }

    template <class T>
    typename Traits<T>::real_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstNormSq(m); } 

    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m,
        const typename ConstUpperTriMatrixView<T>::float_type scale)
    { return DoInstNormSq(m,scale); }

    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstNormF(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstNormF(m); }

    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstMaxAbsElement(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstMaxAbsElement(m); }

    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstMaxAbs2Element(m); }

    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstNorm1(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstNorm1(m); }

    template <class T>
    typename ConstLowerTriMatrixView<T>::float_type InstNorm1(
        const ConstLowerTriMatrixView<T>& m)
    { return DoInstNorm1(m); }

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

    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, 
        const ConstUpperTriMatrixView<T,C>& m)
    {
        if (m.iscm()) InlineWriteCompact(os,m.cmView());
        else if (m.isrm()) InlineWriteCompact(os,m.rmView());
        else InlineWriteCompact(os,m);
    }

    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, 
        const ConstLowerTriMatrixView<T,C>& m)
    {
        if (m.iscm()) InlineWriteCompact(os,m.cmView());
        else if (m.isrm()) InlineWriteCompact(os,m.rmView());
        else InlineWriteCompact(os,m);
    }

    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstUpperTriMatrixView<T,C>& m, 
        typename ConstUpperTriMatrixView<T>::float_type thresh)
    {
        if (m.iscm()) InlineWriteCompact(os,m.cmView(),thresh);
        else if (m.isrm()) InlineWriteCompact(os,m.rmView(),thresh);
        else InlineWriteCompact(os,m,thresh);
    }

    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstLowerTriMatrixView<T,C>& m, 
        typename ConstLowerTriMatrixView<T>::float_type thresh)
    {
        if (m.iscm()) InlineWriteCompact(os,m.cmView(),thresh);
        else if (m.isrm()) InlineWriteCompact(os,m.rmView(),thresh);
        else InlineWriteCompact(os,m,thresh);
    }

    template <class T>
    void InstRead(std::istream& is, UpperTriMatrixView<T> m)
    {
        if (m.iscm()) {
            UpperTriMatrixView<T,ColMajor> mcm = m.cmView();
            InlineRead(is,mcm);
        } else if (m.isrm()) {
            UpperTriMatrixView<T,RowMajor> mrm = m.rmView();
            InlineRead(is,mrm);
        } else 
            InlineRead(is,m);
    }
    template <class T>
    void InstRead(std::istream& is, LowerTriMatrixView<T> m)
    {
        if (m.iscm()) {
            LowerTriMatrixView<T,ColMajor> mcm = m.cmView();
            InlineRead(is,mcm);
        } else if (m.isrm()) {
            LowerTriMatrixView<T,RowMajor> mrm = m.rmView();
            InlineRead(is,mrm);
        } else 
            InlineRead(is,m);
    }

#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


