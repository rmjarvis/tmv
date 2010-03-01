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



#include "TMV_Blas.h"
#include <iostream>
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixIO.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_SwapU.h"
#include "tmv/TMV_NormU.h"

namespace tmv {

    //
    // Copy Matrices
    //

    template <class T1, bool C1, class T2>
    void InstCopy(
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        UpperTriMatrixView<T2,UnknownDiag> m2)
    {
        if (m2.isrm()) {
            UpperTriMatrixView<T2,UnknownDiag,UNKNOWN,1> m2rm = m2;
            if (m1.isrm())
                InlineCopy(m1.rmView(),m2rm);
            else if (m1.iscm())
                InlineCopy(m1.cmView(),m2rm);
            else
                InlineCopy(m1,m2rm);
        } else if (m2.iscm()) {
            UpperTriMatrixView<T2,UnknownDiag,1> m2cm = m2;
            if (m1.isrm())
                InlineCopy(m1.rmView(),m2cm);
            else if (m1.iscm())
                InlineCopy(m1.cmView(),m2cm);
            else
                InlineCopy(m1,m2cm);
        } else {
            if (m1.isrm())
                InlineCopy(m1.rmView(),m2);
            else if (m1.iscm())
                InlineCopy(m1.cmView(),m2);
            else
                InlineCopy(m1,m2);
        }
    }


    //
    // Swap Matrices
    //

    template <class T, bool C> 
    void InstSwap(
        UpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> m1,
        UpperTriMatrixView<T,UnknownDiag> m2)
    {
        if (m2.isrm()) {
            UpperTriMatrixView<T,UnknownDiag,UNKNOWN,1> m2rm = m2;
            if (m1.isrm()) {
                UpperTriMatrixView<T,UnknownDiag,UNKNOWN,1,C> m1rm = m1;
                InlineSwap(m1rm,m2rm);
            } else if (m1.iscm()) {
                UpperTriMatrixView<T,UnknownDiag,1,UNKNOWN,C> m1cm = m1;
                InlineSwap(m1cm,m2rm);
            } else
                InlineSwap(m1,m2rm);
        } else if (m2.iscm()) {
            UpperTriMatrixView<T,UnknownDiag,1> m2cm = m2;
            if (m1.isrm()) {
                UpperTriMatrixView<T,UnknownDiag,UNKNOWN,1,C> m1rm = m1;
                InlineSwap(m1rm,m2cm);
            } else if (m1.iscm()) {
                UpperTriMatrixView<T,UnknownDiag,1,UNKNOWN,C> m1cm = m1;
                InlineSwap(m1cm,m2cm);
            } else
                InlineSwap(m1,m2cm);
        } else {
            if (m1.isrm()) {
                UpperTriMatrixView<T,UnknownDiag,UNKNOWN,1,C> m1rm = m1;
                InlineSwap(m1rm,m2);
            } else if (m1.iscm()) {
                UpperTriMatrixView<T,UnknownDiag,1,UNKNOWN,C> m1cm = m1;
                InlineSwap(m1cm,m2);
            } else
                InlineSwap(m1,m2);
        }
    }


    //
    // Norms, SumElements
    //

    template <class T>
    static T DoInstSumElements(const ConstUpperTriMatrixView<T>& m)
    { 
        if (m.isrm()) return InlineSumElements(m.rmView());
        else if (m.iscm()) return InlineSumElements(m.cmView());
        else return InlineSumElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbsElements(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineSumAbsElements(m.rmView());
        else if (m.iscm()) return InlineSumAbsElements(m.cmView());
        else return InlineSumAbsElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbs2Elements(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineSumAbs2Elements(m.rmView());
        else if (m.iscm()) return InlineSumAbs2Elements(m.cmView());
        else return InlineSumAbs2Elements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineNormSq(m.rmView());
        else if (m.iscm()) return InlineNormSq(m.cmView());
        else return InlineNormSq(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(
        const ConstUpperTriMatrixView<T>& m,
        const typename Traits<T>::real_type scale)
    {
        if (m.isrm()) return InlineNormSq(m.rmView());
        else if (m.iscm()) return InlineNormSq(m.cmView());
        else return InlineNormSq(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormF(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineNormF(m.rmView());
        else if (m.iscm()) return InlineNormF(m.cmView());
        else return InlineNormF(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstMaxAbsElement(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineMaxAbsElement(m.rmView());
        else if (m.iscm()) return InlineMaxAbsElement(m.cmView());
        else return InlineMaxAbsElement(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNorm1(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineNorm1(m.rmView());
        else if (m.iscm()) return InlineNorm1(m.cmView());
        else return InlineNorm1(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormInf(
        const ConstUpperTriMatrixView<T>& m)
    {
        if (m.isrm()) return InlineNormInf(m.rmView());
        else if (m.iscm()) return InlineNormInf(m.cmView());
        else return InlineNormInf(m);
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

    static double DoInstNorm1(const ConstUpperTriMatrixView<double>& m)
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
        char c = m.isrm() ? 'I' : '1';
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

    static double DoInstNormInf(const ConstUpperTriMatrixView<double>& m)
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
        double norm = LAPNAME(dlantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
    static double DoInstNormInf(
        const ConstUpperTriMatrixView<std::complex<double> >& m)
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

    static float DoInstNorm1(const ConstUpperTriMatrixView<float>& m)
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
        char c = m.isrm() ? 'I' : '1';
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

    static float DoInstNormInf(const ConstUpperTriMatrixView<float>& m)
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
        float norm = LAPNAME(slantr) (
            LAPCM LAPV(c),
            m.iscm() ? LAPCH_UP : LAPCH_LO , m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N), LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
    static float DoInstNormInf(
        const ConstUpperTriMatrixView<std::complex<float> >& m)
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
    typename Traits<T>::real_type InstSumAbsElements(
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
    typename Traits<T>::real_type InstNormSq(
        const ConstUpperTriMatrixView<T>& m,
        const typename Traits<T>::real_type scale)
    { return DoInstNormSq(m,scale); }

    template <class T>
    typename Traits<T>::real_type InstNormF(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstNormF(m); }

    template <class T>
    typename Traits<T>::real_type InstMaxAbsElement(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstMaxAbsElement(m); }

    template <class T>
    typename Traits<T>::real_type InstNorm1(
        const ConstUpperTriMatrixView<T>& m)
    { return DoInstNorm1(m); }

    template <class T>
    typename Traits<T>::real_type InstNormInf(
        const ConstUpperTriMatrixView<T>& m)
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

    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os, 
        const ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m)
    {
        if (m.isrm()) InlineWriteCompact(os,m.rmView());
        else if (m.iscm()) InlineWriteCompact(os,m.cmView());
        else InlineWriteCompact(os,m);
    }

    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os, 
        const ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m)
    {
        if (m.isrm()) InlineWriteCompact(os,m.rmView());
        else if (m.iscm()) InlineWriteCompact(os,m.cmView());
        else InlineWriteCompact(os,m);
    }

    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m, 
        typename Traits<T>::real_type thresh)
    {
        if (m.isrm()) InlineWriteCompact(os,m.rmView(),thresh);
        else if (m.iscm()) InlineWriteCompact(os,m.cmView(),thresh);
        else InlineWriteCompact(os,m,thresh);
    }

    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m, 
        typename Traits<T>::real_type thresh)
    {
        if (m.isrm()) InlineWriteCompact(os,m.rmView(),thresh);
        else if (m.iscm()) InlineWriteCompact(os,m.cmView(),thresh);
        else InlineWriteCompact(os,m,thresh);
    }

    template <class T, bool C>
    void InstRead(
        std::istream& is, UpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> m)
    {
        if (m.isrm()) {
            UpperTriMatrixView<T,UnknownDiag,UNKNOWN,1,C> mrm = m.rmView();
            InlineRead(is,mrm);
        } else if (m.iscm()) {
            UpperTriMatrixView<T,UnknownDiag,1,UNKNOWN,C> mcm = m.cmView();
            InlineRead(is,mcm);
        } else 
            InlineRead(is,m);
    }
    template <class T, bool C>
    void InstRead(
        std::istream& is, LowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> m)
    {
        if (m.isrm()) {
            LowerTriMatrixView<T,UnknownDiag,UNKNOWN,1,C> mrm = m.rmView();
            InlineRead(is,mrm);
        } else if (m.iscm()) {
            LowerTriMatrixView<T,UnknownDiag,1,UNKNOWN,C> mcm = m.cmView();
            InlineRead(is,mcm);
        } else 
            InlineRead(is,m);
    }

#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


