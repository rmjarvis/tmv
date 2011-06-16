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
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_SwapM.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_LUD.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_SmallMatrix.h"

namespace tmv {

    //
    // Copy Matrices
    //

    template <class T1, int C1, class T2>
    static void NonLapCopy(
        const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,ColMajor> m2cm = m2;
            if (m1.iscm()) {
                InlineCopy(m1.cmView(),m2cm);
            } else if (m1.isrm()) {
                InlineCopy(m1.rmView(),m2cm);
            } else
                InlineCopy(m1,m2cm);
        } else if (m2.isrm()) {
            NonLapCopy(m1.transpose(),m2.transpose());
        } else {
            if (m1.isrm())
                InlineCopy(m1.rmView(),m2);
            else if (m2.iscm())
                InlineCopy(m1.cmView(),m2);
            else
                InlineCopy(m1,m2);
        }
    }

#ifdef ELAP
    template <class T1, int C1, class T2>
    static void LapCopy(
        const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2)
    { NonLapCopy(m1,m2); }
#ifdef TMV_INST_DOUBLE
    static void LapCopy(
        const ConstMatrixView<double>& m1, MatrixView<double> m2)
    {
        //TMVAssert(m1.cptr() != m2.cptr());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        if (ld1 < 0 || ld2 < 0) NonLapCopy(m1,m2);
        TMVAssert(ld1 >= m);
        TMVAssert(ld2 >= m);
        LAPNAME(dlacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
    static void LapCopy(
        const ConstMatrixView<std::complex<double> >& m1,
        MatrixView<std::complex<double> > m2)
    {
        //TMVAssert(m1.cptr() != m2.cptr());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        if (ld1 < 0 || ld2 < 0) NonLapCopy(m1,m2);
        TMVAssert(ld1 >= m);
        TMVAssert(ld2 >= m);
        LAPNAME(zlacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
#endif
#ifdef TMV_INST_FLOAT
    static void LapCopy(
        const ConstMatrixView<float>& m1, MatrixView<float> m2)
    {
        //TMVAssert(m1.cptr() != m2.cptr());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        if (ld1 < 0 || ld2 < 0) NonLapCopy(m1,m2);
        TMVAssert(ld1 >= m);
        TMVAssert(ld2 >= m);
        LAPNAME(slacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
    static void LapCopy(
        const ConstMatrixView<std::complex<float> >& m1,
        MatrixView<std::complex<float> > m2)
    {
        //TMVAssert(m1.cptr() != m2.cptr());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        if (ld1 < 0 || ld2 < 0) NonLapCopy(m1,m2);
        TMVAssert(ld1 >= m);
        TMVAssert(ld2 >= m);
        LAPNAME(clacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
#endif
#endif

    template <class T1, int C1, class T2>
    void InstCopy(
        const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2)
    {
#ifdef ELAP
        if (m1.iscm() && m2.iscm() && m1.stepj()>0 && m2.stepj()>0)
            LapCopy(m1,m2);
        else
#endif
            NonLapCopy(m1,m2); 
    }

    template <class T1, int C1, class T2>
    void InstAliasCopy(
        const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2)
    { InlineAliasCopy(m1,m2); }

    //
    // Swap Matrices
    //

    template <class T, int C> 
    void InstSwap(MatrixView<T,C> m1, MatrixView<T> m2)
    {
        if (m2.iscm()) {
            MatrixView<T,ColMajor> m2cm = m2.cmView();
            if (m1.iscm()) {
                MatrixView<T,C|ColMajor> m1cm = m1.cmView();
                InlineSwap(m1cm,m2cm);
            } else if (m1.isrm()) {
                MatrixView<T,C|RowMajor> m1rm = m1.rmView();
                InlineSwap(m1rm,m2cm);
            } else
                InlineSwap(m1,m2cm);
        } else if (m2.isrm()) {
            InstSwap(m1.transpose(),m2.transpose());
        } else {
            if (m1.isrm()) {
                MatrixView<T,C|RowMajor> m1rm = m1.rmView();
                InlineSwap(m1rm,m2);
            } else if (m1.iscm()) {
                MatrixView<T,C|ColMajor> m1cm = m1.cmView();
                InlineSwap(m1cm,m2);
            } else
                InlineSwap(m1,m2);
        }
    }

    template <class T, int C> 
    void InstAliasSwap(MatrixView<T,C> m1, MatrixView<T> m2)
    { InlineAliasSwap(m1,m2); }


    // 
    // TransposeSelf
    //

    template <class T>
    void InstTransposeSelf(MatrixView<T> m)
    {
        if (m.iscm()) {
            MatrixView<T,ColMajor> mcm = m;
            InlineTransposeSelf(mcm);
        } else if (m.isrm()) {
            MatrixView<T,ColMajor> mt = m.transpose();
            InlineTransposeSelf(mt);
        } else {
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
    void InstPermuteRows(
        MatrixView<T> m, const int* p, const int i1, const int i2)
    {
        if (m.iscm()) {
            MatrixView<T,ColMajor> mcm = m;
            InlinePermuteRows(mcm,p,i1,i2);
        } else if (m.isrm()) {
            MatrixView<T,RowMajor> mrm = m;
            InlinePermuteRows(mrm,p,i1,i2);
        } else {
            InlinePermuteRows(m,p,i1,i2);
        }
    }

    template <class T>
    void InstReversePermuteRows(
        MatrixView<T> m, const int* p, const int i1, const int i2)
    {
        if (m.iscm()) {
            MatrixView<T,ColMajor> mcm = m;
            InlineReversePermuteRows(mcm,p,i1,i2);
        } else if (m.isrm()) {
            MatrixView<T,RowMajor> mrm = m;
            InlineReversePermuteRows(mrm,p,i1,i2);
        } else {
            InlineReversePermuteRows(m,p,i1,i2);
        }
    }




    //
    // Norms, SumElements
    //

    template <class T>
    static T DoInstSumElements(const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumElements(m.cmView());
        else if (m.isrm()) return InlineSumElements(m.transpose().cmView()); 
        else return InlineSumElements(m);
    }

    template <class T>
    static typename ConstMatrixView<T>::float_type DoInstSumAbsElements(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumAbsElements(m.cmView());
        else if (m.isrm()) return InlineSumAbsElements(m.transpose().cmView()); 
        else return InlineSumAbsElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbs2Elements(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumAbs2Elements(m.cmView());
        else if (m.isrm()) return InlineSumAbs2Elements(
            m.transpose().cmView()); 
        else return InlineSumAbs2Elements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormSq(m.cmView());
        else if (m.isrm()) return InlineNormSq(m.transpose().cmView()); 
        else return InlineNormSq(m);
    }

    template <class T>
    static typename ConstMatrixView<T>::float_type DoInstNormSq(
        const ConstMatrixView<T>& m, 
        const typename ConstMatrixView<T>::float_type scale)
    {
        if (m.iscm()) return InlineNormSq(m.cmView(),scale);
        else if (m.isrm()) return InlineNormSq(m.transpose().cmView(),scale); 
        else return InlineNormSq(m,scale);
    }

    template <class T>
    static typename ConstMatrixView<T>::float_type DoInstNormF(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormF(m.cmView());
        else if (m.isrm()) return InlineNormF(m.transpose().cmView()); 
        else return InlineNormF(m);
    }

    template <class T>
    static typename ConstMatrixView<T>::float_type DoInstMaxAbsElement(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineMaxAbsElement(m.cmView());
        else if (m.isrm()) return InlineMaxAbsElement(m.transpose().cmView()); 
        else return InlineMaxAbsElement(m);
    }
 
    template <class T>
    static typename Traits<T>::real_type DoInstMaxAbs2Element(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineMaxAbs2Element(m.cmView());
        else if (m.isrm()) return InlineMaxAbs2Element(m.transpose().cmView()); 
        else return InlineMaxAbs2Element(m);
    }

    template <class T>
    static typename ConstMatrixView<T>::float_type DoInstNorm1(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNorm1(m.cmView());
        else if (m.isrm()) return InlineNorm1(m.rmView()); 
        else return InlineNorm1(m);
    }

#ifdef XLAP
#ifdef TMV_INST_DOUBLE
    static double DoInstNormF(const ConstMatrixView<double>& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        if (m.isrm()) return DoInstNormF(m.transpose());
        char c = 'F';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(dlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }
    static double DoInstNormF(const ConstMatrixView<std::complex<double> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        if (m.isrm()) return DoInstNormF(m.transpose());
        char c = 'F';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(zlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }

    static double DoInstMaxAbsElement(const ConstMatrixView<double>& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        if (m.isrm()) return DoInstMaxAbsElement(m.transpose());
        char c = 'M';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(dlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }
    static double DoInstMaxAbsElement(
        const ConstMatrixView<std::complex<double> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        if (m.isrm()) return DoInstMaxAbsElement(m.transpose());
        char c = 'M';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        double work;
        double norm = LAPNAME(zlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }

    static double DoInstMaxAbs2Element(const ConstMatrixView<double>& m)
    { return DoInstMaxAbsElement(m); }
    static double DoInstMaxAbs2Element(
        const ConstMatrixView<std::complex<double> >& m)
    { return InlineMaxAbs2Element(m); }

    static double DoInstNorm1(const ConstMatrixView<double>& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? '1' : 'I';
        int M = m.iscm() ? m.colsize() : m.rowsize();
        int N = m.iscm() ? m.rowsize() : m.colsize();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(dlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
        return norm;
    }
    static double DoInstNorm1(const ConstMatrixView<std::complex<double> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? '1' : 'I';
        int M = m.iscm() ? m.colsize() : m.rowsize();
        int N = m.iscm() ? m.rowsize() : m.colsize();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(zlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
        return norm;
    }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static float DoInstNormF(const ConstMatrixView<float>& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        if (m.isrm()) return DoInstNormF(m.transpose());
        char c = 'F';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(slange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }
    static float DoInstNormF(const ConstMatrixView<std::complex<float> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNormF(m.copy());
        if (m.isrm()) return DoInstNormF(m.transpose());
        char c = 'F';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(clange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }

    static float DoInstMaxAbsElement(const ConstMatrixView<float>& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        if (m.isrm()) return DoInstMaxAbsElement(m.transpose());
        char c = 'M';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(slange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }
    static float DoInstMaxAbsElement(
        const ConstMatrixView<std::complex<float> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstMaxAbsElement(m.copy());
        if (m.isrm()) return DoInstMaxAbsElement(m.transpose());
        char c = 'M';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
        float work;
        float norm = LAPNAME(clange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(&work) LAP1);
        return norm;
    }

    static float DoInstMaxAbs2Element(const ConstMatrixView<float>& m)
    { return DoInstMaxAbsElement(m); }
    static float DoInstMaxAbs2Element(
        const ConstMatrixView<std::complex<float> >& m)
    { return InlineMaxAbs2Element(m); }

    static float DoInstNorm1(const ConstMatrixView<float>& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? '1' : 'I';
        int M = m.iscm() ? m.colsize() : m.rowsize();
        int N = m.iscm() ? m.rowsize() : m.colsize();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(slange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
        return norm;
    }
    static float DoInstNorm1(const ConstMatrixView<std::complex<float> >& m)
    {
        if (!m.iscm() && !m.isrm()) return DoInstNorm1(m.copy());
        char c = m.iscm() ? '1' : 'I';
        int M = m.iscm() ? m.colsize() : m.rowsize();
        int N = m.iscm() ? m.rowsize() : m.colsize();
        int lda = m.iscm() ? m.stepj() : m.stepi();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(clange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
        return norm;
    }
#endif // FLOAT
#endif // XLAP

    template <class T>
    T InstSumElements(const ConstMatrixView<T>& m)
    { return DoInstSumElements(m); }

    template <class T>
    typename ConstMatrixView<T>::float_type InstSumAbsElements(
        const ConstMatrixView<T>& m)
    { return DoInstSumAbsElements(m); }

    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstMatrixView<T>& m)
    { return DoInstSumAbs2Elements(m); }

    template <class T>
    typename Traits<T>::real_type InstNormSq(const ConstMatrixView<T>& m)
    { return DoInstNormSq(m); }

    template <class T>
    typename ConstMatrixView<T>::float_type InstNormSq(
        const ConstMatrixView<T>& m, 
        const typename ConstMatrixView<T>::float_type scale)
    { return DoInstNormSq(m,scale); }

    template <class T>
    typename ConstMatrixView<T>::float_type InstNormF(
        const ConstMatrixView<T>& m)
    { return DoInstNormF(m); }

    template <class T>
    typename ConstMatrixView<T>::float_type InstMaxAbsElement(
        const ConstMatrixView<T>& m)
    { return DoInstMaxAbsElement(m); }

    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstMatrixView<T>& m)
    { return DoInstMaxAbs2Element(m); }

    template <class T>
    typename ConstMatrixView<T>::float_type InstNorm1(
        const ConstMatrixView<T>& m)
    { return DoInstNorm1(m); }

#if 0
    template <class T> 
    RT GenMatrix<T>::doNorm2() const
    {
        if (this->colsize() < this->rowsize()) return transpose().doNorm2();
        if (this->rowsize() == 0) return RT(0);
        Matrix<T> m = *this;
        DiagMatrix<RT> S(this->rowsize());
        SV_Decompose(m.view(),S.view(),false);
        return S(0);
    }

    template <class T> 
    RT GenMatrix<T>::DoCondition() const
    {
        if (this->colsize() < this->rowsize()) return transpose().doNorm2();
        if (this->rowsize() == 0) return RT(1);
        Matrix<T> m = *this;
        DiagMatrix<RT> S(this->rowsize());
        SV_Decompose(m.view(),S.view(),false);
        return S(0)/S(S.size()-1);
    }

    template <class T> 
    QuotXM<T,T> GenMatrix<T>::QInverse() const
    { return QuotXM<T,T>(T(1),*this); }
#endif

    //
    // I/O
    //

    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,C>& m)
    {
        if (m.iscm()) InlineWrite(os,m.cmView());
        else if (m.isrm()) InlineWrite(os,m.rmView());
        else InlineWrite(os,m);
    }

    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,C>& m,
        typename ConstMatrixView<T>::float_type thresh)
    {
        if (m.iscm()) InlineWrite(os,m.cmView(),thresh);
        else if (m.isrm()) InlineWrite(os,m.rmView(),thresh);
        else InlineWrite(os,m,thresh);
    }

    template <class T>
    void InstRead(std::istream& is, MatrixView<T> m)
    {
        if (m.iscm()) {
            MatrixView<T,ColMajor> mcm = m.cmView();
            InlineRead(is,mcm);
        } else if (m.isrm()) {
            MatrixView<T,RowMajor> mrm = m.rmView();
            InlineRead(is,mrm);
        } else 
            InlineRead(is,m);
    }

    //
    // MatrixDivHelper
    //

    template <class T>
    MatrixDivHelper2<T>::MatrixDivHelper2() {}

    template <class T>
    MatrixDivHelper2<T>::~MatrixDivHelper2() {}

    template <class T>
    bool MatrixDivHelper2<T>::mIsSquare() const 
    { return getConstView().isSquare(); }

    template <class T>
    Matrix<T> MatrixDivHelper2<T>::getM() const
    { return Matrix<T>(getConstView()); }

    template <class T>
    void MatrixDivHelper2<T>::setDiv() const
    {
        TMVStaticAssert(!Traits<T>::isinteger);
        if (!this->divIsSet()) {
            DivType dt = this->getDivType();
            TMVAssert(dt == tmv::LU || dt == tmv::QR 
                      /*|| dt == tmv::QRP
                       *|| dt == tmv::SV */
            );
            switch (dt) {
              case tmv::LU : 
                   this->divider.reset(
                       new InstLUD<T>(
                           this->getConstView(),this->divIsInPlace()));
                   break;
              case tmv::QR : 
                   this->divider.reset(
                       new InstQRD<T>(
                           this->getConstView(),this->divIsInPlace()));
                   break;
              default :
                   // The above assert should have already failed.
                   // So go ahead and fall through.
                   break;
            }
        }
    }

    template <class T>
    const InstLUD<T>& MatrixDivHelper2<T>::lud() const
    {
        this->divideUsing(LU);
        setDiv();
        TMVAssert(dynamic_cast<const InstLUD<T>*>(this->getDiv()));
        return static_cast<const InstLUD<T>&>(*this->getDiv());
    }

    template <class T>
    const InstQRD<T>& MatrixDivHelper2<T>::qrd() const
    {
        this->divideUsing(QR);
        setDiv();
        TMVAssert(dynamic_cast<const InstQRD<T>*>(this->getDiv()));
        return static_cast<const InstQRD<T>&>(*this->getDiv());
    }

    template <class T>
    void MatrixDivHelper2<T>::qrpd() const {}
    template <class T>
    void MatrixDivHelper2<T>::svd() const {}

#if 0

    template <class T>
    const InstQRPD<T>& MatrixDivHelper2<T>::qrpd() const
    {
        this->divideUsing(QRP);
        setDiv();
        TMVAssert(dynamic_cast<const InstQRPD<T>*>(this->getDiv()));
        return static_cast<const InstQRPD<T>&>(*this->getDiv());
    }

    template <class T>
    const InstSVD<T>& MatrixDivHelper2<T>::svd() const
    {
        this->divideUsing(SV);
        setDiv();
        TMVAssert(dynamic_cast<const InstSVD<T>*>(this->getDiv()));
        return static_cast<const InstSVD<T>&>(*this->getDiv());
    }
#endif

#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


