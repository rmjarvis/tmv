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
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_SwapM.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ProdXM.h"

namespace tmv {

#if 0 // Put this in TMV_DivHelper.cpp
    template <class T> 
    void GenMatrix<T>::NewDivider() const
    {
        switch (this->GetDivType()) {
          case LU : this->SetDiv(
                        new LUDiv<T>(*this, this->IsDivInPlace())); 
                    break;
          case QR : this->SetDiv(
                        new QRDiv<T>(*this, this->IsDivInPlace())); 
                    break;
          case QRP : this->SetDiv(
                         new QRPDiv<T>(*this, this->IsDivInPlace())); 
                     break;
          case SV : this->SetDiv(
                        new SVDiv<T>(*this, this->IsDivInPlace())); 
                    break;
          default : TMVAssert(FALSE);
        }
    }
#endif

    //
    // Copy Matrices
    //

    template <class T1, bool C1, class T2>
    static void NonLapCopy(
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, MatrixView<T2> m2)
    {
        if (m2.isrm() && !m2.iscm()) {
            NonLapCopy(m1.transpose(),m2.transpose());
        } else if (m2.iscm()) {
            MatrixView<T2,1> m2cm = m2;
            if (m1.isrm())
                InlineCopy(m1.rmView(),m2cm);
            else if (m1.iscm()) {
                if (m1.canLinearize() && m2.canLinearize()) {
                    VectorView<T2> m2l = m2.linearView();
                    InstCopy(m1.linearView().xView(),m2l);
                } else
                    InlineCopy(m1.cmView(),m2cm);
            } else
                InlineCopy(m1,m2cm);
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
    template <class T1, bool C1, class T2>
    static void LapCopy(
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, MatrixView<T2> m2)
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

    template <class T1, bool C1, class T2>
    void InstCopy(
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, MatrixView<T2> m2)
    {
#ifdef ELAP
        if (m1.iscm() && m2.iscm() && m1.stepj()>0 && m2.stepj()>0)
            LapCopy(m1,m2);
        else
#endif
            NonLapCopy(m1,m2); 
    }

    //
    // Swap Matrices
    //

    template <class T, bool C> 
    void InstSwap(MatrixView<T,UNKNOWN,UNKNOWN,C> m1, MatrixView<T> m2)
    {
        if (m2.isrm() && !m2.iscm()) {
            InstSwap(m1.transpose(),m2.transpose());
        } else if (m2.iscm()) {
            MatrixView<T,1> m2cm = m2.cmView();
            if (m1.isrm()) {
                MatrixView<T,UNKNOWN,1,C> m1rm = m1.rmView();
                InlineSwap(m1rm,m2cm);
            } else if (m1.iscm()) {
                if (m1.canLinearize() && m2.canLinearize()) {
                    VectorView<T,UNKNOWN,C> m1l = m1.linearView();
                    VectorView<T> m2l = m2.linearView();
                    InstSwap(m1l,m2l);
                } else {
                    MatrixView<T,1,UNKNOWN,C> m1cm = m1.cmView();
                    InlineSwap(m1cm,m2cm);
                }
            } else
                InlineSwap(m1,m2cm);
        } else {
            if (m1.isrm()) {
                MatrixView<T,UNKNOWN,1,C> m1rm = m1.rmView();
                InlineSwap(m1rm,m2);
            } else if (m1.iscm()) {
                MatrixView<T,1,UNKNOWN,C> m1cm = m1.cmView();
                InlineSwap(m1cm,m2);
            } else
                InlineSwap(m1,m2);
        }
    }


    // 
    // TransposeSelf
    //

    template <class T>
    void InstTransposeSelf(MatrixView<T> m)
    {
        if (m.isrm() && !m.iscm()) {
            InstTransposeSelf(m.transpose());
        } else if (m.iscm()) {
            MatrixView<T,1> mcm = m;
            InlineTransposeSelf(mcm);
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
        MatrixView<T> m, const int*const p, const int i1, const int i2)
    {
        if (m.isrm()) {
            MatrixView<T,UNKNOWN,1> mrm = m;
            InlinePermuteRows(mrm,p,i1,i2);
        } else if (m.iscm()) {
            MatrixView<T,1> mcm = m;
            InlinePermuteRows(mcm,p,i1,i2);
        } else {
            InlinePermuteRows(m,p,i1,i2);
        }
    }

    template <class T>
    void InstReversePermuteRows(
        MatrixView<T> m, const int*const p, const int i1, const int i2)
    {
        if (m.isrm()) {
            MatrixView<T,UNKNOWN,1> mrm = m;
            InlineReversePermuteRows(mrm,p,i1,i2);
        } else if (m.iscm()) {
            MatrixView<T,1> mcm = m;
            InlineReversePermuteRows(mcm,p,i1,i2);
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
        if (m.canLinearize()) return m.linearView().sumElements();
        else if (m.iscm()) return InlineSumElements(m.cmView());
        else if (m.isrm()) return InlineSumElements(m.transpose().cmView()); 
        else return InlineSumElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbsElements(
        const ConstMatrixView<T>& m)
    {
        if (m.canLinearize()) return m.linearView().sumAbsElements();
        else if (m.iscm()) return InlineSumAbsElements(m.cmView());
        else if (m.isrm()) return InlineSumAbsElements(m.transpose().cmView()); 
        else return InlineSumAbsElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbs2Elements(
        const ConstMatrixView<T>& m)
    {
        if (m.canLinearize()) return m.linearView().sumAbs2Elements();
        else if (m.iscm()) return InlineSumAbs2Elements(m.cmView());
        else if (m.isrm()) return InlineSumAbs2Elements(
            m.transpose().cmView()); 
        else return InlineSumAbs2Elements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(const ConstMatrixView<T>& m)
    {
        if (m.canLinearize()) return m.linearView().normSq();
        else if (m.iscm()) return InlineNormSq(m.cmView());
        else if (m.isrm()) return InlineNormSq(m.transpose().cmView()); 
        else return InlineNormSq(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(
        const ConstMatrixView<T>& m, const typename Traits<T>::real_type scale)
    {
        if (m.canLinearize()) return m.linearView().normSq(scale);
        else if (m.iscm()) return InlineNormSq(m.cmView(),scale);
        else if (m.isrm()) return InlineNormSq(m.transpose().cmView(),scale); 
        else return InlineNormSq(m,scale);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormF(const ConstMatrixView<T>& m)
    {
        if (m.canLinearize()) return m.linearView().norm2();
        else if (m.iscm()) return InlineNormF(m.cmView());
        else if (m.isrm()) return InlineNormF(m.transpose().cmView()); 
        else return InlineNormF(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstMaxAbsElement(
        const ConstMatrixView<T>& m)
    {
        if (m.canLinearize()) return m.linearView().maxAbsElement();
        else if (m.iscm()) return InlineMaxAbsElement(m.cmView());
        else if (m.isrm()) return InlineMaxAbsElement(m.transpose().cmView()); 
        else return InlineMaxAbsElement(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNorm1(const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNorm1(m.cmView());
        else if (m.isrm()) return InlineNormInf(m.transpose().cmView()); 
        else return InlineNorm1(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormInf(
        const ConstMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormInf(m.cmView());
        else if (m.isrm()) return InlineNorm1(m.transpose().cmView()); 
        else return InlineNormInf(m);
    }

#ifdef XLAP
#ifdef TMV_INST_DOUBLE
    static double DoInstNormF(const ConstMatrixView<double>& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNormF(m.transpose());
        if (!m.iscm()) return InlineNormF(m);
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
        if (m.isrm() && !m.iscm()) return DoInstNormF(m.transpose());
        if (!m.iscm()) return InlineNormF(m);
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
        if (m.isrm() && !m.iscm()) return DoInstMaxAbsElement(m.transpose());
        if (!m.iscm()) return InlineMaxAbsElement(m);
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
        if (m.isrm() && !m.iscm()) return DoInstMaxAbsElement(m.transpose());
        if (!m.iscm()) return InlineMaxAbsElement(m);
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

    static double DoInstNorm1(const ConstMatrixView<double>& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNormInf(m.transpose());
        if (!m.iscm()) return InlineNorm1(m);
        char c = '1';
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
    static double DoInstNorm1(const ConstMatrixView<std::complex<double> >& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNormInf(m.transpose());
        if (!m.iscm()) return InlineNorm1(m);
        char c = '1';
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

    static double DoInstNormInf(const ConstMatrixView<double>& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNorm1(m.transpose());
        if (!m.iscm()) return InlineNormInf(m);
        char c = 'I';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(dlange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
        return norm;
    }
    static double DoInstNormInf(const ConstMatrixView<std::complex<double> >& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNorm1(m.transpose());
        if (!m.iscm()) return InlineNormInf(m);
        char c = 'I';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<double> work(new double[M]);
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
        if (m.isrm() && !m.iscm()) return DoInstNormF(m.transpose());
        if (!m.iscm()) return InlineNormF(m);
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
        if (m.isrm() && !m.iscm()) return DoInstNormF(m.transpose());
        if (!m.iscm()) return InlineNormF(m);
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
        if (m.isrm() && !m.iscm()) return DoInstMaxAbsElement(m.transpose());
        if (!m.iscm()) return InlineMaxAbsElement(m);
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
        if (m.isrm() && !m.iscm()) return DoInstMaxAbsElement(m.transpose());
        if (!m.iscm()) return InlineMaxAbsElement(m);
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

    static float DoInstNorm1(const ConstMatrixView<float>& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNormInf(m.transpose());
        if (!m.iscm()) return InlineNorm1(m);
        char c = '1';
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
    static float DoInstNorm1(const ConstMatrixView<std::complex<float> >& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNormInf(m.transpose());
        if (!m.iscm()) return InlineNorm1(m);
        char c = '1';
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

    static float DoInstNormInf(const ConstMatrixView<float>& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNorm1(m.transpose());
        if (!m.iscm()) return InlineNormInf(m);
        char c = 'I';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(slange) (
            LAPCM LAPV(c),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1);
        return norm;
    }
    static float DoInstNormInf(const ConstMatrixView<std::complex<float> >& m)
    {
        if (m.isrm() && !m.iscm()) return DoInstNorm1(m.transpose());
        if (!m.iscm()) return InlineNormInf(m);
        char c = 'I';
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
        TMVAssert(lda >= m);
#ifndef LAPNOWORK
        auto_array<float> work(new float[M]);
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
    typename Traits<T>::real_type InstSumAbsElements(
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
    typename Traits<T>::real_type InstNormSq(
        const ConstMatrixView<T>& m, const typename Traits<T>::real_type scale)
    { return DoInstNormSq(m,scale); }

    template <class T>
    typename Traits<T>::real_type InstNormF(const ConstMatrixView<T>& m)
    { return DoInstNormF(m); }

    template <class T>
    typename Traits<T>::real_type InstMaxAbsElement(const ConstMatrixView<T>& m)
    { return DoInstMaxAbsElement(m); }

    template <class T>
    typename Traits<T>::real_type InstNorm1(const ConstMatrixView<T>& m)
    { return DoInstNorm1(m); }

    template <class T>
    typename Traits<T>::real_type InstNormInf(const ConstMatrixView<T>& m)
    { return DoInstNormInf(m); }


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

    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,UNKNOWN,UNKNOWN,C>& m)
    {
        if (m.isrm()) InlineWrite(os,m.rmView());
        else if (m.iscm()) InlineWrite(os,m.cmView());
        else InlineWrite(os,m);
    }

    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,UNKNOWN,UNKNOWN,C>& m,
        typename Traits<T>::real_type thresh)
    {
        if (m.isrm()) InlineWrite(os,m.rmView(),thresh);
        else if (m.iscm()) InlineWrite(os,m.cmView(),thresh);
        else InlineWrite(os,m,thresh);
    }

    template <class T, bool C>
    void InstRead(std::istream& is, MatrixView<T,UNKNOWN,UNKNOWN,C> m)
    {
        if (m.isrm()) {
            MatrixView<T,UNKNOWN,1,C> mrm = m.rmView();
            InlineRead(is,mrm);
        } else if (m.iscm()) {
            MatrixView<T,1,UNKNOWN,C> mcm = m.cmView();
            InlineRead(is,mcm);
        } else 
            InlineRead(is,m);
    }

#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


