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
#include "tmv/TMV_LUInverse.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DivMU.h"

// since I use an AliasCopy in LU_Inverse
#include "tmv/TMV_TransposeM.h" 

// since I use normal arithmetic in LU_InverseATA
#include "tmv/TMV_QuotXM.h"
#include "tmv/TMV_QuotMM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_InvertM.h"
#include "tmv/TMV_ScaleU.h"

namespace tmv {

    template <class T2>
    static void NonLapLUInverse(const int* P, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,1> m2cm = m2.cmView();
            InlineLU_Inverse(m2cm.constView(),P,m2cm);
        } else {
            MatrixView<T2,UNKNOWN,1> m2rm = m2.rmView();
            InlineLU_Inverse(m2rm.constView(),P,m2rm);
        }
    }

#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
    template <class T2> 
    static inline void LapLUInverse(const int* P, MatrixView<T2,1> m2)
    { NonLapLUInverse(P,m2); }
#ifdef INST_DOUBLE
    template <>
    static void LapLUInverse(const int* P, MatrixView<double,1> m2)
    {
        int n = m2.colsize();
        int lda = m2.stepj();
#ifdef CLAP
        const int* ipiv = P;
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<double> work(new double[lwork]);
#else
        int lwork = -1;
        auto_array<double> work(new double[1]);
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new double[lwork]);
#endif
#endif
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("dgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"dgetri");
#endif
    }
    template <>
    static void LapLUInverse(
        const int* P, MatrixView<std::complex<double>,1> m2)
    {
        int n = m2.colsize();
        int lda = m2.stepj();
#ifdef CLAP
        const int* ipiv = P;
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<double> > work(new std::complex<double>[1]);
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<double>[lwork]);
#endif
#endif
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("zgetri");
#else
        LAP_Results(int(std::real(work[0])),n,n,lwork,"zgetri");
#endif
    }
#endif
#ifdef INST_FLOAT
    template <>
    static void LapLUInverse(const int* P, MatrixView<float,1> m2)
    {
        int n = m2.colsize();
        int lda = m2.stepj();
#ifdef CLAP
        const int* ipiv = P;
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<float> work(new float[lwork]);
#else
        int lwork = -1;
        auto_array<float> work(new float[1]);
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new float[lwork]);
#endif
#endif
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("sgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"sgetri");
#endif
    }
    template <>
    static void LapLUInverse(
        const int* P, MatrixView<std::complex<float>,1> m2)
    {
        int n = m2.colsize();
        int lda = m2.stepj();
#ifdef CLAP
        const int* ipiv = P;
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<float> > work(new std::complex<float>[1]);
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<float>[lwork]);
#endif
#endif
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m2.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("cgetri");
#else
        LAP_Results(int(std::real(work[0])),n,n,lwork,"cgetri");
#endif
    }
#endif // FLOAT
#endif // ALAP

    template <class T1, class T2, bool C1>
    void InstLU_Inverse(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const int* P,
        MatrixView<T2> m2)
    {
#ifdef ALAP
        if (m2.iscm() && m2.stepj() > 0) {
            InstCopy(m1.xView(),m2);
            LapLUInverse(P,m2.cmView());
        } else {
            Matrix<T2,ColMajor> m2c(m1);
            LapLUInverse(P,m2c.view());
            InstCopy(m2c.xView().constView(),m2);
        }
#else
        if (m2.iscm() || m2.isrm()) {
            InstCopy(m1.xView(),m2);
            NonLapLUInverse(P,m2);
        } else {
            Matrix<T2,ColMajor> m2c(m1);
            NonLapLUInverse(P,m2c.xView());
            InstCopy(m2c.xView().constView(),m2);
        }
#endif
    }

    template <class T1, class T2, bool C1>
    void InstLU_InverseATA(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const int* P,
        const bool trans, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,1> m2cm = m2.cmView();
            InlineLU_InverseATA(m1,P,trans,m2cm);
        } else if (m2.isrm()) {
            MatrixView<T2,UNKNOWN,1> m2rm = m2.rmView();
            InlineLU_InverseATA(m1,P,trans,m2rm);
        } else {
            Matrix<T2,ColMajor> m2c(m2);
            MatrixView<T2,1> m2cm = m2c.view();
            InlineLU_InverseATA(m1,P,trans,m2cm);
            InstCopy(m2c.xView().constView(),m2);
        }
    }

#define InstFile "TMV_LUInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


