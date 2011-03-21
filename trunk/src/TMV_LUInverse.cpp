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

//#define PRINTALGO_LU

#include "TMV_Blas.h"
#include "tmv/TMV_LUInverse.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SimpleMatrix.h"
#include "tmv/TMV_DivMU.h"

namespace tmv {

    template <class T1>
    static inline void NonLapLUInverse(MatrixView<T1> m1, const Permutation& P)
    {
        if (m1.iscm()) {
            MatrixView<T1,1> m1cm = m1.cmView();
            InlineLU_Inverse(m1cm,P);
        } else {
            MatrixView<T1,UNKNOWN,1> m1rm = m1.rmView();
            InlineLU_Inverse(m1rm,P);
        }
    }

#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
    template <class T1> 
    static inline void LapLUInverse(MatrixView<T1,1> m1, const Permutation& P)
    { NonLapLUInverse(m1,P); }
#ifdef INST_DOUBLE
    template <>
    static void LapLUInverse(MatrixView<double,1> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<double> work(new double[lwork]);
#else
        int lwork = -1;
        auto_array<double> work(new double[1]);
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new double[lwork]);
#endif
#endif
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("dgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"dgetri");
#endif
    }
    template <>
    static void LapLUInverse(
        MatrixView<std::complex<double>,1> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<double> > work(new std::complex<double>[1]);
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<double>[lwork]);
#endif
#endif
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
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
    static void LapLUInverse(MatrixView<float,1> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<float> work(new float[lwork]);
#else
        int lwork = -1;
        auto_array<float> work(new float[1]);
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new float[lwork]);
#endif
#endif
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("sgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"sgetri");
#endif
    }
    template <>
    static void LapLUInverse(
        MatrixView<std::complex<float>,1> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<float> > work(new std::complex<float>[1]);
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<float>[lwork]);
#endif
#endif
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("cgetri");
#else
        LAP_Results(int(std::real(work[0])),n,n,lwork,"cgetri");
#endif
    }
#endif // FLOAT
#endif // ALAP

    template <class T1>
    void InstLU_Inverse(MatrixView<T1> m1, const Permutation& P)
    {
#ifdef ALAP
        if (m1.iscm() && m1.stepj() > 0) {
            LapLUInverse(m1.cmView(),P);
        } else {
            SimpleMatrix<T1,ColMajor> m1c(m1);
            LapLUInverse(m1c.view(),P);
            InstCopy(m1c.xView().constView(),m1);
        }
#else
        if (m1.iscm() || m1.isrm()) {
            NonLapLUInverse(m1,P);
        } else {
            SimpleMatrix<T1,ColMajor> m1c(m1);
            NonLapLUInverse(m1c.xView(),P);
            InstCopy(m1c.xView().constView(),m1);
        }
#endif
    }

    template <class T1, class T2, bool C1>
    void InstLU_InverseATA(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const Permutation& P,
        const bool trans, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,1> m2cm = m2.cmView();
            InlineLU_InverseATA(m1,P,trans,m2cm);
        } else if (m2.isrm()) {
            MatrixView<T2,UNKNOWN,1> m2rm = m2.rmView();
            InlineLU_InverseATA(m1,P,trans,m2rm);
        } else {
            SimpleMatrix<T2,ColMajor> m2c(m2);
            MatrixView<T2,1> m2cm = m2c.view();
            InlineLU_InverseATA(m1,P,trans,m2cm);
            InstCopy(m2c.xView().constView(),m2);
        }
    }

#define InstFile "TMV_LUInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


