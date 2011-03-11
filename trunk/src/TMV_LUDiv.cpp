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
#include "tmv/TMV_LUDiv.h"

namespace tmv {

    template <bool trans, class M1, class M2>
    struct NonLapLUSolve // trans == false here
    {
        static void call(const M1& m1, const Permutation& P, M2& m2)
        {
            if (m2.iscm()) {
                typename M2::cmview_type m2cm = m2.cmView();
                InlineLU_SolveInPlace(m1,P,m2cm);
            } else {
                typename M2::rmview_type m2rm = m2.rmView();
                InlineLU_SolveInPlace(m1,P,m2rm);
            }
        }
        static void call2(const M1& m1, const Permutation& P, M2& v2)
        {
            typename M2::unitview_type v2u = v2.unitView();
            InlineLU_SolveInPlace(m1,P,v2u);
        }
    };
    template <class M1, class M2>
    struct NonLapLUSolve<true,M1,M2>
    {
        static void call(const M1& m1, const Permutation& P, M2& m2)
        {
            if (m2.iscm()) {
                typename M2::cmview_type m2cm = m2.cmView();
                InlineLU_SolveTransposeInPlace(m1,P,m2cm);
            } else {
                typename M2::rmview_type m2rm = m2.rmView();
                InlineLU_SolveTransposeInPlace(m1,P,m2rm);
            }
        }
        static void call2(const M1& m1, const Permutation& P, M2& v2)
        {
            typename M2::unitview_type v2u = v2.unitView();
            InlineLU_SolveTransposeInPlace(m1,P,v2u);
        }
    };

#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
    template <bool trans, class M1, class T2, class T1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, MatrixView<T2>& m2, T1)
    { NonLapLUSolve<trans,M1,MatrixView<T2> >::call(m1,P,m2); }
    template <bool trans, class M1, class T2, class T1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, VectorView<T2>& v2, T1)
    { NonLapLUSolve<trans,M1,VectorView<T2> >::call2(m1,P,v2); }
#ifdef INST_DOUBLE
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, MatrixView<double>& m2, double)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv[i] = P.getValues()[i] LAPPLUS1;
#endif
        LAPNAME(dgetrs) (
            LAPCM trans?LAPCH_T:LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m2.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("dgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<double> > m2, std::complex<double> )
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv[i] = P.getValues()[i] LAPPLUS1;
#endif
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAPNAME(zgetrs) (
            LAPCM trans?(m1.isconj()?LAPCH_CT:LAPCH_T):LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAP_Results("zgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<double> > m2, double t1)
    {
        Matrix<double,ColMajor> temp = m2.realPart();
        LapLUSolve(m1,P,temp,t1);
        m2.realPart() = temp;
        temp = m2.imagPart();
        LapLUSolve(m1,P,temp,t1);
        m2.imagPart() = temp;
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, VectorView<double> v2, double t1)
    { LapLUSolve(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, 
        VectorView<std::complex<double> > v2, std::complex<double> t1)
    { LapLUSolve(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, 
        VectorView<std::complex<double> > v2, double t1)
    { LapLUSolve(m1,P,ColVectorViewOf(v2).xView(),t1); }
#endif
#ifdef INST_FLOAT
#ifndef MKL // TODO: Check if this is still required...
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, MatrixView<float>& m2, float)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv[i] = P.getValues()[i] LAPPLUS1;
#endif
        LAPNAME(sgetrs) (
            LAPCM trans?LAPCH_T:LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m2.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("sgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<float> > m2, std::complex<float> )
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv[i] = P.getValues()[i] LAPPLUS1;
#endif
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAPNAME(cgetrs) (
            LAPCM trans?(m1.isconj()?LAPCH_CT:LAPCH_T):LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m.ptr()),LAPV(ldb) LAPINFO LAP1);
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAP_Results("cgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<float> > m2, float t1)
    {
        Matrix<float,ColMajor> temp = m2.realPart();
        LapLUSolve(m1,P,temp,t1);
        m2.realPart() = temp;
        temp = m2.imagPart();
        LapLUSolve(m1,P,temp,t1);
        m2.imagPart() = temp;
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, VectorView<float> v2, float t1)
    { LapLUSolve(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, 
        VectorView<std::complex<float> > v2, std::complex<float> t1)
    { LapLUSolve(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        VectorView<std::complex<float> > v2, float t1)
    { LapLUSolve(m1,P,ColVectorViewOf(v2).xView(),t1); }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <bool trans, class M1, class T2>
    static void DoLUSolve(const M1& m1, const Permutation& P, MatrixView<T2> m2)
    {
#ifdef ALAP
        if (m1.stepj()>0) {
            const typename M1::value_type t1(0);
            if ( (m2.iscm() && m2.stepj()>0) ) {
                LapLUSolve(m1,P,m2,t1);
            } else {
                Matrix<T2,ColMajor> m2c(m2);
                LapLUSolve(m1,P,m2c.constView().xView(),t1);
                InstCopy(m2c.constView().xView(),m2);
            }
        } else {
            Matrix<M1::value_type,ColMajor> m1c(m1);
            DoLUSolve(m1c.xView(),P,m2);
        }
#else
        if (m2.iscm() || m2.isrm()) {
            NonLapLUSolve<trans,M1,MatrixView<T2> >::call(m1,P,m2);
        } else {
            Matrix<T2,ColMajor> m2c(m2);
            MatrixView<T2> m2cv = m2c.xView();
            NonLapLUSolve<trans,M1,MatrixView<T2> >::call(m1,P,m2cv);
            InstCopy(m2cv.constView(),m2);
        }
#endif
    }

    template <bool trans, class M1, class T2>
    static void DoLUSolve(const M1& m1, const Permutation& P, VectorView<T2> v2)
    {
#ifdef ALAP
        if (m1.stepj()>0) {
            const typename M1::value_type t1(0);
            if ( v2.step() == 1 ) {
                LapLUSolve(m1,P,v2,t1);
            } else  {
                Vector<T2> v2c(v2);
                LapLUSolve(m1,P,v2c.constView().xView(),t1);
                InstCopy(v2c.constView().xView(),v2);
            }
        } else {
            Matrix<M1::value_type,ColMajor> m1c(m1);
            DoLUSolve(m1c.xView(),P,v2);
        }
#else
        if ( v2.step() == 1 ) {
            NonLapLUSolve<trans,M1,VectorView<T2> >::call2(m1,P,v2);
        } else {
            Vector<T2> v2c(v2);
            VectorView<T2> v2cv = v2c.xView();
            NonLapLUSolve<trans,M1,VectorView<T2> >::call2(m1,P,v2cv);
            InstCopy(v2cv.constView(),v2);
        }
#endif
    }


    template <class T1, class T2, bool C1>
    void InstLU_SolveInPlace(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const Permutation& P,
        MatrixView<T2> m2)
    { DoLUSolve<false>(m1,P,m2); }
    template <class T1, class T2, bool C1>
    void InstLU_SolveTransposeInPlace(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const Permutation& P,
        MatrixView<T2> m2)
    { DoLUSolve<true>(m1,P,m2); }
    template <class T1, class T2, bool C1>
    void InstLU_SolveInPlace(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const Permutation& P,
        VectorView<T2> v2)
    { DoLUSolve<false>(m1,P,v2); }
    template <class T1, class T2, bool C1>
    void InstLU_SolveTransposeInPlace(
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1, const Permutation& P,
        VectorView<T2> v2)
    { DoLUSolve<true>(m1,P,v2); }

#define InstFile "TMV_LUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


