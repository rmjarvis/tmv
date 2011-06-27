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

//#define PRINTALGO_QR

#include "TMV_Blas.h"
#include "tmv/TMV_QRDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultMM_Block.h"
#include "tmv/TMV_MultMM_Winograd.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_Det.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_MultXU.h"


namespace tmv {

    // 
    // A quick helper to turn the trans template parameter into the 
    // Solve or SolveTranspose call.
    //
    template <bool trans>
    struct QRDivHelper // trans == false here
    {
        template <class M1, class V1, class M2>
        static void callInPlace(
            const M1& QR, const V1& beta, 
            const Permutation* P, int N1, M2& m2)
        { InlineQR_SolveInPlace(QR,beta,P,N1,m2); }
        template <class M1, class V1, class M2, class M3>
        static void call(
            const M1& QR, const V1& beta, 
            const Permutation* P, int N1, const M2& m2, M3& m3)
        { InlineQR_Solve(QR,beta,P,N1,m2,m3); }
    };
    template <>
    struct QRDivHelper<true>
    {
        template <class M1, class V1, class M2>
        static void callInPlace(
            const M1& QR, const V1& beta, 
            const Permutation* P, int N1, M2& m2)
        { InlineQR_SolveTransposeInPlace(QR,beta,P,N1,m2); }
        template <class M1, class V1, class M2, class M3>
        static void call(
            const M1& QR, const V1& beta, 
            const Permutation* P, int N1, const M2& m2, M3& m3)
        { InlineQR_SolveTranspose(QR,beta,P,N1,m2,m3); }
    };

    //
    // The SolveInPlace functions:
    //

    template <bool trans, class M1, class V1, class M2>
    static void DoQR_SolveInPlace(
        const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
    {
        if (beta.step() == 1) {
            if (QR.iscm()) {
                QRDivHelper<trans>::callInPlace(
                    QR.cmView(),beta.unitView(),P,N1,m2);
            } else if (QR.isrm()) {
                QRDivHelper<trans>::callInPlace(
                    QR.rmView(),beta.unitView(),P,N1,m2);
            } else {
                QRDivHelper<trans>::callInPlace(
                    QR.copy().view(),beta.unitView(),P,N1,m2);
            }
        } else {
            DoQR_SolveInPlace<trans>(QR,beta.copy().xView(),P,N1,m2);
        }
    }

    template <bool trans, class M1, class V1, class T2>
    static void DoQR_SolveInPlace1(
        const M1& QR, const V1& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,ColMajor> m2cm = m2;
            DoQR_SolveInPlace<trans>(QR,beta,P,N1,m2cm);
        } else if (m2.isrm()) {
            MatrixView<T2,RowMajor> m2rm = m2;
            DoQR_SolveInPlace<trans>(QR,beta,P,N1,m2rm);
        } else {
            Matrix<T2,ColMajor|NoDivider> m2c(m2);
            MatrixView<T2,ColMajor> m2cm = m2c.view();
            DoQR_SolveInPlace<trans>(QR,beta,P,N1,m2cm);
            InstCopy(m2cm.constView().xView(),m2);
        }
    }

    template <bool trans, class M1, class V1, class T2>
    static void DoQR_SolveInPlace1(
        const M1& QR, const V1& beta,
        const Permutation* P, int N1, VectorView<T2> v2)
    {
        if (v2.step() == 1) {
            VectorView<T2,Unit> v2u = v2;
            DoQR_SolveInPlace<trans>(QR,beta,P,N1,v2u);
        } else {
            Vector<T2> v2c(v2);
            VectorView<T2,Unit> v2u = v2c.view();
            DoQR_SolveInPlace<trans>(QR,beta,P,N1,v2u);
            InstCopy(v2u.constView().xView(),v2);
        }
    }

    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    { DoQR_SolveInPlace1<false>(QR,beta,P,N1,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    { DoQR_SolveInPlace1<true>(QR,beta,P,N1,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, VectorView<T2> v2)
    { DoQR_SolveInPlace1<false>(QR,beta,P,N1,v2); }
    template <class T1, int C1, class RT1, class T2>
    void InstQR_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, VectorView<T2> v2)
    { DoQR_SolveInPlace1<true>(QR,beta,P,N1,v2); }


    //
    // The Solve functions:
    //

    template <bool trans, class M1, class V1, class M2, class M3>
    static void DoQR_Solve(
        const M1& QR, const V1& beta, const Permutation* P, int N1, 
        const M2& m2, M3& m3)
    {
        if (beta.step() == 1) {
            if (QR.iscm()) {
                QRDivHelper<trans>::call(
                    QR.cmView(),beta.unitView(),P,N1,m2,m3);
            } else if (QR.isrm()) {
                QRDivHelper<trans>::call(
                    QR.rmView(),beta.unitView(),P,N1,m2,m3);
            } else {
                QRDivHelper<trans>::call(
                    QR.copy().view(),beta.unitView(),P,N1,m2,m3);
            }
        } else {
            DoQR_Solve<trans>(QR,beta.copy().xView(),P,N1,m2,m3);
        }
    }

    template <bool trans, class M1, class V1, class T2, int C2, class M3>
    static void DoQR_Solve2(
        const M1& QR, const V1& beta, const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, M3& m3)
    {
        if (m2.iscm()) {
            DoQR_Solve<trans>(QR,beta,P,N1,m2.cmView(),m3);
        } else if (m2.isrm()) {
            DoQR_Solve<trans>(QR,beta,P,N1,m2.rmView(),m3);
        } else {
            DoQR_Solve<trans>(QR,beta,P,N1,m2.copy().view(),m3);
        }
    }

    template <bool trans, class M1, class V1, class M2, class T3>
    static void DoQR_Solve1(
        const M1& QR, const V1& beta, const Permutation* P, int N1,
        const M2& m2, MatrixView<T3> m3)
    {
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3;
            DoQR_Solve2<trans>(QR,beta,P,N1,m2,m3cm);
        } else if (m3.isrm()) {
            MatrixView<T3,RowMajor> m3rm = m3;
            DoQR_Solve2<trans>(QR,beta,P,N1,m2,m3rm);
        } else {
            Matrix<T3,ColMajor|NoDivider> m3c(m3);
            MatrixView<T3,ColMajor> m3cm = m3c.view();
            DoQR_Solve2<trans>(QR,beta,P,N1,m2,m3cm);
            InstCopy(m3cm.constView().xView(),m3);
        }
    }

    template <bool trans, class M1, class V1, class T2, int C2, class V3>
    static void DoQR_Solve2(
        const M1& QR, const V1& beta, const Permutation* P, int N1,
        const ConstVectorView<T2,C2>& v2, V3& v3)
    {
        if (v2.step() == 1) {
            DoQR_Solve<trans>(QR,beta,P,N1,v2.unitView(),v3);
        } else {
            DoQR_Solve<trans>(QR,beta,P,N1,v2.copy().view(),v3);
        }
    }

    template <bool trans, class M1, class V1, class V2, class T3>
    static void DoQR_Solve1(
        const M1& QR, const V1& beta, const Permutation* P, int N1,
        const V2& v2, VectorView<T3> v3)
    {
        if (v3.step() == 1) {
            VectorView<T3,Unit> v3u = v3;
            DoQR_Solve2<trans>(QR,beta,P,N1,v2,v3u);
        } else {
            Vector<T3> v3c(v3);
            VectorView<T3,Unit> v3u = v3c.view();
            DoQR_Solve2<trans>(QR,beta,P,N1,v2,v3u);
            InstCopy(v3u.constView().xView(),v3);
        }
    }

    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_Solve(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoQR_Solve1<false>(QR,beta,P,N1,m2,m3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_SolveTranspose(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { DoQR_Solve1<true>(QR,beta,P,N1,m2,m3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_Solve(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { DoQR_Solve1<false>(QR,beta,P,N1,v2,v3); }
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstQR_SolveTranspose(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { DoQR_Solve1<true>(QR,beta,P,N1,v2,v3); }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


