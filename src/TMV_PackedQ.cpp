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
#include "tmv/TMV_PackedQ.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"

namespace tmv {

    // 
    // A quick helper to turn the trans template parameter into MultEq or LDivEq
    //
    template <bool trans>
    struct PQMultEqHelper // trans == false here
    {
        template <class M1, class V1, class M2>
        static void call(const M1& Q, const V1& beta, M2& m2)
        { InlinePackedQ_MultEq(Q,beta,m2); }
    };
    template <>
    struct PQMultEqHelper<true>
    {
        template <class M1, class V1, class M2>
        static void call(const M1& Q, const V1& beta, M2& m2)
        { InlinePackedQ_LDivEq(Q,beta,m2); }
    };

    template <bool trans, class M1, class V1, class M2>
    static void DoPackedQ_MultEq(const M1& Q, const V1& beta, M2& m2)
    {
        if (beta.step() == 1) {
            if (Q.iscm()) {
                PQMultEqHelper<trans>::call(Q.cmView(),beta.unitView(),m2);
            } else if (Q.isrm()) {
                PQMultEqHelper<trans>::call(Q.rmView(),beta.unitView(),m2);
            } else {
                PQMultEqHelper<trans>::call(Q.copy().view(),beta.unitView(),m2);
            }
        } else {
            DoPackedQ_MultEq<trans>(Q,beta.copy().xView(),m2);
        }
    }

    template <bool trans, class M1, class V1, class T2>
    static void DoPackedQ_MultEq1(
        const M1& Q, const V1& beta, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,ColMajor> m2cm = m2;
            DoPackedQ_MultEq<trans>(Q,beta,m2cm);
        } else if (m2.isrm()) {
            MatrixView<T2,RowMajor> m2rm = m2;
            DoPackedQ_MultEq<trans>(Q,beta,m2rm);
        } else {
            Matrix<T2,ColMajor|NoDivider> m2c(m2);
            MatrixView<T2,ColMajor> m2cm = m2c.view();
            DoPackedQ_MultEq<trans>(Q,beta,m2cm);
            InstCopy(m2cm.constView().xView(),m2);
        }
    }

    template <bool trans, class M1, class V1, class T2>
    static void DoPackedQ_MultEq1(
        const M1& Q, const V1& beta, VectorView<T2> v2)
    {
        if (v2.step() == 1) {
            VectorView<T2,Unit> v2u = v2;
            DoPackedQ_MultEq<trans>(Q,beta,v2u);
        } else {
            Vector<T2> v2c(v2);
            VectorView<T2,Unit> v2u = v2c.view();
            DoPackedQ_MultEq<trans>(Q,beta,v2u);
            InstCopy(v2u.constView().xView(),v2);
        }
    }

    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2)
    { DoPackedQ_MultEq1<false>(Q,beta,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        MatrixView<T2> m2)
    { DoPackedQ_MultEq1<true>(Q,beta,m2); }
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_MultEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2)
    { DoPackedQ_MultEq1<false>(Q,beta,v2); }
    template <class T1, int C1, class RT1, class T2>
    void InstPackedQ_LDivEq(
        const ConstMatrixView<T1,C1>& Q, const ConstVectorView<RT1>& beta,
        VectorView<T2> v2)
    { DoPackedQ_MultEq1<true>(Q,beta,v2); }


#define InstFile "TMV_PackedQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


