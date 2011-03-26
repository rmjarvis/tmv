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

#include "tmv/TMV_MultMD.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_MultXD.h"

namespace tmv {

    // 
    //
    // MultMD
    //

    template <class M1, class M2>
    static void DoMultEq(M1& m1, const M2& m2)
    {
        typedef typename M1::real_type RT;
        Scaling<1,RT> one;
        if (m2.step() == 1) {
            InlineMultMM<false>(one,m1.constView(),m2.unitView(),m1);
        } else {
            InlineMultMM<false>(one,m1.constView(),m2.copy().view(),m1);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>  
    void InstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        if (SameStorage(m1,m3)) {
            // Must be exact same storage.  Just check conj.
            Maybe<C1>::conjself(m3);
            InstScale(x,m3);
        } else {
            InstMultXM(x,m1,m3);
        }
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3.cmView();
            DoMultEq(m3cm,m2);
        } else if (m3.isrm()) {
            MatrixView<T3,RowMajor> m3rm = m3.rmView();
            DoMultEq(m3rm,m2);
        } else {
            Matrix<T3,ColMajor|NoDivider> m3c = m3;
            MatrixView<T3,ColMajor> m3cv = m3c.view();
            DoMultEq(m3cv,m2);
            InstCopy(m3c.constView().xView(),m3);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>  
    void InstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        typedef typename Traits2<T1,T2>::type T12;
        Matrix<T12,ColMajor|NoDivider> m1c = m1;
        MatrixView<T12,ColMajor> m1cv = m1c.view();
        DoMultEq(m1cv,m2);
        InstAddMultXM(x,m1c.constView().xView(),m3);
    }

    template <class T1, int C1, class T2, int C2, class T3>  
    void InstAliasMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>  
    void InstAliasAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }


#define InstFile "TMV_MultMD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


