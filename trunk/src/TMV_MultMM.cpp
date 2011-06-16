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

#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"

namespace tmv {

    // Defined in TMV_MultMM_CCC.cpp
    // Defined in TMV_MultMM_CRC.cpp
    // Defined in TMV_MultMM_RCC.cpp
    // Defined in TMV_MultMM_RRC.cpp
    template <bool add, class T1, int C1, class T2, int C2, class T3>
    void DoInstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if TMV_OPT <= 2
        m3.setZero();
        InstAddMultMM(x,m1,m2,m3);
#else
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3.cmView();
            if (m1.iscm()) {
                if (m2.iscm())
                    DoInstMultMM<false>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<false>(x,m1.cmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<false>(x,m1.cmView(),m2c.constView(),m3cm);
                }
            } else if (m1.isrm()) {
                if (m2.iscm())
                    DoInstMultMM<false>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<false>(x,m1.rmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<false>(x,m1.rmView(),m2c.constView(),m3cm);
                }
            } else {
                Matrix<T1,RowMajor|NoDivider> m1c = m1;
                InstMultMM(x,m1c.constView().xView(),m2,m3);
            }
        } else if (m3.isrm()) {
            InstMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
        } else  {
            Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
            InstMultMM(x,m1,m2,m3c.xView());
            InstCopy(m3c.constView().xView(),m3);
        }
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
        if (m3.iscm()) {
            MatrixView<T3,ColMajor> m3cm = m3.cmView();
            if (m1.iscm()) {
                if (m2.iscm())
                    DoInstMultMM<true>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<true>(x,m1.cmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<true>(x,m1.cmView(),m2c.constView(),m3cm);
                }
            } else if (m1.isrm()) {
                if (m2.iscm())
                    DoInstMultMM<true>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<true>(x,m1.rmView(),m2.rmView(),m3cm);
                else {
                    Matrix<T2,ColMajor|NoDivider> m2c = m2;
                    DoInstMultMM<true>(x,m1.rmView(),m2c.constView(),m3cm);
                }
            } else {
                Matrix<T1,RowMajor|NoDivider> m1c = m1;
                InstAddMultMM(x,m1c.constView().xView(),m2,m3);
            }
        } else if (m3.isrm()) {
            InstAddMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
        } else  {
            Matrix<T3,ColMajor|NoDivider> m3c(m3.colsize(),m3.rowsize());
            InstMultMM(T3(1),m1,m2,m3c.xView());
            InstAddMultXM(x,m3c.constView().xView(),m3);
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<false>(Scaling<0,T3>(x),m1,m2,m3); }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    { InlineAliasMultMM<true>(Scaling<0,T3>(x),m1,m2,m3); }

#define InstFile "TMV_MultMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


