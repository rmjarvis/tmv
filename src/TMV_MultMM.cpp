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
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ProdXM.h"

#ifdef BLAS
#include "TMV_MultMM_Blas.h"
#endif

namespace tmv {

    // Defined in TMV_MultMM_CCC.cpp
    template <bool add, class T1, bool C1, class T2, bool C2, class T3>
    void DoInstMultMM(
        const T3 x,
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,1,UNKNOWN,C2>& m2, MatrixView<T3,1> m3);

    // Defined in TMV_MultMM_CRC.cpp
    template <bool add, class T1, bool C1, class T2, bool C2, class T3>
    void DoInstMultMM(
        const T3 x,
        const ConstMatrixView<T1,1,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,1,C2>& m2, MatrixView<T3,1> m3);

    // Defined in TMV_MultMM_RCC.cpp
    template <bool add, class T1, bool C1, class T2, bool C2, class T3>
    void DoInstMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,1,C1>& m1,
        const ConstMatrixView<T2,1,UNKNOWN,C2>& m2, MatrixView<T3,1> m3);

    // Defined in TMV_MultMM_RRC.cpp
    template <bool add, class T1, bool C1, class T2, bool C2, class T3>
    void DoInstMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,1,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,1,C2>& m2, MatrixView<T3,1> m3);

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    {
        if (m3.iscm()) {
            MatrixView<T3,1> m3cm = m3.cmView();
            if (m1.iscm()) 
                if (m2.iscm())
                    DoInstMultMM<false>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<false>(x,m1.cmView(),m2.rmView(),m3cm);
                else
                    InstMultMM(
                        T3(1),m1,
                        Matrix<T3,ColMajor>(x*m2).constView().xView(),m3);
            else if (m1.isrm())
                if (m2.iscm())
                    DoInstMultMM<false>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<false>(x,m1.rmView(),m2.rmView(),m3cm);
                else
                    InstMultMM(
                        T3(1),m1,
                        Matrix<T3,ColMajor>(x*m2).constView().xView(),m3);
            else
                InstMultMM(
                    T3(1),Matrix<T3,RowMajor>(x*m1).constView().xView(),m2,m3);
        } else if (m3.isrm()) {
            InstMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
        } else  {
            Matrix<T3,ColMajor> m3c(m3.colsize(),m3.rowsize());
            InstMultMM(x,m1,m2,m3c.xView());
            InstCopy(m3c.constView().xView(),m3);
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    {
        if (m3.iscm()) {
            MatrixView<T3,1> m3cm = m3.cmView();
            if (m1.iscm()) 
                if (m2.iscm())
                    DoInstMultMM<true>(x,m1.cmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<true>(x,m1.cmView(),m2.rmView(),m3cm);
                else
                    InstAddMultMM(
                        T3(1),m1,
                        Matrix<T3,ColMajor>(x*m2).constView().xView(),m3);
            else if (m1.isrm())
                if (m2.iscm())
                    DoInstMultMM<true>(x,m1.rmView(),m2.cmView(),m3cm);
                else if (m2.isrm())
                    DoInstMultMM<true>(x,m1.rmView(),m2.rmView(),m3cm);
                else
                    InstAddMultMM(
                        T3(1),m1,
                        Matrix<T3,ColMajor>(x*m2).constView().xView(),m3);
            else
                InstAddMultMM(
                    T3(1),Matrix<T3,RowMajor>(x*m1).constView().xView(),m2,m3);
        } else if (m3.isrm()) {
            InstAddMultMM(x,m2.transpose(),m1.transpose(),m3.transpose());
        } else  {
            Matrix<T3,ColMajor> m3c(m3.colsize(),m3.rowsize());
            InstMultMM(T3(1),m1,m2,m3c.xView());
            InstAddMultXM(x,m3c.constView().xView(),m3);
        }
    }

#define InstFile "TMV_MultMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


