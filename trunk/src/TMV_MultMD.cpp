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
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_MultXD.h"
#include "tmv/TMV_ProdXM.h"

namespace tmv {

    // 
    //
    // MultMD
    //

    template <bool add, class M1, class M2, class M3>
    static void DoMultMD(const M1& m1, const M2& m2, M3& m3)
    {
        typedef typename M3::real_type RT;
        Scaling<1,RT> one;
        if (m3.iscm()) {
            typename M3::cmview_type m3cm = m3.cmView();
            if (m1.iscm())
                InlineMultMM<add>(one,m1.cmView(),m2,m3cm);
            else
                InlineMultMM<add>(one,m1.rmView(),m2,m3cm);
        } else {
            typename M3::rmview_type m3rm = m3.rmView();
            if (m1.iscm())
                InlineMultMM<add>(one,m1.cmView(),m2,m3rm);
            else
                InlineMultMM<add>(one,m1.rmView(),m2,m3rm);
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>  
    void InstMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    {
        if ( x == T3(1) && (m1.iscm() || m1.isrm()) && 
             m2.step() == 1 && (m3.iscm() || m3.isrm()) )
            DoMultMD<false>(m1,m2.cmView(),m3);
        else if (m2.step() != 1 || x != T3(1)) 
            InstMultMM(T3(1),m1,(x*m2).calc().constView().xView(),m3);
        else if (!(m1.isrm() || m1.iscm())) 
            InstMultMM(x,m1.copy().constView().xView(),m2,m3);
        else { // !(m3.isrm() || m3.iscm())
            Matrix<T3,ColMajor> m3c(m3.colsize(),m3.rowsize());
            DoMultMD<false>(m1,m2.cmView(),m3c);
            InstCopy(m3c.constView().xView(),m3);
        }
    }
 

    template <class T1, bool C1, class T2, bool C2, class T3>  
    void InstAddMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3)
    {
        if ( x == T3(1) && (m1.iscm() || m1.isrm()) && 
             m2.step() == 1 && (m3.iscm() || m3.isrm()) )
            DoMultMD<true>(m1,m2.cmView(),m3);
        else if (m2.step() != 1 || x != T3(1)) 
            InstAddMultMM(T3(1),m1,(x*m2).calc().constView().xView(),m3);
        else if (!(m1.isrm() || m1.iscm())) 
            InstAddMultMM(x,m1.copy().constView().xView(),m2,m3);
        else { // !(m3.isrm() || m3.iscm())
            Matrix<T3,ColMajor> m3c(m3.colsize(),m3.rowsize());
            DoMultMD<false>(m1,m2.cmView(),m3c);
            InstAddMultXM(T3(1),m3c.constView().xView(),m3);
        }
    }


#define InstFile "TMV_MultMD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


