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

#include "tmv/TMV_MultUD.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_MultXD.h"
#include "tmv/TMV_ProdXM.h"

namespace tmv {

    // 
    //
    // MultUD
    //

    template <bool add, class M1, class M2, class M3>
    static void DoMultUD(const M1& m1, const M2& m2, M3& m3)
    {
        TMVAssert(m1.isrm() || m1.iscm());
        TMVAssert(m3.isrm() || m3.iscm());
        TMVStaticAssert(M2::mdiagstep == 1);
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


    template <class T, class M1, class M2, class M3>  
    static void GenInstMultMM(const T x, const M1& m1, const M2& m2, M3& m3)
    { 
        if ( x == T(1) && (m1.iscm() || m1.isrm()) &&
             m2.step() == 1 && (m3.iscm() || m3.isrm()))
            DoMultUD<false>(m1,m2.cmView(),m3);
        else if (m2.step() != 1 || x != T(1))
            GenInstMultMM(T(1),m1,(x*m2).calc().constView().xView(),m3);
        else if (!(m1.isrm() || m1.iscm()))
            GenInstMultMM(x,m1.copy().constView().xdView(),m2,m3);
        else { // !(m3.isrm() || m3.iscm())
            typename M3::copy_type m3c(m3.size());
            DoMultUD<false>(m1,m2.cmView(),m3c);
            InstCopy(m3c.constView().xdView(),m3.xdView());
        }
    }



    template <class T, class M1, class M2, class M3>  
    static void GenInstAddMultMM(const T x, const M1& m1, const M2& m2, M3& m3)
    {
        if ( x == T(1) && (m1.iscm() || m1.isrm()) &&
             m2.step() == 1 && (m3.iscm() || m3.isrm()))
            DoMultUD<true>(m1,m2.cmView(),m3);
        else if (m2.step() != 1 || x != T(1))
            GenInstAddMultMM(T(1),m1,(x*m2).calc().constView().xView(),m3);
        else if (!(m1.isrm() || m1.iscm()))
            GenInstAddMultMM(x,m1.copy().constView().xdView(),m2,m3);
        else { // !(m3.isrm() || m3.iscm())
            typename M3::copy_type m3c(m3.size());
            DoMultUD<false>(m1,m2.cmView(),m3c);
            InstAddMultXM(T(1),m3c.constView().xdView(),m3.viewAsNonUnitDiag());
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,NonUnitDiag> m3)
    { GenInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,NonUnitDiag> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2,
        LowerTriMatrixView<T3,NonUnitDiag> m3)
    { GenInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2,
        LowerTriMatrixView<T3,NonUnitDiag> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }


#define InstFile "TMV_MultUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


