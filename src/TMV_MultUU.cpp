//////////////////////////////////////////////////////////////////////////////
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

#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_ProdMM.h"

namespace tmv {

    template <class M2, class M3>
    void DoMultEqUU(const M2& m2, M3& m3)
    {
        TMVAssert(m2.isrm() || m2.iscm());
        TMVAssert(m3.isrm() || m3.iscm());
        if (m3.iscm()) {
            typename M3::cmview_type m3cm = m3.cmView();
            if (m2.iscm())
                InlineMultMM<false>(m3cm,m2.cmView(),m3cm);
            else 
                InlineMultMM<false>(m3cm,m2.rmView(),m3cm);
        } else {
            typename M3::rmview_type m3rm = m3.rmView();
            if (m2.iscm())
                InlineMultMM<false>(m3rm,m2.cmView(),m3rm);
            else
                InlineMultMM<false>(m3rm,m2.rmView(),m3rm);
        }
    }

    template <class T, class M1, class M2, class M3>
    void GenInstMultMM(const T x, const M1& m1, const M2& m2, M3& m3)
    {
        if (m2.iscm() || m2.isrm()) {
            if (m3.iscm() || m3.isrm()) {
                if (m3.isunit()) {
                    TMVAssert(x == T(1));
                    TMVAssert(m1.isunit());
                    InstCopy(m1,m3);
                } else {
                    InstMultXM(x,m1,m3.viewAsNonUnitDiag());
                }
                DoMultEqUU(m2,m3);
            } else {
                if (m3.isunit()) {
                    typename M3::copy_type m3c(m3.size());
                    typename M3::copy_type::unitdiag_type::xdview_type m3x = 
                        m3c.viewAsUnitDiag().xdView();
                    TMVAssert(x == T(1));
                    TMVAssert(m1.isunit());
                    InstCopy(m1,m3x);
                    DoMultEqUU(m2,m3x);
                    InstCopy(m3x.constView(),m3);
                } else {
                    typename M3::copy_type m3c(m3.size());
                    typename M3::copy_type::xdview_type m3x = m3c.xdView();
                    InstMultXM(x,m1,m3c.xView());
                    DoMultEqUU(m2,m3x);
                    InstCopy(m3x.constView(),m3);
                }
            }
        } else {
            GenInstMultMM(x,m1,m2.copy().constView().xdView(),m3);
        }
    }


    template <class T, class M1, class M2, class M3>
    void GenInstAddMultMM(const T x, const M1& m1, const M2& m2, M3& m3)
    {
        if (m2.iscm() || m2.isrm()) {
            typename M3::copy_type m3c(m3.size());
            InstMultXM(x,m1,m3c.xView());
            DoMultEqUU(m2,m3c);
            InstAddMultXM(T(1),m3c.constView().xdView(),m3);
        } else {
            GenInstAddMultMM(x,m1,m2.copy().constView().xdView(),m3);
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,UnknownDiag> m3)
    { GenInstMultMM(x,m1,m2,m3); }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,NonUnitDiag> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        LowerTriMatrixView<T3,UnknownDiag> m3)
    { GenInstMultMM(x,m1,m2,m3); }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        LowerTriMatrixView<T3,NonUnitDiag> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }


#define InstFile "TMV_MultUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


