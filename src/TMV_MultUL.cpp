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

#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdXM.h"

namespace tmv {

    template <class M1, class M2, class M3>
    void DoMultUL(const M1& m1, const M2& m2, M3& m3)
    {
        // We are allowed to have m2 in the same storage as m3, so 
        // we use that fact to save a bit on the rm, cm options,
        // by explicitly copying m2 to the corresponding portion of m3.
        TMVAssert(m1.isrm() || m1.iscm());
        TMVAssert(m3.isrm() || m3.iscm());

        typedef typename TypeSelect<M2::mupper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type::xdview_type M2x;
        M2x m2x = Maybe<M2::mupper>::uppertri(m3);

        if (m2.isunit()) {
            typename M2x::unitdiag_type::xdview_type m2u = m2x.viewAsUnitDiag();
            AliasCopy(m2.viewAsUnitDiag(),m2u);
            if (m3.iscm()) {
                typename M3::cmview_type m3cm = m3.cmView();
                if (m1.iscm())
                    InlineMultMM<false>(m1.cmView(),m2u.cmView(),m3cm);
                else
                    InlineMultMM<false>(m1.rmView(),m2u.cmView(),m3cm);
            } else {
                typename M3::rmview_type m3rm = m3.rmView();
                if (m1.iscm())
                    InlineMultMM<false>(m1.cmView(),m2u.rmView(),m3rm);
                else
                    InlineMultMM<false>(m1.rmView(),m2u.rmView(),m3rm);
            }
        } else {
            AliasCopy(m2,m2x);
            if (m3.iscm()) {
                typename M3::cmview_type m3cm = m3.cmView();
                if (m1.iscm())
                    InlineMultMM<false>(m1.cmView(),m2x.cmView(),m3cm);
                else
                    InlineMultMM<false>(m1.rmView(),m2x.cmView(),m3cm);
            } else {
                typename M3::rmview_type m3rm = m3.rmView();
                if (m1.iscm())
                    InlineMultMM<false>(m1.cmView(),m2x.rmView(),m3rm);
                else
                    InlineMultMM<false>(m1.rmView(),m2x.rmView(),m3rm);
            }
        }
    }

    template <class T, class M1, class M2>
    void GenInstMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        if (!(m1.isrm() || m1.iscm())) 
            GenInstMultMM(x,m1.copy().constView().xdView(),m2,m3);
        else if (!(m3.iscm() || m3.isrm())) {
            Matrix<T,ColMajor> m3c(m3.colsize(),m3.rowsize());
            DoMultUL(m1,m2,m3c);
            InstMultXM(x,m3c.constView().xView(),m3);
        } else {
            DoMultUL(m1,m2,m3);
            InstScale(x,m3);
        }
    }

    template <class T, class M1, class M2>
    void GenInstAddMultMM(
        const T x, const M1& m1, const M2& m2, MatrixView<T>& m3)
    {
        if (!(m1.isrm() || m1.iscm()))
            GenInstAddMultMM(x,m1.copy(),m2,m3);
        else {
            Matrix<T,ColMajor> m3c(m3.colsize(),m3.rowsize());
            DoMultUL(m1,m2,m3c);
            InstAddMultXM(x,m3c.constView().xView(),m3);
        }
    }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { GenInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { GenInstMultMM(x,m1,m2,m3); }
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3)
    { GenInstAddMultMM(x,m1,m2,m3); }


#define InstFile "TMV_MultUL.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


