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

#include "tmv/TMV_DivUU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyU.h"

namespace tmv {

    template <class M1, class M2>
    static void DoLDivEq2(M1& m1, const M2& m2)
    {
        TMVAssert(m1.isrm() || m1.iscm());
        TMVAssert(m2.isrm() || m2.iscm());
        if (m2.iscm()) {
            typename M1::cmview_type m1cm = m1.cmView();
            if (m2.iscm())
                InlineLDivEq(m1cm,m2.cmView());
            else
                InlineLDivEq(m1cm,m2.rmView());
        } else {
            typename M1::rmview_type m1rm = m1.rmView();
            if (m2.iscm())
                InlineLDivEq(m1rm,m2.cmView());
            else
                InlineLDivEq(m1rm,m2.rmView());
        }
    }

    template <class M1, class M2>
    static void DoLDivEq(M1& m1, const M2& m2)
    {
        if (m1.iscm() || m1.isrm()) {
            if (m2.iscm() || m2.isrm()) {
                DoLDivEq2(m1,m2);
            } else {
                if (m2.isunit()) 
                    DoLDivEq2(
                        m1,m2.copy().viewAsUnitDiag().constView().xdView());
                else 
                    DoLDivEq2(m1,m2.copy().constView().xdView());
            }
        } else {
            if (m1.isunit()) {
                typename M1::copy_type m1c(m1);
                typename M1::copy_type::unitdiag_type::xdview_type m1cv = 
                    m1c.viewAsUnitDiag().xdView();
                DoLDivEq(m1cv,m2);
                InstCopy(m1cv.constView(),m1);
            } else {
                typename M1::copy_type m1c(m1);
                typename M1::copy_type::xdview_type m1cv = m1c.xdView();
                DoLDivEq(m1cv,m2);
                InstCopy(m1cv.constView(),m1);
            }
        }
    }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        UpperTriMatrixView<T1> m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoLDivEq(m1,m2); }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        LowerTriMatrixView<T1> m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoLDivEq(m1,m2); }

#define InstFile "TMV_DivUU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


