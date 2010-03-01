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

#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_MultMM.h"

#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MultXM.h"

namespace tmv {

    template <class T1, class M2>
    static void DoInstLDivEq(MatrixView<T1> m1, const M2& m2)
    {
        if (m2.iscm() || m2.isrm()) {
            if (m1.iscm()) {
                MatrixView<T1,1> m1cm = m1.cmView();
                if (m2.iscm()) InlineLDivEq(m1cm,m2.cmView());
                else InlineLDivEq(m1cm,m2.rmView());
            } else if (m1.isrm()) {
                MatrixView<T1,UNKNOWN,1> m1rm = m1.rmView();
                if (m2.iscm()) InlineLDivEq(m1rm,m2.cmView());
                else InlineLDivEq(m1rm,m2.rmView());
            } else {
                Matrix<T1,RowMajor> m1c(m1.colsize(),m1.rowsize());
                MatrixView<T1,UNKNOWN,1> m1rm = m1c.rmView();
                InstCopy(m1.constView(),m1c.xView());
                InlineLDivEq(m1rm,m2);
                InstCopy(m1c.constView().xView(),m1);
            }
        } else {
            DoInstLDivEq(m1,m2.copy().xView());
        }
    }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        MatrixView<T1> m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoInstLDivEq(m1,m2); }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        MatrixView<T1> m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoInstLDivEq(m1,m2); }

#define InstFile "TMV_DivMU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


