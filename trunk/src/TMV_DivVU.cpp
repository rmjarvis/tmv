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


#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MultXV.h"

namespace tmv {

    template <class T1, class M2>
    void DoInstLDivEq(VectorView<T1> v1, const M2& m2)
    {
        if (v1.step() == 1) {
            VectorView<T1,1> v1u = v1.unitView();
            if (m2.iscm()) InlineLDivEq(v1u,m2.cmView());
            else if (m2.isrm()) InlineLDivEq(v1u,m2.rmView());
            else InlineLDivEq(v1u,m2);
        } else {
            InlineLDivEq(v1,m2);
        }
    }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        VectorView<T1> v1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoInstLDivEq(v1,m2); }

    template <class T1, class T2, bool C2>
    void InstLDivEq(
        VectorView<T1> v1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2)
    { DoInstLDivEq(v1,m2); }

#define InstFile "TMV_DivVU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


