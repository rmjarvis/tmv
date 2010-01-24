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


#include "tmv/TMV_InvertD.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_MultXD.h"

namespace tmv {

    //
    // Element-wise division
    //

    template <class T1, bool C1, class T2>
    void InstInvert(
        const T2 x,
        const ConstDiagMatrixView<T1,UNKNOWN,C1>& m1, DiagMatrixView<T2> m2)
    {
        typedef typename Traits<T2>::real_type RT;
        Scaling<1,RT> one;
        if (m1.step() == 1 && m2.step() == 1) {
            ConstDiagMatrixView<T1,1,C1> m1unit = m1.cmView();
            DiagMatrixView<T2,1> m2unit = m2.cmView();
            InlineInvert(one,m1unit,m2unit);
        } else 
            InlineInvert(one,m1,m2);
        InstScale(x,m2.diag());
    }

    template <class T1, bool C1, class T2>
    void InstElemInvert(
        const T2 x,
        const ConstVectorView<T1,UNKNOWN,C1>& v1, VectorView<T2> v2)
    {
        typedef typename Traits<T2>::real_type RT;
        Scaling<1,RT> one;
        if (v1.step() == 1 && v2.step() == 1) {
            ConstVectorView<T1,1,C1> v1unit = v1.unitView();
            VectorView<T2,1> v2unit = v2.unitView();
            InlineElemInvert(one,v1unit,v2unit);
        } else 
            InlineElemInvert(one,v1,v2);
        InstScale(x,v2);
    }

#define InstFile "TMV_InvertD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


