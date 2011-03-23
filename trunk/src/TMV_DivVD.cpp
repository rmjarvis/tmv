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

//#define PRINTALGO_DIVVD

#include "tmv/TMV_DivVD.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_ScaleV.h"

namespace tmv {

    //
    // Element-wise division
    //

    template <class T1, int C1, class T2, int C2, class T3>
    void InstElemDivVV(
        const T3 x,
        const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
        typedef typename Traits<T3>::real_type RT;
        const Scaling<1,RT> one;
        if (v1.step() == 1 && v2.step() == 1 && v3.step() == 1) {
            ConstVectorView<T1,C1|Unit> v1unit = v1.unitView();
            ConstVectorView<T2,C2|Unit> v2unit = v2.unitView();
            VectorView<T3,Unit> v3unit = v3.unitView();
            InlineElemDivVV(one,v1unit,v2unit,v3unit);
        } else 
            InlineElemDivVV(one,v1,v2,v3);
        InstScale(x,v3);
    }

#define InstFile "TMV_DivVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


