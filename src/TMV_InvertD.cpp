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
#include "tmv/TMV_Det.h"
#include "tmv/TMV_DiagMatrixIO.h"
#include "tmv/TMV_MatrixIO.h"

namespace tmv {

    //
    // Element-wise division
    //

    template <class T>
    void InstElemInvert(VectorView<T> v)
    {
        if (v.step() == 1) {
            VectorView<T,1> vunit = v.unitView();
            InlineElemInvert(vunit);
        } else {
            InlineElemInvert(v);
        }
    }

#define InstFile "TMV_InvertD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


