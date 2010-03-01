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


#include <iostream>
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixIO.h"

namespace tmv {

    template <class T, bool C>
    void InstWrite(std::ostream& os, const ConstDiagMatrixView<T,UNKNOWN,C>& m)
    {
        if (m.step() == 1)
            InlineWrite(os,m.cmView()); 
        else
            InlineWrite(os,m); 
    }

    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstDiagMatrixView<T,UNKNOWN,C>& m,
        typename Traits<T>::real_type thresh)
    {
        if (m.step() == 1)
            InlineWrite(os,m.cmView(),thresh); 
        else
            InlineWrite(os,m,thresh); 
    }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


