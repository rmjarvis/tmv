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
#ifndef TMV_LUDiv_H
#define TMV_LUDiv_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T, class T1> 
    void LU_LDivEq(
        const GenMatrix<T1>& LUx, const int* P, const MatrixView<T>& m);

    template <class T, class T1> 
    void LU_RDivEq(
        const GenMatrix<T1>& LUx, const int* P, const MatrixView<T>& m);

    template <class T, class T1> 
    void LU_Inverse(
        const GenMatrix<T1>& LUx, const int* P, const MatrixView<T>& m);

    // Specialize disallowed complex combinations:
#define CT std::complex<T>
    template <class T> 
    inline void LU_LDivEq(
        const GenMatrix<CT>& , const int* , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <class T> 
    inline void LU_RDivEq(
        const GenMatrix<CT>& , const int* , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <class T> 
    inline void LU_Inverse(
        const GenMatrix<CT>& , const int* , const MatrixView<T>& )
    { TMVAssert(TMV_FALSE); }
#undef CT

}

#endif
