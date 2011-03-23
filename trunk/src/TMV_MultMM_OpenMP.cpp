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

#include "TMV_Blas.h"
#include "tmv/TMV_MultMM_OpenMP.h"
#include "tmv/TMV_Matrix.h"

namespace tmv {

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM_OpenMP(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && TMV_OPT > 3 && defined(_OPENMP)
        InlineMultMM_OpenMP<false>(Scaling<0,T3>(x),m1,m2,m3); 
#else
        InstMultMM(x,m1,m2,m3);
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM_OpenMP(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3)
    {
#if !defined(BLAS) && TMV_OPT > 0 && defined(_OPENMP)
        InlineMultMM_OpenMP<true>(Scaling<0,T3>(x),m1,m2,m3); 
#else
        InstAddMultMM(x,m1,m2,m3);
#endif
    }

#define InstFile "TMV_MultMM_OpenMP.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


