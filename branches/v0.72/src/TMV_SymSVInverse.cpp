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

#include "TMV_SymSVDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include <ostream>

namespace tmv {


    template <class T, class T1> 
    void HermSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<TMV_RealType(T1)>& SS,
        ptrdiff_t kmax, SymMatrixView<T> sinv)
    {
        TMVAssert(sinv.isherm());
        Matrix<T,ColMajor> SinvUt = U.adjoint().rowRange(0,kmax) /
            SS.subDiagMatrix(0,kmax);
        SymMultMM<false>(T(1),U.colRange(0,kmax),SinvUt,sinv);
    }

    template <class T, class T1> 
    void SymSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<TMV_RealType(T1)>& SS, 
        const GenMatrix<T1>& Vt, ptrdiff_t kmax, SymMatrixView<T> sinv)
    {
        // A = U S Vt
        // A^-1 = V S^-1 Ut
        Matrix<T,ColMajor> SinvUt = U.adjoint().rowRange(0,kmax) /
            SS.subDiagMatrix(0,kmax);
        SymMultMM<false>(T(1),Vt.adjoint().colRange(0,kmax),SinvUt,sinv);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymSVInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


