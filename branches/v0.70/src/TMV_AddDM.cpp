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



#include "tmv/TMV_DiagMatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"

namespace tmv {

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
        // C = alpha*A + beta*B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == B.rowsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == C.rowsize());

        if (A.size() > 0) {
            if (SameStorage(A.diag(),C)) {
                DiagMatrix<Ta> tempA = A;
                C = beta * B;
                AddMM(alpha,tempA,C);
            } else {
                C = beta * B;
                AddMM(alpha,A,C);
            }
        }
    }

#define InstFile "TMV_AddDM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


