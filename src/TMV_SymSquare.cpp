///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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



#include "TMV_SymSquare.h"
#include "TMV_Matrix.h"
#include "TMV_MatrixArith.h"

namespace tmv {

  template <bool herm, class T> void SymSquare(const MatrixView<T>& A)
  {
    const int N = A.colsize();
    if (N == 1) {
      const T A00 = *A.ptr();
#ifdef TMVFLDEBUG
      TMVAssert(A.ptr() >= A.first);
      TMVAssert(A.ptr() < A.last);
#endif
      if (herm)
	*A.ptr() = NORM(REAL(A00));
      else 
	*A.ptr() = SQR(A00);
    } else {
      const int K = N/2;
      MatrixView<T> A00 = A.SubMatrix(0,K,0,K);
      MatrixView<T> A10 = A.SubMatrix(K,N,0,K);
      MatrixView<T> A01 = A.SubMatrix(0,K,K,N);
      MatrixView<T> A11 = A.SubMatrix(K,N,K,N);
      MatrixView<T> A10t = herm ? A10.Adjoint() : A10.Transpose();

      // [ A00 A10t ] [ A00 A10t ] 
      // [ A10 A11  ] [ A10 A11  ]
      // = [ A00^2 + A10t A10    A00 A10t + A10t A11 ]
      //   [ A10 A00 + A11 A10   A10 A10t + A11^2    ]

      // A10 stores the actual data for A10
      // We can therefore write to A01 as a temp matrix.
      
      A01 = A00 * A10t;
      A01 += A10t * A11;

      SymSquare<herm>(A00);
      A00 += A10t*A10;
      SymSquare<herm>(A11);
      A11 += A10*A10t;

      A10t = A01;
    }
  }
  
#define InstFile "TMV_SymSquare.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


