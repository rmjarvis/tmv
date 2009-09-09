///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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


#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <ostream>

namespace tmv {

  //
  // LDivEq
  //

  template <class T1, class T2> void CH_LDivEq(
      const GenSymMatrix<T1>& LL, const MatrixView<T2>& m)
  {
    TMVAssert(LL.size() == m.colsize());
    // m = (LLt)^-1 m
    //   = Lt^-1 L^-1 m
    m /= LL.LowerTri();
    m /= LL.UpperTri();
  }

  //
  // RDivEq Matrix
  //

  template <class T1, class T2> void CH_RDivEq(
      const GenSymMatrix<T1>& LL, const MatrixView<T2>& m)
  {
    TMVAssert(LL.size() == m.rowsize());
    // m = m (LLt)^-1 
    //   = m Lt^-1 L^-1
    m %= LL.UpperTri();
    m %= LL.LowerTri();
  }

#define InstFile "TMV_SymCHDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


