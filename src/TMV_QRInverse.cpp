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


#include "TMV_QRDiv.h"
#include "TMV_Matrix.h"
#include "TMV_Vector.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T, class T1> void QR_Inverse(
      const GenMatrix<T>& QRx, const GenVector<T>& beta,
      const MatrixView<T1>& minv)
  {
    // minv = R^-1 Qt
    TMVAssert(minv.colsize() == QRx.rowsize());
    TMVAssert(minv.rowsize() == QRx.colsize());

    const int N = QRx.rowsize();

    minv.Zero();
    UpperTriMatrixView<T1> R = UpperTriMatrixViewOf(minv.Cols(0,N));
    R = UpperTriMatrixViewOf(QRx);
    R.InvertSelf();
    Q_RDivEq(QRx,beta,minv);
  }

#define InstFile "TMV_QRInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


