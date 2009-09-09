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

#include "TMV_QRDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

  template <class T, class T1> void QR_Inverse(
      const GenMatrix<T>& QRx, const GenVector<T>& beta,
      const MatrixView<T1>& minv)
  {
    //std::cout<<"Start QR_Inverse"<<std::endl;
    //std::cout<<"QRx = "<<TypeText(QRx)<<std::endl;
    //std::cout<<"beta = "<<TypeText(beta)<<std::endl;
    //std::cout<<"minv = "<<TypeText(minv)<<std::endl;
    // minv = R^-1 Qt
    TMVAssert(minv.colsize() == QRx.rowsize());
    TMVAssert(minv.rowsize() == QRx.colsize());

    const int N = QRx.rowsize();

    minv.Zero();
    //std::cout<<"After Zero"<<std::endl;
    UpperTriMatrixView<T1> R = minv.Cols(0,N).UpperTri();
    //std::cout<<"After Make R"<<std::endl;
    R = QRx.UpperTri();
    //std::cout<<"After Set R"<<std::endl;
    R.InvertSelf();
    //std::cout<<"After InvertSelf"<<std::endl;
    //std::cout<<"R = "<<R<<std::endl;
    //std::cout<<"QRx = "<<QRx<<std::endl;
    //std::cout<<"beta = "<<beta<<std::endl;
    //std::cout<<"minv = "<<minv<<std::endl;
    Q_RDivEq(QRx,beta,minv);
    //std::cout<<"End QR_Inverse"<<std::endl;
  }

#define InstFile "TMV_QRInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


