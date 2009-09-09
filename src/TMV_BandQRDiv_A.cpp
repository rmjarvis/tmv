///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_BandMatrix.h"
#include "TMV_BandQRDiv.h"
#include "TMV_Householder.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> void BandQR_Decompose(
      const BandMatrixView<T>& QRx,
      const VectorView<T>& Qbeta, T& det)
  {
    // Decompose A (input as QRx) into A = Q R 
    // where Q is unitary, and R is upper triangular
    // Q and R are stored in the same matrix (QRx)
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(!Qbeta.isconj());
    TMVAssert(Qbeta.step() == 1);

#ifdef XDEBUG
    Matrix<T> A0(QRx);
#endif
    size_t endcol = QRx.nlo()+1;
    size_t endrow = QRx.nhi()+1;
    T* Qbj = Qbeta.ptr();
    for(size_t j=0;j<QRx.rowsize();++j,++Qbj) {
      // Apply the Householder Reflection for this column
      *Qbj = Householder_Reflect(QRx.SubMatrix(j,endcol,j,endrow),det);
      if (endcol < QRx.colsize()) ++endcol;
      if (endrow < QRx.rowsize()) ++endrow;
    }
#ifdef XDEBUG
    Matrix<T> Q(GetQFromBandQR(QRx,Qbeta));
    Matrix<T> R(QRx.Diags(0,QRx.nhi()+1));
    Matrix<T> QR = Q*R;
    if (Norm(QR-A0) > 0.00001*Norm(A0)) {
      cerr<<"BandQR_Decompose: \n";
      cerr<<"A = "<<Type(QRx)<<"  "<<A0<<endl;
      cerr<<"QRx = "<<QRx<<endl;
      cerr<<"Q = "<<Q<<endl;
      cerr<<"R = "<<R<<endl;
      cerr<<"QR = "<<QR<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandQRDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


