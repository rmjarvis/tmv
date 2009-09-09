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



#include "TMV_BandQRDiv.h"
#include "TMV_BandQRD.h"
#include "TMV_BandMatrix.h"
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

  template <class T> void QR_Decompose(
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
    int endcol = QRx.nlo()+1;
    int endrow = QRx.nhi()+1;
    T* Qbj = Qbeta.ptr();
    const int M = QRx.colsize();
    const int N = QRx.rowsize();
    for(int j=0;j<N;++j,++Qbj) {
      // Apply the Householder Reflection for this column
#ifdef TMVFLDEBUG
      TMVAssert(Qbj >= Qbeta.first);
      TMVAssert(Qbj < Qbeta.last);
#endif
      *Qbj = Householder_Reflect(QRx.SubMatrix(j,endcol,j,endrow),det);
      if (endcol < M) ++endcol;
      if (endrow < N) ++endrow;
    }
#ifdef XDEBUG
    Matrix<T> Q(M,N);
    Q = QRx;
    GetQFromBandQR(Q.View(),Qbeta,QRx.nlo());
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

  template <class T> void QR_Decompose(
      const GenBandMatrix<T>& A,
      const MatrixView<T>& Q, const BandMatrixView<T>& R)
  {
    // Decompose A = Q R 
    // where Q is unitary and R is upper banded

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(Q.colsize() == A.colsize());
    TMVAssert(Q.rowsize() == A.rowsize());
    TMVAssert(R.colsize() == A.rowsize());
    TMVAssert(R.rowsize() == A.rowsize());
    TMVAssert(R.nlo() >= 0);
    TMVAssert(R.nhi() == int(R.rowsize())-1 || R.nhi() >= A.nlo() + A.nhi());

    if (Q.isconj()) {
      QR_Decompose(A.Conjugate(),Q.Conjugate(),R.Conjugate());
    } else {
      Vector<T> beta(A.rowsize());
      T d(0);
      Q.Zero();
      BandMatrixView<T>(Q,A.nlo(),A.nhi()) = A;
      QR_Decompose(BandMatrixViewOf(Q,A.nlo(),A.nlo()+A.nhi()),
	  beta.View(),d);
      R = BandMatrixViewOf(Q,0,A.nlo()+A.nhi());
      GetQFromBandQR(Q,beta.View(),A.nlo());
    }
  }

  template <class T> void QR_Decompose(
      const GenBandMatrix<T>& A, const BandMatrixView<T>& R)
  {
    // Decompose A = Q R 
    // where Q is unitary and R is upper banded

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(R.colsize() == A.rowsize());
    TMVAssert(R.rowsize() == A.rowsize());
    TMVAssert(R.nlo() >= 0);
    TMVAssert(R.nhi() == int(R.rowsize())-1 || R.nhi() >= A.nlo() + A.nhi());

    Vector<T> beta(A.rowsize());
    T d(0);
    BandMatrix<T> QR(MIN(A.colsize(),A.rowsize()+A.nlo()),A.rowsize(),
	A.nlo(),A.nlo()+A.nhi(),T(0));
    QR.Zero();
    BandMatrixViewOf(QR,A.nlo(),A.nhi()) = A.Rows(0,QR.colsize());
    QR_Decompose(QR.View(),beta.View(),d);
    R = BandMatrixViewOf(QR,0,A.nlo()+A.nhi());
  }

#define InstFile "TMV_BandQRDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


