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
#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_BaseMatrix.h"
#include "TMV_BaseTriMatrix.h"

namespace tmv {

  // Decompose A (input as QRx) into Q R.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors, and the beta vector.
  // R is the upper triangle part of QRx
  template <class T> void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& beta, T& signdet);

  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, T& signdet);

  template <class T> void GetQFromQR(const MatrixView<T>& Q,
      const GenVector<T>& beta);

  template <class T, class T1> void Q_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m);
  template <class T, class T1> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m);
  template <class T, class T1> inline void Q_LMultEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { Q_RDivEq(Q,beta,m.Adjoint()); }
  template <class T, class T1> inline void Q_RMultEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { Q_LDivEq(Q,beta,m.Adjoint()); }


  template <class T, class T1, class T2> void QR_LDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const int* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, int N1=0);
  template <class T, class T1> void QR_LDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const int* P,
      const MatrixView<T>& m, int N1=0);
  template <class T, class T1, class T2> void QR_RDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const int* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, int N1=0);
  template <class T, class T1> void QR_RDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const int* P,
      const MatrixView<T>& m, int N1=0);
  template <class T, class T1> void QR_Inverse(
      const GenMatrix<T>& QRx, const GenVector<T>& beta,
      const MatrixView<T1>& minv);

}

#endif
