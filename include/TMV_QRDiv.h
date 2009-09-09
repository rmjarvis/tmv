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


//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// QR Decomposition.
//
// The basic idea of an QR decomposition is that any 
// matrix A can be decomposed into a unitary matrix
// times an upper triangle matrix.
//
// A = Q R
//
// We do this using Householder transformations, which can
// be stored in a lower triangle matrix, thus other than the 
// diagonal, which they both need, we can store them in the 
// place of the original matrix. It is more convenient to 
// keep the diagonal of Q in place and take the diagonal of R
// separately.
//
// If R is not singular, the solution to A x = b is found from
// Q R x = b
// R x = Qt b
// which can be solved by back-substitutaion.
//
// If m > n, this does not actually give a solution to A x = b,
// since Q Qt != I  (Q is only column orthogonal if m > n.)
// But it does give the value of x which minimizes the 2-norm
// of (A x - b)
//
// The 2-norm of v is the square root of vt v, so
// |A x - b|^2 =
// (A x - b)t (A x - b) = 
// (xt At - bt) (A x - b) =
// xt At A x - bt A x - xt At b + bt b =
// xt Rt Qt Q R x - bt Q R x - xt Rt Qt b + bt b =
// |R x - Qt b|^2 + |b|^2 - |Qt b|^2
// Clearly the x which minimizes this is the solution of R x = Qt b.
//
// If R is singular, then you need QRP Decomposition (see TMV_QRPDiv.h).
//


#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_Divider.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T> void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& beta, T& det);
  template <class T> inline void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& beta)
  { T d=0; QR_Decompose(QRx,beta,d); }
  // Decompose A (input as QRx) into Q R.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors, and the beta vector.
  // R is the upper triangle part of QRx
 
  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, T& det);
  template <class T> inline void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R)
  { T d=0; QR_Decompose(Q,R,d); }
  // Decompose A (input as Q) into Q R.
 
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

  template <class T> void QR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A);
  // Given that A0 = Q0 R0
  // Find R1, so that [ A0 ] = Q1 R1
  //                  [ A  ] 
  // On input R is R0; on output it is R1.
  template <class T> inline void QR_Update(
      const UpperTriMatrixView<T>& R, const VectorView<T>& z)
  { return QR_Update(R,RowVectorViewOf(z)); }
  template <class T> void QR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A);
  // The opposite of an Update:
  // Given that [ A0 ] = Q1 R1
  //            [ A  ]
  // Find R0, so that A0 = Q0 R0
  // On input R is R1; on output it is R0.
  template <class T> inline void QR_Downdate(
      const UpperTriMatrixView<T>& R, const VectorView<T>& z)
  { QR_Downdate(R,RowVectorViewOf(z)); }

  template <class T, class T1, class T2> void QR_LDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, size_t N1=0);
  template <class T, class T1> void QR_LDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const MatrixView<T>& m, size_t N1=0);
  template <class T, class T1, class T2> void QR_RDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, size_t N1=0);
  template <class T, class T1> void QR_RDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const MatrixView<T>& m, size_t N1=0);

  template <class T> class QRDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      QRDiv(const GenMatrix<T>& A, bool _inplace);
      ~QRDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const;

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool Singular() const;

      //
      // Access Decomposition
      //

      bool IsTrans() const;
      Matrix<T> GetQ() const;
      ConstUpperTriMatrixView<T> GetR() const;
      const GenMatrix<T>& GetQR() const;
      const GenVector<T>& Getbeta() const;

      std::string Type() const;
      DivType GetDivType() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected :

      struct QRDiv_Impl;
      QRDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
  };

} // namespace mv

#endif
