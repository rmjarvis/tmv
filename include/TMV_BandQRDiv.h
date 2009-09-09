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
// This file contains the code for doing division of BandMatrices using 
// QR Decomposition.
//
// I don't implement the QRP method, since the permutations screw up the
// band structure.  It could be implemented using a similar technique 
// as I used for the BandLUDiv class, but if you want to use it you need 
// to cast the band matrix to a regular matrix first.
//


#ifndef TMV_BandQRDiv_H
#define TMV_BandQRDiv_H

#include "TMV_Divider.h"
#include "TMV_BandMatrix.h"
#include "TMV_Matrix.h"

namespace tmv {

  template <class T> void BandQR_Decompose(
      const BandMatrixView<T>& QRx, const VectorView<T>& Qbeta, T& det);
  template <class T> inline void BandQR_Decompose(
      const BandMatrixView<T>& QRx, const VectorView<T>& Qbeta)
  { T d; BandQR_Decompose(QRx,Qbeta,d); }
  // Decompose A (input as QRx) into Q R.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors.  The upper triangle part and Qbeta hold R.

  template <class T1, class T2> void BandQR_LDivEq(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);
  template <class T1, class T2> void BandQR_RDivEq(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);

  template <class T1, class T2, class T3> void BandQR_LDiv(
	const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta,
	const GenMatrix<T2>& m, const MatrixView<T3>& x);
  template <class T1, class T2, class T3> void BandQR_RDiv(
	const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta, 
	const GenMatrix<T2>& m, const MatrixView<T3>& x);

  template <class T1, class T2> void BandQ_LDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);
  template <class T1, class T2> void BandQ_RDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);
  template <class T> Matrix<T> GetQFromBandQR(
      const GenBandMatrix<T>& QRx, const GenVector<T>& Qbeta);

  template <class T> class BandQRDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      BandQRDiv(const GenBandMatrix<T>& A, bool _inplace);
      ~BandQRDiv();

      //
      // Div, DivEq
      //
      
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

      //
      // Determinant, Inverse
      //

      T Det() const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool Singular() const;

#include "TMV_AuxAllDiv.h"

      //
      // Access Decomposition
      //

      bool IsTrans() const;
      Matrix<T> GetQ() const;
      ConstBandMatrixView<T> GetR() const;
      const GenBandMatrix<T>& GetQR() const;
      const GenVector<T>& GetQbeta() const;

      std::string Type() const;
      DivType GetDivType() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

      struct BandQRDiv_Impl;
      BandQRDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
  };

} // namespace mv

#endif
