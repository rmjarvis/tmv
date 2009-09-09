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
// LU Decomposition.
//
// The basics of LU decomposition for band matrices are the same as 
// for regular matrices.  However, there are a few wrinkles about doing
// it efficiently.  
//
// We leave the details to the comments in TMV_BandLUDiv.cpp, but 
// the main difference for the routines in this file is that L can
// be stored in a lower band matrix with m.nlo() subdiagonals.  
// However, U needs m.nlo() + m.nhi() superdiagonals for its storage.
//
//


#ifndef TMV_BandLUDiv_H
#define TMV_BandLUDiv_H

#include "TMV_Divider.h"
#include "TMV_BandMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T> void BandLU_Decompose(
      const BandMatrixView<T>& LUx, size_t* p, T& det, int Anhi);
  // Decompose A (input as LUx) into L * U
  
  template <class T, class Ta> void BandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B, DiagType dt);
  template <class T, class Ta> void BandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& v, DiagType dt);
  // Solve A X = B  where A is upper or lower band triangular

  template <class T, class T1> void BandLU_LDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* p, const MatrixView<T>& m);
  template <class T, class T1> void BandLU_RDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* p, const MatrixView<T>& m);

  template <class T> class BandLUDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      BandLUDiv(const GenBandMatrix<T>& A, bool _inplace);
      BandLUDiv(const AssignableToBandMatrix<T>& A);
      ~BandLUDiv();

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

      inline bool IsTrans() const;
      inline ConstBandMatrixView<T> GetU() const;
      inline LowerTriMatrix<T,UnitDiag> GetL() const;
      inline const GenBandMatrix<T>& GetLU() const;
      const size_t* GetP() const;

      std::string Type() const;
      DivType GetDivType() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

      struct BandLUDiv_Impl;
      BandLUDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
  };

} // namespace mv

#endif
