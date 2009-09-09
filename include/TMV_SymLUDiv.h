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
// This file contains the code for doing division of Symmetric/Hermitian
// matrices using LDLt Decomposition.  (Although it is still referred
// to as LU decomposition for the purposes of the DivideUsing method.)
//
// This Decomposition is appropriate for symmetric matrices which may
// not be positive definite (ie. they may have negative eigenvalues).
// The matrix is decomposed into P * L * D * Lt * Pt.
//
// P is a permutation.
// L is a unit diagonal lower triangle matrix.
// D is a pseudo-diagonal matrix - meaning that there may be 2x2 blocks
//   occasionally along the diagonal.
// Lt is L.Adjoint for Hermitian matrices, or L.Transpose for symmetric.
//
// The determinant of A is just the determinant of D.
//


#ifndef TMV_SYMLUDiv_H
#define TMV_SYMLUDiv_H

#include "TMV_SymDivider.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T> void SymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      size_t* P, T& det);
  // Decompose A into P * L * D * Lt * Pt
  
  template <class T, class T1> void SymLU_LDivEq(
      const GenSymMatrix<T1>& L, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m);
  template <class T, class T1> void SymLU_RDivEq(
      const GenSymMatrix<T1>& L, const GenVector<T1>& xD, const size_t* P, 
      const MatrixView<T>& m);
  template <class T, class T1> void SymLU_Inverse(
      const GenSymMatrix<T>& L, const GenVector<T>& xD, const size_t* P,
      const SymMatrixView<T1>& sinv);

  template <class T> class SymLUDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      SymLUDiv(const GenSymMatrix<T>& A, bool _inplace);
      ~SymLUDiv();

      //
      // Divider Versions of DivEq and Div
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;
      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;

      //
      // Determinant, Inverse
      //

      T Det() const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      template <class T1> void DoInverse(const SymMatrixView<T1>& minv) const;
      inline void Inverse(const SymMatrixView<RealType(T)>& sinv) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(sinv.size() == colsize());
	DoInverse(sinv);
      }
      inline void Inverse(const SymMatrixView<ComplexType(T)>& sinv) const
      {
	TMVAssert(sinv.size() == colsize());
	TMVAssert(sinv.isherm() == isherm());
	TMVAssert(sinv.issym() == issym());
	DoInverse(sinv);
      }
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool Singular() const;

#include "TMV_AuxAllDiv.h"

      //
      // Access Decomposition
      //

      const ConstLowerTriMatrixView<T> GetL() const;
      const Matrix<T> GetD() const;
      const size_t* GetP() const;
      const GenSymMatrix<T>& GetLL() const;
      const GenVector<T>& GetxD() const;

      std::string Type() const;
      DivType GetDivType() const;
      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :
      struct SymLUDiv_Impl;
      SymLUDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
      bool isherm() const;
      bool issym() const;
  };

} // namespace tmv

#endif
