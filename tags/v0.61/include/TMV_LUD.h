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


//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// LU Decomposition.
//
// The name LU Decomposition is traditional, but somewhat
// annoying.  Usually U represents a unitary matrix, not an
// upper traiangular matrix.  The latter are usually represented
// with R.  (for "Right Triangular")  
// For example, in a QR decomposition, the R is an upper
// trianular matrix.  (Q also typically represents unitary matrices.)
// However, I will use U rather than R here, since that is 
// the usual representation in this context.
//
// The basic idea of an LU decomposition is that any 
// square matrix A can be decomposed into a lower triangular
// matrix time an upper triangular matrix.
//
// For stability reasons, we actually decompose a permutation
// of A instead, so:
//
// A = P L U
//
// Only one of L or U needs a non-unit diagonal, so we choose L to 
// have unit diagonal, and U to have the non-unit diagonal.
//
// This means that we can store L and U both in a square matrix
// the same size as A, with L being the elements below the diagonal
// and U being the elements above and including the diagonal.
//
// The determinant of A can be calculated easily from the LU
// decomposition:
//
// det(P) * det(A) = det(L) * det(U)
// +-1 * det(A) = 1 * det(U)
// As we calculate the decomposition, we keep track of whether
// det(P) is +-1 
// The determinant of U is just the product of the diagonal elements.
// So the determinant of A is just det(P) times the diagonal elements 
// of U.
//


#ifndef TMV_LUD_H
#define TMV_LUD_H

#include "TMV_Divider.h"
#include "TMV_BaseMatrix.h"
#include "TMV_BaseTriMatrix.h"

namespace tmv {

  // Decompose A into P * L * U
  // L is returned as LowerTriMatrixViewOf(A,UnitDiag).
  // U is returned as UpperTriMatrixViewOf(A).
  template <class T> void LU_Decompose(
      const MatrixView<T>& A, int* P);

  template <class T> class LUDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      LUDiv(const GenMatrix<T>& A, bool _inplace);
      ~LUDiv();

      //
      // Divider Versions of DivEq and Div
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;
      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      RealType(T) LogDet(T* sign) const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool Singular() const;

      //
      // Access Decomposition
      //

      bool IsTrans() const;
      ConstLowerTriMatrixView<T> GetL() const;
      ConstUpperTriMatrixView<T> GetU() const;
      const GenMatrix<T>& GetLU() const;
      const int* GetP() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

      struct LUDiv_Impl;
      LUDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;

    private :

      inline LUDiv(const LUDiv<T>&) : pimpl(0)
      { TMVAssert(FALSE); }
      inline LUDiv<T>& operator=(const LUDiv<T>&)
      { TMVAssert(FALSE); return *this; }

  };

} // namespace tmv

#endif
