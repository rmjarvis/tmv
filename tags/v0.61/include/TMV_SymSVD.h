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
#ifndef TMV_SymSVD_H
#define TMV_SymSVD_H

#include "TMV_SymDivider.h"
#include "TMV_BaseDiagMatrix.h"

namespace tmv {

  // Find EigenValues and EigenVectors of hermitian matrix, A.
  // For each lambda(i), A V.col(i) = lambda(i) V.col(i).
  // In other words, A * V = V * DiagMatrixViewOf(lambda)
  // Or, A = V * DiagMatrixViewOf(lambda) * V.Inverse()
  // Furthermore, since A is hermitian, V.Inverse() = V.Adjoint().
  // On input, lambda and V must have the same size as A.
  // On output, the lambda are sorted to be increasing in value.
  template <class T> void Eigen(
      const GenSymMatrix<T>& A,
      const MatrixView<T>& V, const VectorView<RealType(T)>& lambda);

  // The same, but don't return V
  // For this one, A is junk on output.
  template <class T> void Eigen(
      const SymMatrixView<T>& A, const VectorView<RealType(T)>& lambda);

  // Decompose A into U S V
  template <class T> void SV_Decompose(
      const GenSymMatrix<T>& A, const MatrixView<T>& U,
      const DiagMatrixView<RealType(T)>& S, const MatrixView<T>& V);

  // The same, but don't return U and/or V
  template <class T> void SV_Decompose(
      const GenSymMatrix<T>& A, 
      const MatrixView<T>& U, const DiagMatrixView<RealType(T)>& S);
  template <class T> void SV_Decompose(
      const GenSymMatrix<T>& A,
      const DiagMatrixView<RealType(T)>& S, const MatrixView<T>& V);
  // For this last one, A is junk on output.
  template <class T> void SV_Decompose(
      const SymMatrixView<T>& A, const DiagMatrixView<RealType(T)>& S);


  template <class T> class HermSVDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      HermSVDiv(const GenSymMatrix<T>& A);
      ~HermSVDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const;

      //
      // Determinant, Inverse
      //

      T Det() const;
      RealType(T) LogDet(T* sign) const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      template <class T1> void DoInverse(const SymMatrixView<T1>& sinv) const;
      inline void Inverse(const SymMatrixView<RealType(T)>& sinv) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(sinv.size() == colsize());
	DoInverse(sinv);
      }
      inline void Inverse(const SymMatrixView<ComplexType(T)>& sinv) const
      {
	TMVAssert(sinv.size() == colsize());
	TMVAssert(sinv.isherm());
	DoInverse(sinv);
      }
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool Singular() const;
      RealType(T) Norm2() const;
      RealType(T) Condition() const;

#include "TMV_AuxAllDiv.h"

      //
      // Determine which (if any) S values to zero out
      //

      void Thresh(RealType(T) toler, std::ostream* debugout=0) const;
      void Top(int neigen, std::ostream* debugout=0) const;
      int GetKMax() const;

      //
      // Access Decomposition
      //

      ConstMatrixView<T> GetU() const;
      DiagMatrix<RealType(T)> GetS() const;
      Matrix<T> GetV() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected:

      struct HermSVDiv_Impl;
      HermSVDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;

    private :

      inline HermSVDiv(const HermSVDiv<T>&) : pimpl(0)
      { TMVAssert(FALSE); }
      inline HermSVDiv<T>& operator=(const HermSVDiv<T>&)
      { TMVAssert(FALSE); return *this; }

  }; // HermSVDiv

  template <class T> class SymSVDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      SymSVDiv(const GenSymMatrix<T>& A);
      ~SymSVDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      RealType(T) LogDet(T* sign) const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      template <class T1> void DoInverse(const SymMatrixView<T1>& sinv) const;
      inline void Inverse(const SymMatrixView<RealType(T)>& sinv) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(sinv.size() == colsize());
	DoInverse(sinv);
      }
      inline void Inverse(const SymMatrixView<ComplexType(T)>& sinv) const
      {
	TMVAssert(sinv.size() == colsize());
	TMVAssert(sinv.issym());
	DoInverse(sinv);
      }
      void DoInverseATA(const MatrixView<T>& minv) const;

      bool Singular() const;
      RealType(T) Norm2() const;
      RealType(T) Condition() const;

      //
      // Determine which (if any) S values to zero out
      //

      void Thresh(RealType(T) toler, std::ostream* debugout=0) const;
      void Top(int neigen, std::ostream* debugout=0) const;
      int GetKMax() const;

      //
      // Access Decomposition
      //

      ConstMatrixView<T> GetU() const;
      ConstDiagMatrixView<RealType(T)> GetS() const;
      ConstMatrixView<T> GetV() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected:

      struct SymSVDiv_Impl;
      SymSVDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;


    private :

      inline SymSVDiv(const SymSVDiv<T>&) : pimpl(0)
      { TMVAssert(FALSE); }
      inline SymSVDiv<T>& operator=(const SymSVDiv<T>&)
      { TMVAssert(FALSE); return *this; }

  }; // SymSVDiv

} // namespace tmv;

#endif
