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

#include "TMV_SymMatrix.h"
#include "TMV_SymDivider.h"

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

#define CT complex<T>
  template <class T> inline void SymLU_LDivEq(
      const GenSymMatrix<CT>& , const GenVector<CT>& , const size_t* , 
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SymLU_RDivEq(
      const GenSymMatrix<CT>& , const GenVector<CT>& , const size_t* , 
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> class SymLUDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      SymLUDiv(const GenSymMatrix<T>& A, bool _inplace);
      ~SymLUDiv() {}

      //
      // Divider Versions of DivEq and Div
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      { 
	TMVAssert(LLx.size() == m.colsize());
	SymLU_LDivEq(LLx,xD,P.get(),m); 
      } 

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(LLx.size() == m.rowsize());
	SymLU_RDivEq(LLx,xD,P.get(),m); 
      } 

      template <class T1, class T2> inline void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const 
      { 
	TMVAssert(m1.colsize() == m0.colsize()); 
	TMVAssert(m1.rowsize() == m0.rowsize());
	TMVAssert(LLx.size() == m1.colsize());
	SymLU_LDivEq(LLx,xD,P.get(),m0=m1); 
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const 
      { 
	TMVAssert(m1.colsize() == m0.colsize()); 
	TMVAssert(m1.rowsize() == m0.rowsize());
	TMVAssert(LLx.size() == m1.rowsize());
	SymLU_RDivEq(LLx,xD,P.get(),m0=m1); 
      } 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const { return det; }
      void Inverse(const SymMatrixView<T>& sinv) const;
      void Inverse(const MatrixView<T>& minv) const;
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      inline const ConstLowerTriMatrixView<T> GetL() const 
      { return LLx.LowerTri().MakeUnitDiag(); }
      const Matrix<T> GetD() const;
      inline const size_t* GetP() const { return P.get(); }
      inline const GenSymMatrix<T>& GetLL() const { return LLx; }
      inline const GenVector<T>& GetxD() const { return xD; }

      inline string Type() const
      { return string("SymLUDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return LU; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    private :
      bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      SymMatrixView<T> LLx;
      Vector<T> xD;
      auto_array<size_t> P;
      mutable T det;
  };

} // namespace tmv

#endif
