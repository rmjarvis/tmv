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


#ifndef TMV_LUDiv_H
#define TMV_LUDiv_H

#include "TMV_Matrix.h"
#include "TMV_Divider.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T> void LU_Decompose(const MatrixView<T>& A, size_t* P, T& det);
  // Decompose A into P * L * U
  
  template <class T, class T1> void LU_LDivEq(
      const GenMatrix<T1>& LUx, const size_t* P, const MatrixView<T>& m);
  template <class T, class T1> void LU_RDivEq(
      const GenMatrix<T1>& LUx, const size_t* P, const MatrixView<T>& m);

#define CT complex<T>
  template <class T> inline void LU_LDivEq(
      const GenMatrix<CT>& , const size_t* , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void LU_RDivEq(
      const GenMatrix<CT>& , const size_t* , const MatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> class LUDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      LUDiv(const GenMatrix<T>& A, bool _inplace);
      ~LUDiv() {}

      //
      // Divider Versions of DivEq and Div
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      { 
	TMVAssert(LUx.colsize() == m.colsize());
	if (istrans) LU_RDivEq(LUx,P.get(),m.Transpose()); 
	else LU_LDivEq(LUx,P.get(),m); 
      } 

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(LUx.colsize() == m.rowsize());
	if (istrans) LU_LDivEq(LUx,P.get(),m.Transpose()); 
	else LU_RDivEq(LUx,P.get(),m); 
      } 

      template <class T1, class T2> inline void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const 
      { 
	TMVAssert(m1.colsize() == m0.colsize()); 
	TMVAssert(m1.rowsize() == m0.rowsize());
	TMVAssert(LUx.colsize() == m1.colsize());
	DoLDivEq(m0=m1); 
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const 
      { 
	TMVAssert(m1.colsize() == m0.colsize()); 
	TMVAssert(m1.rowsize() == m0.rowsize());
	TMVAssert(LUx.colsize() == m1.rowsize());
	DoRDivEq(m0=m1); 
      } 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const MatrixView<T>& minv) const;
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      inline bool IsTrans() const { return istrans; }
      inline ConstLowerTriMatrixView<T> GetL() const
      { return LowerTriMatrixViewOf(LUx,UnitDiag); }
      inline ConstUpperTriMatrixView<T> GetU() const
      { return UpperTriMatrixViewOf(LUx,NonUnitDiag); }
      const GenMatrix<T>& GetLU() const { return LUx; }
      const size_t* GetP() const { return P.get(); }

      inline string Type() const
      { return string("LUDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return LU; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    private :
      const bool istrans;
      const bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      MatrixView<T> LUx;
      auto_array<size_t> P;
      mutable T det;
      mutable bool donedet;
  };

} // namespace tmv

#endif
