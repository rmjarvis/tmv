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
// det(P) is +-1 (actually the Permutation class does this for us).
// The determinant of U is just the product of the diagonal elements.
// So the determinant of A is just det(P) times the diagonal elements 
// of U.
//
// Finally, since we actually want P.Transpose() in the calculations
// for solving A x = b and xt A = bt, we store the transpose as Pt.
//
//


#ifndef TMV_LUDiv_H
#define TMV_LUDiv_H

#include "TMV_Matrix.h"
#include "TMV_Permutation.h"
#include "TMV_Divider.h"

namespace tmv {

  template <class T> void LU_Decompose(const MatrixView<T>& A, Permutation& P);
  // Decompose A into P * L * U
  
  template <class T1, class T2> void LU_LDivEq(
      const GenMatrix<T1>& LUx, const Permutation& P, 
      const MatrixView<T2>& m);
  template <class T1, class T2> void LU_RDivEq(
      const GenMatrix<T1>& LUx, const Permutation& P, 
      const MatrixView<T2>& m);

#define CT complex<T>
  template <class T> inline void LU_LDivEq(
      const GenMatrix<CT>& LUx, const Permutation& P, 
      const MatrixView<T>& m)
  { TMVAssert(false); }
  template <class T> inline void LU_RDivEq(
      const GenMatrix<CT>& LUx, const Permutation& P, 
      const MatrixView<T>& m)
  { TMVAssert(false); }
#undef CT

  template <class T> class LUDiv : 
    virtual public BaseLUDiv<T> 
  {

    public :

      //
      // Constructors
      //

      LUDiv(const BaseMatrix<T>& A) :
	LUx(A),P(A.colsize()),det(T(0)),donedet(false)
      { 
	TMVAssert(A.IsSquare());
	LU_Decompose(LUx.QuickView(),P); 
      }

      ~LUDiv() {}

      //
      // Divider Versions of DivEq and Div
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      { 
	TMVAssert(LUx.colsize() == m.colsize());
	LU_LDivEq(LUx,P,m); 
      } 

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(LUx.colsize() == m.rowsize());
	LU_RDivEq(LUx,P,m); 
      } 

      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m1,
	    const MatrixView<T2>& m0) const 
	{ 
	  TMVAssert(m1.colsize() == m0.colsize()); 
	  TMVAssert(m1.rowsize() == m0.rowsize());
	  TMVAssert(LUx.colsize() == m1.colsize());
	  LU_LDivEq(LUx,P,m0=m1); 
	}

      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m1,
	    const MatrixView<T2>& m0) const 
	{ 
	  TMVAssert(m1.colsize() == m0.colsize()); 
	  TMVAssert(m1.rowsize() == m0.rowsize());
	  TMVAssert(LUx.colsize() == m1.rowsize());
	  LU_RDivEq(LUx,P,m0=m1); 
	} 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> InverseATA() const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      Matrix<T,ColMajor> LU_GetL() const;
      Matrix<T,ColMajor> LU_GetU() const;
      inline Permutation LU_GetP() const
      { return P; }
      inline Permutation LU_GetQ() const
      { return Permutation(LUx.rowsize()); }
      const Matrix<T,ColMajor>& GetLU() const { return LUx; }
      const Permutation& GetP() const { return P; }

      inline std::string Type() const
      { return std::string("LUDiv<") + tmv::Type(T()) + ">"; }

    private :
      Matrix<T,ColMajor> LUx;
      Permutation P;
      mutable T det;
      mutable bool donedet;
  };

} // namespace tmv

#endif
