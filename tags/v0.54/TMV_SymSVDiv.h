//---------------------------------------------------------------------------
#ifndef TMV_SymSVDiv_H
#define TMV_SymSVDiv_H

#include "TMV_SymMatrix.h"
#include "TMV_SymDivider.h"

namespace tmv {

  template <class T> void HermSV_Decompose_From_Tridiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, bool SetU);
  // Reduce the offdiagonal elements, E(i), to 0.
  // U is multiplied to preserve the constancy of U T Ut (or U T UT)

  template <class T> void HermSV_Decompose(const MatrixView<T>& U,
      const VectorView<RealType(T)>& S, bool StoreU=true);
  // The Hermitian Matrix A is input as the LowerTri portion of U.
  // If StoreU = true, U is output as a unitary matrix with A = U S Ut.

  template <class T> void SymSV_Decompose(const MatrixView<T>& U,
      const VectorView<RealType(T)>& S, const MatrixView<T>* V, T& det);
  // The Symmetric Matrix A is input as the LowerTri portion of U.
  // If V != 0, U,V are output as a unitary matrices with A = U S V.

  template <class T> class HermSVDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      HermSVDiv(const GenSymMatrix<T>& A, bool StoreU);
      ~HermSVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get());
	TMVAssert(m.colsize() == U->colsize());
	SV_LDiv(*U,S,U->Adjoint(),kmax,m,m); 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get());
	TMVAssert(m.rowsize() == U->rowsize());
	SV_RDiv(*U,S,U->Adjoint(),kmax,m,m); 
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get());
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == U->colsize());
	TMVAssert(x.colsize() == U->rowsize());
	SV_LDiv(*U,S,U->Adjoint(),kmax,m,x); 
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get());
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == U->rowsize());
	TMVAssert(x.rowsize() == U->colsize());
	SV_RDiv(*U,S,U->Adjoint(),kmax,m,x); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const SymMatrixView<T>& sinv) const;
      void Inverse(const MatrixView<T>& minv) const;
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return abs(S(0)); }
      inline RealType(T) Condition() const { return abs(S(0)/S(S.size()-1)); }

      //
      // Determine which (if any) S values to zero out
      //

      void Thresh(RealType(T) toler,ostream* debugout=0) const;
      void Top(size_t neigen,ostream* debugout=0) const;
      inline size_t GetKMax() const { return kmax; }

      //
      // Access Decomposition
      //

      inline ConstMatrixView<T> GetU() const 
      { TMVAssert(U.get()); return U->View(); }
      inline ConstVectorView<RealType(T)> GetS() const 
      { return S.View(); }
      inline ConstMatrixView<T> GetV() const 
      { TMVAssert(U.get()); return U->Adjoint(); }

      inline string Type() const
      { return string("HermSVDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return SV; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected:

      auto_ptr<Matrix<T,ColMajor> > U;
      Vector<RealType(T)> S;
      mutable RealType(T) det;
      mutable bool calcdet;
      mutable size_t kmax;
  }; // HermSVDiv

  template <class T> class SymSVDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      SymSVDiv(const GenSymMatrix<T>& A, bool StoreU, bool StoreV);
      ~SymSVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.colsize() == U->colsize());
	SV_LDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.rowsize() == U->rowsize());
	SV_RDiv(*U,S,*V,kmax,m,m); 
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == U->colsize());
	TMVAssert(x.colsize() == U->rowsize());
	SV_LDiv(*U,S,*V,kmax,m,x); 
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == U->rowsize());
	TMVAssert(x.rowsize() == U->colsize());
	SV_RDiv(*U,S,*V,kmax,m,x); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const { return det; }
      void Inverse(const SymMatrixView<T>& sinv) const;
      void Inverse(const MatrixView<T>& minv) const;
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0); }
      inline RealType(T) Condition() const { return S(0)/S(S.size()-1); }

      //
      // Determine which (if any) S values to zero out
      //

      void Thresh(RealType(T) toler,ostream* debugout=0) const;
      void Top(size_t neigen,ostream* debugout=0) const;
      inline size_t GetKMax() const { return kmax; }

      //
      // Access Decomposition
      //

      inline ConstMatrixView<T> GetU() const 
      { TMVAssert(U.get()); return U->View(); }
      inline ConstVectorView<RealType(T)> GetS() const 
      { return S.View(); }
      inline ConstMatrixView<T> GetV() const 
      { TMVAssert(V.get()); return V->View(); }

      inline string Type() const
      { return string("SymSVDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return SV; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected:

      auto_ptr<Matrix<T,ColMajor> > U;
      Vector<RealType(T)> S;
      auto_ptr<Matrix<T,ColMajor> > V;
      T det;
      mutable size_t kmax;
  }; // SymSVDiv

} // namespace tmv;

#endif
