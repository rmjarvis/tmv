//---------------------------------------------------------------------------
#ifndef TMV_SVDiv_H
#define TMV_SVDiv_H

#include "TMV_Matrix.h"
#include "TMV_Vector.h"
#include "TMV_Divider.h"

namespace tmv {

  template <class T, class T1, class T2> void SV_LDiv(
	const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
	const GenMatrix<T1>& V, size_t kmax,
	const GenMatrix<T2>& m, const MatrixView<T>& x);
  template <class T, class T1, class T2> void SV_RDiv(
	const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
	const GenMatrix<T1>& V, size_t kmax,
	const GenMatrix<T2>& m, const MatrixView<T>& x);

#define CT complex<T>
  template <class T> inline void SV_LDiv(
      const GenMatrix<CT>& , const GenVector<T>& , 
      const GenMatrix<CT>& , size_t ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SV_RDiv(
      const GenMatrix<CT>& , const GenVector<T>& , 
      const GenMatrix<CT>& , size_t ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV=false);
  // If SetUV, then (U,V must be !=0) U,V are Set to the appropriate matrices.
  // If !SetUV (default), then (if U,V!=0) U,V are multiplied to preserve
  // constancy of U B V.

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det, bool StoreU=true);
  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      T& det, bool StoreU=true);

  // Without det:
  template <class T> inline void SV_Decompose(
      const MatrixView<T>& UV, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& S)
  { const T det=0; SV_Decompose(UV,S,det); }
  template <class T> inline void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, bool StoreU = true)
  { const T det=0; SV_Decompose(U,S,V,det,StoreU); }
  template <class T> inline void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      bool StoreU = true)
  { const T det=0; SV_Decompose(U,S,det,StoreU); }

  template <class T> class SVDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      SVDiv(const GenMatrix<T>& A, bool _inplace, bool StoreU, bool StoreV);
      ~SVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.colsize() == U->colsize());
	TMVAssert(U->IsSquare());
	if (istrans) SV_RDiv(*U,S,*V,kmax,m.Transpose(),m.Transpose()); 
	else SV_LDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.rowsize() == U->rowsize());
	TMVAssert(U->IsSquare());
	if (istrans) SV_LDiv(*U,S,*V,kmax,m.Transpose(),m.Transpose()); 
	else SV_RDiv(*U,S,*V,kmax,m,m); 
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? U->rowsize() : U->colsize());
	TMVAssert(x.colsize() == istrans ? U->colsize() : U->rowsize());
	if (istrans) SV_RDiv(*U,S,*V,kmax,m.Transpose(),x.Transpose()); 
	else SV_LDiv(*U,S,*V,kmax,m,x); 
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? U->colsize() : U->rowsize());
	TMVAssert(x.rowsize() == istrans ? U->rowsize() : U->colsize());
	if (istrans) SV_LDiv(*U,S,*V,kmax,m.Transpose(),x.Transpose()); 
	else SV_RDiv(*U,S,*V,kmax,m,x); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const { return det; }
      void Inverse(const MatrixView<T>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0); }
      inline RealType(T) Condition() const { return S(0)/S(S.size()-1); }
      inline void InverseATA(const MatrixView<T>& minv) const
      {
	TMVAssert(V.get());
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning: InverseATA called for short matrix in SVDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	DoInverseATA(minv);
      }

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
      { 
	if (istrans) { TMVAssert(V.get()); return V->Transpose(); }
	else { TMVAssert(U.get()); return U->View(); }
      }
      inline ConstVectorView<RealType(T)> GetS() const 
      { return S.View(); }
      inline ConstMatrixView<T> GetV() const 
      { 
	if (istrans) { TMVAssert(U.get()); return U->Transpose(); }
	else { TMVAssert(V.get()); return V->View(); }
      }

      inline string Type() const
      { return string("SVDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const 
      { return U.get() ? V.get() ? SV : SVU : V.get() ? SVV : SVS; }
      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected:

      bool istrans;
      bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      auto_ptr<MatrixView<T> > U;
      Vector<RealType(T)> S;
      auto_ptr<Matrix<T,ColMajor> > V;
      T det;
      mutable size_t kmax;
  };

} // namespace tmv;

#endif
