//---------------------------------------------------------------------------
#ifndef TMV_BandSVDiv_H
#define TMV_BandSVDiv_H

#include "TMV_BandMatrix.h"
#include "TMV_Divider.h"
#include "TMV_SVDiv.h"

namespace tmv {

  template <class T> void BandSV_Decompose(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det);
  // Decompose A into U S V
  // where S is a diagonal real matrix, and U,V are unitary matrices.
  // U,S,V are N x N

  template <class T> class BandSVDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      BandSVDiv(const GenBandMatrix<T>& A, bool StoreU, bool StoreV);
      ~BandSVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.colsize() == U->colsize());
	TMVAssert(U->IsSquare());
	TMVAssert(!istrans);
	if (istrans) SV_RDiv(*U,S,*V,kmax,m.Transpose(),m.Transpose());
	else SV_LDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.rowsize() == U->colsize());
	TMVAssert(U->IsSquare());
	if (istrans) SV_LDiv(*U,S,*V,kmax,m.Transpose(),m.Transpose());
	else SV_RDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(U.get() && V.get());
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? U->rowsize() : U->colsize());
	TMVAssert(x.colsize() == istrans ? U->colsize() : U->rowsize());
	if (istrans) SV_RDiv(*U,S,*V,kmax,m.Transpose(),x.Transpose());
	else SV_LDiv(*U,S,*V,kmax,m,x);
      } 

      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m,
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
      inline void InverseATA(const MatrixView<T>& minv) const
      {
	TMVAssert(V.get());
#ifdef TMVDEBUG
	if (!V->IsSquare()) {
	  cout<<"Warning: InverseATA called for short matrix in BandSVDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	DoInverseATA(minv);
      }

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

      inline std::string Type() const
      { return std::string("BandSVDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return SV; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected:

      bool istrans;
      auto_ptr<Matrix<T,ColMajor> > U;
      Vector<RealType(T)> S;
      auto_ptr<Matrix<T,ColMajor> > V;
      T det;
      mutable size_t kmax;
  }; // BandSVDiv

} // namespace tmv;

#endif
