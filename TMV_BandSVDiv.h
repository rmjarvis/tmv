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
    virtual public BaseSVDiv<T> 
  {

    public :

      //
      // Constructors
      //

      explicit BandSVDiv(const GenBandMatrix<T>& A, bool StoreU,
	  bool StoreV) :
	istrans(A.colsize() < A.rowsize()),
	U(0), S(min(A.rowsize(),A.colsize())), V(0), det(T(1))
      {
	size_t M = max(A.colsize(),A.rowsize());
	size_t N = min(A.colsize(),A.rowsize());

	if (istrans) swap(StoreU,StoreV);

	if (StoreU) U = new Matrix<T,ColMajor>(M,N);
        if (StoreV) V = new Matrix<T,ColMajor>(N,N);
	MatrixView<T>* Uv = U ? new MatrixView<T>(U->View()) : 0;
	MatrixView<T>* Vv = V ? new MatrixView<T>(V->View()) : 0;

	if (istrans)
	  BandSV_Decompose(A.QuickTranspose(),Uv,S.View(),Vv,det);
	else
	  BandSV_Decompose(A,Uv,S.View(),Vv,det);

	delete Uv;
	delete Vv;

	Thresh(Epsilon<T>());
      }

      ~BandSVDiv() { if (U) delete U; if (V) delete V; }

      //
      // Div, DivEq
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U && V);
	TMVAssert(m.colsize() == U->colsize());
	TMVAssert(U->IsSquare());
	TMVAssert(!istrans);
	SV_LDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U && V);
	TMVAssert(m.rowsize() == U->colsize());
	TMVAssert(U->IsSquare());
	TMVAssert(!istrans);
	SV_RDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(U && V);
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? U->rowsize() : U->colsize());
	TMVAssert(x.colsize() == istrans ? U->colsize() : U->rowsize());
	if (istrans)
	  SV_RDiv(*U,S,*V,kmax,m.QuickTranspose(),x.QuickTranspose());
	else
	  SV_LDiv(*U,S,*V,kmax,m,x);
      } 

      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(U && V);
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? U->colsize() : U->rowsize());
	TMVAssert(x.rowsize() == istrans ? U->rowsize() : U->colsize());
	if (istrans)
	  SV_LDiv(*U,S,*V,kmax,m.QuickTranspose(),x.QuickTranspose());
	else
	  SV_RDiv(*U,S,*V,kmax,m,x);
      } 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const { return det; }
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> DoInverseATA() const;
      inline Matrix<T,ColMajor> InverseATA() const
      {
	TMVAssert(V);
#ifdef TMVDEBUG
	if (!V->IsSquare()) {
	  cout<<"Warning: InverseATA called for short matrix in BandSVDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	return DoInverseATA();
      }

      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0); }

      //
      // Determine which (if any) S values to zero out
      //

      void Thresh(RealType(T) toler,ostream* debugout=0) const;
      void Top(size_t neigen,ostream* debugout=0) const;
      inline size_t GetKMax() const { return kmax; }

      //
      // Access Decomposition
      //

      inline Matrix<T,ColMajor> SV_GetU() const 
      {
	if (istrans) { TMVAssert(V); return V->QuickTranspose(); }
	else { TMVAssert(U); return *U; }
      }
      inline Vector<RealType(T)> SV_GetS() const { return S; }
      inline Matrix<T,ColMajor> SV_GetV() const 
      {
	if (istrans) { TMVAssert(U); return U->QuickTranspose(); }
	else { TMVAssert(V); return *V; }
      }
      inline const Matrix<T,ColMajor>* GetU() const { return U; }
      inline const Vector<RealType(T)>& GetS() const { return S; }
      inline const Matrix<T,ColMajor>* GetV() const { return V; }
      inline bool IsTrans() const { return istrans; }

      inline std::string Type() const
      { return std::string("BandSVDiv<") + tmv::Type(T()) + ">"; }

    protected:

      bool istrans;
      Matrix<T,ColMajor>* U;
      Vector<RealType(T)> S;
      Matrix<T,ColMajor>* V;
      T det;
      mutable size_t kmax;
  }; // BandSVDiv

}; // namespace tmv;

#endif
