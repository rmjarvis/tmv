//---------------------------------------------------------------------------
#ifndef TMV_BandSVDiv_H
#define TMV_BandSVDiv_H

#include "TMV_BandMatrix.h"
#include "TMV_Divider.h"
#include "TMV_SVDiv.h"

namespace tmv {

  template <class T> void BandSV_Decompose(
      const GenBandMatrix<T>& A,
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det);
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

      explicit BandSVDiv(const GenBandMatrix<T>& A) :
	istrans(A.colsize() < A.rowsize()),
	U(max(A.colsize(),A.rowsize()),min(A.colsize(),A.rowsize())), 
	S(U.rowsize()), V(U.rowsize(),U.rowsize()), det(T(1))
      {
	if (istrans)
	  BandSV_Decompose(A.QuickTranspose(),U.QuickView(),S.View(),
	      V.QuickView(),det);
	else
	  BandSV_Decompose(A,U.QuickView(),S.View(),V.QuickView(),det);
	Thresh(Epsilon<T>());
      }

      ~BandSVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.colsize() == U.colsize());
	TMVAssert(U.IsSquare());
	TMVAssert(!istrans);
	SV_LDiv(U,S,V,kmax,m,m); 
      } 

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.rowsize() == U.colsize());
	TMVAssert(U.IsSquare());
	TMVAssert(!istrans);
	SV_RDiv(U,S,V,kmax,m,m); 
      } 

      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? U.rowsize() : U.colsize());
	TMVAssert(x.colsize() == istrans ? U.colsize() : U.rowsize());
	if (istrans)
	  SV_RDiv(U,S,V,kmax,m.QuickTranspose(),x.QuickTranspose());
	else
	  SV_LDiv(U,S,V,kmax,m,x);
      } 

      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? U.colsize() : U.rowsize());
	TMVAssert(x.rowsize() == istrans ? U.rowsize() : U.colsize());
	if (istrans)
	  SV_LDiv(U,S,V,kmax,m.QuickTranspose(),x.QuickTranspose());
	else
	  SV_RDiv(U,S,V,kmax,m,x);
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
#ifdef TMVDEBUG
	if (!V.IsSquare()) {
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
      { return istrans ? V.QuickTranspose() : U.QuickView(); }
      inline Vector<RealType(T)> SV_GetS() const 
      { return S; }
      inline Matrix<T,ColMajor> SV_GetV() const 
      { return istrans ? U.QuickTranspose() : V.QuickView(); }
      inline const Matrix<T,ColMajor>& GetU() const { return U; }
      inline const Vector<RealType(T)>& GetS() const { return S; }
      inline const Matrix<T,ColMajor>& GetV() const { return V; }
      inline bool IsTrans() const { return istrans; }

      inline std::string Type() const
      { return std::string("BandSVDiv<") + tmv::Type(T()) + ">"; }

    protected:

      bool istrans;
      Matrix<T,ColMajor> U;
      Vector<RealType(T)> S;
      Matrix<T,ColMajor> V;
      T det;
      mutable size_t kmax;
  };

}; // namespace tmv;

#endif
