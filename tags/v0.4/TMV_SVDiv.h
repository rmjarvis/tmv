//---------------------------------------------------------------------------
#ifndef TMV_SVDiv_H
#define TMV_SVDiv_H

#include "TMV_Matrix.h"
#include "TMV_Vector.h"
#include "TMV_Divider.h"

namespace tmv {

  template <class T1, class T2, class T3> void SV_LDiv(
	const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
	const GenMatrix<T1>& V, size_t kmax,
	const GenMatrix<T2>& m, const MatrixView<T3>& x);
  template <class T1, class T2, class T3> void SV_RDiv(
	const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
	const GenMatrix<T1>& V, size_t kmax,
	const GenMatrix<T2>& m, const MatrixView<T3>& x);

  template <class T1, class T2, class T3> void SV_LDiv(
	const GenMatrix<T1>& UV,
	const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta,
	const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
	const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
	const GenVector<RealType(T1)>& S, size_t kmax,
	const GenMatrix<T2>& m, const MatrixView<T3>& x);
  template <class T1, class T2, class T3> void SV_RDiv(
	const GenMatrix<T1>& UV,
	const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta,
	const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
	const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
	const GenVector<RealType(T1)>& S, size_t kmax,
	const GenMatrix<T2>& m, const MatrixView<T3>& x);

#define CT complex<T>
  template <class T> void SV_LDiv(
      const GenMatrix<CT>& U, const GenVector<T>& S, 
      const GenMatrix<CT>& V, size_t kmax,
      const GenMatrix<T>& m, const MatrixView<T>& x)
  { TMVAssert(false); }
  template <class T> void SV_RDiv(
      const GenMatrix<CT>& U, const GenVector<T>& S, 
      const GenMatrix<CT>& V, size_t kmax,
      const GenMatrix<T>& m, const MatrixView<T>& x)
  { TMVAssert(false); }
  template <class T> void SV_LDiv(
	const GenMatrix<CT>& UV,
	const GenVector<CT>& Ubeta, const GenVector<CT>& Vbeta,
	const GenMatrix<CT>* UV2, const GenVector<CT>* Qbeta,
	const GenMatrix<T>& U1, const GenMatrix<T>& V1,
	const GenVector<T>& S, size_t kmax,
	const GenMatrix<T>& m, const MatrixView<T>& x)
  { TMVAssert(false); }
  template <class T> void SV_RDiv(
	const GenMatrix<CT>& UV,
	const GenVector<CT>& Ubeta, const GenVector<CT>& Vbeta,
	const GenMatrix<CT>* UV2, const GenVector<CT>* Qbeta,
	const GenMatrix<T>& U1, const GenMatrix<T>& V1,
	const GenVector<T>& S, size_t kmax,
	const GenMatrix<T>& m, const MatrixView<T>& x)
  { TMVAssert(false); }
#undef CT

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV=false);
  // If SetUV, then (U,V must be !=0) U,V are Set to the appropriate matrices.
  // If !SetUV (default), then (if U,V!=0) U,V are multiplied to preserve
  // constancy of U B V.

  template <class T> void SV_Decompose(
      const MatrixView<T>& UV, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& S, T& det);
  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det, bool StoreU=true);
  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      T& det, bool StoreU=true);
  template <class T> void SV_Decompose(
      const MatrixView<T>& UV, 
      const VectorView<T>& Ubeta, const VectorView<T>& Vbeta,
      Matrix<T,ColMajor>*& UV2, Vector<T>*& Qbeta,
      const MatrixView<RealType(T)>& U1, const MatrixView<RealType(T)>& V1,
      const VectorView<RealType(T)>& S, T& det);

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

  template <class T> class SVFDiv : 
    virtual public BaseSVDiv<T> 
  {

    public :

      //
      // Constructors
      //

      SVFDiv(const BaseMatrix<T>& A, bool StoreU, bool StoreV) :
	istrans(A.colsize() < A.rowsize()), 
	U(0), S(min(A.colsize(),A.rowsize())), V(0), det(T(1))
      {
	size_t M = max(A.colsize(),A.rowsize());
	size_t N = min(A.colsize(),A.rowsize());

	U = new Matrix<T,ColMajor>(M,N);
	if (istrans) U->QuickTranspose() = A;
	else *U = A;

	if (istrans) swap(StoreU,StoreV);

	if (StoreV) {
	  V = new Matrix<T,ColMajor>(N,N);
	  SV_Decompose(U->QuickView(),S.View(),V->QuickView(),det,StoreU);
	}
	else SV_Decompose(U->QuickView(),S.View(),det,StoreU);

	if (!StoreU) { delete U; U = 0; }

	// Set kmax for actual 0 elements (to within machine precision).
	// Any further cut in the number of singular values to use
	// should be done by the user.
	Thresh(Epsilon<T>());
      }

      ~SVFDiv() { if (U) delete U; if (V) delete V; }

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U && V);
	TMVAssert(m.colsize() == U->colsize());
	TMVAssert(U->IsSquare());
	TMVAssert(!istrans);
	SV_LDiv(*U,S,*V,kmax,m,m); 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(U && V);
	TMVAssert(m.rowsize() == U->rowsize());
	TMVAssert(U->IsSquare());
	TMVAssert(!istrans);
        SV_RDiv(*U,S,*V,kmax,m,m); 
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
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

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
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
      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0); }
      inline Matrix<T,ColMajor> InverseATA() const
      {
	TMVAssert(V);
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning: InverseATA called for short matrix in SVFDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	return DoInverseATA();
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
      inline bool IsTrans() const { return istrans; }
      inline const Matrix<T,ColMajor>* GetU() const { return U; }
      inline const Vector<RealType(T)>& GetS() const { return S; }
      inline const Matrix<T,ColMajor>* GetV() const { return V; }

      inline std::string Type() const
      { return std::string("SVFDiv<") + tmv::Type(T()) + ">"; }

    protected:

      bool istrans;
      Matrix<T,ColMajor>* U;
      Vector<RealType(T)> S;
      Matrix<T,ColMajor>* V;
      T det;
      mutable size_t kmax;
  };

  template <class T> class SVDiv : 
    virtual public BaseSVDiv<T> 
  {

    public :

      //
      // Constructors
      //

      // Normally A = U U1 S V1 V
      // where U,V are stored in UV, Ubeta, Vbeta
      //
      // If UV2 and Qbeta are set (ie. pointers are not 0), then
      // A = Q U U1 S V1 V
      // where U,V are stored in *UV2, Ubeta, Vbeta
      // and Q is stored in UV, *Qbeta
      //
      // If istrans = true, then all of this is the storage for A.Transpose().
      //
      explicit SVDiv(const BaseMatrix<T>& A) :
	istrans(A.colsize() < A.rowsize()), 
	UV(max(A.colsize(),A.rowsize()),min(A.colsize(),A.rowsize())), 
	Ubeta(UV.rowsize()), Vbeta(UV.rowsize()-1),
	UV2(0), Qbeta(0), U1(UV.rowsize(),UV.rowsize()), 
	V1(UV.rowsize(),UV.rowsize()), S(UV.rowsize()), det(T(1))
      {
	if (istrans) UV.QuickTranspose() = A;
	else UV = A;

	SV_Decompose(UV.QuickView(),Ubeta.View(),Vbeta.View(),UV2,Qbeta,
	    U1.QuickView(),V1.QuickView(),S.View(),det);

	// Set kmax for actual 0 elements (to within machine precision).
	// Any further cut in the number of singular values to use
	// should be done by the user.
	Thresh(Epsilon<T>());
      }

      ~SVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.colsize() == UV.colsize());
	TMVAssert(UV.IsSquare());
	TMVAssert(!istrans);
	SV_LDiv(UV,Ubeta,Vbeta,UV2,Qbeta,U1,V1,S,kmax,m,m); 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.rowsize() == UV.rowsize());
	TMVAssert(UV.IsSquare());
	TMVAssert(!istrans);
        SV_RDiv(UV,Ubeta,Vbeta,UV2,Qbeta,U1,V1,S,kmax,m,m); 
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? UV.rowsize() : UV.colsize());
	TMVAssert(x.colsize() == istrans ? UV.colsize() : UV.rowsize());
	if (istrans)
	  SV_RDiv(UV,Ubeta,Vbeta,UV2,Qbeta,U1,V1,S,kmax,
	      m.QuickTranspose(),x.QuickTranspose()); 
	else
	  SV_LDiv(UV,Ubeta,Vbeta,UV2,Qbeta,U1,V1,S,kmax,m,x); 
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? UV.colsize() : UV.rowsize());
	TMVAssert(x.rowsize() == istrans ? UV.rowsize() : UV.colsize());
	if (istrans)
	  SV_LDiv(UV,Ubeta,Vbeta,UV2,Qbeta,U1,V1,S,kmax,
	      m.QuickTranspose(),x.QuickTranspose()); 
	else
	  SV_RDiv(UV,Ubeta,Vbeta,UV2,Qbeta,U1,V1,S,kmax,m,x); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const { return det; }
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> DoInverseATA() const;
      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0); }
      inline Matrix<T,ColMajor> InverseATA() const
      {
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning: InverseATA called for short matrix in SVDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	return DoInverseATA();
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

      inline Matrix<T,ColMajor> SV_GetU() const 
      { if (istrans) return GetV().QuickTranspose(); else return GetU(); }
      inline Vector<RealType(T)> SV_GetS() const { return S; }
      inline Matrix<T,ColMajor> SV_GetV() const 
      { if (istrans) return GetU().QuickTranspose(); else return GetV(); }

      Matrix<T,ColMajor> GetU() const;
      Matrix<T,ColMajor> GetV() const;

      inline const Vector<RealType(T)>& GetS() const { return S; }
      inline bool IsTrans() const { return istrans; }
      inline const Matrix<T,ColMajor>& GetUV() const { return UV; }
      inline const Vector<T>& GetUbeta() const { return Ubeta; }
      inline const Vector<T>& GetVbeta() const { return Vbeta; }
      inline const Matrix<T,ColMajor>* GetUV2() const { return UV2; }
      inline const Vector<T>* GetQbeta() const { return Qbeta; }
      inline const Matrix<RealType(T),ColMajor>& GetU1() const { return U1; }
      inline const Matrix<RealType(T),ColMajor>& GetV1() const { return V1; }

      inline std::string Type() const
      { return std::string("SVDiv<") + tmv::Type(T()) + ">"; }

    protected:

      bool istrans;
      Matrix<T,ColMajor> UV;
      Vector<T> Ubeta;
      Vector<T> Vbeta;
      Matrix<T,ColMajor>* UV2;
      Vector<T>* Qbeta;
      Matrix<RealType(T),ColMajor> U1;
      Matrix<RealType(T),ColMajor> V1;
      Vector<RealType(T)> S;
      T det;
      mutable size_t kmax;
  };

}; // namespace tmv;

#endif
