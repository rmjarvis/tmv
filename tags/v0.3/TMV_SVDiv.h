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
#undef CT

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>& V);

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det);


  // MJ: Write a version of SVDiv which keeps U,V in compact form
  template <class T> class SVDiv : 
    virtual public BaseSVDiv<T> 
  {

    public :

      //
      // Constructors
      //

      explicit SVDiv(const BaseMatrix<T>& A) :
	istrans(A.colsize() < A.rowsize()), 
	U(max(A.colsize(),A.rowsize()),min(A.colsize(),A.rowsize())), 
	S(U.rowsize()), V(U.rowsize(),U.rowsize()), det(T(1))
      {
	if (istrans) U.QuickTranspose() = A;
	else U = A;

	SV_Decompose(U.QuickView(),S.View(),V.QuickView(),det);

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
	TMVAssert(m.colsize() == U.colsize());
	TMVAssert(U.IsSquare());
	TMVAssert(!istrans);
	SV_LDiv(U,S,V,kmax,m,m); 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.rowsize() == U.rowsize());
	TMVAssert(U.IsSquare());
	TMVAssert(!istrans);
        SV_RDiv(U,S,V,kmax,m,m); 
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
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

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
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
      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0); }
      inline Matrix<T,ColMajor> InverseATA() const
      {
#ifdef TMVDEBUG
	if (!V.IsSquare()) {
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
      { if (istrans) return V.QuickTranspose(); else return U; }
      inline Vector<RealType(T)> SV_GetS() const { return S; }
      inline Matrix<T,ColMajor> SV_GetV() const 
      { if (istrans) return U.QuickTranspose(); else return V; }
      inline bool IsTrans() const { return istrans; }
      inline const Matrix<T,ColMajor>& GetU() const { return U; }
      inline const Vector<RealType(T)>& GetS() const { return S; }
      inline const Matrix<T,ColMajor>& GetV() const { return V; }

      inline std::string Type() const
      { return std::string("SVDiv<") + tmv::Type(T()) + ">"; }

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
