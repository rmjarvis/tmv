//---------------------------------------------------------------------------
//
// This file defines the Divider helper for the DiagMatrix class.
//
#ifndef TMV_DiagSVDiv_H
#define TMV_DiagSVDiv_H

#include "TMV_DiagMatrix.h"
#include "TMV_Divider.h"
#include "TMV_DiagLUDiv.h"
#include "TMV_Permutation.h"

namespace tmv {

  template <class T> void DiagSV_Decompose(
      const GenDiagMatrix<T>& m, Permutation& P,
      const DiagMatrixView<RealType(T)>& S,
      const DiagMatrixView<T>& V1, T& det);

  template <class T1, class T2> void SV_LDivEq(
      const Permutation& P, const GenDiagMatrix<RealType(T1)>& S,
      const GenDiagMatrix<T1>& V1, const size_t kmax,
      const MatrixView<T2>& m);

#define CT complex<T>
  template <class T> void SV_LDivEq(
      const Permutation& P, const GenDiagMatrix<T>& S, 
      const GenDiagMatrix<CT>& V1,
      const size_t kmax, const MatrixView<T>& m)
  { TMVAssert(false); }
#undef CT

  template <class T> class DiagSVDiv : 
    virtual public BaseSVDiv<T>, 
    virtual public DiagDivider<T>
  {

    // A = U S V
    // For a diagonal matrix, 
    // U is a permutation = P.Transpose()
    // S is diagonal real
    // V = V1 * P, where V1 are the exp(it) components of A
    //
    
    public:

      DiagSVDiv(const GenDiagMatrix<T>& m) :
	P(m.size()), S(m.size()), V1(m.size()),
	det(T(1)), calcdet(false), kmax(m.size())
      { 
	DiagSV_Decompose(m,P,S.View(),V1.QuickView(),det); 
	Thresh(Epsilon<T>());
      }

      ~DiagSVDiv() {}

      T Det() const;
      DiagMatrix<T> DInverse() const;
      DiagMatrix<T> DInverseATA() const;

      Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(DInverse()); }
      Matrix<T,ColMajor> InverseATA() const 
      { return Matrix<T,ColMajor>(DInverseATA()); }

      inline bool Singular() const { return kmax < S.size(); }
      inline RealType(T) Norm2() const { return S(0,0); }

      template <class T2> inline void DoLDivEq(const MatrixView<T2>& m) const 
      {
	TMVAssert(m.colsize() == S.size());
	SV_LDivEq(P,S,V1,kmax,m);
      }

      template <class T2> inline void DoRDivEq(const MatrixView<T2>& m) const 
      {
	TMVAssert(m.rowsize() == S.size());
	SV_LDivEq(P,S,V1,kmax,m.QuickTranspose());
      }

      template <class T2, class T3> inline void DoLDiv(
	  const GenMatrix<T2>& m1, const MatrixView<T3>& m0) const 
      { 
	TMVAssert(m1.colsize() == m0.colsize()); 
	TMVAssert(m1.rowsize() == m0.rowsize()); 
	TMVAssert(S.size() == m0.colsize()); 
	SV_LDivEq(P,S,V1,kmax,m0=m1); 
      } 

      template <class T2, class T3> inline void DoRDiv(
	  const GenMatrix<T2>& m1, const MatrixView<T3>& m0) const 
      { 
	TMVAssert(m1.colsize() == m0.colsize()); 
	TMVAssert(m1.rowsize() == m0.rowsize()); 
	TMVAssert(S.size() == m0.rowsize()); 
	SV_LDivEq(P,S,V1,kmax,(m0=m1).QuickTranspose()); 
      } 

#include "TMV_AuxAllDiv.h"

      void Thresh(RealType(T) toler, ostream* debugout=0) const;
      void Top(size_t neigen, ostream* debugout=0) const;
      inline size_t GetKMax() const { return kmax; }

      inline Matrix<T,ColMajor> SV_GetU() const { return P.Transpose(); }
      inline Vector<RealType(T)> SV_GetS() const { return S.diag(); }
      inline Matrix<T,ColMajor> SV_GetV() const 
      { return P * Matrix<T,ColMajor>(V1); }

      inline const Permutation& GetP() const { return P; }
      inline const DiagMatrix<RealType(T)>& GetS() const { return S; }
      inline const DiagMatrix<T>& GetV1() const { return V1; }

      inline std::string Type() const
      { return std::string("DiagSVDiv<") + tmv::Type(T()) + ">"; }

    protected:

      Permutation P;
      DiagMatrix<RealType(T)> S;
      DiagMatrix<T> V1;
      mutable T det;
      mutable bool calcdet;
      mutable size_t kmax;

  };

} // namespace tmv

#endif
