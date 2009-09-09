//---------------------------------------------------------------------------
//
// This file defines the Divider helper for the DiagMatrix class.
//
#ifndef TMV_DiagDiv_H
#define TMV_DiagDiv_H

#include "TMV_DiagMatrix.h"
#include "TMV_Divider.h"

namespace tmv {

  template <class T1, class T2> void DiagLDivEq(
      const GenDiagMatrix<T1>& d, const MatrixView<T2>& m);

#define CT complex<T>
  template <class T> inline void DiagLDivEq(
      const GenDiagMatrix<CT>& d, const MatrixView<T>& m)
  { TMVAssert(false); }
#undef CT

  template <class T> class DiagDiv : 
    virtual public Divider<T>
  {

    public:

      DiagDiv(const GenDiagMatrix<T>& m) :
	itsm(&m), det(T(1)), calcdet(false) {}
      ~DiagDiv() {}

      template <class T2> inline void DoLDivEq(const MatrixView<T2>& m) const 
      { 
	TMVAssert(m.colsize() == itsm->size());
	DiagLDivEq(*itsm,m); 
      } 

      template <class T2> inline void DoRDivEq(const MatrixView<T2>& m) const 
      { 
	TMVAssert(m.rowsize() == itsm->size());
	DiagLDivEq(*itsm,m.Transpose()); 
      } 

      template <class T2, class T3> inline void DoLDiv(
	  const GenMatrix<T2>& m1, const MatrixView<T3>& m0) const 
      {  
	TMVAssert(m1.rowsize() == m0.rowsize()); 
	TMVAssert(m1.colsize() == m0.colsize());
	TMVAssert(m0.colsize() == itsm->size()); 
	DiagLDivEq(*itsm,m0=m1); 
      } 

      template <class T2, class T3> inline void DoRDiv(
	  const GenMatrix<T2>& m1, const MatrixView<T3>& m0) const 
      {  
	TMVAssert(m1.rowsize() == m0.rowsize()); 
	TMVAssert(m1.colsize() == m0.colsize());
	TMVAssert(m0.rowsize() == itsm->size());
	DiagLDivEq(*itsm,(m0=m1).Transpose()); 
      } 

#include "TMV_AuxAllDiv.h"

#define RT RealType(T)
#define CT ComplexType(T)
#define DefDivEq(T1) \
      inline void DivEq(const DiagMatrixView<T1>& m) const\
      { DoLDivEq(ColVectorViewOf(m.diag())); }

      DefDivEq(RT);
      DefDivEq(CT);

#undef DefDivEq

#define DefDiv(T1,T2) \
      inline void Div(const GenDiagMatrix<T1>& b, \
	  const DiagMatrixView<T2>& x) const \
      { DoLDiv(ColVectorViewOf(b.diag()),ColVectorViewOf(x.diag())); }

      DefDiv(RT,T);
      DefDiv(CT,CT);

#undef DefDiv
#undef RT
#undef CT

      T Det() const;
      void DInverse(const DiagMatrixView<T>& minv) const;
      void DInverseATA(const DiagMatrixView<T>& minv) const;
      inline void Inverse(const MatrixView<T>& minv) const 
      { 
	minv.Zero();
	DInverse(DiagMatrixViewOf(minv.diag()));
      }
      inline void InverseATA(const MatrixView<T>& minv) const 
      { 
	minv.Zero();
	DInverseATA(DiagMatrixViewOf(minv.diag()));
      }
      inline bool Singular() const { return Det() == T(0); }

      inline std::string Type() const
      { return std::string("DiagDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return LU; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected:

      const GenDiagMatrix<T>*const itsm;
      mutable T det;
      mutable bool calcdet;
  };

} // namespace tmv

#endif
