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
    virtual public BaseLUDiv<T>, 
    virtual public DiagDivider<T>
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
	DiagLDivEq(*itsm,m.QuickTranspose()); 
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
	DiagLDivEq(*itsm,(m0=m1).QuickTranspose()); 
      } 

#include "TMV_AuxAllDiv.h"

      T Det() const;
      DiagMatrix<T> DInverse() const;
      DiagMatrix<T> DInverseATA() const;
      inline Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(DInverse()); }
      inline Matrix<T,ColMajor> InverseATA() const 
      { return Matrix<T,ColMajor>(DInverseATA()); }

      inline bool Singular() const { return Det() == T(0); }

      inline Matrix<T,ColMajor> LU_GetL() const 
      { return Eye<T,ColMajor>(itsm->size(),T(1)); }
      inline Matrix<T,ColMajor> LU_GetU() const 
      { return Matrix<T,ColMajor>(*itsm); }
      inline Permutation LU_GetP() const { return Permutation(itsm->size()); }
      inline Permutation LU_GetQ() const { return Permutation(itsm->size()); }

      inline std::string Type() const
      { 
	return std::string("DiagDiv<") + tmv::Type(T()) + ">"; 
      }

    protected:

      const GenDiagMatrix<T> *const itsm;
      mutable T det;
      mutable bool calcdet;
  };

} // namespace tmv

#endif
