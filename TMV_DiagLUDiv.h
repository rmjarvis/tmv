//---------------------------------------------------------------------------
//
// This file defines the Divider helper for the DiagMatrix class.
//
#ifndef TMV_DiagLUDiv_H
#define TMV_DiagLUDiv_H

#include "TMV_DiagMatrix.h"
#include "TMV_Divider.h"

namespace tmv {

  template <class T1, class T2> void LU_LDivEq(
      const GenDiagMatrix<T1>& d, const MatrixView<T2>& m);

#define CT complex<T>
  template <class T> inline void LU_LDivEq(
      const GenDiagMatrix<CT>& d, const MatrixView<T>& m)
  { TMVAssert(false); }
#undef CT

  template <class T> class DiagLUDiv : 
    virtual public BaseLUDiv<T>, 
    virtual public BaseQRDiv<T>, 
    virtual public DiagDivider<T>
  {

    public:

      DiagLUDiv(const GenDiagMatrix<T>& m) :
	itsm(&m), det(T(1)), calcdet(false) {}
      ~DiagLUDiv() {}

      template <class T2> inline void DoLDivEq(const MatrixView<T2>& m) const 
      { 
	TMVAssert(m.colsize() == itsm->size());
	LU_LDivEq(*itsm,m); 
      } 

      template <class T2> inline void DoRDivEq(const MatrixView<T2>& m) const 
      { 
	TMVAssert(m.rowsize() == itsm->size());
	LU_LDivEq(*itsm,m.QuickTranspose()); 
      } 

      template <class T2, class T3> inline void DoLDiv(
	  const GenMatrix<T2>& m1, const MatrixView<T3>& m0) const 
      {  
	TMVAssert(m1.rowsize() == m0.rowsize()); 
	TMVAssert(m1.colsize() == m0.colsize());
	TMVAssert(m0.colsize() == itsm->size()); 
	LU_LDivEq(*itsm,m0=m1); 
      } 

      template <class T2, class T3> inline void DoRDiv(
	  const GenMatrix<T2>& m1, const MatrixView<T3>& m0) const 
      {  
	TMVAssert(m1.rowsize() == m0.rowsize()); 
	TMVAssert(m1.colsize() == m0.colsize());
	TMVAssert(m0.rowsize() == itsm->size());
	LU_LDivEq(*itsm,(m0=m1).QuickTranspose()); 
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

      inline Matrix<T,ColMajor> QR_GetQ() const 
      { return Eye<T,ColMajor>(itsm->size(),T(1)); }
      inline Matrix<T,ColMajor> QR_GetR() const 
      { return Matrix<T,ColMajor>(*itsm); }
      inline Permutation QR_GetP() const { return Permutation(itsm->size()); }
      inline bool QR_IsTrans() const { return false; }

      inline std::string Type() const
      { 
	return std::string("DiagLUDiv<") + tmv::Type(T()) + ">"; 
      }

    protected:

      const GenDiagMatrix<T> *const itsm;
      mutable T det;
      mutable bool calcdet;
  };

} // namespace tmv

#endif
