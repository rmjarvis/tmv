//---------------------------------------------------------------------------
//
// This file contains the code for doing division of Triangular matrices.
//
// This is done using back or forward substitution.


#ifndef TMV_TriDiv_H
#define TMV_TriDiv_H

#include "TMV_TriMatrix.h"
#include "TMV_TriDivider.h"

namespace tmv {

  template <class T1, class T2> void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T2>& m);
  template <class T1, class T2> void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T2>& m);
  template <class T1, class T2> void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const UpperTriMatrixView<T2>& m);
  template <class T1, class T2> void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const LowerTriMatrixView<T2>& m);

#define CT complex<T>
  template <class T> inline void TriLDivEq(
      const GenUpperTriMatrix<CT>& A, const MatrixView<T>& m)
  { TMVAssert(false); }
  template <class T> inline void TriLDivEq(
      const GenLowerTriMatrix<CT>& A, const MatrixView<T>& m)
  { TMVAssert(false); }
  template <class T> inline void TriLDivEq(
      const GenUpperTriMatrix<CT>& A, const UpperTriMatrixView<T>& m)
  { TMVAssert(false); }
  template <class T> inline void TriLDivEq(
      const GenLowerTriMatrix<CT>& A, const LowerTriMatrixView<T>& m)
  { TMVAssert(false); }
#undef CT

  template <class T> class UpperTriDiv : 
    virtual public BaseLUDiv<T>,
    virtual public UpperTriDivider<T>
  {

    public :

      //
      // Constructors
      //

      UpperTriDiv(const GenUpperTriMatrix<T>& A) : 
	itsm(&A), det(T(1)), donedet(false) {}
      ~UpperTriDiv() {}

      //
      // Divider Versions of DivEq and Div
      //

      template <class T2> inline void DoLDivEq(const MatrixView<T2>& m) const 
      { TriLDivEq(*itsm,m); } 

      template <class T2> inline void DoRDivEq(const MatrixView<T2>& m) const 
      { TriLDivEq(itsm->QuickTranspose(),m.QuickTranspose()); } 

      template <class T2> inline void DoLDivEq(
	  const UpperTriMatrixView<T2>& m) const 
      { TriLDivEq(*itsm,m); } 

      template <class T2> inline void DoRDivEq(
	  const UpperTriMatrixView<T2>& m) const 
      { TriLDivEq(itsm->QuickTranspose(),m.QuickTranspose()); } 

      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m1, 
	  const MatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.colsize() == m2.colsize()); 
	TMVAssert(m1.rowsize() == m2.rowsize());
	TriLDivEq(*itsm,m2=m1); 
      } 

      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m1, 
	  const MatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.colsize() == m2.colsize()); 
	TMVAssert(m1.rowsize() == m2.rowsize()); 
	TriLDivEq(itsm->QuickTranspose(),(m2=m1).QuickTranspose()); 
      } 

      template <class T1, class T2> inline void DoLDiv(
	  const GenUpperTriMatrix<T1>& m1, 
	  const UpperTriMatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.size() == m2.size()); 
	TriLDivEq(*itsm,m2=m1); 
      } 

      template <class T1, class T2> inline void DoRDiv(
	  const GenUpperTriMatrix<T1>& m1, 
	  const UpperTriMatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.size() == m2.size()); 
	TriLDivEq(itsm->QuickTranspose(),(m2=m1).QuickTranspose()); 
      } 

#include "TMV_AuxAllDiv.h"
#include "TMV_AuxUpperTriAllDiv.h"
      
      //
      // Determinant, Inverse
      //

      T Det() const;
      UpperTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const;
      Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(TInverse()); }
      Matrix<T,ColMajor> InverseATA() const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      inline Matrix<T,ColMajor> LU_GetL() const 
      { return Eye<T,ColMajor>(itsm->size()); }
      inline Matrix<T,ColMajor> LU_GetU() const
      { return Matrix<T,ColMajor>(*itsm); }
      inline Permutation LU_GetP() const { return Permutation(itsm->size()); }
      inline Permutation LU_GetQ() const { return Permutation(itsm->size()); }

      inline std::string Type() const
      { return std::string("UpperTriDiv<") + tmv::Type(T()) + ">"; }

    private :
      const GenUpperTriMatrix<T>*const itsm;
      mutable T det;
      mutable bool donedet;

  }; // UpperTriDiv

  template <class T> class LowerTriDiv : 
    virtual public BaseLUDiv<T>,
    virtual public LowerTriDivider<T>
  {

    public :

      //
      // Constructors
      //

      LowerTriDiv(const GenLowerTriMatrix<T>& A) : 
	itsm(&A), det(T(1)), donedet(false) {}
      ~LowerTriDiv() {}

      //
      // Divider Versions of DivEq and Div
      //

      template <class T2> inline void DoLDivEq(const MatrixView<T2>& m) const 
      { TriLDivEq(*itsm,m); } 

      template <class T2> inline void DoRDivEq(const MatrixView<T2>& m) const 
      { TriLDivEq(itsm->QuickTranspose(),m.QuickTranspose()); } 

      template <class T2> inline void DoLDivEq(
	  const LowerTriMatrixView<T2>& m) const 
      { TriLDivEq(*itsm,m); } 

      template <class T2> inline void DoRDivEq(
	  const LowerTriMatrixView<T2>& m) const 
      { TriLDivEq(itsm->QuickTranspose(),m.QuickTranspose()); } 

      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m1, 
	  const MatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.colsize() == m2.colsize()); 
	TMVAssert(m1.rowsize() == m2.rowsize());
	TriLDivEq(*itsm,m2=m1); 
      } 

      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m1, 
	  const MatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.colsize() == m2.colsize()); 
	TMVAssert(m1.rowsize() == m2.rowsize()); 
	TriLDivEq(itsm->QuickTranspose(),(m2=m1).QuickTranspose()); 
      } 

      template <class T1, class T2> inline void DoLDiv(
	  const GenLowerTriMatrix<T1>& m1, 
	  const LowerTriMatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.size() == m2.size()); 
	TriLDivEq(*itsm,m2=m1); 
      } 

      template <class T1, class T2> inline void DoRDiv(
	  const GenLowerTriMatrix<T1>& m1, 
	  const LowerTriMatrixView<T2>& m2) const 
      { 
	TMVAssert(m1.size() == m2.size()); 
	TriLDivEq(itsm->QuickTranspose(),(m2=m1).QuickTranspose()); 
      } 

#include "TMV_AuxAllDiv.h"
#include "TMV_AuxLowerTriAllDiv.h"
      
      //
      // Determinant, Inverse
      //

      T Det() const;
      LowerTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const;
      Matrix<T,ColMajor> Inverse() const 
      { return Matrix<T,ColMajor>(TInverse()); }
      Matrix<T,ColMajor> InverseATA() const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      inline Matrix<T,ColMajor> LU_GetL() const 
      { return Matrix<T,ColMajor>(*itsm); }
      inline Matrix<T,ColMajor> LU_GetU() const
      { return Eye<T,ColMajor>(itsm->size()); }
      inline Permutation LU_GetP() const { return Permutation(itsm->size()); }
      inline Permutation LU_GetQ() const { return Permutation(itsm->size()); }

      inline std::string Type() const
      { return std::string("LowerTriDiv<") + tmv::Type(T()) + ">"; }

    private :
      const GenLowerTriMatrix<T>*const itsm;
      mutable T det;
      mutable bool donedet;

  }; // LowerTriDiv

} // namespace mv

#endif
