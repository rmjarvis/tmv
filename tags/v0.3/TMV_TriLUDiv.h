//---------------------------------------------------------------------------
//
// This file contains the code for doing division of Triangular matrices
// using LU Decomposition.
//
// This is obviously pretty trivial, since a triangular matrix is already
// in LU form (either IU or LI).
// But these files (this and the corresponding .cpp file) have the code
// for doing the actual division of the L and U matrices that the full LU
// decomposition depends on.
//


#ifndef TMV_TriLUDiv_H
#define TMV_TriLUDiv_H

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

  template <class T> class UpperTriLUDiv : 
    virtual public BaseLUDiv<T>, 
    virtual public BaseQRDiv<T>, 
    virtual public UpperTriDivider<T>
  {

    public :

      //
      // Constructors
      //

      UpperTriLUDiv(const GenUpperTriMatrix<T>& A) : 
	itsm(&A), det(T(1)), donedet(false) {}
      ~UpperTriLUDiv() {}

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

      inline Matrix<T,ColMajor> QR_GetQ() const 
      { return Eye<T,ColMajor>(itsm->size()); }
      inline Matrix<T,ColMajor> QR_GetR() const
      { return Matrix<T,ColMajor>(*itsm); }
      inline Permutation QR_GetP() const 
      { return Permutation(itsm->size()); }
      inline bool QR_IsTrans() const { return false; }

      inline std::string Type() const
      { return std::string("UpperTriLUDiv<") + tmv::Type(T()) + ">"; }

    private :
      const GenUpperTriMatrix<T>*const itsm;
      mutable T det;
      mutable bool donedet;

  }; // UpperTriLUDiv

  template <class T> class LowerTriLUDiv : 
    virtual public BaseLUDiv<T>, 
    virtual public BaseQRDiv<T>, 
    virtual public LowerTriDivider<T>
  {

    public :

      //
      // Constructors
      //

      LowerTriLUDiv(const GenLowerTriMatrix<T>& A) : 
	itsm(&A), det(T(1)), donedet(false) {}
      ~LowerTriLUDiv() {}

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

      inline Matrix<T,ColMajor> QR_GetQ() const 
      { return Matrix<T,ColMajor>(ReversePermutation(itsm->size())); }
      inline Matrix<T,ColMajor> QR_GetR() const
      { 
	Permutation rp = ReversePermutation(itsm->size());
	Matrix<T,ColMajor> temp(*itsm);
	return rp * temp * rp;
      }
      inline Permutation QR_GetP() const 
      { return ReversePermutation(itsm->size()); }
      inline bool QR_IsTrans() const { return false; }

      inline std::string Type() const
      { return std::string("LowerTriLUDiv<") + tmv::Type(T()) + ">"; }

    private :
      const GenLowerTriMatrix<T>*const itsm;
      mutable T det;
      mutable bool donedet;

  }; // LowerTriLUDiv

} // namespace mv

#endif
