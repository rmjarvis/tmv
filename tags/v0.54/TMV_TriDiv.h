//---------------------------------------------------------------------------
//
// This file contains the code for doing division of Triangular matrices.
//
// This is done using back or forward substitution.


#ifndef TMV_TriDiv_H
#define TMV_TriDiv_H

#include "TMV_TriMatrix.h"
#include "TMV_Divider.h"

namespace tmv {

  template <class T, class T1> void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const VectorView<T>& v);
  template <class T, class T1> void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const VectorView<T>& v);
  template <class T, class T1> void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T>& m);
  template <class T, class T1> void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T>& m);
  template <class T, class T1> void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const UpperTriMatrixView<T>& m);
  template <class T, class T1> void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const LowerTriMatrixView<T>& m);

#define CT complex<T>
  template <class T> inline void TriLDivEq(
      const GenUpperTriMatrix<CT>& , const VectorView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void TriLDivEq(
      const GenLowerTriMatrix<CT>& , const VectorView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void TriLDivEq(
      const GenUpperTriMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void TriLDivEq(
      const GenLowerTriMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void TriLDivEq(
      const GenUpperTriMatrix<CT>& , const UpperTriMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void TriLDivEq(
      const GenLowerTriMatrix<CT>& , const LowerTriMatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

   template <class T> class SingularUpperTriMatrix :
     public Singular
     {
       public:
	 const GenUpperTriMatrix<T>& A;

	 SingularUpperTriMatrix(const GenUpperTriMatrix<T>& _A) :
	   Singular("UpperTriMatrix"), A(_A) {}
	 ~SingularUpperTriMatrix() throw() {}
	 void Write(ostream& os) const throw()
	 {
	   Singular::Write(os);
	   os<<A<<endl;
	 }
     };

   template <class T> class SingularLowerTriMatrix :
     public Singular
     {
       public:
	 const GenLowerTriMatrix<T>& A;

	 SingularLowerTriMatrix(const GenLowerTriMatrix<T>& _A) :
	   Singular("LowerTriMatrix"), A(_A) {}
	 ~SingularLowerTriMatrix() throw() {}
	 void Write(ostream& os) const throw()
	 {
	   Singular::Write(os);
	   os<<A<<endl;
	 }
     };

  template <class T> class UpperTriDiv : 
    public Divider<T>
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
      { TriLDivEq(itsm->Transpose(),m.Transpose()); } 

      template <class T2> inline void DoLDivEq(
	  const UpperTriMatrixView<T2>& m) const 
      { TriLDivEq(*itsm,m); } 

      template <class T2> inline void DoRDivEq(
	  const UpperTriMatrixView<T2>& m) const 
      { TriLDivEq(itsm->Transpose(),m.Transpose()); } 

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
	TriLDivEq(itsm->Transpose(),(m2=m1).Transpose()); 
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
	TriLDivEq(itsm->Transpose(),(m2=m1).Transpose()); 
      } 

#include "TMV_AuxAllDiv.h"

#define RT RealType(T)
#define CT ComplexType(T)
#define DefDivEq(T) \
      inline void LDivEq(const UpperTriMatrixView<T>& m) const \
      { DoLDivEq(m); } \
      inline void RDivEq(const UpperTriMatrixView<T>& m) const \
      { DoRDivEq(m); } \

      DefDivEq(RT);
      DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
      inline void LDiv(const GenUpperTriMatrix<T1>& m1, \
	  const UpperTriMatrixView<T2>& m0) const \
      { DoLDiv(m1,m0); } \
      inline void RDiv(const GenUpperTriMatrix<T1>& m1, \
	  const UpperTriMatrixView<T2>& m0) const \
      { DoRDiv(m1,m0); } \

      DefDiv(RT,RT);
      DefDiv(RT,CT);
      DefDiv(CT,CT);
#undef DefDiv
#undef RT
#undef CT
      
      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const UpperTriMatrixView<T>& minv) const;
      inline void Inverse(const MatrixView<T>& minv) const 
      {
	TMVAssert(minv.IsSquare());
	minv.Zero();
	Inverse(UpperTriMatrixViewOf(minv));
      }
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }

      inline string Type() const
      { return string("UpperTriDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return LU; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    private :
      const GenUpperTriMatrix<T>*const itsm;
      mutable T det;
      mutable bool donedet;

  }; // UpperTriDiv

  template <class T> class LowerTriDiv : 
    public Divider<T>
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
      { TriLDivEq(itsm->Transpose(),m.Transpose()); } 

      template <class T2> inline void DoLDivEq(
	  const LowerTriMatrixView<T2>& m) const 
      { TriLDivEq(*itsm,m); } 

      template <class T2> inline void DoRDivEq(
	  const LowerTriMatrixView<T2>& m) const 
      { TriLDivEq(itsm->Transpose(),m.Transpose()); } 

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
	TriLDivEq(itsm->Transpose(),(m2=m1).Transpose()); 
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
	TriLDivEq(itsm->Transpose(),(m2=m1).Transpose()); 
      } 

#include "TMV_AuxAllDiv.h"
      
#define RT RealType(T)
#define CT ComplexType(T)
#define DefDivEq(T) \
      inline void LDivEq(const LowerTriMatrixView<T>& m) const \
      { DoLDivEq(m); } \
      inline void RDivEq(const LowerTriMatrixView<T>& m) const \
      { DoRDivEq(m); } \

      DefDivEq(RT);
      DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
      inline void LDiv(const GenLowerTriMatrix<T1>& m1, \
	  const LowerTriMatrixView<T2>& m0) const \
      { DoLDiv(m1,m0); } \
      inline void RDiv(const GenLowerTriMatrix<T1>& m1, \
	  const LowerTriMatrixView<T2>& m0) const \
      { DoRDiv(m1,m0); } \

      DefDiv(RT,RT);
      DefDiv(RT,CT);
      DefDiv(CT,CT);
#undef DefDiv
#undef RT
#undef CT

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const LowerTriMatrixView<T>& minv) const;
      inline void Inverse(const MatrixView<T>& minv) const 
      {
	TMVAssert(minv.IsSquare());
	Inverse(LowerTriMatrixViewOf(minv));
	if (minv.colsize() > 0)
	  UpperTriMatrixViewOf(minv,UnitDiag).OffDiag().Zero(); 
      }
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }

      inline string Type() const
      { return string("LowerTriDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return LU; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    private :
      const GenLowerTriMatrix<T>*const itsm;
      mutable T det;
      mutable bool donedet;

  }; // LowerTriDiv

} // namespace mv

#endif
