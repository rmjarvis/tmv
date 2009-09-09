//---------------------------------------------------------------------------
//
// This file defines the TMV BaseMatrix class.
//
// This base class defines some of the things that all 
// matrices need to be able to do, as well as some of the
// arithmetic operations (those that return a Vector).
// This should be used as the base class for generic
// matrices as well as any special ones (eg. sparse,
// symmetric, etc.)
//
//

#ifndef TMV_BaseMatrix_H
#define TMV_BaseMatrix_H

#include "TMV_Vector.h"

namespace tmv {

  template <class T> class Divider;

  template <class T> class BaseMatrix
  {

    public:

      //
      // Constructors
      //

      BaseMatrix<T>(const DivType dt=XXX) : 
	itsdiv(0), itsdt(dt), divinplace(false), cachediv(false) {}
      BaseMatrix(const BaseMatrix<T>& rhs) : 
	itsdiv(0), itsdt(rhs.itsdt), divinplace(false), cachediv(false) {}
      virtual ~BaseMatrix() { if (itsdiv) delete itsdiv; }

      //
      // Access Functions
      //

      virtual size_t colsize() const = 0;
      virtual size_t rowsize() const = 0;
      inline size_t ncols() const { return rowsize(); }
      inline size_t nrows() const { return colsize(); }
      inline virtual bool IsSquare() const { return colsize() == rowsize(); }

      // Cannot be virtual because of FortranStyle indexing.
      // Plus, no advantage to virtual here.
      inline T operator()(size_t i, size_t j) const 
      {
	TMVAssert(i<colsize() && j<rowsize());
	return cref(i,j); 
      }

      virtual void CopyToMatrix(
	  const MatrixView<RealType(T)>& m) const = 0;
      virtual void CopyToMatrix(
	  const MatrixView<ComplexType(T)>& m) const = 0;

      //
      // Functions of Matrix
      //

      virtual inline T Det() const
      { SetDiv(); T det = itsdiv->Det(); DoneDiv(); return det; }

      virtual T Trace() const = 0;

      virtual inline void Inverse(const MatrixView<T>& minv) const
      { SetDiv(); itsdiv->Inverse(minv); DoneDiv(); }

      virtual inline void InverseATA(const MatrixView<T>& minv) const
      { SetDiv(); itsdiv->InverseATA(minv); DoneDiv(); }

      virtual inline bool Singular() const
      { SetDiv(); bool s = itsdiv->Singular(); DoneDiv(); return s; }

      inline RealType(T) Norm() const { return NormF(); }
      virtual RealType(T) NormSq() const = 0;
      inline RealType(T) NormF() const { return SQRT(NormSq()); }
      virtual RealType(T) Norm1() const = 0;
      virtual RealType(T) Norm2() const 
      {
	if (colsize() == 0 || rowsize() == 0) return RealType(T)(0);
	else {
	  if (SV_Type(itsdt)) SetDiv();
	  else {
	    if (!itsdiv || !itsdiv->IsSV()) {
#ifdef TMVDEBUG
	      cout<<"Warning: calling Norm2 without previously calling ";
	      cout<<"DivideUsing(SV)\n";
#endif
	      DivideUsing(SVS);
	    }
	    SetDiv();
	  }
          RealType(T) norm2 = itsdiv->Norm2();
	  DoneDiv();
	  return norm2;
	}
      }
      virtual RealType(T) Condition() const 
      {
	if (colsize() == 0 || rowsize() == 0) return RealType(T)(0);
	else {
	  if (SV_Type(itsdt)) SetDiv();
	  else {
	    if (!itsdiv || !itsdiv->IsSV()) {
#ifdef TMVDEBUG
	      cout<<"Warning: calling Condition without previously calling ";
	      cout<<"DivideUsing(SV)\n";
#endif
	      DivideUsing(SVS);
	    }
	    SetDiv();
	  }
	  RealType(T) cond = itsdiv->Condition();
	  DoneDiv();
	  return cond;
	}
      }
      virtual RealType(T) NormInf() const = 0;

      virtual BaseMatrix<T>* NewTranspose() const = 0;
      virtual BaseMatrix<T>* NewConjugate() const = 0;
      virtual BaseMatrix<T>* NewAdjoint() const = 0;
      virtual inline BaseMatrix<T>* NewInverse() const 
      { 
	Matrix<T,ColMajor>* minv = new Matrix<T,ColMajor>(rowsize(),colsize());
	Inverse(minv->View());
	return minv;
      }
      virtual BaseMatrix<T>* NewView() const = 0;
      virtual BaseMatrix<T>* NewCopy() const = 0;

      // 
      // I/O: Write
      //

      virtual void Write(ostream& os) const
      {
	os << colsize() <<"  "<<rowsize()<<endl;
	for(size_t i=0;i<colsize();++i) {
	  os << "( ";
	  for(size_t j=0;j<rowsize();++j) os << ' ' << cref(i,j) << ' ';
	  os << " )"<<endl;
	}
      }

      virtual void Write(ostream& os, RealType(T) thresh) const
      {
	os << colsize() <<"  "<<rowsize()<<endl;
	for(size_t i=0;i<colsize();++i) {
	  os << "( ";
	  for(size_t j=0;j<rowsize();++j) {
	    const T mij = cref(i,j);
	    os << ' ' <<(abs(mij) < thresh ? T(0) : mij) << ' ';
	  }
	  os << " )"<<endl;
	}
      }

      //
      // Arithmetic Helpers
      //

      // m^-1 * v -> v
      template <class T1> inline void LDivEq(const VectorView<T1>& v) const 
      { SetDiv(); itsdiv->LDivEq(ColVectorViewOf(v)); DoneDiv(); }

      template <class T1> inline void LDivEq(const MatrixView<T1>& m) const 
      { SetDiv(); itsdiv->LDivEq(m); DoneDiv(); }

      // v * m^-1 -> v
      template <class T1> inline void RDivEq(const VectorView<T1>& v) const 
      { SetDiv(); itsdiv->RDivEq(RowVectorViewOf(v)); DoneDiv(); }

      template <class T1> inline void RDivEq(const MatrixView<T1>& m) const 
      { SetDiv(); itsdiv->RDivEq(m); DoneDiv(); }

      // m^-1 * v1 -> v0
      template <class T1, class T0> inline void LDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { SetDiv(); itsdiv->LDiv(ColVectorViewOf(v1),ColVectorViewOf(v0)); DoneDiv(); }
     
      template <class T1, class T0> inline void LDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { SetDiv(); itsdiv->LDiv(m1,m0); DoneDiv(); }

      // v1 * m^-1 -> v0
      template <class T1, class T0> inline void RDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { SetDiv(); itsdiv->RDiv(RowVectorViewOf(v1),RowVectorViewOf(v0)); DoneDiv(); }

      template <class T1, class T0> inline void RDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { SetDiv(); itsdiv->RDiv(m1,m0); DoneDiv(); }

      //
      // Division Control
      //

      inline void DivideInPlace() const 
      { divinplace = true; }

      inline void SaveDiv() const 
      { cachediv = true; }

      inline void DivideUsing(DivType dt) const 
      {
	if (dt != itsdt) UnSetDiv();
	itsdt = dt;
      }

      inline void SetDiv() const 
      { 
	if (!itsdiv) {
	  if (itsdt == XXX) itsdt = IsSquare() ? LU : QR;
	  NewDivider();
	} 
      }

      inline void UnSetDiv() const 
      { if (itsdiv) { delete itsdiv; itsdiv = 0; } }

      inline void ReSetDiv() const 
      { UnSetDiv(); SetDiv(); }

      inline void DoneDiv() const
      { if (!cachediv) UnSetDiv(); }

      inline const Divider<T>* GetDiv() const { return itsdiv; }
      inline DivType GetDivType() const { return itsdt; }
      inline bool CheckDecomp(ostream* fout=0) const 
      { 
	TMVAssert(itsdiv);
	return itsdiv->CheckDecomp(*this,fout);
      }
      inline bool CheckDecomp(const BaseMatrix<T>& m2, ostream* fout=0) const 
      { 
	TMVAssert(itsdiv);
	return itsdiv->CheckDecomp(m2,fout);
      }

    protected :

      mutable const Divider<T>* itsdiv;
      mutable DivType itsdt;
      mutable bool divinplace;
      mutable bool cachediv;

      virtual T cref(size_t i, size_t j) const = 0;
      virtual void NewDivider() const = 0;

    private :

      void operator=(const BaseMatrix<T>&) { TMVAssert(false); }

  }; // BaseMatrix

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const BaseMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const BaseMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const BaseMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormSq(const BaseMatrix<T>& m)
  { return m.NormSq(); }
  
  template <class T> inline RealType(T) NormF(const BaseMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm1(const BaseMatrix<T>& m)
  { return m.Norm1(); }

  template <class T> inline RealType(T) Norm2(const BaseMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const BaseMatrix<T>& m)
  { return m.NormInf(); }

  //
  // I/O
  //

  template <class T> ostream& operator<<(ostream& os, const BaseMatrix<T>& m)
  { m.Write(os); return os; }

}; // namespace tmv

#endif
