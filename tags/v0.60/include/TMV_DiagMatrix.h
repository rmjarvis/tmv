///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// This file defines the TMV DiagMatrix class.
//
// The DiagMatrix class is provided for efficient storage of a diagonal
// matrix.  You can do most of the things that you can do with a 
// regular Matrix, but it will do them more efficiently.
//
// Constructors:
//
//    DiagMatrix<T>(size_t size)
//        Makes a DiagMatrix with column size and row size = size
//        with _uninitialized_ values
//
//    DiagMatrix<T>(size_t size, T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    DiagMatrix<T>(const Vector<T>& vv)
//        Make a DiagMatrix which copies the elements of vv.
//
//    ConstDiagMatrixView<T>(const Vector<T>& v)
//        Make a constant DiagMatrix view with v as the diagonal.
//        While this view cannon be modified, changing the original v or m
//        will cause corresponding changes in this view.
//
//    DiagMatrixView<T>(Vector<T>& v)
//        Make a mutable DiagMatrix view with v as the diagonal.
//
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the DiagMatrix
//
//    T& operator()(size_t i)
//    T operator()(size_t i) const
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the DiagMatrix
//        For the single paramter version, j=i
//
//    VectorView& diag()
//    ConstVectorView& diag() const
//        Return the diagonal of the DiagMatrix as a VectorView
//
//
// Modifying Functions - The same as the regular Matrix counterparts
//
//    DiagMatrix& Zero()
//    DiagMatrix& SetAllTo(T x)
//    DiagMatrix<T>& TransposeSelf() 
//        (Does nothing.)
//    DiagMatrix& ConjugateSelf()
//    DiagMatrix& SetToIdentity(x = 1)
//    void Swap(DiagMatrix& m1, DiagMatrix& m2)
//
//
// SubDiagMatrix:
//
//    SubDiagMatrix(int i1, int i2, int istep=1)
//        Returns a Sub-DiagMatrix which extends from i1 to i2 (step istep)
//        which refers to the same physical elements as the original.
//        As usual, i2 is the "one past the end" element.
//
//
// Functions of DiagMatrices - Same as for regular Matrices:
//
//    Det(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//    MaxAbsElement(m) 
//        (Note - for diagonal matrices, 
//        Norm1 = Norm2 = NormInf = MaxAbsElement.)
//    Transpose(m)
//        (same as the original).
//    Conjugate(m)
//    Adjoint(m)
//
//    m.Inverse()
//    Inverse(m)
//    m.InvertSelf()
//    m.Inverse(minv) (takes either Matrix or DiagMatrix argument)
//    m.InverseATA(invata) (takes either Matrix or DiagMatrix argument)
//
// Operators:
//       By default with TMV_Diag.h, we define the operators for
//       two DiagMatrix objects and for a DiagMatrix and either a 
//       Matrix or a Vector.  For combinations with other Sparse Types
//       you need to include TMV_Diag.h before them in sequence.  
//       ie. 
//       #include "TMV.h"
//       #include "TMV_Diag.h"
//       #include "TMV_Band.h"
//       will include the necessary files to do arithmetic with DiagMatrix 
//       and BandMatrix objects.
//
//
// I/O: 
//
//    os << d 
//        Writes d to ostream os as a full matrix
//
//    d.WriteCompact(os)
//        Writes only the diagonal Vector to os
//
//    is >> d
//        Reads in d in the compact format
//
//

#ifndef TMV_DiagMatrix_H
#define TMV_DiagMatrix_H

#include "TMV_BaseDiagMatrix.h"
#include "TMV_Vector.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T> class GenDiagMatrix : 
    virtual public AssignableToDiagMatrix<T>,
    virtual public AssignableToUpperTriMatrix<T>,
    virtual public AssignableToLowerTriMatrix<T>,
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline GenDiagMatrix<T>() {}
      inline GenDiagMatrix(const GenDiagMatrix<T>&) {}
      virtual inline ~GenDiagMatrix() {}

      //
      // Access Functions
      //

      inline size_t size() const { return diag().size(); }
      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }
      inline DiagType dt() const { return NonUnitDiag; }

      inline T operator()(size_t i, size_t j) const 
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	if (i==j) return diag()(i); 
	else return T(0);
      }

      inline T operator()(size_t i) const 
      {
	TMVAssert(i<size());
	return diag()(i); 
      }

      inline ConstVectorView<T> diag() const { return cdiag(); }

      template <class T2> inline bool SameAs(const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameAs(const GenDiagMatrix<T>& m2) const
      { 
	if (this == &m2) return true;
	else return (diag().SameAs(m2.diag())); 
      }

      void DoAssignToM(const MatrixView<RealType(T)>& m2) const;
      void DoAssignToM(const MatrixView<ComplexType(T)>& m2) const;
      inline void AssignToM(const MatrixView<RealType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	TMVAssert(IsReal(T()));
	DoAssignToM(m2);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	DoAssignToM(m2);
      }
      
      void DoAssignToU(const UpperTriMatrixView<RealType(T)>& m2) const;
      void DoAssignToU(const UpperTriMatrixView<ComplexType(T)>& m2) const;
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m2) const
      {
	TMVAssert(m2.size() == size());
	TMVAssert(IsReal(T()));
	DoAssignToU(m2);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.size() == size());
	DoAssignToU(m2);
      }

      void DoAssignToL(const LowerTriMatrixView<RealType(T)>& m2) const;
      void DoAssignToL(const LowerTriMatrixView<ComplexType(T)>& m2) const;
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m2) const
      {
	TMVAssert(m2.size() == size());
	TMVAssert(IsReal(T()));
	DoAssignToL(m2);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.size() == size());
	DoAssignToL(m2);
      }

      inline void AssignToD(const DiagMatrixView<RealType(T)>& m2) const
      { 
	TMVAssert(m2.size() == size());
	TMVAssert(IsReal(T()));
	if (!SameAs(m2)) m2.diag() = diag(); 
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.size() == size());
	if (!SameAs(m2)) m2.diag() = diag(); 
      }

      //
      // SubDiagMatrix
      //

      inline ConstDiagMatrixView<T> SubDiagMatrix(int i1, int i2) const
      { 
	TMVAssert(diag().OKSubVector(i1,i2,1));
	return ConstDiagMatrixView<T>(diag().SubVector(i1,i2)); 
      }

      inline ConstDiagMatrixView<T> SubDiagMatrix(int i1, int i2,
	  int istep) const
      { 
	TMVAssert(diag().OKSubVector(i1,i2,istep));
	return ConstDiagMatrixView<T>(diag().SubVector(i1,i2,istep)); 
      }

      inline ConstDiagMatrixView<RealType(T)> Real() const
      { return ConstDiagMatrixView<RealType(T)>(diag().Real()); }

      inline ConstDiagMatrixView<RealType(T)> Imag() const
      { return ConstDiagMatrixView<RealType(T)>(diag().Imag()); }

      inline ConstDiagMatrixView<T> View() const
      { 
	return ConstDiagMatrixView<T>(diag());
      }

      inline ConstDiagMatrixView<T> Transpose() const
      { return View(); }

      inline ConstDiagMatrixView<T> Conjugate() const
      { 
	return ConstDiagMatrixView<T>(diag().Conjugate());
      }

      inline ConstDiagMatrixView<T> Adjoint() const
      { return Conjugate(); }


      //
      // Functions of DiagMatrix
      //

      T Det() const;

      inline T Trace() const
      { return diag().SumElements(); }

      inline RealType(T) Norm() const 
      { return NormF(); }

      inline RealType(T) NormF() const 
      { return SQRT(diag().NormSq()); }

      inline RealType(T) NormSq() const 
      { return diag().NormSq(); }

      inline RealType(T) Norm1() const 
      { return diag().MaxAbsElement(); }

      inline RealType(T) Norm2() const 
      { return diag().MaxAbsElement(); }

      inline RealType(T) Condition() const 
      { return diag().MaxAbsElement()/diag().MinAbsElement(); }

      inline RealType(T) NormInf() const
      { return diag().MaxAbsElement(); }

      inline RealType(T) MaxAbsElement() const
      { return diag().MaxAbsElement(); }

      inline bool Singular() const 
      { return diag().MinAbsElement() == RealType(T)(0); }

      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      template <class T1> void DoInverse(const DiagMatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& ata) const;
      void DoInverseATA(const DiagMatrixView<T>& ata) const;

      inline void Inverse(const MatrixView<T>& minv) const
      { 
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	DoInverse(minv);
      }
      template <class T1> inline void Inverse(const MatrixView<T1>& minv) const
      { 
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	DoInverse(minv);
      }
      template <class T1> inline void Inverse(
	  const DiagMatrixView<T1>& minv) const
      { 
	TMVAssert(minv.size() == size());
	DoInverse(minv);
      }

      QuotXD<T,T> QInverse() const;
      inline QuotXD<T,T> Inverse() const
      { return QInverse(); }

      inline void InverseATA(const MatrixView<T>& ata) const
      { 
	TMVAssert(ata.colsize() == size());
	TMVAssert(ata.rowsize() == size());
	DoInverseATA(ata);
      }
      inline void InverseATA(const DiagMatrixView<T>& ata) const
      { 
	TMVAssert(ata.size() == size());
	DoInverseATA(ata);
      }

      template <class T1, IndexStyle I> inline void Inverse(
	  DiagMatrix<T1,I>& minv) const
      { Inverse(minv.View()); }

      template <class T1, StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T1,S,I>& minv) const
      { Inverse(minv.View()); }

      template <IndexStyle I> inline void InverseATA(
	  DiagMatrix<T,I>& minv) const
      { InverseATA(minv.View()); }

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& minv) const
      { InverseATA(minv.View()); }

      auto_ptr<BaseMatrix<T> > NewTranspose() const;
      auto_ptr<BaseMatrix<T> > NewConjugate() const;
      auto_ptr<BaseMatrix<T> > NewAdjoint() const;
      auto_ptr<BaseMatrix<T> > NewInverse() const;
      auto_ptr<BaseMatrix<T> > NewView() const;
      auto_ptr<BaseMatrix<T> > NewCopy() const;

      // 
      // I/O
      //

      void Write(std::ostream& os) const;
      void Write(std::ostream& os, RealType(T) thresh) const;

      inline void WriteCompact(std::ostream& os) const
      { os << "D " << diag() << std::endl; }
      inline void WriteCompact(std::ostream& os, RealType(T) thresh) const
      { os << "D "; diag().Write(os,thresh); os << std::endl; }

      // 
      // Arithmetic Helpers
      //

      template <class T1> void DoLDivEq(const VectorView<T1>& v) const;
      template <class T1, class T0> void DoLDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const;
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T0> void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;

      template <class T1> inline void LDivEq(const VectorView<T1>& v) const
      {
	TMVAssert(v.size() == size());
	DoLDivEq(v);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      {
	TMVAssert(v1.size() == size());
	TMVAssert(v0.size() == size());
	DoLDiv(v1,v0);
      }
      template <class T1> inline void RDivEq(const VectorView<T1>& v) const
      { LDivEq(v); }
      template <class T1, class T0> inline void RDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { LDiv(v1,v0); }

      template <class T1> inline void LDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(m.colsize() == size());
	DoLDivEq(m);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      {
	TMVAssert(m1.colsize() == size());
	TMVAssert(m0.colsize() == size());
	TMVAssert(m1.rowsize() == m0.rowsize());
	DoLDiv(m1,m0);
      }
      template <class T1> inline void RDivEq(const MatrixView<T1>& m) const
      { LDivEq(m.Transpose()); }
      template <class T1, class T0> void RDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { LDiv(m1.Transpose(),m0.Transpose()); }

      template <class T1> void DivEq(const DiagMatrixView<T1>& m) const
      { LDivEq(m.diag()); }
      template <class T1, class T0> void Div(
	  const GenDiagMatrix<T1>& m1, const DiagMatrixView<T0>& m0) const
      { LDiv(m1.diag(),m0.diag()); }

      // For easier compatibility with regular matrices:
      inline void DivideInPlace() const {}
      inline void SaveDiv() const {}
      inline void SetDiv() const {}
      inline void UnSetDiv() const {}
      inline void ReSetDiv() const {}
      inline void DivideUsing(DivType DEBUGPARAM(dt)) const 
      { TMVAssert(dt == LU); }

    protected :

      virtual ConstVectorView<T> cdiag() const = 0;
      inline T cref(size_t i, size_t j) const
      { return i==j ? cdiag()(i) : 0; }

    private :

      inline void operator=(const GenDiagMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenDiagMatrix

  template <class T, IndexStyle I> class ConstDiagMatrixView : 
    public GenDiagMatrix<T>
  {
    public :

      inline ConstDiagMatrixView(const ConstDiagMatrixView<T,I>& rhs) :
	itsdiag(rhs.cdiag()) {}

      inline ConstDiagMatrixView(const GenDiagMatrix<T>& rhs) :
	itsdiag(rhs.diag()) {}

      explicit inline ConstDiagMatrixView(const GenVector<T>& v) :
	itsdiag(v) {}

      virtual inline ~ConstDiagMatrixView() {}

    protected :

      ConstVectorView<T> itsdiag;

      inline ConstVectorView<T> cdiag() const { return itsdiag; }

    private :

      inline void operator=(const ConstDiagMatrixView<T,I>&) 
      { TMVAssert(FALSE); }

  }; // ConstDiagMatrixView

  template <class T> class ConstDiagMatrixView<T,FortranStyle> : 
    public ConstDiagMatrixView<T,CStyle>
  {
    public :

      inline ConstDiagMatrixView(
	  const ConstDiagMatrixView<T,FortranStyle>& rhs) :
	ConstDiagMatrixView<T,CStyle>(rhs) {}

      inline ConstDiagMatrixView(const GenDiagMatrix<T>& rhs) :
	ConstDiagMatrixView<T,CStyle>(rhs) {}

      explicit inline ConstDiagMatrixView(const GenVector<T>& v) :
	ConstDiagMatrixView<T,CStyle>(v) {}

      virtual inline ~ConstDiagMatrixView() {}

      //
      // Access Functions
      //
      
      inline T operator()(size_t i, size_t j) const 
      { 
	TMVAssert(i>0 && i<=size() && j>0 && j<=size());
	if (i==j) return diag()(i);
	else return T(0);
      }

      inline T operator()(size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return diag()(i);
      }

      inline ConstVectorView<T,FortranStyle> diag() const 
      { 
	return ConstVectorView<T,FortranStyle>(ConstDiagMatrixView<T>::diag());
      }

      //
      // SubDiagMatrix
      //

      inline ConstDiagMatrixView<T,FortranStyle> SubDiagMatrix(
	  int i1, int i2) const
      {
	TMVAssert(diag().OKSubVector(i1-1,i2,1));
	return ConstDiagMatrixView<T,FortranStyle>(diag().SubVector(i1,i2)); 
      }

      inline ConstDiagMatrixView<T,FortranStyle> SubDiagMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(diag().OKSubVector(i1-1,i2-1+istep,istep));
	return ConstDiagMatrixView<T,FortranStyle>(
	    diag().SubVector(i1,i2,istep)); 
      }

      inline ConstDiagMatrixView<RealType(T),FortranStyle> Real() const
      { return ConstDiagMatrixView<RealType(T),FortranStyle>(diag().Real()); }

      inline ConstDiagMatrixView<RealType(T),FortranStyle> Imag() const
      { return ConstDiagMatrixView<RealType(T),FortranStyle>(diag().Imag()); }

      inline ConstDiagMatrixView<T,FortranStyle> View() const
      { 
	return ConstDiagMatrixView<T,FortranStyle>(diag());
      }

      inline ConstDiagMatrixView<T,FortranStyle> Transpose() const
      { return View(); }

      inline ConstDiagMatrixView<T,FortranStyle> Conjugate() const
      { 
	return ConstDiagMatrixView<T,FortranStyle>(diag().Conjugate());
      }

      inline ConstDiagMatrixView<T,FortranStyle> Adjoint() const
      { return Conjugate(); }

      using GenDiagMatrix<T>::size;

    private :

      inline void operator=(const ConstDiagMatrixView<T,FortranStyle>&) 
      { TMVAssert(FALSE); }

  }; // ConstDiagMatrixView - FortranStyle

  template <class T, IndexStyle I> class DiagMatrixView : 
    public GenDiagMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline DiagMatrixView(const DiagMatrixView<T,I>& rhs) :
	itsdiag(rhs.diag()) {}

      explicit inline DiagMatrixView(const VectorView<T>& _diag) :
	itsdiag(_diag) {}

      virtual inline ~DiagMatrixView() {} 

      //
      // Op=
      //

      inline const DiagMatrixView<T,I>& operator=(
	  const DiagMatrixView<T,I>& m2) const
      { m2.AssignToD(*this); return *this; }

      inline const DiagMatrixView<T,I>& operator=(
	  const GenDiagMatrix<RealType(T)>& m2) const
      { m2.AssignToD(*this); return *this; }

      inline const DiagMatrixView<T,I>& operator=(
	  const GenDiagMatrix<ComplexType(T)>& m2) const
      { m2.AssignToD(*this); return *this; }

      template <class T2> inline const DiagMatrixView<T,I>& operator=(
	  const GenDiagMatrix<T2>& m2) const
      { itsdiag = m2.diag(); return *this; }

      inline const DiagMatrixView<T,I>& operator=(T x) const 
      { return SetToIdentity(x); }

      inline const DiagMatrixView<T,I>& operator=(
	  const AssignableToDiagMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	m2.AssignToD(*this);
	return *this;
      }

      inline const DiagMatrixView<T,I>& operator=(
	  const AssignableToDiagMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	m2.AssignToD(*this);
	return *this;
      }


      //
      // Access
      //
 
      inline RefType(T) operator()(size_t i) const 
      { 
	TMVAssert(i<size());
	return diag()(i); 
      }
      inline RefType(T) operator()(size_t i, size_t DEBUGPARAM(j)) const 
      { 
	TMVAssert(i<size());
	TMVAssert(i==j); 
	return diag()(i); 
      }

      inline VectorView<T,I> diag() const { return itsdiag; }

      //
      // Modifying Functions
      //

      inline const DiagMatrixView<T,I>& Zero() const 
      { return SetAllTo(0); }

      inline const DiagMatrixView<T,I>& Clip(RealType(T) thresh) const
      { diag().Clip(thresh); return *this; }

      inline const DiagMatrixView<T,I>& SetAllTo(T x) const
      { diag().SetAllTo(x); return *this; }

      inline const DiagMatrixView<T,I>& TransposeSelf() const
      { return *this; }

      inline const DiagMatrixView<T,I>& ConjugateSelf() const
      { diag().ConjugateSelf(); return *this; }

      const DiagMatrixView<T,I>& InvertSelf() const;

      inline const DiagMatrixView<T,I>& SetToIdentity(T x=T(1)) const 
      { return SetAllTo(x); }

      //
      // SubDiagMatrix
      //

      inline DiagMatrixView<T,I> SubDiagMatrix(int i1, int i2) const
      { return DiagMatrixView<T,I>(itsdiag.SubVector(i1,i2)); }

      inline DiagMatrixView<T,I> SubDiagMatrix(int i1, int i2,
	  int istep) const
      { return DiagMatrixView<T,I>(itsdiag.SubVector(i1,i2,istep)); }

      inline DiagMatrixView<RealType(T),I> Real() const
      { return DiagMatrixView<RealType(T),I>(diag().Real()); }

      inline DiagMatrixView<RealType(T),I> Imag() const
      { return DiagMatrixView<RealType(T),I>(diag().Imag()); }

      inline DiagMatrixView<T,I> View() const
      { return *this; }

      inline DiagMatrixView<T,I> Transpose() const
      { return *this; }

      inline DiagMatrixView<T,I> Conjugate() const
      { return DiagMatrixView<T,I>(diag().Conjugate()); }

      inline DiagMatrixView<T,I> Adjoint() const
      { return Conjugate(); }

#if 0
      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> inline void operator=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator+=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator-=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator*=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator/=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator%=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
#endif

      using GenDiagMatrix<T>::size;

    protected:

      VectorView<T> itsdiag;
      inline ConstVectorView<T> cdiag() const { return itsdiag; }

  }; // DiagMatrixView

  template <class T> class DiagMatrixView<T,FortranStyle> : 
    public DiagMatrixView<T,CStyle>
  {

    public:

      //
      // Constructors
      //

      inline DiagMatrixView(const DiagMatrixView<T,FortranStyle>& rhs) :
	DiagMatrixView<T,CStyle>(rhs) {}

      inline DiagMatrixView(const DiagMatrixView<T,CStyle>& rhs) :
	DiagMatrixView<T,CStyle>(rhs) {}

      explicit inline DiagMatrixView(const VectorView<T>& _diag) :
	DiagMatrixView<T,CStyle>(_diag) {}

      virtual inline ~DiagMatrixView() {} 

      //
      // Op=
      //

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const DiagMatrixView<T,FortranStyle>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const GenDiagMatrix<RealType(T)>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const GenDiagMatrix<ComplexType(T)>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      template <class T2> 
	inline const DiagMatrixView<T,FortranStyle>& operator=(
	    const GenDiagMatrix<T2>& m2) const
	{ DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(T x) const 
      { DiagMatrixView<T,CStyle>::operator=(x); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const AssignableToDiagMatrix<RealType(T)>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const AssignableToDiagMatrix<ComplexType(T)>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      //
      // Access
      //
 
      inline RefType(T) operator()(size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return diag()(i); 
      }
      inline RefType(T) operator()(size_t i, size_t DEBUGPARAM(j)) const 
      { 
	TMVAssert(i==j); 
	TMVAssert(i>0 && i<=size());
	return diag()(i); 
      }

      inline VectorView<T,FortranStyle> diag() const 
      { return DiagMatrixView<T,CStyle>::diag(); }

      //
      // Modifying Functions
      //

      inline const DiagMatrixView<T,FortranStyle>& Zero() const 
      { diag().Zero(); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& Clip(
	  RealType(T) thresh) const
      { diag().Clip(thresh); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& SetAllTo(T x) const
      { diag().SetAllTo(x); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& TransposeSelf() const
      { return *this; }

      inline const DiagMatrixView<T,FortranStyle>& ConjugateSelf() const
      { diag().ConjugateSelf(); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& InvertSelf() const
      { return DiagMatrixView<T,CStyle>::InvertSelf(); }

      inline const DiagMatrixView<T,FortranStyle>& SetToIdentity(
	  T x=T(1)) const 
      { diag().SetAllTo(x); return *this; }

      //
      // SubDiagMatrix
      //

      inline DiagMatrixView<T,FortranStyle> SubDiagMatrix(int i1, int i2) const
      { 
	TMVAssert(diag().OKSubVector(i1,i2,1));
        return DiagMatrixView<T,FortranStyle>(diag().SubVector(i1,i2)); 
      }

      inline DiagMatrixView<T,FortranStyle> SubDiagMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(diag().OKSubVector(i1,i2,istep));
        return DiagMatrixView<T,FortranStyle>(diag().SubVector(i1,i2,istep)); 
      }

      inline DiagMatrixView<RealType(T),FortranStyle> Real() const
      { return DiagMatrixView<RealType(T),FortranStyle>(diag().Real()); }

      inline DiagMatrixView<RealType(T),FortranStyle> Imag() const
      { return DiagMatrixView<RealType(T),FortranStyle>(diag().Imag()); }

      inline DiagMatrixView<T,FortranStyle> View() const
      { return *this; }

      inline DiagMatrixView<T,FortranStyle> Transpose() const
      { return *this; }

      inline DiagMatrixView<T,FortranStyle> Conjugate() const
      { return DiagMatrixView<T,FortranStyle>(diag().Conjugate()); }

      inline DiagMatrixView<T,FortranStyle> Adjoint() const
      { return Conjugate(); }

      template <class T2> inline void operator=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }

      using GenDiagMatrix<T>::size;

  }; // FortranStyle DiagMatrixView


  template <class T, IndexStyle I> class DiagMatrix : 
    public GenDiagMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline explicit DiagMatrix(size_t _size) : itsdiag(_size) {}

      inline DiagMatrix(size_t _size, T x) : itsdiag(_size,x)  {}

      inline explicit DiagMatrix(const GenVector<T>& rhs) : itsdiag(rhs) {}

      inline explicit DiagMatrix(const GenMatrix<T>& m) : itsdiag(m.diag()) {}

      inline DiagMatrix(const DiagMatrix<T,I>& rhs) : itsdiag(rhs.diag()) {}

      template <IndexStyle I2> inline DiagMatrix(const DiagMatrix<T,I2>& rhs) :
	itsdiag(rhs.diag()) {}

      inline DiagMatrix(const GenDiagMatrix<RealType(T)>& rhs) :
	itsdiag(rhs.size()) 
      { rhs.AssignToD(View()); }

      inline DiagMatrix(const GenDiagMatrix<ComplexType(T)>& rhs) :
	itsdiag(rhs.size()) 
      { 
	TMVAssert(IsComplex(T()));
	rhs.AssignToD(View()); 
      }

      template <class T2> inline DiagMatrix(const GenDiagMatrix<T2>& rhs) :
	itsdiag(rhs.diag()) {}

      inline DiagMatrix(const AssignableToDiagMatrix<RealType(T)>& m2) :
	itsdiag(m2.colsize())
      {
	TMVAssert(m2.colsize() == m2.rowsize());
	m2.AssignToD(View()); 
      }

      inline DiagMatrix(const AssignableToDiagMatrix<ComplexType(T)>& m2) :
	itsdiag(m2.colsize())
      {
	TMVAssert(m2.colsize() == m2.rowsize());
	TMVAssert(IsComplex(T()));
	m2.AssignToD(View()); 
      }

      virtual inline ~DiagMatrix() {}


      //
      // Op=
      //

      inline DiagMatrix<T,I>& operator=(const DiagMatrix<T,I>& m2)
      {
	TMVAssert(m2.size() == size());
	m2.AssignToD(View());
	return *this; 
      }

      inline DiagMatrix<T,I>& operator=(const GenDiagMatrix<RealType(T)>& m2)
      {
	TMVAssert(m2.size() == size());
	m2.AssignToD(View()); 
	return *this; 
      }

      inline DiagMatrix<T,I>& operator=(const GenDiagMatrix<ComplexType(T)>& m2)
      {
	TMVAssert(m2.size() == size());
	TMVAssert(IsComplex(T()));
	m2.AssignToD(View()); 
	return *this; 
      }

      template <class T2> inline DiagMatrix<T,I>& operator=(
	  const GenDiagMatrix<T2>& m2)
      { 
	TMVAssert(m2.size() == size());
	View() = m2; 
	return *this; 
      }

      inline DiagMatrix<T,I>& operator=(T x) { View() = x; return *this; }

      inline DiagMatrix<T,I>& operator=(
	  const AssignableToDiagMatrix<RealType(T)>& m2)
      { 
	TMVAssert(m2.size() == size());
	m2.AssignToD(View());
	return *this; 
      }

      inline DiagMatrix<T,I>& operator=(
	  const AssignableToDiagMatrix<ComplexType(T)>& m2)
      { 
	TMVAssert(m2.size() == size());
	TMVAssert(IsComplex(T()));
	m2.AssignToD(View());
	return *this; 
      }


      //
      // Access
      //

      inline T& operator()(size_t i) 
      { 
	if (I==CStyle) { TMVAssert(i<size()); return itsdiag(i); }
	else { TMVAssert(i>0 && i<=size()); return itsdiag(i-1); }
      }

      inline T& operator()(size_t i, size_t DEBUGPARAM(j)) 
      { 
	if (I==CStyle) { 
	  TMVAssert(i<size()); 
	  TMVAssert(j<size()); 
	} else { 
	  TMVAssert(i>0 && i<=size()); 
	  TMVAssert(j>0 && j<=size()); 
	}
	TMVAssert(i==j);
	return operator()(i);
      }

      inline T operator()(size_t i) const 
      {
	if (I==CStyle) { TMVAssert(i<size()); return itsdiag(i); }
	else { TMVAssert(i>0 && i<=size()); return itsdiag(i-1); }
      }

      inline T operator()(size_t i,size_t j) const 
      {
	if (I==CStyle) { 
	  TMVAssert(i<size()); 
	  TMVAssert(j<size()); 
	} else { 
	  TMVAssert(i>0 && i<=size()); 
	  TMVAssert(j>0 && j<=size()); 
	}
	if (i==j) return operator()(i);
	else return T(0);
      }

      inline VectorView<T,I> diag() { return itsdiag.View(); }

      inline ConstVectorView<T,I> diag() const { return itsdiag.View(); }

      //
      // Modifying Functions
      //

      inline DiagMatrix<T,I>& Zero() { return SetAllTo(0); }

      inline DiagMatrix<T,I>& Clip(RealType(T) thresh)
      { diag().Clip(thresh); return *this; }

      inline DiagMatrix<T,I>& SetAllTo(T x) 
      { itsdiag.SetAllTo(x); return *this; }

      inline DiagMatrix<T,I>& TransposeSelf() 
      { return *this; }

      inline DiagMatrix<T,I>& ConjugateSelf() 
      { itsdiag.ConjugateSelf(); return *this; }

      inline DiagMatrix<T,I>& InvertSelf()
      { View().InvertSelf(); return *this; }

      inline DiagMatrix<T,I>& SetToIdentity(T x=T(1)) 
      { itsdiag.SetAllTo(x); return *this; }

      //
      // SubDiagMatrix
      //

      inline DiagMatrixView<T,I> SubDiagMatrix(int i1, int i2) 
      { return DiagMatrixView<T,I>(itsdiag.SubVector(i1,i2)); }

      inline DiagMatrixView<T,I> SubDiagMatrix(int i1, int i2, int istep) 
      { return DiagMatrixView<T,I>(itsdiag.SubVector(i1,i2,istep)); }

      inline DiagMatrixView<RealType(T),I> Real()
      { return DiagMatrixView<RealType(T),I>(diag().Real()); }

      inline DiagMatrixView<RealType(T),I> Imag()
      { return DiagMatrixView<RealType(T),I>(diag().Imag()); }

      inline ConstDiagMatrixView<T,I> SubDiagMatrix(int i1, int i2) const
      { return ConstDiagMatrixView<T,I>(itsdiag.SubVector(i1,i2)); }

      inline ConstDiagMatrixView<T,I> SubDiagMatrix(int i1, int i2,
	  int istep) const
      { return ConstDiagMatrixView<T,I>(itsdiag.SubVector(i1,i2,istep)); }

      inline ConstDiagMatrixView<RealType(T),I> Real() const
      { return ConstDiagMatrixView<RealType(T),I>(diag().Real()); }

      inline ConstDiagMatrixView<RealType(T),I> Imag() const
      { return ConstDiagMatrixView<RealType(T),I>(diag().Imag()); }

      inline ConstDiagMatrixView<T,I> View() const
      { return ConstDiagMatrixView<T,I>(itsdiag); }

      inline ConstDiagMatrixView<T,I> Transpose() const
      { return View(); }

      inline ConstDiagMatrixView<T,I> Conjugate() const
      { return ConstDiagMatrixView<T,I>(itsdiag.Conjugate()); }

      inline ConstDiagMatrixView<T,I> Adjoint() const
      { return Conjugate(); }

      inline DiagMatrixView<T,I> View() 
      { return DiagMatrixView<T,I>(itsdiag.View()); }

      inline DiagMatrixView<T,I> Transpose() 
      { return View(); }

      inline DiagMatrixView<T,I> Conjugate()
      { return DiagMatrixView<T,I>(itsdiag.Conjugate()); }

      inline DiagMatrixView<T,I> Adjoint() 
      { return Conjugate(); }

      using GenDiagMatrix<T>::size;

    protected :

      Vector<T> itsdiag;

      inline ConstVectorView<T> cdiag() const { return itsdiag.View(); }

  }; // DiagMatrix

//---------------------------------------------------------------------------

  //
  // Special Creators:
  //   DiagMatrixViewOf(m)
  //   DiagMatrixViewOf(v)
  //

  template <class T> inline ConstDiagMatrixView<T> DiagMatrixViewOf(
      const GenMatrix<T>& m)
  { return ConstDiagMatrixView<T>(m.diag()); }

  template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> DiagMatrixViewOf(
	const ConstMatrixView<T,I>& m)
    { return ConstDiagMatrixView<T,I>(m.diag()); }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> DiagMatrixViewOf(const Matrix<T,S,I>& m)
    { return ConstDiagMatrixView<T,I>(m.diag()); }

  template <class T, StorageType S, IndexStyle I> 
    inline DiagMatrixView<T,I> DiagMatrixViewOf(Matrix<T,S,I>& m)
    { return DiagMatrixView<T,I>(m.diag()); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> DiagMatrixViewOf(
      const MatrixView<T,I>& m)
  { return DiagMatrixView<T,I>(m.diag()); }

  template <class T> inline ConstDiagMatrixView<T> DiagMatrixViewOf(
      const GenVector<T>& v)
  { return ConstDiagMatrixView<T>(v); }

  template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> DiagMatrixViewOf(
	const ConstVectorView<T,I>& v)
    { return ConstDiagMatrixView<T,I>(v); }

  template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> DiagMatrixViewOf(
	const Vector<T,I>& v)
    { return ConstDiagMatrixView<T,I>(v.View()); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> DiagMatrixViewOf(
      const VectorView<T,I>& v)
  { return DiagMatrixView<T,I>(v); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> DiagMatrixViewOf(
      Vector<T,I>& v)
  { return DiagMatrixView<T,I>(v.View()); }

  template <class T> inline ConstDiagMatrixView<T> DiagMatrixViewOf(
      const T* v, size_t size)
  { return ConstDiagMatrixView<T>(VectorViewOf(v,size)); }

  template <class T> inline DiagMatrixView<T> DiagMatrixViewOf(
      T* v, size_t size)
  { return DiagMatrixView<T>(VectorViewOf(v,size)); }

  //
  // Swap Matrices
  //

  template <class T> inline void Swap(
      const DiagMatrixView<T>& m1, const DiagMatrixView<T>& m2)
  { Swap(m1.diag(),m2.diag()); }

  template <class T, IndexStyle I> inline void Swap(
      const DiagMatrix<T,I>& m1, const DiagMatrixView<T>& m2)
  { Swap(m1.diag(),m2.diag()); }

  template <class T, IndexStyle I> inline void Swap(
      const DiagMatrixView<T>& m1, const DiagMatrix<T,I>& m2)
  { Swap(m1.diag(),m2.diag()); }

  template <class T, IndexStyle I1, IndexStyle I2> inline void Swap(
      const DiagMatrix<T,I1>& m1, const DiagMatrix<T,I2>& m2)
  { Swap(m1.diag(),m2.diag()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline ConstDiagMatrixView<T> Transpose(
      const GenDiagMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> inline ConstDiagMatrixView<T,I> Transpose(
      const ConstDiagMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> inline ConstDiagMatrixView<T,I> Transpose(
      const DiagMatrix<T,I>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> Transpose(
      const DiagMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> Transpose(
      DiagMatrix<T,I>& m)
  { return m.Transpose(); }

  template <class T> inline ConstDiagMatrixView<T> Conjugate(
      const GenDiagMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> inline ConstDiagMatrixView<T,I> Conjugate(
      const ConstDiagMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> inline ConstDiagMatrixView<T,I> Conjugate(
      const DiagMatrix<T,I>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> Conjugate(
      const DiagMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> Conjugate(
      DiagMatrix<T,I>& m)
  { return m.Conjugate(); }

  template <class T> inline ConstDiagMatrixView<T> Adjoint(
      const GenDiagMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline ConstDiagMatrixView<T,I> Adjoint(
      const ConstDiagMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline ConstDiagMatrixView<T,I> Adjoint(
      const DiagMatrix<T,I>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> Adjoint(
      const DiagMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline DiagMatrixView<T,I> Adjoint(
      DiagMatrix<T,I>& m)
  { return m.Adjoint(); }

  template <class T> inline QuotXD<T,T> Inverse(const GenDiagMatrix<T>& m)
  { return m.Inverse(); }

  //
  // DiagMatrix ==, != DiagMatrix
  //

  template <class T1, class T2> inline bool operator==(
      const GenDiagMatrix<T1>& m1, const GenDiagMatrix<T2>& m2)
  { return m1.diag() == m2.diag(); }

  template <class T1, class T2> inline bool operator!=(
      const GenDiagMatrix<T1>& m1, const GenDiagMatrix<T2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, IndexStyle I> std::istream& operator>>(std::istream& fin, 
      auto_ptr<DiagMatrix<T,I> >& m);

  template <class T> std::istream& operator>>(std::istream& fin,
      const DiagMatrixView<T>& m);

  template <class T, IndexStyle I> inline std::istream& operator>>(
      std::istream& fin, DiagMatrix<T,I>& m)
  { return fin >> m.View(); }

} // namespace tmv

#endif
