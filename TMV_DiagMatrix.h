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
//    DiagMatrix<T>(size_t colsize, size_t rowsize)
//        colsize must = rowsize
//        Included for easier conversion of code written for Matrix.
//
//    DiagMatrix<T>(size_t size, T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    DiagMatrix<T>(size_t size, const Vector<T>& vv)
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
//    bool IsSquare()
//        (Always true.)
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
//        (Note - for diagonal matrices, Norm1 = Norm2 = NormInf.)
//    Transpose(m)
//        (same as the original).
//    Conjugate(m)
//    Adjoint(m)
//    Inverse(m)
//    InverseATA(m)
//    DInverse(m)   (This and DInverseATA return DiagMatrices.)
//    DInverseATA(m)
//
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

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenDiagMatrix;
  template <class T> class ConstDiagMatrixView;
  template <class T> class DiagMatrixView;
  template <class T> class DiagMatrix;
  template <class T> class DiagMatrixComposite;

  template <class T> class DiagDiv;

}

#include "TMV_DiagDivider.h"

namespace tmv {

  template <class T> class GenDiagMatrix : public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      GenDiagMatrix<T>() : BaseMatrix<T>(LU) {}
      GenDiagMatrix<T>(const MetaDiagDivider<T>* mdiv) : 
	BaseMatrix<T>(mdiv) {}
      GenDiagMatrix(const GenDiagMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs) {}
      ~GenDiagMatrix() {}

      //
      // Access Functions
      //

      inline virtual size_t size() const { return diag().size(); }
      inline virtual size_t colsize() const { return size(); }
      inline virtual size_t rowsize() const { return size(); }

      inline ConstVectorView<T> diag() const { return cdiag(); }

      using BaseMatrix<T>::operator();
      inline T operator()(size_t i) const { return diag()(i); }

      inline bool IsSquare() const { return true; }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (diag().cptr() == m2.cptr()); }

      inline bool SameStorageAs(const GenDiagMatrix<T>& m2) const
      { return (diag().cptr() == m2.diag().cptr()); }

      template <class T2> inline bool SameAs(
	  const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameAs(const GenDiagMatrix<T>& m2) const
      { 
	if (this == &m2) return true;
	else return (diag().SameAs(m2.diag())); 
      }

      inline void CopyToMatrix(const MatrixView<RealType(T)>& m) const
      { 
	TMVAssert(IsReal(T())); 
	TMVAssert(m.colsize() == size());
	TMVAssert(m.rowsize() == size());
	m.Zero();
	m.diag() = diag();
      }
      inline void CopyToMatrix(const MatrixView<ComplexType(T)>& m) const
      {
	TMVAssert(m.colsize() == size());
	TMVAssert(m.rowsize() == size());
	m.Zero();
	m.diag() = diag();
      }

      //
      // SubDiagMatrix
      //

      inline ConstDiagMatrixView<T> SubDiagMatrix(int i1, int i2) const
      { return ConstDiagMatrixView<T>(diag().SubVector(i1,i2)); }

      inline ConstDiagMatrixView<T> SubDiagMatrix(int i1, int i2,
	  int istep) const
      { return ConstDiagMatrixView<T>(diag().SubVector(i1,i2,istep)); }

      inline ConstDiagMatrixView<RealType(T)> Real() const
      { return ConstDiagMatrixView<RealType(T)>(diag().Real()); }

      inline ConstDiagMatrixView<RealType(T)> Imag() const
      { return ConstDiagMatrixView<RealType(T)>(diag().Imag()); }

      inline ConstDiagMatrixView<T> View() const
      { 
	return ConstDiagMatrixView<T>(diag(),
	    MakeMetaDiagDivider(false,false,this));
      }

      inline ConstDiagMatrixView<T> Transpose() const
      { return View(); }

      inline ConstDiagMatrixView<T> Conjugate() const
      { 
	return ConstDiagMatrixView<T>(diag().Conjugate(),
	    MakeMetaDiagDivider(false,true,this));
      }

      inline ConstDiagMatrixView<T> Adjoint() const
      { return Conjugate(); }

      inline ConstDiagMatrixView<T> QuickView() const
      { 
	return ConstDiagMatrixView<T>(diag());
      }

      inline ConstDiagMatrixView<T> QuickTranspose() const
      { return QuickView(); }

      inline ConstDiagMatrixView<T> QuickConjugate() const
      { 
	return ConstDiagMatrixView<T>(diag().Conjugate());
      }

      inline ConstDiagMatrixView<T> QuickAdjoint() const
      { return QuickConjugate(); }


      //
      // Functions of DiagMatrix
      //

      inline T Trace() const
      { return diag().SumElements(); }

      inline DiagMatrix<T> DInverse() const
      { 
	this->SetDiv(); 
	const DiagDivider<T>* ddiv =
	  dynamic_cast<const DiagDivider<T>*>(this->GetDiv());
	TMVAssert(ddiv);
	return ddiv->DInverse(); 
      }

      inline DiagMatrix<T> DInverseATA() const
      {
	this->SetDiv(); 
	const DiagDivider<T>* ddiv =
	  dynamic_cast<const DiagDivider<T>*>(this->GetDiv());
	TMVAssert(ddiv);
	return ddiv->DInverseATA(); 
      }

      // NormF()^2
      inline RealType(T) NormSq() const 
      { return diag().NormSq(); }

      // 1-Norm = max_j (sum_i |a_ij|)
      inline RealType(T) Norm1() const 
      { return diag().MaxAbsElement(); }

      // 2-Norm = sqrt(max_eigenvalue(A At))
      inline RealType(T) Norm2() const 
      { return diag().MaxAbsElement(); }

      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return diag().MaxAbsElement(); }

      inline BaseMatrix<T>* NewTranspose() const
      { return new ConstDiagMatrixView<T>(Transpose()); }

      inline BaseMatrix<T>* NewConjugate() const
      { return new ConstDiagMatrixView<T>(Conjugate()); }

      inline BaseMatrix<T>* NewAdjoint() const
      { return new ConstDiagMatrixView<T>(Adjoint()); }

      inline BaseMatrix<T>* NewView() const
      { return new ConstDiagMatrixView<T>(View()); }

      inline BaseMatrix<T>* NewCopy() const
      { return new DiagMatrix<T>(*this); }

      // 
      // I/O
      //

      void WriteCompact(ostream& os) const
      { os << diag(); }

      // 
      // Arithmetic Helpers
      //

      using BaseMatrix<T>::LDivEq;
      using BaseMatrix<T>::RDivEq;
      using BaseMatrix<T>::LDiv;
      using BaseMatrix<T>::RDiv;

      template <class T1> inline void DivEq(const DiagMatrixView<T1>& m) const
      {
	this->SetDiv();
	const DiagDivider<T>* ddiv =
	  dynamic_cast<const DiagDivider<T>*>(this->GetDiv());
	TMVAssert(ddiv);
	ddiv->DivEq(m);
      }

      template <class T1, class T0> inline void Div(
	  const GenDiagMatrix<T1>& m1, const DiagMatrixView<T0>& m0) const
      {
	this->SetDiv();
	const DiagDivider<T>* ddiv =
	  dynamic_cast<const DiagDivider<T>*>(this->GetDiv());
	TMVAssert(ddiv);
	ddiv->Div(m1,m0);
      }

      //
      // Division Control 
      //

      inline void DivideUsing(DivType dt) const
      { TMVAssert(dt == LU); }

    protected :

      virtual ConstVectorView<T> cdiag() const = 0;
      inline T cref(size_t i, size_t j) const
      { return i==j ? cdiag()(i) : 0; }

      void NewDivider() const;

    private :

      void operator=(const GenDiagMatrix<T>&) { TMVAssert(false); }

  }; // GenDiagMatrix

  template <class T> class ConstDiagMatrixView : 
    public GenDiagMatrix<T>
  {
    public :

      ConstDiagMatrixView(const ConstDiagMatrixView<T>& rhs) :
	GenDiagMatrix<T>(rhs), itsdiag(rhs.cdiag()) {}

      ConstDiagMatrixView(const GenDiagMatrix<T>& rhs) :
	GenDiagMatrix<T>(rhs), itsdiag(rhs.cdiag()) {}

      explicit ConstDiagMatrixView(const GenVector<T>& v) :
	itsdiag(v) {}

      ConstDiagMatrixView(const ConstVectorView<T>& _diag, 
	  const MetaDiagDivider<T>* mdiv) :
	GenDiagMatrix<T>(mdiv), itsdiag(_diag) {}

      ~ConstDiagMatrixView() {}

    protected :

      ConstVectorView<T> itsdiag;

      inline ConstVectorView<T> cdiag() const { return itsdiag; }

    private :

      void operator=(const ConstDiagMatrixView<T>&) { TMVAssert(false); }

  }; // ConstDiagMatrixView

  template <class T> class DiagMatrixComposite;
  // Defined in TMV_DiagMatrixArtith.h

  template <class T> class DiagMatrixView : 
    public GenDiagMatrix<T>
  {

    public:

      //
      // Constructors
      //

      DiagMatrixView(const DiagMatrixView<T>& rhs) :
	GenDiagMatrix<T>(rhs), itsdiag(rhs.diag()) {}

      explicit DiagMatrixView(const VectorView<T>& _diag) :
	itsdiag(_diag) {}

      explicit DiagMatrixView(Vector<T>& _diag) : itsdiag(_diag.View()) {}

      DiagMatrixView(const VectorView<T>& _diag, 
	  const MetaDiagDivider<T>* mdiv) :
	GenDiagMatrix<T>(mdiv), itsdiag(_diag) {}

      ~DiagMatrixView() {} 

      //
      // Op=
      //

      inline const DiagMatrixView<T>& operator=(
	  const DiagMatrixView<T>& m2) const
      { itsdiag = m2.diag(); return *this; }

      inline const DiagMatrixView<T>& operator=(
	  const GenDiagMatrix<T>& m2) const
      { itsdiag = m2.diag(); return *this; }

      template <class T2> inline const DiagMatrixView<T>& operator=(
	  const GenDiagMatrix<T2>& m2) const
      { itsdiag = m2.diag(); return *this; }

      inline const DiagMatrixView<T>& operator=(T x) const 
      { return SetToIdentity(x); }

      inline const DiagMatrixView<T>& operator=(
	  const DiagMatrixComposite<T>& dcomp) const
      { 
	TMVAssert(this->size() == dcomp.size());
	dcomp.AssignTo(*this);
	return *this;
      }


      //
      // Access
      //
 
      inline VectorView<T> diag() const { return itsdiag; }

      inline RefType(T) operator()(size_t i) const { return diag()(i); }
      inline RefType(T) operator()(size_t i, size_t j) const 
      { TMVAssert(i==j); return diag()(i); }

      //
      // Modifying Functions
      //

      inline const DiagMatrixView<T>& Zero() const { return SetAllTo(0); }

      inline const DiagMatrixView<T>& SetAllTo(T x) const
      { diag().SetAllTo(x); return *this; }

      inline const DiagMatrixView<T>& TransposeSelf() const
      { return *this; }

      inline const DiagMatrixView<T>& ConjugateSelf() const
      { diag().ConjugateSelf(); return *this; }

      inline const DiagMatrixView<T>& SetToIdentity(T x=T(1)) const 
      { return SetAllTo(x); }

      //
      // SubDiagMatrix
      //

      inline DiagMatrixView<T> SubDiagMatrix(int i1, int i2) const
      { return DiagMatrixView<T>(itsdiag.SubVector(i1,i2)); }

      inline DiagMatrixView<T> SubDiagMatrix(int i1, int i2,
	  int istep) const
      { return DiagMatrixView<T>(itsdiag.SubVector(i1,i2,istep)); }

      inline DiagMatrixView<RealType(T)> Real() const
      { return DiagMatrixView<RealType(T)>(diag().Real()); }

      inline DiagMatrixView<RealType(T)> Imag() const
      { return DiagMatrixView<RealType(T)>(diag().Imag()); }

      inline DiagMatrixView<T> View() const
      { return *this; }

      inline DiagMatrixView<T> Transpose() const
      { return *this; }

      inline DiagMatrixView<T> Conjugate() const
      { 
	return DiagMatrixView<T>(diag().Conjugate(),
	    MakeMetaDiagDivider(false,true,this));
      }

      inline DiagMatrixView<T> Adjoint() const
      { return Conjugate(); }

      inline DiagMatrixView<T> QuickView() const
      { return DiagMatrixView<T>(diag()); }

      inline DiagMatrixView<T> QuickTranspose() const
      { return QuickView(); }

      inline DiagMatrixView<T> QuickConjugate() const
      { 
	return DiagMatrixView<T>(diag().Conjugate());
      }

      inline DiagMatrixView<T> QuickAdjoint() const
      { return QuickConjugate(); }

      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> void operator=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator+=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator-=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator*=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator/=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator%=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }

    protected:

      VectorView<T> itsdiag;
      inline ConstVectorView<T> cdiag() const { return itsdiag; }

  }; // DiagMatrixView


  template <class T> class DiagMatrix : 
    public GenDiagMatrix<T>
  {

    public:

      //
      // Constructors
      //

      explicit DiagMatrix(size_t _size) : itsdiag(_size) {}

      DiagMatrix(size_t colsize, size_t rowsize) : itsdiag(colsize) 
      { TMVAssert(colsize == rowsize); }

      DiagMatrix(size_t _size, T x) : itsdiag(_size,x)  {}

      explicit DiagMatrix(const GenVector<T>& rhs) : itsdiag(rhs) {}

      explicit DiagMatrix(const GenMatrix<T>& m) : itsdiag(m.diag()) {}

      DiagMatrix(const DiagMatrix<T>& rhs) : itsdiag(rhs.diag()) {}

      template <class T2> DiagMatrix(const GenDiagMatrix<T2>& rhs) :
	itsdiag(rhs.diag()) {}

      DiagMatrix(const DiagMatrixComposite<T>& dcomp) : itsdiag(dcomp.size())
      { dcomp.AssignTo(QuickView()); }

      ~DiagMatrix() {}


      //
      // Op=
      //

      inline DiagMatrix<T>& operator=(const DiagMatrix<T>& m2)
      {
	if (&m2 != this) QuickView() = m2.QuickView(); 
	return *this; 
      }

      template <class T2> inline DiagMatrix<T>& operator=(
	  const GenDiagMatrix<T2>& m2)
      { QuickView() = m2; return *this; }

      inline DiagMatrix<T>& operator=(T x) { QuickView() = x; return *this; }

      inline DiagMatrix<T>& operator=(const DiagMatrixComposite<T>& dcomp)
      { QuickView() = dcomp; return *this; }


      //
      // Access
      //

      inline VectorView<T> diag() { return itsdiag.View(); }

      inline T& operator()(size_t i) { return itsdiag(i); }
      inline T& operator()(size_t i, size_t j) { return itsdiag(i); }

      inline ConstVectorView<T> diag() const { return itsdiag.View(); }

      inline T operator()(size_t i) const { return itsdiag(i); }
      inline T operator()(size_t i,size_t j) const 
      { if (i==j) return itsdiag(i); else return T(0); }


      //
      // Modifying Functions
      //

      inline DiagMatrix<T>& Zero() { return SetAllTo(0); }

      inline DiagMatrix<T>& SetAllTo(T x) 
      { itsdiag.SetAllTo(x); return *this; }

      inline DiagMatrix<T>& TransposeSelf() 
      { return *this; }

      inline DiagMatrix<T>& ConjugateSelf() 
      { itsdiag.ConjugateSelf(); return *this; }

      inline DiagMatrix<T>& SetToIdentity(T x=T(1)) 
      { itsdiag.SetAllTo(x); return *this; }

      //
      // SubDiagMatrix
      //

      inline DiagMatrixView<T> SubDiagMatrix(int i1, int i2) 
      { return DiagMatrixView<T>(itsdiag.SubVector(i1,i2)); }

      inline DiagMatrixView<T> SubDiagMatrix(int i1, int i2, int istep) 
      { return DiagMatrixView<T>(itsdiag.SubVector(i1,i2,istep)); }

      inline DiagMatrixView<RealType(T)> Real()
      { return DiagMatrixView<RealType(T)>(diag().Real()); }

      inline DiagMatrixView<RealType(T)> Imag()
      { return DiagMatrixView<RealType(T)>(diag().Imag()); }

      inline ConstDiagMatrixView<T> SubDiagMatrix(int i1, int i2) const
      { return ConstDiagMatrixView<T>(itsdiag.SubVector(i1,i2)); }

      inline ConstDiagMatrixView<T> SubDiagMatrix(int i1, int i2,
	  int istep) const
      { return ConstDiagMatrixView<T>(itsdiag.SubVector(i1,i2,istep)); }

      inline ConstDiagMatrixView<RealType(T)> Real() const
      { return ConstDiagMatrixView<RealType(T)>(diag().Real()); }

      inline ConstDiagMatrixView<RealType(T)> Imag() const
      { return ConstDiagMatrixView<RealType(T)>(diag().Imag()); }

      inline ConstDiagMatrixView<T> View() const
      { 
	return ConstDiagMatrixView<T>(itsdiag.View(),
	    new MetaDiagDivider<T>(false,false,this));
      }

      inline ConstDiagMatrixView<T> Transpose() const
      { return View(); }

      inline ConstDiagMatrixView<T> Conjugate() const
      { 
	return ConstDiagMatrixView<T>(
	    itsdiag.Conjugate(), new MetaDiagDivider<T>(false,true,this));
      }

      inline ConstDiagMatrixView<T> Adjoint() const
      { return Conjugate(); }

      inline DiagMatrixView<T> View() 
      {
	return DiagMatrixView<T>(itsdiag.View(),
	    new MetaDiagDivider<T>(false,false,this));
      }

      inline DiagMatrixView<T> Transpose() 
      { return View(); }

      inline DiagMatrixView<T> Conjugate()
      { 
	return DiagMatrixView<T>(itsdiag.Conjugate(),
	    new MetaDiagDivider<T>(false,true,this));
      }

      inline DiagMatrixView<T> Adjoint() 
      { return Conjugate(); }

      inline ConstDiagMatrixView<T> QuickView() const
      { return ConstDiagMatrixView<T>(itsdiag); }

      inline ConstDiagMatrixView<T> QuickTranspose() const
      { return QuickView(); }

      inline ConstDiagMatrixView<T> QuickConjugate() const
      { return ConstDiagMatrixView<T>(itsdiag.Conjugate()); }

      inline ConstDiagMatrixView<T> QuickAdjoint() const
      { return QuickConjugate(); }

      inline DiagMatrixView<T> QuickView() 
      { return DiagMatrixView<T>(itsdiag); }

      inline DiagMatrixView<T> QuickTranspose() 
      { return QuickView(); }

      inline DiagMatrixView<T> QuickConjugate()
      { return DiagMatrixView<T>(itsdiag.Conjugate()); }

      inline DiagMatrixView<T> QuickAdjoint() 
      { return QuickConjugate(); }

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

  template <class T> inline DiagMatrixView<T> DiagMatrixViewOf(
      const MatrixView<T>& m)
  { return DiagMatrixView<T>(m.diag()); }

  template <class T, StorageType S> inline DiagMatrixView<T> DiagMatrixViewOf(
      Matrix<T,S>& m)
  { return DiagMatrixView<T>(m.diag()); }

  template <class T> inline ConstDiagMatrixView<T> DiagMatrixViewOf(
      const GenVector<T>& v)
  { return ConstDiagMatrixView<T>(v); }

  template <class T> inline DiagMatrixView<T> DiagMatrixViewOf(
      const VectorView<T>& v)
  { return DiagMatrixView<T>(v); }

  template <class T> inline DiagMatrixView<T> DiagMatrixViewOf(
      Vector<T>& v)
  { return DiagMatrixView<T>(v.View()); }


  //
  // Swap Matrices
  //

  template <class T> inline void Swap(
      const DiagMatrixView<T>& m1, const DiagMatrixView<T>& m2)
  { Swap(m1.diag(),m2.diag()); }

  template <class T> inline void Swap(
      const DiagMatrix<T>& m1, const DiagMatrixView<T>& m2)
  { Swap(m1.diag(),m2.diag()); }

  template <class T> inline void Swap(
      const DiagMatrixView<T>& m1, const DiagMatrix<T>& m2)
  { Swap(m1.diag(),m2.diag()); }

  template <class T> inline void Swap(
      const DiagMatrix<T>& m1, const DiagMatrix<T>& m2)
  { Swap(m1.diag(),m2.diag()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline ConstDiagMatrixView<T> Transpose(
      const GenDiagMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline DiagMatrixView<T> Transpose(
      const DiagMatrixView<T>& m)
  { return m.Transpose(); }

  template <class T> inline DiagMatrixView<T> Transpose(DiagMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline ConstDiagMatrixView<T> Conjugate(
      const GenDiagMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T> inline DiagMatrixView<T> Conjugate(
      const DiagMatrixView<T>& m)
  { return m.Conjugate(); }

  template <class T> inline DiagMatrixView<T> Conjugate(DiagMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T> inline ConstDiagMatrixView<T> Adjoint(
      const GenDiagMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T> inline DiagMatrixView<T> Adjoint(
      const DiagMatrixView<T>& m)
  { return m.Adjoint(); }

  template <class T> inline DiagMatrixView<T> Adjoint(DiagMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T> inline DiagMatrix<T> Inverse(const GenDiagMatrix<T>& m)
  { return m.DInverse(); }


  //
  // DiagMatrix ==, != DiagMatrix
  //

  template <class T> inline bool operator==(
      const GenDiagMatrix<T>& m1, const GenDiagMatrix<T>& m2)
  { return m1.diag() == m2.diag(); }

  template <class T> inline bool operator!=(
      const GenDiagMatrix<T>& m1, const GenDiagMatrix<T>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T> istream& operator>>(istream& fin, DiagMatrix<T>* m);

  template <class T> istream& operator>>(istream& fin,
      const DiagMatrixView<T>& m);

  template <class T> inline istream& operator>>(istream& fin, DiagMatrix<T>& m)
  { return fin >> m.View(); }

  template <class T> inline std::string Type(const DiagMatrix<T>& m)
  { return std::string("DiagMatrix<")+Type(T())+">"; }
  template <class T> inline std::string Type(const GenDiagMatrix<T>& m)
  { return std::string("GenDiagMatrix<")+Type(T())+","+Text(m.diag().ct())+">"; }
  template <class T> inline std::string Type(const ConstDiagMatrixView<T>& m)
  { return std::string("ConstDiagMatrixView<")+Type(T())+","+Text(m.diag().ct())+">"; }
  template <class T> inline std::string Type(const DiagMatrixView<T>& m)
  { return std::string("DiagMatrixView<")+Type(T())+","+Text(m.diag().ct())+">"; }

}; // namespace tmv

#endif
