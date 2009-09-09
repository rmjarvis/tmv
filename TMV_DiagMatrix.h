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
//
//    m.Inverse()
//    Inverse(m)
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

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenDiagMatrix;
  template <class T, IndexStyle I=CStyle> class ConstDiagMatrixView;
  template <class T, IndexStyle I=CStyle> class DiagMatrixView;
  template <class T, IndexStyle I=CStyle> class DiagMatrix;
  template <class T> class DiagMatrixComposite;

  template <class T> class DiagDiv;
  template <class T, class Tm> class QuotXD;

  template <class T> class DiagMatrixReadError :
    public ReadError
  {
    public :
      size_t i;
      mutable auto_ptr<DiagMatrix<T> > m;
      char exp,got;
      size_t s;
      bool is, iseof, isbad;

      inline DiagMatrixReadError(istream& _is) throw() :
	ReadError("DiagMatrix"),
	i(0), m(0), exp(0), got(0), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline DiagMatrixReadError(size_t _i, const GenDiagMatrix<T>& _m,
	  char _e, char _g, size_t _s,
	  bool _is, bool _iseof, bool _isbad) throw() :
	ReadError("DiagMatrix"),
	i(_i), m(new DiagMatrix<T>(_m)), exp(_e), got(_g), s(_s),
	is(_is), iseof(_iseof), isbad(_isbad) {}
      inline DiagMatrixReadError(const GenDiagMatrix<T>& _m,
	  istream& _is, size_t _s) throw() :
	ReadError("DiagMatrix"),
	i(0), m(new DiagMatrix<T>(_m)), s(_s),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline DiagMatrixReadError(istream& _is, char _e, char _g) throw() :
	ReadError("DiagMatrix"),
	i(0), m(0), exp(_e), got(_g), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline DiagMatrixReadError(const DiagMatrixReadError<T>& rhs) :
	i(rhs.i), m(rhs.m), exp(rhs.exp), got(rhs.got), s(rhs.s),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~DiagMatrixReadError() throw() {}

      virtual void Write(ostream& os) const throw();
  };


  template <class T> class GenDiagMatrix : public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline GenDiagMatrix<T>() : BaseMatrix<T>(LU) {}
      inline GenDiagMatrix(const GenDiagMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs) {}
      inline ~GenDiagMatrix() {}

      //
      // Access Functions
      //

      inline virtual size_t size() const { return diag().size(); }
      inline virtual size_t colsize() const { return size(); }
      inline virtual size_t rowsize() const { return size(); }

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

      inline bool IsSquare() const { return true; }

      inline ConstVectorView<T> diag() const { return cdiag(); }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (diag().cptr() == m2.cptr()); }

      inline bool SameStorageAs(const GenDiagMatrix<T>& m2) const
      { return (diag().cptr() == m2.diag().cptr()); }

      template <class T2> inline bool SameAs(const BaseMatrix<T2>& ) const
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

      inline T Trace() const
      { return diag().SumElements(); }

      inline QuotXD<T,T> Inverse() const
      { return QuotXD<T,T>(T(1),*this); }

      void Inverse(const DiagMatrixView<T>& minv) const;

      void InverseATA(const DiagMatrixView<T>& minv) const;

      inline void Inverse(const MatrixView<T>& minv) const
      {
	minv.Zero();
	Inverse(DiagMatrixViewOf(minv.diag()));
      }

      inline void InverseATA(const MatrixView<T>& minv) const
      { 
	minv.Zero();
	InverseATA(DiagMatrixViewOf(minv.diag()));
      }

      template <IndexStyle I> inline void Inverse(DiagMatrix<T,I>& minv) const
      { Inverse(minv.View()); }

      template <IndexStyle I> inline void InverseATA(
	  DiagMatrix<T,I>& minv) const
      { InverseATA(minv.View()); }

      template <StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T,S,I>& minv) const
      { Inverse(minv.View()); }

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& minv) const
      { InverseATA(minv.View()); }

      // NormF()^2
      inline RealType(T) NormSq() const 
      { return diag().NormSq(); }

      // 1-Norm = max_j (sum_i |a_ij|)
      inline RealType(T) Norm1() const 
      { return diag().MaxAbsElement(); }

      // 2-Norm = sqrt(max_eigenvalue(A At))
      inline RealType(T) Norm2() const 
      { return diag().MaxAbsElement(); }

      inline RealType(T) Condition() const 
      { return diag().MaxAbsElement()/diag().MinAbsElement(); }

      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return diag().MaxAbsElement(); }

      inline auto_ptr<BaseMatrix<T> > NewTranspose() const
      {
	auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(Transpose())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewConjugate() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(Conjugate())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewAdjoint() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(Adjoint())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewView() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(View())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewCopy() const
      {
	auto_ptr<BaseMatrix<T> > a(new DiagMatrix<T>(*this)); 
	return a;
      }

      // 
      // I/O
      //

      inline void WriteCompact(ostream& os) const
      { os << "D " << diag(); }
      inline void WriteCompact(ostream& os, RealType(T) thresh) const
      { os << "D "; diag().Write(os,thresh); }

      // 
      // Arithmetic Helpers
      //

      using BaseMatrix<T>::LDivEq;
      using BaseMatrix<T>::RDivEq;
      using BaseMatrix<T>::LDiv;
      using BaseMatrix<T>::RDiv;

      template <class T1> void DivEq(const DiagMatrixView<T1>& m) const;
      template <class T1, class T0> void Div(
	  const GenDiagMatrix<T1>& m1, const DiagMatrixView<T0>& m0) const;

      //
      // Division Control 
      //

      inline void DivideUsing(DivType 
#ifdef TMVDEBUG
	  dt
#endif
	  ) const
      { TMVAssert(dt == LU); }

    protected :

      virtual ConstVectorView<T> cdiag() const = 0;
      inline T cref(size_t i, size_t j) const
      { return i==j ? cdiag()(i) : 0; }

      void NewDivider() const;

    private :

      inline void operator=(const GenDiagMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenDiagMatrix

  template <class T, IndexStyle I> class ConstDiagMatrixView : 
    public GenDiagMatrix<T>
  {
    public :

      inline ConstDiagMatrixView(const ConstDiagMatrixView<T,I>& rhs) :
	GenDiagMatrix<T>(rhs), itsdiag(rhs.cdiag()) {}

      inline ConstDiagMatrixView(const GenDiagMatrix<T>& rhs) :
	GenDiagMatrix<T>(rhs), itsdiag(rhs.diag()) {}

      explicit inline ConstDiagMatrixView(const GenVector<T>& v) :
	itsdiag(v) {}

      inline ~ConstDiagMatrixView() {}

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

      inline ~ConstDiagMatrixView() {}

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

  template <class T> class DiagMatrixComposite;
  // Defined in TMV_DiagMatrixArtith.h

  template <class T, IndexStyle I> class DiagMatrixView : 
    public GenDiagMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline DiagMatrixView(const DiagMatrixView<T,I>& rhs) :
	GenDiagMatrix<T>(rhs), itsdiag(rhs.diag()) {}

      explicit inline DiagMatrixView(const VectorView<T>& _diag) :
	itsdiag(_diag) {}

      inline ~DiagMatrixView() {} 

      //
      // Op=
      //

      inline const DiagMatrixView<T,I>& operator=(
	  const DiagMatrixView<T,I>& m2) const
      { itsdiag = m2.diag(); return *this; }

      inline const DiagMatrixView<T,I>& operator=(
	  const GenDiagMatrix<T>& m2) const
      { itsdiag = m2.diag(); return *this; }

      template <class T2> inline const DiagMatrixView<T,I>& operator=(
	  const GenDiagMatrix<T2>& m2) const
      { itsdiag = m2.diag(); return *this; }

      inline const DiagMatrixView<T,I>& operator=(T x) const 
      { return SetToIdentity(x); }

      inline const DiagMatrixView<T,I>& operator=(
	  const DiagMatrixComposite<T>& mcomp) const
      { 
	TMVAssert(size() == mcomp.size());
	mcomp.AssignTo(*this);
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
      inline RefType(T) operator()(size_t i, size_t 
#ifdef TMVDEBUG
	  j
#endif
	  ) const 
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

      inline ~DiagMatrixView() {} 

      //
      // Op=
      //

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const DiagMatrixView<T,FortranStyle>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const GenDiagMatrix<T>& m2) const
      { DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      template <class T2> 
	inline const DiagMatrixView<T,FortranStyle>& operator=(
	    const GenDiagMatrix<T2>& m2) const
	{ DiagMatrixView<T,CStyle>::operator=(m2); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(T x) const 
      { DiagMatrixView<T,CStyle>::operator=(x); return *this; }

      inline const DiagMatrixView<T,FortranStyle>& operator=(
	  const DiagMatrixComposite<T>& mcomp) const
      { DiagMatrixView<T,CStyle>::operator=(mcomp); return *this; }

      //
      // Access
      //
 
      inline RefType(T) operator()(size_t i) const 
      { 
	TMVAssert(i>0 && i<=size());
	return diag()(i); 
      }
      inline RefType(T) operator()(size_t i, size_t 
#ifdef TMVDEBUG
	  j
#endif
	  ) const 
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

      explicit inline DiagMatrix(size_t _size) : itsdiag(_size) {}

      inline DiagMatrix(size_t _size, T x) : itsdiag(_size,x)  {}

      explicit inline DiagMatrix(const GenVector<T>& rhs) : itsdiag(rhs) {}

      explicit inline DiagMatrix(const GenMatrix<T>& m) : itsdiag(m.diag()) {}

      inline DiagMatrix(const DiagMatrix<T,I>& rhs) : itsdiag(rhs.diag()) {}

      inline DiagMatrix(const GenDiagMatrix<T>& rhs) : itsdiag(rhs.diag()) {}

      template <class T2> inline DiagMatrix(const GenDiagMatrix<T2>& rhs) :
	itsdiag(rhs.diag()) {}

      inline DiagMatrix(const DiagMatrixComposite<T>& mcomp) :
	itsdiag(mcomp.size())
      { mcomp.AssignTo(View()); }

      inline ~DiagMatrix() {}


      //
      // Op=
      //

      inline DiagMatrix<T,I>& operator=(const DiagMatrix<T,I>& m2)
      {
	if (&m2 != this) View() = m2.View(); 
	return *this; 
      }

      inline DiagMatrix<T,I>& operator=(const GenDiagMatrix<T>& m2)
      { View() = m2; return *this; }

      template <class T2> inline DiagMatrix<T,I>& operator=(
	  const GenDiagMatrix<T2>& m2)
      { View() = m2; return *this; }

      inline DiagMatrix<T,I>& operator=(T x) { View() = x; return *this; }

      inline DiagMatrix<T,I>& operator=(const DiagMatrixComposite<T>& mcomp)
      { View() = mcomp; return *this; }


      //
      // Access
      //

      inline T& operator()(size_t i) 
      { 
	if (I==CStyle) { TMVAssert(i<size()); return itsdiag(i); }
	else { TMVAssert(i>0 && i<=size()); return itsdiag(i-1); }
      }

      inline T& operator()(size_t i, size_t 
#ifdef TMVDEBUG
	  j
#endif
	  ) 
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
 
  template <class T, IndexStyle I> istream& operator>>(istream& fin, 
      auto_ptr<DiagMatrix<T,I> >& m);

  template <class T> istream& operator>>(istream& fin,
      const DiagMatrixView<T>& m);

  template <class T, IndexStyle I> inline istream& operator>>(istream& fin, 
      DiagMatrix<T,I>& m)
  { return fin >> m.View(); }

  template <class T, IndexStyle I> inline string Type(const DiagMatrix<T,I>& )
  { return string("DiagMatrix<")+Type(T())+","+Text(I)+">"; }
  template <class T> inline string Type(const GenDiagMatrix<T>& m)
  {
    return string("GenDiagMatrix<")+Type(T())+","+Text(m.diag().ct())+">"; 
  }
  template <class T, IndexStyle I> inline string Type(
      const ConstDiagMatrixView<T,I>& m)
  { 
    return string("ConstDiagMatrixView<")+Type(T())+","+Text(I)+","+
      Text(m.diag().ct())+">"; 
  }
  template <class T, IndexStyle I> inline string Type(
      const DiagMatrixView<T,I>& m)
  { 
    return string("DiagMatrixView<")+Type(T())+","+Text(I)+","+
      Text(m.diag().ct())+">"; 
  }

} // namespace tmv

#endif
