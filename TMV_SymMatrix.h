//---------------------------------------------------------------------------
//
// This file defines the TMV SymMatrix class.
//
// The SymMatrix class and all associated functions are contained
// in the namespace tmv.  Also, the SymMatrix class is a template, 
// so for a SymMatrix of long doubles, one would write 
// tmv::SymMatrix<long double>.  
//
// Note: This implementation of Symmetric Matrices does not have any 
// improvement in storage space.  The storage requirements are still
// an NxN array.  However, there are improvements in calculation time,
// typically a factor of 2 for most operations.  Use PackSymMatrix
// if you require the factor of 2 savings in storage as well.
//
// Caveat: Symmetric Matrices are such that At = A, which implies
// that their diagonal elements are real.  Most routines involving
// Symmetric Matrices assume the reality of the diagonal.  However,
// it is possible in some cases to assign a non-real value to the 
// diagonal.  If the user does this, incorrect answers are likely
// to result.
//
// Constructors:
//
//    SymMatrix<T>(size_t size)
//        Makes a SymMatrix with column size = row size = size 
//        with _uninitialized_ values
//
//    SymMatrix<T>(const Matrix<T>&, bool upper)
//        Turns a regular Matrix into a SymMatrix.
//        If upper = true, this takes the values from the upper triangle.
//        If upper = false, this takes the values from the lower triangle.
//
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//        Return the dimensions of the SymMatrix
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row_a(size_t i, size_t j1, size_t j2)
//    const Vector& row_a(size_t i) const
//        Return the ith row of the SymMatrix up to (not including)
//        the diagonal element
//        With optional j1,j2, give the Subvector with the row from j1-j2.
//        j1,j2 must be < j
//
//    Vector& row_b(size_t i)
//    const Vector& row_b(size_t i, size_t j1, size_t j2) const
//        Return the ith row of the SymMatrix from (including)
//        the diagonal element to the end of the row
//        With optional j1,j2, give the Subvector with the row from j1-j2.
//        j1,j2 must be >= j
//
//    Vector& col_a(size_t j)
//    const Vector& col_a(size_t j) const
//    Vector& col_b(size_t j)
//    const Vector& col_b(size_t j) const
//        Likewise for the columns
//
//    Vector& diag()
//    const Vector& diag() const
//        Return the diagonal of the SymMatrix as a Vector
//        The returned Vector is explicitly real.
//
//    Vector& superdiag(int i)
//    const Vector& superdiag(int i) const
//        i must be > 0
//        Return the super-diagonal i starting at m_0i
//
//    Vector& subdiag(int i)
//    const Vector& subdiag(int i) const
//        i must be < 0
//        Return the sub-diagonal i starting at m_|i|0
//
// Modifying Functions
//
//    SymMatrix& Zero()
//        Sets all elements to 0
//
//    SymMatrix<T>& TransposeSelf() 
//        Transposes the elements of a square SymMatrix or SubSymMatrix
//        This is equivalent to ConjugateSelf.
//
//    SymMatrix& ConjugateSelf()
//        Sets all elements to its conjugate
//
//    SymMatrix& SetToIdentity(x = 1)
//        Set to Identity SymMatrix, or 
//        with a parameter, set to x times Identity SymMatrix
//
//    void Swap(SymMatrix& m1, SymMatrix& m2)
//        Swap the values of two Matrices
//        The Matrices must be the same size
//
//
// SubSymMatrix:
//
//    SubSymMatrix(int i1, int i2)
//        This member function will return a SubSymMatrix
//        from i1 to i2 along the main diagonal, and with the same
//        values for nhi and nlo.
//        The submatrix refers to the same physical elements as the original.
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        Just like the regular matrix version, but all entries in the 
//        submatrix must be within the upper triangle
//
//
// Functions of SymMatrices:
//        (These are all both member functions and functions of a SymMatrix,
//         so Norm(m) and m.Norm() for example are equivalent.)
//
//    Det(m)
//        Returns the determinant of a SymMatrix.
//
//    Trace(m)
//        Returns the trace of a SymMatrix.
//        = sum_i ( a_ii )
//
//    Norm(m) or NormF(m)
//        Return the Frobenius norm of a SymMatrix.
//        = sqrt( sum_ij |a_ij|^2 )
//
//    NormSq(m)
//        Returns the square of Norm().
//
//    Norm1(m) 
//        Returns the 1-norm of a SymMatrix.
//        = max_j (sum_i |a_ij|)
//
//    Norm2(m) 
//        Returns the 2-norm of a SymMatrix.
//        = sqrt( Max Eigenvalue of (A.Adjoint * A) )
//        = Maximum singular value
//        Note: This norm is costly to calculate if one is not 
//              otherwise doing a singular value decomposition
//              of the SymMatrix.
//
//    NormInf(m) 
//        Returns the infinity-norm of a SymMatrix.
//        = max_i (sum_j |a_ij|)
//
//    Transpose(m)
//        Returns the transpose of a SymMatrix.
//
//    Conjugate(m)
//        Returns the conjugate of a SymMatrix
//
//    Adjoint(m)
//        Returns the conjugate of the transpose of a SymMatrix
//        Note: Some people define the adjoint as the cofactor matrix.
//              This is not the same as our definition of the Adjoint.
//
//    Inverse(m)
//        Returns the inverse of m
//
//    m.NewTranspose()
//    m.NewConjugate()
//    m.NewAdjoint()
//    m.NewInverse()
//    m.NewCopy()
//    m.NewView()
//        These all return pointers to new BaseMatrix's equal to the 
//        Transpose, etc. created with new.
//
// Operators:
//
//    SymMatrices have the same operators as a regular Matrix plus mixing
//    of SymMatrices with regular Matrices.
//
//
// I/O: 
//
//    os << m 
//        Writes m to ostream os in the usual Matrix format
//
//    m.WriteCompact(os)
//        Writes m to ostream os in the following compact format:
//          size 
//          m(0,0) m(0,1) m(0,2) ... m(0,size) 
//          m(1,1) ... m(1,size)
//          ...
//          m(size,size)
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the SymMatrix to be read, you can
//        use this form where mptr is a pointer to an undefined SymMatrix.
//        Then *mptr will be created with the right size using new,
//        so you should subsequently delete it.
//
//
// Division Control Functions:
//
//    In addition to the usual division methods, SymMatrices
//    can use Cholesky decomposition if they are positive definite.
//    (This is often the case - especially if their diagonal elements
//    are all positive - see Golub & van Loan, Section 4.2 for more
//    info about positive definite matrices.)  The DivType code for
//    Cholesky decomposition is tmv::CH.
//
//


#ifndef TMV_SymMatrix_H
#define TMV_SymMatrix_H

#include "TMV_Matrix.h"

namespace tmv {

  template <class T, class IT=It<T> > class GenSymMatrix;
  template <class T, class IT=It<T> > class ConstSymMatrixView;
  template <class T, class IT=It<T> > class SymMatrixView;
  template <class T> class SymMatrix;
  template <class T> class SymMatrixComposite;

  template <class T> class SymLUDiv;
  template <class T> class SymSVDiv;
  template <class T> class SymQRDiv;
  template <class T> class SymQRPDiv;

  template <class T> inline const T* ConvertToTCPtr(const void* ptr)
  { return static_cast<const T*>(ptr); }
  template <class T> inline T* ConvertToTCPtr(void* ptr)
  { return static_cast<T*>(ptr); }
  template <class T> inline size_t ConvertToRealStep(size_t step)
  { return step; }
  template <class T> inline size_t ConvertToRealStep<complex<T> >(size_t step)
  { return 2*step; }

  template <class T, class IT> class GenSymMatrix : public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      GenSymMatrix() : BaseMatrix<T>(tmv::LU) {}
      GenSymMatrix(const GenSymMatrix<T,IT>& rhs) : BaseMatrix<T>(rhs.itsdt) {}

      virtual ~GenSymMatrix() {}

      //
      // Access Functions
      //

      size_t colsize() const { return size(); }
      size_t rowsize() const { return size(); }
      virtual size_t size() const =0;

      inline ConstSubVector<T,IT> row_a(const size_t i, 
	  const size_t j1=0, const size_t j2=i) const
      { 
	return ConstSubVector<T,IT>(
	    cptr()+i*stepi()+j1*stepj(), j2-j1, stepj()); 
      }

      inline ConstSubVector<T,typename IT::ConjIT> row_b(const size_t i,
	  const size_t j1=i, const size_t j2=size()) const
      { 
	return ConstSubVector<T,typename IT::ConjIT>(
	    cptr()+j1*stepi()+i*stepj(), j2-j1, stepi()); 
      }

      inline ConstSubVector<T,typename IT::ConjIT> col_a(const size_t j,
	  const size_t i1=0, const size_t i2=j) const
      { 
	return ConstSubVector<T,typename IT::ConjIT>(
	    cptr()+j*stepi()+i1*stepj(), i2-i1, stepj()); 
      }

      inline ConstSubVector<T,IT> col_b(const size_t j,
	  const size_t i1=j, const size_t i2=size()) const
      {
	return ConstSubVector<T,IT>(
	    cptr()+i1*stepi()+j*stepj(), i2-i1, stepi()); 
      }

      inline ConstSubVector<RealType(T)> > diag() const
      {
	return ConstSubVector<RealType(T)>(realcptr(),size(),realdiagstep());
      }

      inline ConstSubVector<T,typename IT::ConjIT> superdiag(int i) const
      {
	TMVAssert(i>=0);
	TMVAssert(i<int(size()));
	return ConstSubVector<T,typename IT::ConjIT>(cptr()+i*stepi(),size()-i,stepi()+stepj());
      }

      inline ConstSubVector<T,IT> subdiag(int i) const
      {
	TMVAssert(i<=0);
	i = -i;
	TMVAssert(i<int(size));
	return ConstSubVector<T,IT>(cptr()+i*stepi(),size()-i,stepi()+stepj());
      }

      inline bool IsSquare() const { return true; }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& m2) const
      { return false; }

      template <class IT2> inline bool SameStorageAs(
	  const GenSymMatrix<T,IT2>& m2) const
      { return (cptr()==m2.cptr()); }

      template <class T2, class IT2> inline bool SameAs(
	  const GenSymMatrix<T2,IT2>& m2) const
      { return false; }

      inline bool SameAs(const GenSymMatrix<T,IT>& m2) const
      { 
	if (this == &m2) return true;
	else return (cptr()==m2.cptr() && size()==m2.size() && 
	    stepi()==m2.stepi() && stepj()==m2.stepj());
      }

      inline operator Matrix<T>() const 
      { 
	Matrix<T> temp(size(),size());
	CopyToMatrix(*this,temp.View());
	return temp;
      }

      //
      // SubSymMatrix
      //

      bool OKSubSymMatrix(int i1, int i2, int istep) const;

      inline ConstSymMatrixView<T,IT> SubSymMatrix(int i1, int i2, int istep=1) const 
      {
	TMVAssert(OKSubSymMatrix(i1,i2,istep));
	int len = (i2-i1)/istep;
	return ConstSymMatrixView<T,IT>(cptr()+i1*(stepi()+stepj()),
	    len, istep*stepi(), istep*stepj());
      }

      bool OKSubMatrix(int i1, int i2, int j1, int j2,
          int istep, int jstep) const;

      inline ConstSubMatrix<T,IT> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep=1, int jstep=1) const 
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstSubMatrix<T,IT>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj());
      }

      inline ConstSymMatrixView<T,IT> View() const
      { 
	return ConstSymMatrixView<T,IT>(cptr(),size(),nlo(),nhi(),
	    stepi(),stepj()); 
      }

      //
      // Functions of SymMatrix
      //

      inline RealType(T) Trace() const
      { return diag().SumElements(); }

      // Default Norm = NormF()
      inline RealType(T) Norm() const
      { return NormF(); }

      // NormF()^2
      RealType(T) NormSq() const;

      // Frobenius norm = sqrt(sum_ij |a_ij|^2 )
      inline RealType(T) NormF() const 
      { return SQRT(NormSq()); }

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const;

      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return Norm1(); }

      inline ConstSymMatrixView<T,typename IT::ConjIT> Transpose() const
      { return ConstSymMatrixView<T,typename IT::ConjIT>(cptr(),size(),stepi(),stepj()); }

      inline ConstSymMatrixView<T,typename IT::ConjIT> Conjugate() const
      { return ConstSymMatrixView<T,typename IT::ConjIT>(cptr(),size(),stepi(),stepj()); }

      inline ConstSymMatrixView<T,IT> Adjoint() const
      { return View(); }

      inline BaseMatrix<T>* NewTranspose() const 
      { return new ConstSymMatrixView<T,typename IT::ConjIT>(Transpose()); }

      inline BaseMatrix<T>* NewConjugate() const
      { return new ConstSymMatrixView<T,typename IT::ConjIT>(Conjugate()); }

      inline BaseMatrix<T>* NewAdjoint() const
      { return new ConstSymMatrixView<T,IT>(Adjoint()); }

      inline BaseMatrix<T>* NewView() const
      { return new ConstSymMatrixView<T,IT>(View()); }

      inline BaseMatrix<T>* NewCopy() const
      { return new SymMatrix<T>(*this); }

      void WriteCompact(ostream& fout) const;

      //
      // Division Control
      //

      virtual inline const T* cptr() const = 0;
      virtual inline const RealType(T)* realcptr() const
      { return ConvertToTCPtr<RealType(T)>(cptr()); }
      virtual inline int stepi() const = 0;
      virtual inline int stepj() const = 0;
      virtual inline int realdiagstep() const
      { return ConvertToRealStep<T>(stepi()+stepj()); }

    protected :

      virtual void setunchanged() const = 0;

      inline T cref(size_t i, size_t j) const 
      {
	TMVAssert(i<colsize() && j<rowsize());
	if (i < j) return IT::ConjIT::VAL(cptr()+j*stepi()+i*stepj());
	else return IT::VAL(cptr()+i*stepi()+j*stepj());
      }

      inline void NewDivider() const 
      {
	switch (this->itsdt) {
	  case tmv::LU : this->itsdiv = new SymLUDiv<T>(*this); break;
	  case tmv::QR : this->itsdiv = new SymQRDiv<T>(*this); break;
	  case tmv::QRP : this->itsdiv = new SymQRPDiv<T>(*this); break;
	  case tmv::SV : this->itsdiv = new SymSVDiv<T>(*this); break;
	  case tmv::CH : TMVAssert(IsHerm());
			 this->itsdiv = new SymCHDiv<T>(*this); break;
	  default : TMVAssert(false);
	}
      }

    private :

      inline void operator=(const GenSymMatrix<T,IT>&) { TMVAssert(false); }

  }; // GenSymMatrix

  template <class T, class IT> class ConstSymMatrixView : public GenSymMatrix<T,IT>
  {
    public :

      ConstSymMatrixView(const ConstSymMatrixView<T,IT>& rhs) :
	itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj) {}

      ConstSymMatrixView(const T* _m, size_t _s, int _si, int _sj) : 
	itsm(_m), itss(_s), itssi(_si), itssj(_sj) {}

      virtual ~ConstSymMatrixView() {}

      inline size_t size() const { return itss; }

      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }

    protected :

      const T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;

      inline bool ischanged() const { return false; }
      inline void setunchanged() const {}

    private :

      inline void operator=(const ConstSymMatrixView<T,IT>&) { TMVAssert(false); }

  }; // ConstSymMatrixView

  template <class T> class SymMatrixComposite;
  // Defined in TMV_SymMatrixArith.h

  template <class T, class IT> class SymMatrixView : public GenSymMatrix<T,IT>
  {

    public:

      //
      // Constructors
      //

      SymMatrixView(const SymMatrixView<T,IT>& rhs) : 
	itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
	changed(false) DEFFIRSTLAST(rhs.first,rhs.last) {}

#ifdef TMVFLDEBUG
      SymMatrixView(T* _m, size_t _s, int _si, int _sj,
	  const T* _first, const T* _last) :
	itsm(_m), itss(_s), itssi(_si), itssj(_sj), 
	first(_first), last(_last) {}
#else
      SymMatrixView(T* _m, size_t _s, int _si, int _sj) :
	itsm(_m), itss(_s), itssi(_si), itssj(_sj) {}
#endif

      virtual ~SymMatrixView() {}

      //
      // Op=
      //

    private:
      template <class T2, class IT2> inline void Copy(const GenSymMatrix<T2,IT2>& m2) const
      {
	changed = true;
	TMVAssert(size() == m2.size());
	diag() = m2.diag();
	for(int i=0; i<size(); ++i) row_a(i) = m2.row_a(i);
      }
    public:

      inline const SymMatrixView<T,IT>& operator=(const SymMatrixView<T,IT>& m2) const
      { if (!SameAs(m2)) Copy(m2); return *this; }

      inline const SymMatrixView<T,IT>& operator=(const GenSymMatrix<T,IT>& m2) const
      { if (!SameAs(m2)) Copy(m2); return *this; }

      template <class T2, class IT2> inline const SymMatrixView<T,IT>& operator=(
	  const GenSymMatrix<T2,IT2>& m2) const
      { Copy(m2); return *this; }

      inline const SymMatrixView<T,IT>& operator=(T x) const 
      { return SetToIdentity(x); }

      inline const SymMatrixView<T,IT>& operator=(const SymMatrixComposite<T>& mcomp) const
      { 
	changed = true;
	TMVAssert(size() == mcomp.size());
	TMVAssert(nlo() == mcomp.nlo());
	TMVAssert(nhi() == mcomp.nhi());
	mcomp.AssignTo(*this);
	return *this;
      }

      //
      // Access
      //
 
      inline size_t size() const { return itss; }

      inline VectorView<T,IT> row_a(size_t i) const
      { return VectorView<T,IT>(ptr()+i*stepi(),i,stepj()); }

      inline VectorView<T,typename IT::ConjIT> row_b(size_t i) const
      { return VectorView<T,typename IT::ConjIT>(ptr()+i*(stepi()+stepj()),size()-i,stepj()); }

      inline VectorView<T,typename IT::ConjIT> col_a(size_t j) const
      { return VectorView<T,typename IT::ConjIT>(ptr()+j*stepi(),j,stepj()); }

      inline VectorView<T,IT> col_b(size_t j) const
      { return VectorView<T,IT>(ptr()+j*(stepi()-stepj()),size()-j,stepj()); }

      inline VectorView<RealType(T)> > diag() const
      {
	return VectorView<RealType(T)>(realptr(),size(),realdiagstep());
      }

      inline VectorView<T,typename IT::ConjIT> superdiag(int i) const
      {
	TMVAssert(i>=0);
	TMVAssert(i<int(size()));
	return VectorView<T,typename IT::ConjIT>(ptr()+i*stepi(),size()-i,stepi()+stepj());
      }

      inline VectorView<T,IT> subdiag(int i) const
      {
	TMVAssert(i<=0);
	i = -i;
	TMVAssert(i<int(size));
	return VectorView<T,IT>(ptr()+i*stepi(),size()-i,stepi()+stepj());
      }

      inline typename IT::REF operator()(size_t i,size_t j) const 
      { return ref(i,j); }


      //
      // Modifying Functions
      //

      inline const SymMatrixView<T,IT>& Zero() const 
      { 
	for(size_t i=0;i<size();++i) row_a(i).Zero();
	diag().Zero();
	return *this; 
      }

      inline const SymMatrixView<T,IT>& TransposeSelf() const
      {
	for(size_t i=0;i<size();++i) row_a(i).ConjugateSelf();
	return *this; 
      }

      inline const SymMatrixView<T,IT>& ConjugateSelf() const
      { 
	for(size_t i=0;i<size();++i) row_a(i).ConjugateSelf();
	return *this; 
      }

      inline const SymMatrixView<T,IT>& SetToIdentity(
	  RealType(T) x=RealType(T)(1)) const 
      { Zero(); diag().SetAllTo(x); return *this; }

      //
      // SubSymMatrix
      //

      inline SymMatrixView<T,IT> SubSymMatrix(int i1, int i2, int istep=1) const
      {
	TMVAssert(this->OKSubSymMatrix(i1,i2,istep));
	const int len = (i2-i1)/istep;
	return SymMatrixView<T,IT>(ptr()+i1*(stepi()+stepj()),
	    len, istep*stepi(), istep*stepj() 
	    FIRSTLAST );
      }

      inline MatrixView<T,IT> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep=1, int jstep=1) const 
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T,IT>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj() 
	    FIRSTLAST );
      }

      inline SymMatrixView<T,IT> View() const
      { return *this; }

      inline SymMatrixView<T,IT> Conjugate() const
      { 
	return SymMatrixView<T,typename IT::ConjIT>(ptr(),size(),
	    stepi(),stepj() FIRSTLAST ); 
      }

      inline SymMatrixView<T,typename IT::ConjIT> Transpose() const
      {
	return SymMatrixView<T,typename IT::ConjIT>(ptr(),size(),
	    stepi(),stepj() FIRSTLAST ); 
      }

      inline SymMatrixView<T,typename IT::ConjIT> Adjoint() const
      { return View(); }


      //
      // I/O
      //

      void Read(istream& fin) const;

      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }

    protected:

      T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;

      mutable bool changed; 

      inline bool ischanged() const { return changed; }
      inline void setunchanged() const { changed = false; }
      inline typename IT::REF ref(size_t i, size_t j) const 
      {
	TMVAssert(i<colsize() && j<rowsize());
	if (i < j) {
	  T* r = ptr()+j*stepi()+i*stepj();
#ifdef TMVFLDEBUG
	  TMVAssert(r >= first);
	  TMVAssert(r < last);
#endif
	  return typename IT::ConjIT::REF(*r);
	} else {
	  T* r = ptr()+i*stepi()+j*stepj();
#ifdef TMVFLDEBUG
	  TMVAssert(r >= first);
	  TMVAssert(r < last);
#endif
	  return typename IT::REF(*r);
	}
      }

#ifdef TMVFLDEBUG
      const T*const first;
      const T*const last;
#endif

  }; // SymMatrixView

  template <class T> class SymMatrix : public GenSymMatrix<T,It<T> > 
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(s) itsm(new T[(s)*(s)]), itss(s), itssi(s), itssj(1), changed(false) DEFFIRSTLAST(itsm,itsm+(s)*(s))

      SymMatrix(size_t _size) : NEW_SIZE(_size) {}

      template <class IT> explicit SymMatrix(
	  const GenMatrix<T,IT>& rhs) : NEW_SIZE(rhs.colsize())
      { 
	TMVAssert(rhs.IsSquare());
	typename It<T>::VIT diagit=subdiag(i).begin();
	typename IT::CVIT rhsdiagit=rhs.diag().begin();
	for(size_t i=0;i<size();++i,++diagit,++rhsdiagit) {
	  row_a(i) = rhs.row(i).SubVector(0,i);
	  *diagit = *rhsdiagit;
	}
      }

      SymMatrix(const SymMatrix<T>& rhs) : NEW_SIZE(rhs.size())
      {
	for(size_t i=0;i<size();++i) row_a(i) = rhs.row_a(i);
	diag() = rhs.diag():
      }

      template <class T2, class IT2> SymMatrix(
	  const GenSymMatrix<T2,IT2>& rhs) : NEW_SIZE(rhs.size())
      {
	for(size_t i=0;i<size();++i) row_a(i) = rhs.row_a(i);
	diag() = rhs.diag():
      }

      SymMatrix<T>(const SymMatrixComposite<T>& mcomp) : NEW_SIZE(mcomp.size())
      { mcomp.AssignTo(View()); }

#undef NEW_SIZE

      virtual ~SymMatrix() { TMVAssert(itsm); delete[] itsm; }


      //
      // Op=
      //

      inline SymMatrix<T>& operator=(const SymMatrix<T>& m2)
      {
	if (&m2 != this) View() = m2.View();
	return *this; 
      }

      template <class IT>
      inline SymMatrix<T>& operator=(const GenSymMatrix<T,IT>& m2)
      { View() = m2; return *this; }

      template <class T2, class IT2> inline SymMatrix<T>& operator=(
	  const GenSymMatrix<T2>& m2)
      { View() = m2; return *this; }

      inline SymMatrix<T>& operator=(T x) { return SetToIdentity(x); }

      inline SymMatrix<T>& operator=(const SymMatrixComposite<T>& mcomp)
      { View() = mcomp; return *this; }


      //
      // Access
      //

      inline size_t size() const { return itss; }

      inline VectorView<T> row_a(size_t i) 
      { return View().row_a(i); }
      inline VectorView<T,typename It<T>::ConjIT> row_b(size_t i) 
      { return View().row_b(i); }
      inline VectorView<T,typename It<T>::ConjIT> col_a(size_t j) 
      { return View().col_a(j); }
      inline VectorView<T> col_b(size_t j) 
      { return View().col_b(j); }
      inline VectorView<RealType(T)> diag()
      { return View().diag(); }
      inline VectorView<T,typename It<T>::ConjIT> superdiag(int i)
      { return View().superdiag(i); }
      inline VectorView<T> subdiag(int i)
      { return View().subdiag(i); }

      inline T& operator()(size_t i,size_t j) 
      { return View()(i,j); }

      inline ConjSubVector<T> row_a(size_t i) const
      { return GenSymMatrix<T,It<T> >::row_a(i); }
      inline ConjSubVector<T,typename It<T>::ConjIT> row_b(size_t i) const
      { return GenSymMatrix<T,It<T> >::row_b(i); }
      inline ConjSubVector<T,typename It<T>::ConjIT> col_a(size_t j) const
      { return GenSymMatrix<T,It<T> >::col_a(j); }
      inline ConjSubVector<T> col_b(size_t j) const
      { return GenSymMatrix<T,It<T> >::col_b(j); }
      inline ConjSubVector<RealType(T)> diag() const
      { return GenSymMatrix<T,It<T> >::diag(); }
      inline ConjSubVector<T,typename It<T>::ConjIT> superdiag(int i) const
      { return GenSymMatrix<T,It<T> >::superdiag(i); }
      inline ConjSubVector<T> subdiag(int i) const
      { return GenSymMatrix<T,It<T> >::subdiag(i); }

      inline T operator()(size_t i,size_t j) const 
      { return GenSymMatrix<T,It<T> >::operator()(i,j); }


      //
      // Modifying Functions
      //

      inline SymMatrix<T>& Zero() { View().Zero(); return *this; }

      inline SymMatrix<T>& TransposeSelf() 
      { View().TransposeSelf(); return *this; }

      inline SymMatrix<T>& ConjugateSelf() 
      { View().ConjugateSelf(); return *this; }

      inline SymMatrix<T>& SetToIdentity(RealType(T) x=RealType(T)(1)) 
      { View().SetToIdentity(x);  return *this; }

      //
      // SubSymMatrix
      //

      inline SymMatrixView<T> SubSymMatrix(size_t i1, size_t i2)
      { return View().SubSymMatrix(i1,i2); }

      inline ConstSymMatrixView<T> SubSymMatrix(size_t i1, size_t i2) const
      { return GenSymMatrix<T>::SubSymMatrix(i1,i2); }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep=1, int jstep=1)
      { return View().SubMatrix(i1,i2,j1,j2,istep,jstep); }

      inline ConstSubMatrix<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep=1, int jstep=1) const
      { return GenSymMatrix<T>::SubMatrix(i1,i2,j1,j2,istep,jstep); }

      inline ConstSymMatrixView<T> View() const
      { return GenSymMatrix<T>::View(); }

      inline SymMatrixView<T> View() 
      { 
	changed = true;
	return SymMatrixView<T>(ptr(),size(),nlo(),nhi(),stepi(),stepj() 
	  FIRSTLAST );
      }

      inline SymMatrixView<T,typename IT<T>::ConjIT> Transpose() 
      { return View().Transpose(); }

      inline ConstSymMatrixView<T,typename IT<T>::ConjIT> Transpose() const
      { return GenSymMatrix<T>::Transpose(); }

      inline SymMatrixView<T,typename It<T>::ConjIT> Conjugate() 
      { return View().Conjugate(); }

      inline ConstSymMatrixView<T,typename It<T>::ConjIT> Conjugate() const
      { return GenSymMatrix<T,It<T> >::Conjugate(); }

      inline SymMatrixView<T> Adjoint()
      { return View().Adjoint(); }

      inline ConstSymMatrixView<T> Adjoint() const
      { return GenSymMatrix<T,It<T> >::Adjoint(); }

      inline const T* cptr() const { return itsm; }
      inline T* ptr() { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }

    protected :

      T*const itsm;
      const size_t itss;
      const size_t itssi;
      const size_t itssj;

      mutable bool changed; 

      inline bool ischanged() const { return changed; }
      inline void setunchanged() const { changed = false; }

#ifdef TMVFLDEBUG
      const T*const first;
      const T*const last;
#endif

  }; // SymMatrix

//---------------------------------------------------------------------------

  //
  // Op= Matrix
  //

  template <class T1, class IT1, class T2, class IT2> 
    inline void CopyToMatrix(const GenSymMatrix<T1,IT1>& m1,
	const MatrixView<T2,IT2>& m2)
    {
      TMVAssert(m2.colsize() == m1.colsize());
      TMVAssert(m2.rowsize() == m1.rowsize());
      temp.Zero();
      for(size_t i=0;i<size();++i) {
	ConstSubVector<T,IT> rowai = row_a(i);
	temp.row(i).SubVector(0,i) = rowai;
	temp.col(i).SubVector(0,i) = rowai.Conjugate();
      }
      temp.diag() = diag();
    }

  //
  // Swap Matrices
  //

  template <class T, class IT1, class IT2> inline void Swap(
      const SymMatrixView<T,IT1>& m1, const SymMatrixView<T,IT2>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    for(size_t i=0;i<m1.size();++i) Swap(m1.row_a(i),m2.row_a(i));
    Swap(m1.diag(),m2.diag());
  }

  template <class T, class IT1> inline void Swap(
      const SymMatrixView<T,IT1>& m1, SymMatrix<T>& m2)
  { Swap(m1,m2.View()); }

  template <class T, class IT2> inline void Swap(
      SymMatrix<T>& m1, const SymMatrixView<T,IT2>& m2)
  { Swap(m1.View(),m2); }

  template <class T> inline void Swap(SymMatrix<T>& m1, SymMatrix<T>& m2)
  { Swap(m1.View(),m2.View()); }

  //
  // Functions of Matrices:
  //

  template <class T, class IT> inline T Det(const GenSymMatrix<T,IT>& m)
  { return m.Det(); }

  template <class T, class IT> inline T Trace(const GenSymMatrix<T,IT>& m)
  { return m.Trace(); }

  template <class T, class IT> inline RealType(T) Norm(const GenSymMatrix<T,IT>& m)
  { return m.Norm(); }

  template <class T, class IT> inline RealType(T) NormSq(const GenSymMatrix<T,IT>& m)
  { return m.NormSq(); }
  
  template <class T, class IT> inline RealType(T) NormF(const GenSymMatrix<T,IT>& m)
  { return m.NormF(); }

  template <class T, class IT> inline RealType(T) Norm1(const GenSymMatrix<T,IT>& m)
  { return m.Norm1(); }

  template <class T, class IT> inline RealType(T) Norm2(const GenSymMatrix<T,IT>& m)
  { return m.Norm2(); }

  template <class T, class IT> inline RealType(T) NormInf(const GenSymMatrix<T,IT>& m)
  { return m.NormInf(); }

  template <class T, class IT> inline ConstSymMatrixView<T,typename IT::ConjIT> Transpose(const GenSymMatrix<T,IT>& m)
  { return m.Transpose(); }

  template <class T, class IT> inline SymMatrixView<T,typename IT::ConjIT> Transpose(const SymMatrixView<T,IT>& m)
  { return m.Transpose(); }

  template <class T> inline SymMatrixView<T,typename It<T>::ConjIT> Transpose(SymMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, class IT> inline ConstSymMatrixView<T,typename IT::ConjIT> Conjugate(const GenSymMatrix<T,IT>& m)
  { return m.Conjugate(); }

  template <class T, class IT> inline SymMatrixView<T,typename IT::ConjIT> Conjugate(const SymMatrixView<T,IT>& m)
  { return m.Conjugate(); }

  template <class T> inline SymMatrixView<T,typename It<T>::ConjIT> Conjugate(SymMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, class IT> inline ConstSymMatrixView<T,IT> Adjoint(const GenSymMatrix<T,IT>& m)
  { return m.Adjoint(); }

  template <class T, class IT> inline SymMatrixView<T,IT> Adjoint(const SymMatrixView<T,IT>& m)
  { return m.Adjoint(); }

  template <class T> inline SymMatrixView<T,It<T> > Adjoint(SymMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, class IT> inline SymMatrix<T> Inverse(const GenSymMatrix<T,IT>& m)
  { return m.Inverse(); }

  template <class T, class IT> inline SymMatrix<T> InverseATA(const GenSymMatrix<T,IT>& m)
  { return m.InverseATA(); }

  //
  // SymMatrix ==, != SymMatrix
  //

  template <class T, class IT1, class IT2> inline bool operator==(
      const GenSymMatrix<T,IT1>& m1, const GenSymMatrix<T,IT2>& m2)
  { return false; }
  
  template <class T, class IT> bool operator==(
      const GenSymMatrix<T,IT>& m1, const GenSymMatrix<T,IT>& m2);

  template <class T, class IT1, class IT2> inline bool operator!=(
      const GenSymMatrix<T,IT1>& m1, const GenSymMatrix<T,IT2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T> istream& operator>>(istream& fin, SymMatrix<T>& m);

  template <class T> istream& operator>>(istream& fin, SymMatrix<T>* m);

  template <class T, class IT> istream& operator>>(
      istream& fin, const SymMatrixView<T,IT>& m);

  template <class T> inline std::string Type(const SymMatrix<T>&)
  { return std::string("SymMatrix<") + Type(T()) + ">"; }
  template <class T,class IT> inline std::string Type(const GenSymMatrix<T,IT>&)
  { return std::string("GenSymMatrix<") + Type(T()) + "," + Type(IT()) + ">"; }
  template <class T,class IT> inline std::string Type(const ConstSymMatrixView<T,IT>&)
  { return std::string("ConstSymMatrixView<") + Type(T()) + "," + Type(IT()) + ">"; }
  template <class T,class IT> inline std::string Type(const SymMatrixView<T,IT>&)
  { return std::string("SymMatrixView<") + Type(T()) + "," + Type(IT()) + ">"; }

}; // namespace tmv

#endif
