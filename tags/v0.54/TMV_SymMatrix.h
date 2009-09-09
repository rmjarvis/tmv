//---------------------------------------------------------------------------
//
// This file defines the TMV SymMatrix and HermMatrix classes.
//
// Constructors:
//
//    SymMatrix is used for symmetric matrices, and HermMatrix is used
//    for Hermitian matrices.  For real matrices, these are the same thing:
//    A = A.Transpose().
//    But for complex, they are different:
//    A_sym = A_sym.Transpose()
//    A_herm = A_herm.Adjoint()
//
//    For these notes, I will always write SymMatrix, but (except where
//    otherwise indicated) everything applies the same for Sym and Herm.
//    Also, the Views keep track of sym/herm difference with a parameter,
//    so it is always a GenSymMatrix, ConstSymMatrixView, or 
//    SymMatrixView - never Herm in any of these.
//
//    Caveat: Complex Hermitian matrices are such that A = At, which 
//    implies that their diagonal elements are real.  Many routines
//    involving HermMatrixes assume the reality of the diagonal.
//    However, it is possible to assign a non-real value to a diagonal
//    element.  If the user does this, the results are undefined.
//
//    In addition to the type template parameter (T), SymMatrixes have two
//    additional template parameters:
//        UpLoType uplo = Upper || Lower
//        StorageType stor = RowMajor || ColMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {Upper, RowMajor}
//
//        The first, uplo, refers to which triangular half stores the actual 
//        data.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    SymMatrix<T,uplo,stor>(size_t n)
//        Makes a Symmetric Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SymMatrix<T,uplo,stor>(size_t n, T x)
//        Makes a Symmetric Matrix with values all initialized to x
//        For Hermitian matrixces, x must be real.
//
//    SymMatrix<T,uplo,stor>(size_t n, const valarray<T>& m)
//    SymMatrix<T,uplo,stor>(size_t n, const T* m)
//    SymMatrix<T,uplo,stor>(size_t n, const vector<T>& m)
//        Make a Symmetric Matrix with copies the elements of m which
//        fall in tha appropriate upper or lower triangle.
//        The lengths of the arrays must be n*n, but only about half the 
//        elements are used.
//
//    SymMatrix<T,uplo,stor>(const Matrix<T>& m)
//        Makes a SymMatrix which copies the corresponding elements of m.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(const Matrix<T>& m, uplo)
//    ConstSymMatrixView<T> HermMatrixViewOf(const Matrix<T>& m, uplo)
//        Makes a constant SymMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//
//    SymMatrixView<T> SymMatrixViewOf(Matrix<T>& m, uplo)
//    SymMatrixView<T> HermMatrixViewOf(Matrix<T>& m, uplo)
//        Makes a modifiable SymMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the original Matrix.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(const T* m, size, uplo, stor)
//    ConstSymMatrixView<T> HermMatrixViewOf(const T* m, size, uplo, stor)
//    SymMatrixView<T> SymMatrixViewOf(T* m, size, uplo, stor)
//    SymMatrixView<T> HermMatrixViewOf(T* m, size, uplo, stor)
//        View the actual memory pointed to by m as a SymMatrix/HermMatrix
//        with the given size and storage.
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the SymMatrix
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row(size_t i, size_t j1, size_t j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row
//        and be entirely in the upper or lower triangle.
//
//    Vector& col(size_t j, size_t i1, size_t i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column
//        and be entirely in the upper or lower triangle.
//
//    Vector& diag()
//        Return the main diagonal
//
//    Vector& diag(int i, size_t j1, size_t j2)
//    Vector& diag(int i)
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal Vector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//
// Modifying Functions:
//
//    Zero()
//    SetAllTo(T x) 
//        For HermMatrix, x must be real.
//    Clip(RT thresh)
//    ConjugateSelf()
//    TransposeSelf()
//    SetToIdentity(x = 1)
//    SwapRowsCols(i1,i2)
//        Equivalent to swapping rows i1,i2 then swapping cols i1,i2.
//    PermuteRowsCols(const size_t* p)
//    ReversePermuteRowsCols(const size_t* p)
//        Perform a series of row/col swaps.
//    void Swap(SymMatrix& m1, SymMatrix& m2)
//        The SymMatrices must be the same size and Hermitianity.
//
// Views of a SymMatrix:
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the upper 
//        or lower triangle.
//
//    SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a VectorView which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.  The subvector must be completely with the upper or
//        lower triangle.
//
//    SubSymMatrix(int i1, int i2, int istep)
//        Returns the SymMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with a SymMatrix of size 10, the x's below
//        are the original data, the O's are the SubSymMatrix returned
//        with the command SubSymMatrix(3,11,2), and the #'s are the 
//        SubSymMatrix returned with SubSymMatrix(0,3)
//
//        ###xxxxxxx
//        ###xxxxxxx
//        ###xxxxxxx
//        xxxOxOxOxO
//        xxxxxxxxxx
//        xxxOxOxOxO
//        xxxxxxxxxx
//        xxxOxOxOxO
//        xxxxxxxxxx
//        xxxOxOxOxO
//
//    Transpose(m)
//    Adjoint(m)
//    Conjugate(m)
//
//
// Functions of Matrixs:
//
//    Det(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//
//    m.Inverse()
//    Inverse(m)
//    m.Inverse(minv) // Takes either a SymMatrix or Matrix argument
//    m.InverseATA(invata)
//
//    m.NewTranspose()
//    m.NewConjugate()
//    m.NewAdjoint()
//    m.NewInverse()
//    m.NewCopy()
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
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(size-1,0) ... m(size-1,size-1) )
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the SymMatrix to be read, you can
//        use this form where mptr is an auto_ptr to an undefined SymMatrix.
//
//
// Division Control Functions:
//
//    m.DivideUsing(dt)
//    where dt is LU, CH, SV, SVS, SVU, SVV
//     
//    LUD(), CHD(), SymSVD(), and HermSVD() return the
//        corresponding Divider classes.  Since Symmetric and Hermitian
//        LUD and SVD classes are different types, they need to have 
//        different names for the methods that return them.
//
//    For SymMatrixes, LU actually does an LDLT decomposition.
//        ie. L(Lower Triangle) * Diagonal * L.Transpose()
//        For non-positive definite matrices, D may have 2x2 blocks along
//        the diagonal, which can then be further decomposed via normal LU
//        decomposition.
//        
//    The option unique to hermitian matrixes is CH - Cholskey decomposition.
//        This is only appropriate if you know that your HermMatrix is 
//        positive definite (ie. all eigenvalues are positive).
//        (This is guaranteed, for example, if all the square 
//        submatrices have positive determinant.)
//        In this case, the SymMatrix can be decomposed into L*L.Adjoint().
//
//    Finally, a difference for the SV Decomposition for Hermitian matrices
//        is that the decomposition can be done with U = Vt.  In this case,
//        the "singular values" are really the eigenvalues of A.
//        The only caveat is that they may be negative, whereas the usual
//        definition of singular values is that they are all positive.
//        Requiring positivity would destroy U = Vt, so it seemed more
//        useful to leave them as the actual eigenvalues.  Just keep that
//        in mind if you use the singular values for anything that expects
//        them to be positive.


#ifndef TMV_SymMatrix_H
#define TMV_SymMatrix_H

#include "TMV_BaseMatrix.h"
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  template <class T> class GenSymMatrix;
  template <class T, IndexStyle I=CStyle> class ConstSymMatrixView;
  template <class T, IndexStyle I=CStyle> class SymMatrixView;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle> class SymMatrix;
  template <class T, UpLoType U=Upper, StorageType S=RowMajor, IndexStyle I=CStyle> class HermMatrix;
  template <class T> class SymMatrixComposite; 
  template <class T, class Tm> class SumSX;
  template <class T, class Tm> class ProdXS;
  template <class T, class T1, class T2> class SumSS;
  template <class T> class SymDivider;

  template <class T> class SymLUDiv;
  template <class T> class HermCHDiv;
  template <class T> class HermSVDiv;
  template <class T> class SymSVDiv;
  template <class T, class Tm> class QuotXS;

#define UTransOf(U) (U==Upper ? Lower : Upper)

  template <class T> class SymMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<SymMatrix<T> > m;
      char exp,got;
      size_t s;
      bool is, iseof, isbad;

      inline SymMatrixReadError(
	  size_t _i, size_t _j, const GenSymMatrix<T>& _m,
	  istream& _is) throw() :
	ReadError("SymMatrix"), 
	i(_i), j(_j), m(new SymMatrix<T>(_m)), exp(0), got(0), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymMatrixReadError(istream& _is) throw() :
	ReadError("SymMatrix"), 
	i(0), j(0), m(0), exp(0), got(0), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymMatrixReadError(
	  size_t _i, size_t _j, const GenSymMatrix<T>& _m,
	  istream& _is, char _e, char _g) throw() :
	ReadError("SymMatrix"), 
	i(_i), j(_j), m(new SymMatrix<T>(_m)), exp(_e), got(_g), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymMatrixReadError(istream& _is, char _e, char _g) throw() :
	ReadError("SymMatrix"), 
	i(0), j(0), m(0), exp(_e), got(_g), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymMatrixReadError(
	  const GenSymMatrix<T>& _m, istream& _is, size_t _s) throw() :
	ReadError("SymMatrix"), 
	i(0), j(0), m(new SymMatrix<T>(_m)), exp(0), got(0), s(_s),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline SymMatrixReadError(const SymMatrixReadError<T>& rhs) :
	ReadError("SymMatrix"), 
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), s(rhs.s),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~SymMatrixReadError() throw() {}

      virtual void Write(ostream& os) const throw();
  };

  template <class T> class HermMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<HermMatrix<T> > m;
      char exp,got;
      size_t s;
      T dv;
      bool is, iseof, isbad;

      inline HermMatrixReadError(
	  size_t _i, size_t _j, const GenSymMatrix<T>& _m,
	  istream& _is) throw() :
	ReadError("HermMatrix"), 
	i(_i), j(_j), m(new HermMatrix<T>(_m)), exp(0), got(0), s(_m.size()),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermMatrixReadError(istream& _is) throw() :
	ReadError("HermMatrix"), 
	i(0), j(0), m(0), exp(0), got(0), s(0),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermMatrixReadError(
	  size_t _i, size_t _j, const GenSymMatrix<T>& _m,
	  istream& _is, char _e, char _g) throw() :
	ReadError("HermMatrix"), 
	i(_i), j(_j), m(new HermMatrix<T>(_m)), exp(_e), got(_g), s(_m.size()),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermMatrixReadError(
	  size_t _i, size_t _j, const GenSymMatrix<T>& _m,
	  istream& _is, T _dv) throw() :
	ReadError("HermMatrix"), 
	i(_i), j(_j), m(new HermMatrix<T>(_m)), exp(0), got(0), s(_m.size()),
	dv(_dv), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermMatrixReadError(istream& _is, char _e, char _g) throw() :
	ReadError("HermMatrix"), 
	i(0), j(0), m(0), exp(_e), got(_g), s(0),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermMatrixReadError(
	  const GenSymMatrix<T>& _m, istream& _is, size_t _s) throw() :
	ReadError("HermMatrix"), 
	i(0), j(0), m(new HermMatrix<T>(_m)), exp(0), got(0), s(_s),
        dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline HermMatrixReadError(const HermMatrixReadError<T>& rhs) :
	ReadError("HermMatrix"), 
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), s(rhs.s),
	dv(rhs.dv), is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~HermMatrixReadError() throw() {}

      virtual void Write(ostream& os) const throw();
  };


  template <class T> class GenSymMatrix : 
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline GenSymMatrix(SymType h, UpLoType u, StorageType s, ConjItType c) :
	BaseMatrix<T>(LU), itssym(IsComplex(T())?h:Sym),
        itsuplo(u), itsstor(s), itsct(c) {}
      inline GenSymMatrix(const GenSymMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs), itssym(rhs.itssym),
        itsuplo(rhs.itsuplo), itsstor(rhs.itsstor), itsct(rhs.itsct)
      { TMVAssert(IsComplex(T()) || itssym==Sym); }
      inline ~GenSymMatrix() {}

      //
      // Access Functions
      //

      inline T operator()(size_t i, size_t j) const
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	return cref(i,j);
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert( i<=j1 || j2<=i+1 );
	if ((i<=j1 && itsuplo==Upper) || (j2<=i+1 && itsuplo==Lower))
	  return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	      itsct); 
	else
	  return ConstVectorView<T>(cptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
	      issym()?itsct:ConjOf(T,itsct));
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert( j<=i1 || i2<=j+1 );
	if ((i2<=j+1 && itsuplo==Upper) || (j<=i1 && itsuplo==Lower))
	  return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	      itsct); 
	else 
	  return ConstVectorView<T>(cptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
	      issym()?itsct:ConjOf(T,itsct));
      }

      inline ConstVectorView<T> diag() const
      { return ConstVectorView<T>(cptr(),size(),stepi()+stepj(),itsct); }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(-i<=int(size())); 
	if (i>=0)
	  if (itsuplo==Upper)
	    return ConstVectorView<T>(cptr()+i*stepj(),size()-i,
		stepi()+stepj(),itsct); 
	  else
	    return ConstVectorView<T>(cptr()+i*stepi(),size()-i,
		stepi()+stepj(),issym()?itsct:ConjOf(T,itsct));
	else
	  if (itsuplo==Upper)
	    return ConstVectorView<T>(cptr()-i*stepj(),size()+i,
		stepi()+stepj(),issym()?itsct:ConjOf(T,itsct));
	  else
	    return ConstVectorView<T>(cptr()-i*stepi(),size()+i,
		stepi()+stepj(),itsct); 
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(-i<=int(size())); 
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()-abs(i));
	const int ds = stepi()+stepj();
	if (i>=0)
	  if (itsuplo==Upper)
	    return ConstVectorView<T>(cptr()+i*stepj()+j1*ds,j2-j1,ds,itsct); 
	  else
	    return ConstVectorView<T>(cptr()+i*stepi()+j1*ds,j2-j1,ds,
		issym()?itsct:ConjOf(T,itsct));
	else
	  if (itsuplo==Upper)
	    return ConstVectorView<T>(cptr()-i*stepj()+j1*ds,j2-j1,ds,
		issym()?itsct:ConjOf(T,itsct));
	  else
	    return ConstVectorView<T>(cptr()-i*stepi()+j1*ds,j2-j1,ds,itsct); 
      }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      inline bool SameStorageAs(const GenSymMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      template <class T2> inline bool SameAs(const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameAs(const GenSymMatrix<T>& m2) const
      { 
	if (this == &m2) return true;
	else if (cptr()==m2.cptr() && size()==m2.size() && 
	      itssym == m2.itssym) {
	  if (itsuplo == m2.itsuplo)
	    return (stepi() == m2.stepi() && stepj() == m2.stepj()
		&& itsct == m2.itsct);
	  else
	    return (stepi() == m2.stepj() && stepj() == m2.stepi()
		&& issym() == (itsct==m2.itsct));
	} else return false;
      }

      inline void CopyToMatrix(const MatrixView<RealType(T)>& m2) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	UpperTriMatrixViewOf(m2) = UpperTri();
	if (size() > 0)
	  if (issym())
	    LowerTriMatrixViewOf(m2).OffDiag() =
	      UpperTriMatrixViewOf(m2).OffDiag().Transpose();
	  else
	    LowerTriMatrixViewOf(m2).OffDiag() =
	      UpperTriMatrixViewOf(m2).OffDiag().Adjoint();
      }

      inline void CopyToMatrix(const MatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	UpperTriMatrixViewOf(m2) = UpperTri();
	if (size() > 0)
	  if (issym())
	    LowerTriMatrixViewOf(m2).OffDiag() =
	      UpperTriMatrixViewOf(m2).OffDiag().Transpose();
	  else
	    LowerTriMatrixViewOf(m2).OffDiag() =
	      UpperTriMatrixViewOf(m2).OffDiag().Adjoint();
      }

      //
      // SubMatrix
      //

      bool OKSubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const;

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	TMVAssert(i2<=j1 || j2<=i1);
	if ((i2<=j1 && itsuplo==Upper) || (j2<=i1 && itsuplo==Lower))
	  return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	      i2-i1,j2-j1,stepi(),stepj(),itsstor,itsct);
	else
	  return ConstMatrixView<T>(cptr()+i1*stepj()+j1*stepi(),
	      i2-i1,j2-j1,stepj(),stepi(),TransOf(itsstor),
	      issym()?itsct:ConjOf(T,itsct));
      }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	TMVAssert(i2<=j1 || j2<=i1);
	if ((i2<=j1 && itsuplo==Upper) || (j2<=i1 && itsuplo==Lower))
	  return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	      (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	      newstor,itsct);
	else
	  return ConstMatrixView<T>(cptr()+i1*stepj()+j1*stepi(),
	      (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
	      TransOf(newstor),issym()?itsct:ConjOf(T,itsct));
      }

      bool OKSubVector(
	  size_t i, size_t j, int istep, int jstep, size_t n) const;

      inline ConstVectorView<T> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t n) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,n));
	if ((i<=j && itsuplo==Upper) || (j<=i && itsuplo==Lower))
	  return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),n,
	      istep*stepi()+jstep*stepj(),itsct);
	else
	  return ConstVectorView<T>(cptr()+i*stepj()+j*stepi(),n,
	      istep*stepj()+jstep*stepi(),issym()?itsct:ConjOf(T,itsct));
      }

      bool OKSubSymMatrix(int i1, int i2, int istep) const;

      inline ConstSymMatrixView<T> SubSymMatrix(int i1, int i2) const
      {
	TMVAssert(OKSubSymMatrix(i1,i2,1));
	return ConstSymMatrixView<T>(cptr()+i1*(stepi()+stepj()),i2-i1,
	    stepi(),stepj(),itssym,itsuplo,itsstor,itsct);
      }

      inline ConstSymMatrixView<T> SubSymMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(OKSubSymMatrix(i1,i2,istep));
	return ConstSymMatrixView<T>(cptr()+i1*(stepi()+stepj()),
	    (i2-i1)/istep,istep*stepi(),istep*stepj(),itssym,itsuplo,
	    istep==1 ? itsstor : NoMajor,itsct);
      }

      inline ConstUpperTriMatrixView<T> UpperTri(
	  DiagType dt = NonUnitDiag) const
      {
	if (itsuplo == Upper)
	  return ConstUpperTriMatrixView<T>(cptr(),size(),
	      stepi(),stepj(),dt,itsstor,itsct);
	else
	  return ConstUpperTriMatrixView<T>(cptr(),size(),
	      stepj(),stepi(),dt,TransOf(itsstor),
	      issym()?itsct:ConjOf(T,itsct));
      }

      inline ConstLowerTriMatrixView<T> LowerTri(
	  DiagType dt = NonUnitDiag) const
      {
	if (itsuplo == Lower)
	  return ConstLowerTriMatrixView<T>(cptr(),size(),
	      stepi(),stepj(),dt,itsstor,itsct);
	else
	  return ConstLowerTriMatrixView<T>(cptr(),size(),
	      stepj(),stepi(),dt,TransOf(itsstor),
	      issym()?itsct:ConjOf(T,itsct));
      }

      inline ConstSymMatrixView<RealType(T)> Real() const
      {
	return ConstSymMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    Sym, itsuplo, IsReal(T()) ? itsstor : NoMajor,NonConj);
      }

      inline ConstSymMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(issym());
	// The imaginary part of a Hermitian matrix is anti-symmetric
	return ConstSymMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    size(),2*stepi(),2*stepj(), Sym,itsuplo,NoMajor,NonConj);
      }

      // 
      // Views
      //

      inline ConstSymMatrixView<T> View() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
	    itssym,itsuplo,itsstor,itsct);
      }

      inline ConstSymMatrixView<T> Transpose() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
	    itssym,UTransOf(itsuplo),TransOf(itsstor),itsct);
      }

      inline ConstSymMatrixView<T> Conjugate() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
	    itssym,itsuplo,itsstor,ConjOf(T,itsct));
      }

      inline ConstSymMatrixView<T> Adjoint() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
	    itssym,UTransOf(itsuplo),TransOf(itsstor),ConjOf(T,itsct));
      }

      inline SymMatrixView<T> NonConst() const
      {
	return SymMatrixView<T>(const_cast<T*>(cptr()),size(),
	    stepi(),stepj(),itssym,itsuplo,itsstor,itsct
	    FIRSTLAST1(cptr(),row(colsize()-1,0,colsize()).end().GetP()));
      }

      //
      // Functions of Matrix
      //

      inline T Trace() const
      { return diag().SumElements(); }

      inline RealType(T) Norm() const 
      { return NormF(); }

      RealType(T) NormF() const;

      // NormF()^2
      RealType(T) NormSq() const;

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const;

      // 2-Norm defined in BaseMatrix.h

      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return Norm1(); }

      // = max_i,j (|a_ij|)
      RealType(T) MaxAbsElement() const;

      inline QuotXS<T,T> Inverse() const
      { return QuotXS<T,T>(T(1),*this); }

      void Inverse(const SymMatrixView<T>& sinv) const;

      inline void Inverse(const MatrixView<T>& minv) const
      { BaseMatrix<T>::Inverse(minv); }

      template <UpLoType U, StorageType S, IndexStyle I> inline void Inverse(
	  SymMatrix<T,U,S,I>& sinv) const
      {
	TMVAssert(issym());
	Inverse(sinv.View()); 
      }

      template <UpLoType U, StorageType S, IndexStyle I> inline void Inverse(
	  HermMatrix<T,U,S,I>& sinv) const
      {
	TMVAssert(isherm());
	Inverse(sinv.View()); 
      }

      template <StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T,S,I>& minv) const
      { BaseMatrix<T>::Inverse(minv.View()); }

      using BaseMatrix<T>::InverseATA;

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& minv) const
      { Inverse(minv.View()); }

      inline auto_ptr<BaseMatrix<T> > NewTranspose() const 
      {
	auto_ptr<BaseMatrix<T> > a(new ConstSymMatrixView<T>(Transpose())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewConjugate() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstSymMatrixView<T>(Conjugate())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewAdjoint() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstSymMatrixView<T>(Adjoint())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewView() const
      {
	auto_ptr<BaseMatrix<T> > a(new ConstSymMatrixView<T>(View())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewCopy() const
      { 
	auto_ptr<BaseMatrix<T> > a;
	if (issym()) {
	  if (itsuplo == Upper) {
	    if (isrm()) a.reset(new SymMatrix<T,Upper,RowMajor>(*this)); 
	    else a.reset(new SymMatrix<T,Upper,ColMajor>(*this)); 
	  } else {
	    if (isrm()) a.reset(new SymMatrix<T,Lower,RowMajor>(*this)); 
	    else a.reset(new SymMatrix<T,Lower,ColMajor>(*this)); 
	  }
	} else {
	  if (itsuplo == Upper) {
	    if (isrm()) a.reset(new HermMatrix<T,Upper,RowMajor>(*this)); 
	    else a.reset(new HermMatrix<T,Upper,ColMajor>(*this)); 
	  } else {
	    if (isrm()) a.reset(new HermMatrix<T,Lower,RowMajor>(*this)); 
	    else a.reset(new HermMatrix<T,Lower,ColMajor>(*this)); 
	  }
	}
	return a;
      }

      // 
      // Division Control
      //

      using BaseMatrix<T>::GetDiv;

      inline void DivideUsing(DivType dt) const
      {
	TMVAssert(dt == CH || dt == LU || SV_Type(dt));
	TMVAssert(isherm() || dt != CH);
	BaseMatrix<T>::DivideUsing(dt); 
      }

      inline const SymLUDiv<T>& LUD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const SymLUDiv<T>*>(GetDiv()));
	return *dynamic_cast<const SymLUDiv<T>*>(GetDiv());
      }

      inline const HermCHDiv<T>& CHD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const HermCHDiv<T>*>(GetDiv()));
	return *dynamic_cast<const HermCHDiv<T>*>(GetDiv());
      }

      inline const HermSVDiv<T>& HermSVD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const HermSVDiv<T>*>(GetDiv()));
	return *dynamic_cast<const HermSVDiv<T>*>(GetDiv());
      }

      inline const SymSVDiv<T>& SymSVD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const SymSVDiv<T>*>(GetDiv()));
	return *dynamic_cast<const SymSVDiv<T>*>(GetDiv());
      }

      //
      // I/O
      //

      void Write(ostream& os) const;
      void WriteCompact(ostream& os) const;
      void Write(ostream& os, RealType(T) thresh) const;
      void WriteCompact(ostream& os, RealType(T) thresh) const;

      //
      // Arithmetic Helpers
      //

      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }
      virtual size_t size() const = 0;
      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual inline SymType sym() const { return itssym; }
      virtual inline UpLoType uplo() const { return itsuplo; }
      virtual inline StorageType stor() const { return itsstor; }
      virtual inline ConjItType ct() const { return itsct; }
      virtual inline bool isrm() const { return itsstor == RowMajor; }
      virtual inline bool iscm() const { return itsstor == ColMajor; }
      inline bool isconj() const
      { 
	TMVAssert(IsComplex(T()) || itsct==NonConj);
	return IsComplex(T()) && itsct == Conj; 
      }
      inline bool isherm() const { return IsReal(T()) || itssym == Herm; }
      inline bool issym() const { return IsReal(T()) || itssym == Sym; }

    protected :

      virtual T cref(size_t i, size_t j) const;

      void NewDivider() const;

    private :

      const SymType itssym;
      const UpLoType itsuplo;
      const StorageType itsstor;
      const ConjItType itsct;

      inline void operator=(const GenSymMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenSymMatrix

  template <class T, IndexStyle I> class ConstSymMatrixView : 
    public GenSymMatrix<T>
  {
    public :

      inline ConstSymMatrixView(const ConstSymMatrixView<T,I>& rhs) :
	GenSymMatrix<T>(rhs), itsm(rhs.itsm),
	itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj) {}

      inline ConstSymMatrixView(const GenSymMatrix<T>& rhs) :
	GenSymMatrix<T>(rhs), itsm(rhs.cptr()),
	itss(rhs.size()), itssi(rhs.stepi()), itssj(rhs.stepj()) {}

      inline ConstSymMatrixView(
	  const T* _m, size_t _s, int _si, int _sj,
	  SymType insym, UpLoType inuplo, StorageType instor,
	  ConjItType inct) : 
	GenSymMatrix<T>(insym,inuplo,instor,inct), itsm(_m), itss(_s),
	itssi(_si), itssj(_sj)
      { 
	TMVAssert(instor==RowMajor ? _sj == 1 : instor==ColMajor ?
	    _si==1 : true);
      }

      inline ~ConstSymMatrixView() {}

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }

    protected :

      const T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;

    private :

      inline void operator=(const ConstSymMatrixView<T,I>&) 
      { TMVAssert(FALSE); }

  }; // ConstSymMatrixView

  template <class T> class ConstSymMatrixView<T,FortranStyle> : 
    public ConstSymMatrixView<T,CStyle>
  {
    public :

      inline ConstSymMatrixView(const ConstSymMatrixView<T,FortranStyle>& rhs) :
	ConstSymMatrixView<T,CStyle>(rhs) {}

      inline ConstSymMatrixView(const ConstSymMatrixView<T,CStyle>& rhs) :
	ConstSymMatrixView<T,CStyle>(rhs) {}

      inline ConstSymMatrixView(const GenSymMatrix<T>& rhs) :
	ConstSymMatrixView<T,CStyle>(rhs) {}

      inline ConstSymMatrixView(const T* _m, size_t _s, int _si, int _sj,
	  SymType insym, UpLoType inuplo, StorageType instor,
	  ConjItType inct) : 
	ConstSymMatrixView<T,CStyle>(_m,_s,_si,_sj,insym,inuplo,instor,inct) {}

      inline ~ConstSymMatrixView() {}

      //
      // Access Functions
      //

      inline T operator()(size_t i, size_t j) const
      {
	TMVAssert(i>0 && i<=size());
	TMVAssert(j>0 && j<=size());
	return cref(i-1,j-1);
      }

      inline ConstVectorView<T,FortranStyle> row(
	  size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i>0 && i<=size());
	TMVAssert(j1>0 && j1<=j2);
	TMVAssert(j2<=size());
	return GenSymMatrix<T>::row(i-1,j1-1,j2);
      }

      inline ConstVectorView<T,FortranStyle> col(
	  size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j>0 && j<=size());
	TMVAssert(i1>0 && i1<=i2);
	TMVAssert(i2<=size());
	return GenSymMatrix<T>::col(j-1,i1-1,i2);
      }

      inline ConstVectorView<T,FortranStyle> diag() const
      { return GenSymMatrix<T>::diag(); }

      inline ConstVectorView<T,FortranStyle> diag(int i) const
      { return GenSymMatrix<T>::diag(i); }

      inline ConstVectorView<T,FortranStyle> diag(
	  int i, size_t j1, size_t j2) const
      {
	TMVAssert(j1>0);
	return GenSymMatrix<T>::diag(i,j1-1,j2);
      }

      //
      // SubMatrix
      //

      bool OKSubMatrix(int i1, int i2, int j1, int j2, 
	  int istep, int jstep) const;

      bool OKSubVector(size_t i, size_t j, int istep, int jstep, 
	  size_t n) const;

      bool OKSubSymMatrix(int i1, int i2, int istep) const;

      inline ConstMatrixView<T,FortranStyle> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return GenSymMatrix<T>::SubMatrix(i1-1,i2,j1-1,j2);
      }

      inline ConstMatrixView<T,FortranStyle> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return GenSymMatrix<T>::SubMatrix(i1-1,i2-1+istep,j1-1,j2-1+jstep,
	    istep,jstep);
      }

      inline ConstVectorView<T,FortranStyle> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t n) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,n));
	return GenSymMatrix<T>::SubVector(i-1,j-1,istep,jstep,n);
      }

      inline ConstSymMatrixView<T,FortranStyle> SubSymMatrix(
	  int i1, int i2) const
      {
	TMVAssert(OKSubSymMatrix(i1,i2,1));
	return GenSymMatrix<T>::SubSymMatrix(i1-1,i2);
      }

      inline ConstSymMatrixView<T,FortranStyle> SubSymMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(OKSubSymMatrix(i1,i2,istep));
	return GenSymMatrix<T>::SubSymMatrix(i1-1,i2-1+istep,istep);
      }

      inline ConstUpperTriMatrixView<T,FortranStyle> UpperTri(
	  DiagType dt = NonUnitDiag) const
      { return GenSymMatrix<T>::UpperTri(dt); }

      inline ConstLowerTriMatrixView<T,FortranStyle> LowerTri(
	  DiagType dt = NonUnitDiag) const
      { return GenSymMatrix<T>::LowerTri(dt); }

      inline ConstSymMatrixView<RealType(T),FortranStyle> Real() const
      { return GenSymMatrix<T>::Real(); }

      inline ConstSymMatrixView<RealType(T),FortranStyle> Imag() const
      { return GenSymMatrix<T>::Imag(); }

      // 
      // Views
      //

      inline ConstSymMatrixView<T,FortranStyle> View() const
      { return GenSymMatrix<T>::View(); }

      inline ConstSymMatrixView<T,FortranStyle> Transpose() const
      { return GenSymMatrix<T>::Transpose(); }

      inline ConstSymMatrixView<T,FortranStyle> Conjugate() const
      { return GenSymMatrix<T>::Conjugate(); }

      inline ConstSymMatrixView<T,FortranStyle> Adjoint() const
      { return GenSymMatrix<T>::Adjoint(); }

      inline SymMatrixView<T,FortranStyle> NonConst() const
      { return GenSymMatrix<T>::NonConst(); }

      using ConstSymMatrixView<T,CStyle>::size;

    protected :

      using GenSymMatrix<T>::cref;

    private :

      inline void operator=(const ConstSymMatrixView<T,FortranStyle>&) 
      { TMVAssert(FALSE); }

  }; // FortranStyle ConstSymMatrixView
  
  template <class T, IndexStyle I> class SymMatrixView : 
    public GenSymMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline SymMatrixView(const SymMatrixView<T,I>& rhs) : 
	GenSymMatrix<T>(rhs), itsm(rhs.itsm), itss(rhs.itss), 
	itssi(rhs.itssi), itssj(rhs.itssj)
	  DEFFIRSTLAST(rhs.first,rhs.last) {}

      inline SymMatrixView(
	  T* _m, size_t _s, int _si, int _sj, 
	  SymType insym, UpLoType inuplo, StorageType instor, ConjItType inct 
	  PARAMFIRSTLAST(T) ) :
	GenSymMatrix<T>(insym,inuplo,instor,inct), itsm(_m), itss(_s),
	itssi(_si), itssj(_sj) DEFFIRSTLAST(_first,_last)
      {
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ?
	  _si==1 : true); 
      }

      inline ~SymMatrixView() {} 

      //
      // Op=
      //

      template <class T2> inline void DoCopy(
	  const GenSymMatrix<T2>& m2) const
      {
	TMVAssert(size() == m2.size());
	TMVAssert(IsComplex(T()) || IsReal(T2()));
	TMVAssert(IsReal(T2()) || m2.sym() == sym());
	UpperTri() = m2.UpperTri(); 
      }

      inline const SymMatrixView<T,I>& operator=(
	  const SymMatrixView<T,I>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsReal(T()) || m2.sym() == sym());
	if (!SameAs(m2)) DoCopy(m2); return *this; 
      }

      inline const SymMatrixView<T,I>& operator=(
	  const GenSymMatrix<T>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsReal(T()) || m2.sym() == sym());
	if (!SameAs(m2)) DoCopy(m2); return *this; 
      }

      template <class T2> inline const SymMatrixView<T,I>& operator=(
	  const GenSymMatrix<T2>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsComplex(T()) || IsReal(T2()));
	TMVAssert(IsReal(T2()) || m2.sym() == sym());
	DoCopy(m2); return *this; 
      }

      inline const SymMatrixView<T,I>& operator=(T x) const 
      { 
	TMVAssert(issym() || IMAG(x) == RealType(T)(0));
	return SetToIdentity(x); 
      }

      inline const SymMatrixView<T,I>& operator=(
	  const SymMatrixComposite<T>& mcomp) const
      { 
	TMVAssert(size() == mcomp.size());
	TMVAssert(IsReal(T()) || mcomp.sym() == sym());
	mcomp.AssignTo(View());
	return *this;
      }

      template <class Tm> inline const SymMatrixView<T,I>& operator=(
	  const SumSX<T,Tm>& mcomp) const
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(IsReal(Tm()) || issym() == mcomp.GetM().issym());
	TMVAssert(issym() || IMAG(mcomp.GetX1()) == RealType(T)(0));
	TMVAssert(issym() || IMAG(mcomp.GetX2()) == RealType(T)(0));
	mcomp.AssignTo(View());
	return *this;
      }

      template <class Tm> inline const SymMatrixView<T,I>& operator=(
	  const ProdXS<T,Tm>& mcomp) const
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(IsReal(Tm()) || issym() == mcomp.GetM().issym());
	TMVAssert(issym() || IMAG(mcomp.GetX()) == RealType(T)(0));
	mcomp.AssignTo(View());
	return *this;
      }

      template <class T1, class T2> inline const SymMatrixView<T,I>& operator=(
	  const SumSS<T,T1,T2>& mcomp) const
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(IsReal(T1()) || issym() == mcomp.GetM1().issym());
	TMVAssert(IsReal(T2()) || issym() == mcomp.GetM2().issym());
	TMVAssert(issym() || IMAG(mcomp.GetX1()) == RealType(T)(0));
	TMVAssert(issym() || IMAG(mcomp.GetX2()) == RealType(T)(0));
	mcomp.AssignTo(View());
	return *this;
      }

      template <class T2> inline const SymMatrixView<T,I>& operator=(
	  const OProdVV<T,T2,T2>& opvv) const
      {
	TMVAssert(size() == opvv.colsize());
	TMVAssert(size() == opvv.rowsize());
	TMVAssert(opvv.GetV1().SameAs(
	      isherm() ? opvv.GetV2().Conjugate() : opvv.GetV2().View()));
	Rank1Update(opvv.GetX(),opvv.GetV1(),0,*this);
	return *this;
      }

      template <class T2> inline const SymMatrixView<T,I>& operator=(
	  const ProdMM<T,T2,T2>& pmm) const
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(
	      isherm() ? pmm.GetM2().Adjoint() : pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,*this);
	return *this;
      }

      template <class T2> inline const SymMatrixView<T,I>& operator=(
	  const ProdLU<T,T2,T2>& pmm) const
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(
	      isherm() ? pmm.GetM2().Adjoint() : pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,*this);
	return *this;
      }

      template <class T2> inline const SymMatrixView<T,I>& operator=(
	  const ProdUL<T,T2,T2>& pmm) const
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(
	      isherm() ? pmm.GetM2().Adjoint() : pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,*this);
	return *this;
      }

      //
      // Access
      //

      inline RefType(T) operator()(size_t i,size_t j) const 
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	return ref(i,j); 
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert( i<=j1 || j2<=i+1 );
	if ((i<=j1 && uplo()==Upper) || (j2<=i+1 && uplo()==Lower))
	  return VectorView<T>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	      ct() FIRSTLAST); 
	else
	  return VectorView<T>(ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert( j<=i1 || i2<=j+1 );
	if ((i2<=j+1 && uplo()==Upper) || (j<=i1 && uplo()==Lower))
	  return VectorView<T>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	      ct() FIRSTLAST); 
	else 
	  return VectorView<T>(ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }

      inline VectorView<T> diag() const
      { return VectorView<T>(ptr(),size(),stepi()+stepj(),ct() FIRSTLAST); }

      inline VectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(size())); 
	if (i>=0)
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()+i*stepj(),size()-i,
		stepi()+stepj(),ct() FIRSTLAST); 
	  else
	    return VectorView<T>(ptr()+i*stepi(),size()-i,
		stepi()+stepj(),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
	else
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()-i*stepj(),size()+i,
		stepi()+stepj(),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
	  else
	    return VectorView<T>(ptr()-i*stepi(),size()+i,
		stepi()+stepj(),ct() FIRSTLAST); 
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()-abs(i));
	const int ds = stepi()+stepj();
	if (i>=0)
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()+i*stepj()+j1*ds,j2-j1,ds,ct() 
		FIRSTLAST); 
	  else
	    return VectorView<T>(ptr()+i*stepi()+j1*ds,j2-j1,ds, 
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
	else
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()-i*stepj()+j1*ds,j2-j1,ds, 
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
	  else
	    return VectorView<T>(ptr()-i*stepi()+j1*ds,j2-j1,ds,ct() 
		FIRSTLAST); 
      }

      //
      // Modifying Functions
      //

      inline const SymMatrixView<T,I>& Zero() const 
      { return SetAllTo(T(0)); }

      inline const SymMatrixView<T,I>& SetAllTo(T x) const
      { 
	TMVAssert(IsReal(x) || issym());
	UpperTri().SetAllTo(x); return *this; 
      }

      inline const SymMatrixView<T,I>& Clip(RealType(T) thresh) const
      { UpperTri().Clip(thresh); return *this; }

      inline const SymMatrixView<T,I>& ConjugateSelf() const
      { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

      inline const SymMatrixView<T,I>& TransposeSelf() const
      { if (!issym()) UpperTri().ConjugateSelf(); return *this; }

      inline const SymMatrixView<T,I>& SetToIdentity(T x=T(1)) const
      { 
	TMVAssert(IsReal(T()) || issym());
	Zero(); diag().SetAllTo(x); return *this; 
      }

      const SymMatrixView<T,I>& SwapRowsCols(size_t i1, size_t i2) const;

      const SymMatrixView<T,I>& PermuteRowsCols(const size_t* p,
	  size_t i1, size_t i2) const;

      const SymMatrixView<T,I>& ReversePermuteRowsCols(const size_t* p,
	  size_t i1, size_t i2) const;

      inline const SymMatrixView<T,I>& PermuteRowsCols(const size_t* p) const
      { return PermuteRowsCols(p,0,size()); }

      inline const SymMatrixView<T,I>& ReversePermuteRowsCols(
	  const size_t* p) const
      { return ReversePermuteRowsCols(p,0,size()); }


      //
      // SubMatrix
      //

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(GenSymMatrix<T>::OKSubMatrix(i1,i2,j1,j2,1,1));
	TMVAssert(i2<=j1 || j2<=i1);
	if ((i2<=j1 && uplo()==Upper) || (j2<=i1 && uplo()==Lower))
	  return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	      i2-i1,j2-j1,stepi(),stepj(),stor(),ct() FIRSTLAST);
	else
	  return MatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
	      i2-i1,j2-j1,stepj(),stepi(),TransOf(stor()),
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(GenSymMatrix<T>::OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	TMVAssert(i2<=j1 || j2<=i1);
	if ((i2<=j1 && uplo()==Upper) || (j2<=i1 && uplo()==Lower))
	  return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	      (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	      newstor,ct() FIRSTLAST);
	else
	  return MatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
	      (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
	      TransOf(newstor),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t n) const
      {
	TMVAssert(GenSymMatrix<T>::OKSubVector(i,j,istep,jstep,n));
	if ((i<=j && uplo()==Upper) || (j<=i && uplo()==Lower))
	  return VectorView<T>(ptr()+i*stepi()+j*stepj(),n,
	      istep*stepi()+jstep*stepj(),ct() FIRSTLAST);
	else
	  return VectorView<T>(ptr()+i*stepj()+j*stepi(),n,
	      istep*stepj()+jstep*stepi(),issym()?ct():ConjOf(T,ct()) 
	      FIRSTLAST);
      }

      inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
      {
	TMVAssert(GenSymMatrix<T>::OKSubSymMatrix(i1,i2,1));
	return SymMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),i2-i1,
	    stepi(),stepj(),sym(),uplo(),stor(),ct() FIRSTLAST);
      }

      inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(GenSymMatrix<T>::OKSubSymMatrix(i1,i2,istep));
	return SymMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
	    istep*stepi(),istep*stepj(),sym(),uplo(),
	    istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T> UpperTri(DiagType dt = NonUnitDiag) const
      {
	if (uplo() == Upper)
	  return UpperTriMatrixView<T>(ptr(),size(),
	      stepi(),stepj(),dt,stor(),ct() FIRSTLAST);
	else
	  return UpperTriMatrixView<T>(ptr(),size(),
	      stepj(),stepi(),dt,TransOf(stor()),
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> LowerTri(DiagType dt = NonUnitDiag) const
      {
	if (uplo() == Lower)
	  return LowerTriMatrixView<T>(ptr(),size(),
	      stepi(),stepj(),dt,stor(),ct() FIRSTLAST);
	else
	  return LowerTriMatrixView<T>(ptr(),size(),
	      stepj(),stepi(),dt,TransOf(stor()),
	      issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }

      inline SymMatrixView<T,I> View() const
      { return *this; }

      inline SymMatrixView<T,I> Transpose() const
      {
	return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	    sym(),UTransOf(uplo()),TransOf(stor()),ct() FIRSTLAST);
      }

      inline SymMatrixView<T,I> Conjugate() const
      {
	return SymMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	    sym(),uplo(),stor(),ConjOf(T,ct()) FIRSTLAST);
      }

      inline SymMatrixView<T,I> Adjoint() const
      {
	return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	    sym(),UTransOf(uplo()),TransOf(stor()),ConjOf(T,ct()) 
	    FIRSTLAST);
      }

      inline SymMatrixView<RealType(T),I> Real() const
      {
	return SymMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    sym(),uplo(), IsReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline SymMatrixView<RealType(T),I> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(issym());
	return SymMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    size(),2*stepi(),2*stepj(),sym(),uplo(),NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // I/O
      //

      void Read(istream& is) const;

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      using GenSymMatrix<T>::sym;
      using GenSymMatrix<T>::uplo;
      using GenSymMatrix<T>::stor;
      using GenSymMatrix<T>::ct;
      using GenSymMatrix<T>::isrm;
      using GenSymMatrix<T>::iscm;
      using GenSymMatrix<T>::isconj;
      using GenSymMatrix<T>::isherm;
      using GenSymMatrix<T>::issym;

      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> inline void operator=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }
      template <class T2> inline void operator+=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }
      template <class T2> inline void operator-=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }
      template <class T2> inline void operator*=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }
      template <class T2> inline void operator/=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }
      template <class T2> inline void operator%=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }

    protected :

      T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      RefType(T) ref(size_t i, size_t j) const;

  }; // SymMatrixView

  template <class T> class SymMatrixView<T,FortranStyle> : 
    public SymMatrixView<T,CStyle>
    {

      public:

	//
	// Constructors
	//

	inline SymMatrixView(const SymMatrixView<T,FortranStyle>& rhs) : 
	  SymMatrixView<T,CStyle>(rhs) {}

	inline SymMatrixView(const SymMatrixView<T,CStyle>& rhs) : 
	  SymMatrixView<T,CStyle>(rhs) {}

	inline SymMatrixView(
	    T* _m, size_t _s, int _si, int _sj, 
	    SymType insym, UpLoType inuplo, StorageType instor, ConjItType inct 
	    PARAMFIRSTLAST(T) ) :
	  SymMatrixView<T,CStyle>(_m,_s,_si,_sj,insym,inuplo,instor,inct
	      FIRSTLAST1(_first,_last) ) {}

	inline ~SymMatrixView() {} 

	//
	// Op=
	//

	inline const SymMatrixView<T,FortranStyle>& operator=(
	    const SymMatrixView<T,FortranStyle>& m2) const
	{ SymMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const SymMatrixView<T,FortranStyle>& operator=(
	    const GenSymMatrix<T>& m2) const
	{ SymMatrixView<T,CStyle>::operator=(m2); return *this; }

	template <class T2> 
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const GenSymMatrix<T2>& m2) const
	  { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const SymMatrixView<T,FortranStyle>& operator=(T x) const 
	{ SymMatrixView<T,CStyle>::operator=(x); return *this; }

	inline const SymMatrixView<T,FortranStyle>& operator=(
	    const SymMatrixComposite<T>& mcomp) const
	{ SymMatrixView<T,CStyle>::operator=(mcomp); return *this; }

	template <class Tm>
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const SumSX<T,Tm>& mcomp) const
	  { SymMatrixView<T,CStyle>::operator=(mcomp); return *this; }

	template <class Tm>
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const ProdXS<T,Tm>& mcomp) const
	  { SymMatrixView<T,CStyle>::operator=(mcomp); return *this; }

	template <class T1, class T2>
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const SumSS<T,T1,T2>& mcomp) const
	  { SymMatrixView<T,CStyle>::operator=(mcomp); return *this; }

	template <class T2> 
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const OProdVV<T,T2,T2>& opvv) const
	  { SymMatrixView<T,CStyle>::operator=(opvv); return *this; }

	template <class T2> 
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const ProdMM<T,T2,T2>& pmm) const
	  { SymMatrixView<T,CStyle>::operator=(pmm); return *this; }

	template <class T2> 
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const ProdLU<T,T2,T2>& pmm) const
	  { SymMatrixView<T,CStyle>::operator=(pmm); return *this; }

	template <class T2> 
	  inline const SymMatrixView<T,FortranStyle>& operator=(
	      const ProdUL<T,T2,T2>& pmm) const
	  { SymMatrixView<T,CStyle>::operator=(pmm); return *this; }

	//
	// Access
	//

	inline RefType(T) operator()(size_t i,size_t j) const 
	{
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  return ref(i-1,j-1); 
	}

	inline VectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const 
	{ 
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j1>0 && j1<=j2);
	  TMVAssert(j2<=size());
	  return SymMatrixView<T,CStyle>::row(i-1,j1-1,j2);
	}

	inline VectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=size());
	  TMVAssert(i1>0 && i1<=i2 && i2<=size());
	  return SymMatrixView<T,CStyle>::col(j-1,i1-1,i2);
	}

	inline VectorView<T,FortranStyle> diag() const
	{ return SymMatrixView<T,CStyle>::diag(); }

	inline VectorView<T,FortranStyle> diag(int i) const
	{ return SymMatrixView<T,CStyle>::diag(i); }

	inline VectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(j1>0);
	  return SymMatrixView<T,CStyle>::diag(i,j1-1,j2); 
	}

	//
	// Modifying Functions
	//

	inline const SymMatrixView<T,FortranStyle>& Zero() const 
	{ SymMatrixView<T,CStyle>::Zero(); return *this; }

	inline const SymMatrixView<T,FortranStyle>& SetAllTo(T x) const
	{ SymMatrixView<T,CStyle>::SetAllTo(x); return *this; }

	inline const SymMatrixView<T,FortranStyle>& Clip(
	    RealType(T) thresh) const
	{ SymMatrixView<T,CStyle>::Clip(thresh); return *this; }

	inline const SymMatrixView<T,FortranStyle>& ConjugateSelf() const
	{ SymMatrixView<T,CStyle>::ConjugateSelf(); return *this; }

	inline const SymMatrixView<T,FortranStyle>& TransposeSelf() const
	{ SymMatrixView<T,CStyle>::TransposeSelf(); return *this; }

	inline const SymMatrixView<T,FortranStyle>& SetToIdentity(
	    T x=T(1)) const
	{ SymMatrixView<T,CStyle>::SetToIdentity(x); return *this; }

	inline const SymMatrixView<T,FortranStyle>& SwapRowsCols(
	    size_t i1, size_t i2) const
	{
	  TMVAssert(i1>0 && i1<=size());
	  TMVAssert(i2>0 && i2<=size());
	  SymMatrixView<T,CStyle>::SwapRowsCols(i1-1,i2-1); 
	  return *this; 
	}

	inline const SymMatrixView<T,FortranStyle>& PermuteRowsCols(
	    const size_t* p, size_t i1, size_t i2) const
	{
	  TMVAssert(i1>0 && i1<=i2 && i2<=size());
	  SymMatrixView<T,CStyle>::PermuteRowsCols(p,i1-1,i2); 
	  return *this; 
	}

	inline const SymMatrixView<T,FortranStyle>& ReversePermuteRowsCols(
	    const size_t* p, size_t i1, size_t i2) const
	{
	  TMVAssert(i1>0 && i1<=i2 && i2<=size());
	  SymMatrixView<T,CStyle>::ReversePermuteRowsCols(p,i1-1,i2); 
	  return *this; 
	}

	inline const SymMatrixView<T,FortranStyle>& PermuteRowsCols(
	    const size_t* p) const
	{ SymMatrixView<T,CStyle>::PermuteRowsCols(p); return *this; }

	inline const SymMatrixView<T,FortranStyle>& ReversePermuteRowsCols(
	    const size_t* p) const
	{ SymMatrixView<T,CStyle>::ReversePermuteRowsCols(p); return *this; }

	//
	// SubMatrix
	//

	inline bool OKSubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  return ConstSymMatrixView<T,FortranStyle>(*this).OKSubMatrix(
	      i1,i2,j1,j2,istep,jstep); 
	}

	inline bool OKSubVector(
	    size_t i, size_t j, int istep, int jstep, size_t n) const
	{
	  return ConstSymMatrixView<T,FortranStyle>(*this).OKSubVector(
	      i,j,istep,jstep,n); 
	}


	inline bool OKSubSymMatrix(int i1, int i2, int istep) const
	{
	  return ConstSymMatrixView<T,FortranStyle>(*this).OKSubSymMatrix(
	      i1,i2,istep); 
	}

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	  return SymMatrixView<T,CStyle>::SubMatrix(i1-1,i2,j1-1,j2);
	}

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return SymMatrixView<T,CStyle>::SubMatrix(i1-1,i2-1+istep,
	      j1-1,j2-1+jstep,istep,jstep);
	}

	inline VectorView<T,FortranStyle> SubVector(size_t i, size_t j,
	    int istep, int jstep, size_t n) const
	{
	  TMVAssert(OKSubVector(i,j,istep,jstep,n));
	  return SymMatrixView<T,CStyle>::SubVector(i-1,j-1,istep,jstep,n);
	}

	inline SymMatrixView<T,FortranStyle> SubSymMatrix(int i1, int i2) const
	{
	  TMVAssert(OKSubSymMatrix(i1,i2,1));
	  return SymMatrixView<T,CStyle>::SubSymMatrix(i1-1,i2);
	}

	inline SymMatrixView<T,FortranStyle> SubSymMatrix(int i1, int i2,
	    int istep) const
	{
	  TMVAssert(OKSubSymMatrix(i1,i2,istep));
	  return SymMatrixView<T,CStyle>::SubSymMatrix(i1-1,i2-1+istep,istep);
	}

	inline UpperTriMatrixView<T,FortranStyle> UpperTri(
	    DiagType dt = NonUnitDiag) const
	{ return SymMatrixView<T,CStyle>::UpperTri(dt); }

	inline LowerTriMatrixView<T,FortranStyle> LowerTri(
	    DiagType dt = NonUnitDiag) const
	{ return SymMatrixView<T,CStyle>::LowerTri(dt); }

	inline SymMatrixView<T,FortranStyle> View() const
	{ return SymMatrixView<T,CStyle>::View(); }

	inline SymMatrixView<T,FortranStyle> Transpose() const
	{ return SymMatrixView<T,CStyle>::Transpose(); }

	inline SymMatrixView<T,FortranStyle> Conjugate() const
	{ return SymMatrixView<T,CStyle>::Conjugate(); }

	inline SymMatrixView<T,FortranStyle> Adjoint() const
	{ return SymMatrixView<T,CStyle>::Adjoint(); }

	inline SymMatrixView<RealType(T),FortranStyle> Real() const
	{ return SymMatrixView<T,CStyle>::Real(); }

	inline SymMatrixView<RealType(T),FortranStyle> Imag() const
	{ return SymMatrixView<T,CStyle>::Imag(); }

	using SymMatrixView<T,CStyle>::size;

	template <class T2> inline void operator=(const BaseMatrix<T2>&) const 
	{ TMVAssert(FALSE); }

      protected :

	using SymMatrixView<T,CStyle>::ref;

    }; // FortranStyle SymMatrixView

  template <class T, UpLoType U, StorageType S, IndexStyle I> class SymMatrix : 
    public GenSymMatrix<T>
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(s) GenSymMatrix<T>(Sym,U,S,NonConj), itslen((s)*(s)), \
      itsm(new T[itslen]), itss(s)  \
      DEFFIRSTLAST(itsm.get(),itsm.get()+itslen) 

      explicit inline SymMatrix(size_t _size) : NEW_SIZE(_size) 
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      inline SymMatrix(size_t _size, const T x) : NEW_SIZE(_size) 
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	SetAllTo(x);
      }

      inline SymMatrix(size_t _size, const valarray<T>& vv) : NEW_SIZE(_size)
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(vv.size() == itslen);
	T* vi = itsm.get();
	for(size_t i=0;i<itslen;++i) *vi=vv[i];
      }

      inline SymMatrix(size_t _size, const T* vv) : NEW_SIZE(_size)
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	memmove(itsm.get(),vv,itslen*sizeof(T));
      }

      inline SymMatrix(size_t _size, const vector<T>& vv) : NEW_SIZE(_size)
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(vv.size() == itslen);
	T* vi = itsm.get();
	typename vector<T>::const_iterator vvi = vv.begin();
	for(size_t i=itslen;i>0;--i) *vi = *vvi;
      }

      inline SymMatrix(const SymMatrix<T,U,S,I>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	memmove(itsm.get(),rhs.itsm.get(),itslen*sizeof(T));
      }

      template <UpLoType U2, StorageType S2, IndexStyle I2> inline SymMatrix(
	  const SymMatrix<T,U2,S2,I2>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	if ((U==U2 && S==S2) ||
	    ((S2==RowMajor || S2==ColMajor) && S!=S2 && U!=U2))
	  memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	else
	  if (U == Upper)
	    UpperTri() = rhs.UpperTri();
	  else
	    LowerTri() = rhs.LowerTri();
      }

      inline SymMatrix(const GenSymMatrix<T>& rhs) :
	NEW_SIZE(rhs.size())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.issym());
	if (U == Upper)
	  UpperTri() = rhs.UpperTri();
	else
	  LowerTri() = rhs.LowerTri();
      }

      template <class T2> inline SymMatrix(const GenSymMatrix<T2>& rhs) :
	NEW_SIZE(rhs.size())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.issym());
	if (U == Upper)
	  UpperTri() = rhs.UpperTri();
	else
	  LowerTri() = rhs.LowerTri();
      }

      template <StorageType S2, IndexStyle I2> inline SymMatrix(
	  const Matrix<T,S2,I2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.IsSquare());
	if (S2==S)
	  memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	else 
	  if (U == Upper)
	    UpperTri() = UpperTriMatrixViewOf(rhs);
	  else
	    LowerTri() = LowerTriMatrixViewOf(rhs);
      }

      inline SymMatrix(const GenMatrix<T>& rhs) :
	NEW_SIZE(rhs.rowsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.IsSquare());
	if (U == Upper)
	  UpperTri() = UpperTriMatrixViewOf(rhs);
	else
	  LowerTri() = LowerTriMatrixViewOf(rhs);
      }

      template <class T2> inline SymMatrix(const GenMatrix<T2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.IsSquare());
	if (U == Upper)
	  UpperTri() = UpperTriMatrixViewOf(rhs);
	else
	  LowerTri() = LowerTriMatrixViewOf(rhs);
      }

      inline SymMatrix(const SymMatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.size())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.issym());
	mcomp.AssignTo(View());
      }

      template <class Tm> inline SymMatrix(const SumSX<T,Tm>& mcomp) :
	NEW_SIZE(mcomp.colsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.GetM().issym());
	mcomp.AssignTo(View());
      }

      template <class Tm> inline SymMatrix(const ProdXS<T,Tm>& mcomp) :
	NEW_SIZE(mcomp.colsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.GetM().issym());
	mcomp.AssignTo(View());
      }

      template <class T1, class T2> inline SymMatrix(
	  const SumSS<T,T1,T2>& mcomp) :
	NEW_SIZE(mcomp.colsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.GetM1().issym());
	TMVAssert(mcomp.GetM2().issym());
	mcomp.AssignTo(View());
      }

      template <class T2> inline SymMatrix(const OProdVV<T,T2,T2>& opvv) :
	NEW_SIZE(opvv.colsize())
      {
	TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().View()));
	Rank1Update(opvv.GetX(),opvv.GetV1(),0,View());
      }

      template <class T2> inline SymMatrix(const ProdMM<T,T2,T2>& pmm) :
	NEW_SIZE(pmm.colsize())
      {
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
      }

      template <class T2> inline SymMatrix(const ProdLU<T,T2,T2>& pmm) :
	NEW_SIZE(pmm.colsize())
      {
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
      }

      template <class T2> inline SymMatrix(const ProdUL<T,T2,T2>& pmm) :
	NEW_SIZE(pmm.colsize())
      {
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
      }

#undef NEW_SIZE

      inline ~SymMatrix() {}


      //
      // Op=
      //

      inline SymMatrix<T,U,S,I>& operator=(const SymMatrix<T,U,S,I>& m2)
      { 
	TMVAssert(size() == m2.size());
	if (&m2 != this) 
	  memmove(itsm.get(),m2.itsm.get(),size()*size()*sizeof(T));
	return *this;
      }

      inline SymMatrix<T,U,S,I>& operator=(
	  const GenSymMatrix<T>& m2)
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(m2.issym());
	UpperTri() = m2.UpperTri(); 
	return *this; 
      }

      template <class T2> inline SymMatrix<T,U,S,I>& operator=(
	  const GenSymMatrix<T2>& m2)
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsReal(T2()) || m2.issym());
	UpperTri() = m2.UpperTri(); 
	return *this; 
      }

      inline SymMatrix<T,U,S,I>& operator=(T x) 
      { return SetToIdentity(x); }

      inline SymMatrix<T,U,S,I>& operator=(const SymMatrixComposite<T>& mcomp)
      { 
	TMVAssert(size() == mcomp.size());
	TMVAssert(mcomp.issym());
	mcomp.AssignTo(View());
	return *this;
      }

      template <class Tm> inline SymMatrix<T,U,S,I>& operator=(
	  const SumSX<T,Tm>& mcomp) 
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(mcomp.GetM().issym());
	mcomp.AssignTo(View());
	return *this;
      }

      template <class Tm> inline SymMatrix<T,U,S,I>& operator=(
	  const ProdXS<T,Tm>& mcomp) 
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(mcomp.GetM().issym());
	mcomp.AssignTo(View());
	return *this;
      }

      template <class T1, class T2> inline SymMatrix<T,U,S,I>& operator=(
	  const SumSS<T,T1,T2>& mcomp) 
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(mcomp.GetM1().issym());
	TMVAssert(mcomp.GetM2().issym());
	mcomp.AssignTo(View());
	return *this;
      }

      template <class T2> inline SymMatrix<T,U,S,I>& operator=(
	  const OProdVV<T,T2,T2>& opvv)
      {
	TMVAssert(size() == opvv.colsize());
	TMVAssert(size() == opvv.rowsize());
	TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().View()));
	Rank1Update(opvv.GetX(),opvv.GetV1(),0,View());
	return *this;
      }

      template <class T2> inline SymMatrix<T,U,S,I>& operator=(
	  const ProdMM<T,T2,T2>& pmm)
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
	return *this;
      }

      template <class T2> inline SymMatrix<T,U,S,I>& operator=(
	  const ProdLU<T,T2,T2>& pmm)
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
	return *this;
      }

      template <class T2> inline SymMatrix<T,U,S,I>& operator=(
	  const ProdUL<T,T2,T2>& pmm)
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
	return *this;
      }

      //
      // Access
      //

      inline T operator()(size_t i, size_t j) const
      { 
	if (I==CStyle) {
	  TMVAssert(i<size());
	  TMVAssert(j<size());
	  return cref(i,j); 
	} else {
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  return cref(i-1,j-1); 
	}
      }

      inline T& operator()(size_t i, size_t j) 
      { 
	if (I==CStyle) {
	  TMVAssert(i<size());
	  TMVAssert(j<size());
	  return ref(i,j); 
	} else {
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  return ref(i-1,j-1); 
	}
      }

      inline ConstVectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
      { 
	if (I==FortranStyle) { TMVAssert(i>0 && j1>0); }
	const size_t ix = (I == CStyle ? i : i-1);
	const size_t j1x = (I == CStyle ? j1 : j1-1);
	TMVAssert(ix<size());
	TMVAssert(j1x<=j2);
	TMVAssert(j2<=size());
	TMVAssert( ix<=j1x || j2<=ix+1 );
	if ((U==Upper && ix<=j1x) || (U==Lower && j2<=ix+1))
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj); 
	else
	  return ConstVectorView<T,I>(cptr()+ix*stepj()+j1x*stepi(),
	      j2-j1x,stepi(),NonConj); 
      }

      inline ConstVectorView<T,I> col(size_t j, size_t i1, size_t i2) const
      {
	if (I==FortranStyle) { TMVAssert(j>0 && i1>0); }
	const size_t jx = (I == CStyle ? j : j-1);
	const size_t i1x = (I == CStyle ? i1 : i1-1);
	TMVAssert(jx<size());
	TMVAssert(i1x<=i2);
	TMVAssert(i2<=size());
	TMVAssert( jx<=i1x || i2<=jx+1 );
	if ((U==Upper && i2<=jx+1) || (U==Lower && jx<=i1x))
	  return ConstVectorView<T,I>(cptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj); 
	else 
	  return ConstVectorView<T,I>(cptr()+i1x*stepj()+jx*stepi(),
	      i2-i1x,stepj(),NonConj); 
      }

      inline ConstVectorView<T,I> diag() const
      { return ConstVectorView<T,I>(cptr(),size(),stepi()+stepj(),NonConj); }

      inline ConstVectorView<T,I> diag(int i) const
      {
	i = abs(i);
	TMVAssert(i<=int(size())); 
	return ConstVectorView<T,I>(cptr()+i*(U==Upper?stepj():stepi()),
	    size()-i,stepi()+stepj(),NonConj); 
      }

      inline ConstVectorView<T,I> diag(int i, size_t j1, size_t j2) const
      {
	if (I==FortranStyle) { TMVAssert(j1>0); }
	const size_t j1x = (I == CStyle ? j1 : j1-1);
	i = abs(i);
	TMVAssert(j1x <= j2);
	TMVAssert(j2 <= size()-i);
	const int ds = stepi()+stepj();
	return ConstVectorView<T,I>(cptr()+i*(U==Upper?stepj():stepi())+j1x*ds,
	    j2-j1x,ds,NonConj);
      }

      inline VectorView<T,I> row(size_t i, size_t j1, size_t j2)
      { 
	if (I==FortranStyle) { TMVAssert(i>0 && j1>0); }
	const size_t ix = (I == CStyle ? i : i-1);
	const size_t j1x = (I == CStyle ? j1 : j1-1);
	TMVAssert(ix<size());
	TMVAssert(j1x<=j2);
	TMVAssert(j2<=size());
	TMVAssert( ix<=j1x || j2<=ix+1 );
	if ((U==Upper && ix<=j1x) || (U==Lower && j2<=ix+1))
	  return VectorView<T,I>(ptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj FIRSTLAST); 
	else
	  return VectorView<T,I>(ptr()+ix*stepj()+j1x*stepi(),
	      j2-j1x,stepi(),NonConj FIRSTLAST); 
      }

      inline VectorView<T,I> col(size_t j, size_t i1, size_t i2)
      {
	if (I==FortranStyle) { TMVAssert(j>0 && i1>0); }
	const size_t jx = (I == CStyle ? j : j-1);
	const size_t i1x = (I == CStyle ? i1 : i1-1);
	TMVAssert(jx<size());
	TMVAssert(i1x<=i2);
	TMVAssert(i2<=size());
	TMVAssert( jx<=i1x || i2<=jx+1 );
	if ((U==Upper && i2<=jx+1) || (U==Lower && jx<=i1x))
	  return VectorView<T,I>(ptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj FIRSTLAST); 
	else 
	  return VectorView<T,I>(ptr()+i1x*stepj()+jx*stepi(),
	      i2-i1x,stepj(),NonConj FIRSTLAST); 
      }

      inline VectorView<T,I> diag()
      {
	return VectorView<T,I>(ptr(),size(),stepi()+stepj(),NonConj 
	    FIRSTLAST);
      }

      inline VectorView<T,I> diag(int i) 
      {
	i = abs(i);
	TMVAssert(i<=int(size())); 
	return VectorView<T,I>(ptr()+i*(U==Upper?stepj():stepi()),
	    size()-i,stepi()+stepj(),NonConj FIRSTLAST); 
      }

      inline VectorView<T,I> diag(int i, size_t j1, size_t j2) 
      {
	if (I==FortranStyle) { TMVAssert(j1>0); }
	const size_t j1x = (I == CStyle ? j1 : j1-1);
	i = abs(i);
	TMVAssert(j1x <= j2);
	TMVAssert(j2 <= size()-abs(i));
	const int ds = stepi()+stepj();
	return VectorView<T,I>(ptr()+i*(U==Upper?stepj():stepi())+j1x*ds,
	    j2-j1x,ds,NonConj FIRSTLAST);
      }


      //
      // Modifying Functions
      //

      inline SymMatrix<T,U,S,I>& Zero() 
      { return SetAllTo(0); }

      inline SymMatrix<T,U,S,I>& SetAllTo(T x) 
      { UpperTri().SetAllTo(x); return *this; }

      inline SymMatrix<T,U,S,I>& Clip(RealType(T) thresh) 
      { UpperTri().Clip(thresh); return *this; }

      inline SymMatrix<T,U,S,I>& ConjugateSelf() 
      { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

      inline SymMatrix<T,U,S,I>& TransposeSelf() 
      { return *this; }

      inline SymMatrix<T,U,S,I>& SetToIdentity(T x=T(1)) 
      { Zero(); diag().SetAllTo(x); return *this; }

      inline SymMatrix<T,U,S,I>& SwapRowsCols(size_t i1, size_t i2)
      { View().SwapRowsCols(i1,i2); return *this; }

      inline SymMatrix<T,U,S,I>& PermuteRowsCols(const size_t* p,
	  size_t i1, size_t i2)
      { View().PermuteRowsCols(p,i1,i2); return *this; }

      inline SymMatrix<T,U,S,I>& ReversePermuteRowsCols(const size_t* p,
	  size_t i1, size_t i2)
      { View().ReversePermuteRowsCols(p,i1,i2); return *this; }

      inline SymMatrix<T,U,S,I>& PermuteRowsCols(const size_t* p)
      { View().PermuteRowsCols(p); return *this; }

      inline SymMatrix<T,U,S,I>& ReversePermuteRowsCols(const size_t* p)
      { View().ReversePermuteRowsCols(p); return *this; }


      //
      // SubMatrix
      //

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	TMVAssert(i2<=j1x || j2<=i1x);
	if ((U==Upper && i2<=j1x) || (U==Lower && j2<=i1x))
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x,j2-j1x,stepi(),stepj(),S,NonConj);
	else
	  return ConstMatrixView<T,I>(cptr()+i1x*stepj()+j1x*stepi(),
	      i2-i1x,j2-j1x,stepj(),stepi(),TransOf(S),NonConj);
      }

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	const int j2x = (I==CStyle ? j2 : j2-1+jstep); 
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  istep == 1 ? ColMajor : NoMajor;
	TMVAssert(i2x<=j1x || j2x<=i1x);
	if ((U==Upper && i2x<=j1x) || (U==Lower && j2x<=i1x))
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	      newstor,NonConj);
	else
	  return ConstMatrixView<T,I>(cptr()+i1x*stepj()+j1x*stepi(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepj(),jstep*stepi(),
	      TransOf(newstor),NonConj);
      }

      inline ConstVectorView<T,I> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t n) const
      {
	TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
	const int ix = (I==CStyle ? i : i-1); 
	const int jx = (I==CStyle ? j : j-1); 
	if ((U==Upper && ix<=jx) || (U==Lower && jx<=ix))
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+jx*stepj(),n,
	      istep*stepi()+jstep*stepj(),NonConj);
	else
	  return ConstVectorView<T,I>(cptr()+ix*stepj()+jx*stepi(),n,
	      istep*stepj()+jstep*stepi(),NonConj);
      }

      inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	return ConstSymMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),i2-i1x,
	    stepi(),stepj(),Sym,U,S,NonConj);
      }

      inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	return ConstSymMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	    (i2x-i1x)/istep,istep*stepi(),istep*stepj(),Sym,U,
	    istep==1 ? S : NoMajor, NonConj);
      }

      inline ConstUpperTriMatrixView<T,I> UpperTri(
	  DiagType dt = NonUnitDiag) const
      {
	return U==Upper ? 
	  ConstUpperTriMatrixView<T,I>(cptr(),size(),
	    stepi(),stepj(),dt,S,NonConj) :
	  ConstUpperTriMatrixView<T,I>(cptr(),size(),
	    stepj(),stepi(),dt,TransOf(S),NonConj);
      }

      inline ConstLowerTriMatrixView<T,I> LowerTri(
	  DiagType dt = NonUnitDiag) const
      {
	return U==Lower ? 
	  ConstLowerTriMatrixView<T,I>(cptr(),size(),
	    stepi(),stepj(),dt,S,NonConj) :
	  ConstLowerTriMatrixView<T,I>(cptr(),size(),
	    stepj(),stepi(),dt,TransOf(S),NonConj);
      }

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	TMVAssert(i2<=j1x || j2<=i1x);
	if ((U==Upper && i2<=j1x) || (U==Lower && j2<=i1x))
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x,j2-j1x,stepi(),stepj(),S,NonConj FIRSTLAST);
	else
	  return MatrixView<T,I>(ptr()+i1x*stepj()+j1x*stepi(),
	      i2-i1x,j2-j1x,stepj(),stepi(),TransOf(S),NonConj FIRSTLAST);
      }

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep)
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  istep == 1 ? ColMajor : NoMajor;
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	const int j2x = (I==CStyle ? j2 : j2-1+jstep); 
	TMVAssert(i2x<=j1x || j2x<=i1x);
	if ((U==Upper && i2x<=j1x) || (U==Lower && j2x<=i1x))
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	      newstor,NonConj FIRSTLAST);
	else
	  return MatrixView<T,I>(ptr()+i1x*stepj()+j1x*stepi(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepj(),jstep*stepi(),
	      TransOf(newstor),NonConj FIRSTLAST);
      }

      inline VectorView<T,I> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t n)
      {
	TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
	const size_t ix = (I==CStyle ? i : i-1); 
	const size_t jx = (I==CStyle ? j : j-1); 
	if ((U==Upper && ix<=jx) || (U==Lower && jx<=ix))
	  return VectorView<T,I>(ptr()+ix*stepi()+jx*stepj(),n,
	      istep*stepi()+jstep*stepj(),NonConj FIRSTLAST);
	else
	  return VectorView<T,I>(ptr()+ix*stepj()+jx*stepi(),n,
	      istep*stepj()+jstep*stepi(),NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2)
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	return SymMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),i2-i1x,
	    stepi(),stepj(),Sym,U,S,NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) 
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	return SymMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),(i2x-i1x)/istep,
	    istep*stepi(),istep*stepj(),Sym,U,
	    istep==1 ? S : NoMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> UpperTri(DiagType dt = NonUnitDiag)
      {
	return U==Upper ? 
	  UpperTriMatrixView<T,I>(ptr(),size(),
	    stepi(),stepj(),dt,S,NonConj FIRSTLAST) :
	  UpperTriMatrixView<T,I>(ptr(),size(),
	    stepj(),stepi(),dt,TransOf(S),NonConj FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> LowerTri(DiagType dt = NonUnitDiag)
      {
	return U==Lower ? 
	  LowerTriMatrixView<T,I>(ptr(),size(),
	    stepi(),stepj(),dt,S,NonConj FIRSTLAST) :
	  LowerTriMatrixView<T,I>(ptr(),size(),
	    stepj(),stepi(),dt,TransOf(S),NonConj FIRSTLAST);
      }

      inline ConstSymMatrixView<T,I> View() const
      { 
	return ConstSymMatrixView<T,I>(cptr(),size(),stepi(),stepj(),
	    Sym,U,S,NonConj);
      }

      inline ConstSymMatrixView<T,I> Transpose() const
      {
	return ConstSymMatrixView<T,I>(cptr(),size(),stepj(),stepi(),
	    Sym,UTransOf(U),TransOf(S),NonConj);
      }

      inline ConstSymMatrixView<T,I> Conjugate() const
      { 
	return ConstSymMatrixView<T,I>(cptr(),size(),stepi(),stepj(),
	    Sym,U,S,ConjOf(T,NonConj));
      }

      inline ConstSymMatrixView<T,I> Adjoint() const
      {
	return ConstSymMatrixView<T,I>(cptr(),size(),stepj(),stepi(),
	    Sym,UTransOf(U),TransOf(S),ConjOf(T,NonConj));
      }

      inline ConstSymMatrixView<RealType(T),I> Real() const
      {
	return ConstSymMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    Sym,U,IsReal(T())?S:NoMajor,NonConj);
      }

      inline ConstSymMatrixView<RealType(T),I> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return ConstSymMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,size(),
	    2*stepi(),2*stepj(),Sym,U,NoMajor,NonConj);
      } 

      inline SymMatrixView<T,I> View() 
      { 
	return SymMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	    Sym,U,S,NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> Transpose() 
      {
	return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	    Sym,UTransOf(U),TransOf(S),NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> Conjugate() 
      { 
	return SymMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	    Sym,U,S,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline SymMatrixView<T,I> Adjoint() 
      {
	return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	    Sym,UTransOf(U),TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
      }

      inline SymMatrixView<RealType(T),I> Real()
      {
	return SymMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    Sym,U,IsReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline SymMatrixView<RealType(T),I> Imag()
      {
	TMVAssert(IsComplex(T()));
	return SymMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,size(),
	    2*stepi(),2*stepj(),Sym,U,NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      } 

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm.get(); }
      inline T* ptr() { return itsm.get(); }
      inline int stepi() const { return S==RowMajor ? itss : 1; }
      inline int stepj() const { return S==RowMajor ? 1 : itss; }
      inline SymType sym() const { return Sym; }
      inline UpLoType uplo() const { return U; }
      inline StorageType stor() const { return S; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isrm() const { return S==RowMajor; }
      inline bool iscm() const { return S==ColMajor; }
      inline bool isconj() const { return false; }
      inline bool isherm() const { return IsReal(T()); }
      inline bool issym() const { return true; }

    protected :

      const size_t itslen;
      auto_array<T> itsm;
      const size_t itss;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
#endif

      inline T& ref(size_t i, size_t j)
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
	  return *(ptr() + i*stepi() + j*stepj());
	else 
	  return *(ptr() + j*stepi() + i*stepj());
      }

      inline T cref(size_t i, size_t j) const 
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
	  return *(cptr() + i*stepi() + j*stepj());
	else 
	  return *(cptr() + j*stepi() + i*stepj());
      }

  }; // SymMatrix

  template <class T, UpLoType U, StorageType S, IndexStyle I> class HermMatrix : 
    public GenSymMatrix<T>
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(s) GenSymMatrix<T>(Herm,U,S,NonConj), itslen((s)*(s)), \
      itsm(new T[itslen]), itss(s)  \
      DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

      explicit inline HermMatrix(size_t _size) : NEW_SIZE(_size) 
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
	SetAllTo(RealType(T)(888));
#endif
      }

      inline HermMatrix(size_t _size, const RealType(T) x) : NEW_SIZE(_size) 
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	SetAllTo(x);
      }

      inline HermMatrix(size_t _size, const valarray<T>& vv) : NEW_SIZE(_size)
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(vv.size() == itslen);
	T* vi = itsm.get();
	for(size_t i=0;i<itslen;++i) *vi=vv[i];
      }

      inline HermMatrix(size_t _size, const T* vv) : NEW_SIZE(_size)
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	memmove(itsm.get(),vv,itslen*sizeof(T));
      }

      inline HermMatrix(size_t _size, const vector<T>& vv) : NEW_SIZE(_size)
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(vv.size() == itslen);
	T* vi = itsm.get();
	typename vector<T>::const_iterator vvi = vv.begin();
	for(size_t i=itslen;i>0;--i) *vi = *vvi;
      }

      inline HermMatrix(const HermMatrix<T,U,S,I>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	memmove(itsm.get(),rhs.itsm.get(),itslen*sizeof(T));
      }

      template <UpLoType U2, StorageType S2, IndexStyle I2> inline HermMatrix(
	  const HermMatrix<T,U2,S2,I2>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	if ((U==U2 && S==S2) ||
	    (IsReal(T()) && (S2==RowMajor || S2==ColMajor) && S!=S2 && U!=U2))
	  memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	else
	  if (U==Upper)
	    UpperTri() = rhs.UpperTri();
	  else
	    LowerTri() = rhs.LowerTri();
      }

      inline HermMatrix(const GenSymMatrix<T>& rhs) :
	NEW_SIZE(rhs.size())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	if (U == Upper)
	  UpperTri() = rhs.UpperTri();
	else
	  LowerTri() = rhs.LowerTri();
	if (IsComplex(T()) && rhs.issym()) 
	  diag().Imag().Zero();
      }

      template <class T2> inline HermMatrix(const GenSymMatrix<T2>& rhs) :
	NEW_SIZE(rhs.size())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	if (U == Upper)
	  UpperTri() = rhs.UpperTri();
	else
	  LowerTri() = rhs.LowerTri();
	if (IsComplex(T()) && IsComplex(T2()) && rhs.issym()) 
	  diag().Imag().Zero();
      }

      template <StorageType S2, IndexStyle I2> inline HermMatrix(
	  const Matrix<T,S2,I2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.IsSquare());
	if (S2==S)
	  memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	else if (U == Upper)
	  UpperTri() = UpperTriMatrixViewOf(rhs);
	else 
	  LowerTri() = LowerTriMatrixViewOf(rhs);
	if (IsComplex(T())) diag().Imag().Zero();
      }

      inline HermMatrix(const GenMatrix<T>& rhs) :
	NEW_SIZE(rhs.rowsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.IsSquare());
	if (U == Upper)
	  UpperTri() = UpperTriMatrixViewOf(rhs);
	else
	  LowerTri() = LowerTriMatrixViewOf(rhs);
	if (IsComplex(T())) diag().Imag().Zero();
      }

      template <class T2> inline HermMatrix(const GenMatrix<T2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(rhs.IsSquare());
	if (U == Upper)
	  UpperTri() = UpperTriMatrixViewOf(rhs);
	else
	  LowerTri() = LowerTriMatrixViewOf(rhs);
	if (IsComplex(T())) diag().Imag().Zero();
      }

      inline HermMatrix(const SymMatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.size())
      {
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.isherm());
	mcomp.AssignTo(View());
      }

      template <class Tm> inline HermMatrix(const SumSX<T,Tm>& mcomp) :
	NEW_SIZE(mcomp.colsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.GetM().isherm());
	TMVAssert(IMAG(mcomp.GetX1()) == RealType(T)(0));
	TMVAssert(IMAG(mcomp.GetX2()) == RealType(T)(0));
	mcomp.AssignTo(View());
      }

      template <class Tm> inline HermMatrix(const ProdXS<T,Tm>& mcomp) :
	NEW_SIZE(mcomp.colsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.GetM().isherm());
	TMVAssert(IMAG(mcomp.GetX()) == RealType(T)(0));
	mcomp.AssignTo(View());
      }

      template <class T1, class T2> inline HermMatrix(
	  const SumSS<T,T1,T2>& mcomp) :
	NEW_SIZE(mcomp.colsize())
      { 
	TMVAssert(S==RowMajor || S==ColMajor);
	TMVAssert(mcomp.GetM1().isherm());
	TMVAssert(mcomp.GetM2().isherm());
	TMVAssert(IMAG(mcomp.GetX1()) == RealType(T)(0));
	TMVAssert(IMAG(mcomp.GetX2()) == RealType(T)(0));
	mcomp.AssignTo(View());
      }

      template <class T2> inline HermMatrix(const OProdVV<T,T2,T2>& opvv) :
	NEW_SIZE(opvv.colsize())
      {
	TMVAssert(size() == opvv.colsize());
	TMVAssert(size() == opvv.rowsize());
	TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().Conjugate()));
	Rank1Update(opvv.GetX(),opvv.GetV1(),0,View());
      }

      template <class T2> inline HermMatrix(const ProdMM<T,T2,T2>& pmm) :
	NEW_SIZE(pmm.colsize())
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
      }

      template <class T2> inline HermMatrix(const ProdLU<T,T2,T2>& pmm) :
	NEW_SIZE(pmm.colsize())
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
      }

      template <class T2> inline HermMatrix(const ProdUL<T,T2,T2>& pmm) :
	NEW_SIZE(pmm.colsize())
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
      }

#undef NEW_SIZE

      inline ~HermMatrix() {}


      //
      // Op=
      //

      inline HermMatrix<T,U,S,I>& operator=(const HermMatrix<T,U,S,I>& m2)
      { 
	TMVAssert(m2.size() == size());
	if (&m2 != this) 
	  memmove(itsm.get(),m2.itsm.get(),size()*size()*sizeof(T));
	return *this;
      }

      inline HermMatrix<T,U,S,I>& operator=(const GenSymMatrix<T>& m2)
      {
	TMVAssert(m2.size() == size());
	TMVAssert(m2.isherm());
	if (U == Upper)
	  UpperTri() = m2.UpperTri(); 
	else
	  LowerTri() = m2.LowerTri(); 
	return *this; 
      }

      template <class T2> inline HermMatrix<T,U,S,I>& operator=(
	  const GenSymMatrix<T2>& m2)
      {
	TMVAssert(m2.size() == size());
	TMVAssert(IsReal(T2()) || m2.isherm());
	if (U == Upper)
	  UpperTri() = m2.UpperTri(); 
	else
	  LowerTri() = m2.LowerTri(); 
	return *this; 
      }

      inline HermMatrix<T,U,S,I>& operator=(T x) 
      {
	TMVAssert(IMAG(x) == RealType(T)(0));
	return SetToIdentity(x); 
      }

      inline HermMatrix<T,U,S,I>& operator=(const SymMatrixComposite<T>& mcomp)
      { 
	TMVAssert(size() == mcomp.size());
	TMVAssert(mcomp.isherm());
	mcomp.AssignTo(View());
	return *this;
      }

      template <class Tm> inline HermMatrix<T,U,S,I>& operator=(
	  const SumSX<T,Tm>& mcomp) 
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(mcomp.GetM().isherm());
	TMVAssert(IMAG(mcomp.GetX1()) == RealType(T)(0));
	TMVAssert(IMAG(mcomp.GetX2()) == RealType(T)(0));
	mcomp.AssignTo(View());
	return *this;
      }

      template <class Tm> inline HermMatrix<T,U,S,I>& operator=(
	  const ProdXS<T,Tm>& mcomp) 
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(mcomp.GetM().isherm());
	TMVAssert(IMAG(mcomp.GetX()) == RealType(T)(0));
	mcomp.AssignTo(View());
	return *this;
      }

      template <class T1, class T2> inline HermMatrix<T,U,S,I>& operator=(
	  const SumSS<T,T1,T2>& mcomp) 
      { 
	TMVAssert(size() == mcomp.colsize());
	TMVAssert(mcomp.GetM1().isherm());
	TMVAssert(mcomp.GetM2().isherm());
	TMVAssert(IMAG(mcomp.GetX1()) == RealType(T)(0));
	TMVAssert(IMAG(mcomp.GetX2()) == RealType(T)(0));
	mcomp.AssignTo(View());
	return *this;
      }

      template <class T2> inline HermMatrix<T,U,S,I>& operator=(
	  const OProdVV<T,T2,T2>& opvv)
      {
	TMVAssert(size() == opvv.colsize());
	TMVAssert(size() == opvv.rowsize());
	TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().Conjugate()));
	Rank1Update(opvv.GetX(),opvv.GetV1(),0,View());
	return *this;
      }

      template <class T2> inline HermMatrix<T,U,S,I>& operator=(
	  const ProdMM<T,T2,T2>& pmm)
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
	return *this;
      }

      template <class T2> inline HermMatrix<T,U,S,I>& operator=(
	  const ProdLU<T,T2,T2>& pmm)
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
	return *this;
      }

      template <class T2> inline HermMatrix<T,U,S,I>& operator=(
	  const ProdUL<T,T2,T2>& pmm)
      {
	TMVAssert(size() == pmm.colsize());
	TMVAssert(size() == pmm.rowsize());
	TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
	RankKUpdate(pmm.GetX(),pmm.GetM1(),0,View());
	return *this;
      }

      //
      // Access
      //

      inline T operator()(size_t i, size_t j) const
      {
	if (I==CStyle) {
	  TMVAssert(i<size());
	  TMVAssert(j<size());
	  return cref(i,j); 
	} else {
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  return cref(i-1,j-1); 
	}
      }

      inline RefType(T) operator()(size_t i, size_t j) 
      { 
	if (I==CStyle) {
	  TMVAssert(i<size());
	  TMVAssert(j<size());
	  return ref(i,j); 
	} else {
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  return ref(i-1,j-1); 
	}
      }

      inline ConstVectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
      { 
	if (I == FortranStyle) { TMVAssert(i>0 && j1>0); }
	const size_t ix = (I==CStyle ? i : i-1);
	const size_t j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(ix<size());
	TMVAssert(j1x<=j2 && j2<=size());
	TMVAssert( ix<=j1x || j2<=ix+1 );
	if ((U==Upper && ix<=j1x) || (U==Lower && j2<=ix+1))
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj); 
	else
	  return ConstVectorView<T,I>(cptr()+ix*stepj()+j1x*stepi(),
	      j2-j1x,stepi(),ConjOf(T,NonConj)); 
      }

      inline ConstVectorView<T,I> col(size_t j, size_t i1, size_t i2) const
      {
	if (I == FortranStyle) { TMVAssert(j>0 && i1>0); }
	const size_t jx = (I==CStyle ? j : j-1);
	const size_t i1x = (I==CStyle ? i1 : i1-1);
	TMVAssert(jx<size());
	TMVAssert(i1x<=i2 && i2<=size());
	TMVAssert( jx<=i1x || i2<=jx+1 );
	if ((U==Upper && i2<=jx+1) || (U==Lower && jx<=i1x))
	  return ConstVectorView<T,I>(cptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj); 
	else 
	  return ConstVectorView<T,I>(cptr()+i1x*stepj()+jx*stepi(),
	      i2-i1x,stepj(),ConjOf(T,NonConj)); 
      }

      inline ConstVectorView<T,I> diag() const
      { return ConstVectorView<T,I>(cptr(),size(),stepi()+stepj(),NonConj); }

      inline ConstVectorView<T,I> diag(int i) const
      {
	ConjItType newct = 
	  ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
	i = abs(i);
	TMVAssert(i<=int(size())); 
	return ConstVectorView<T,I>(cptr()+i*(U==Upper?stepj():stepi()),
	    size()-i,stepi()+stepj(),newct);
      }

      inline ConstVectorView<T,I> diag(int i, size_t j1, size_t j2) const
      {
	if (I == FortranStyle) { TMVAssert(j1>0); }
	const size_t j1x = (I==CStyle ? j1 : j1-1); 
	ConjItType newct = 
	  ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
	i = abs(i);
	TMVAssert(j1x <= j2);
	TMVAssert(j2 <= size()-i);
	const int ds = stepi()+stepj();
	return ConstVectorView<T,I>(cptr()+i*(U==Upper?stepj():stepi())+j1x*ds,
	    j2-j1x,ds,newct);
      }

      inline VectorView<T,I> row(size_t i, size_t j1, size_t j2)
      { 
	if (I == FortranStyle) { TMVAssert(i>0 && j1>0); }
	const size_t ix = (I==CStyle ? i : i-1);
	const size_t j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(ix<size());
	TMVAssert(j1x<=j2 && j2<=size());
	TMVAssert( ix<=j1x || j2<=ix+1 );
	if ((U==Upper && ix<=j1x) || (U==Lower && j2<=ix+1))
	  return VectorView<T,I>(ptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj FIRSTLAST); 
	else
	  return VectorView<T,I>(ptr()+ix*stepj()+j1x*stepi(),
	      j2-j1x,stepi(),ConjOf(T,NonConj) FIRSTLAST); 
      }

      inline VectorView<T,I> col(size_t j, size_t i1, size_t i2)
      {
	if (I == FortranStyle) { TMVAssert(j>0 && i1>0); }
	const size_t jx = (I==CStyle ? j : j-1);
	const size_t i1x = (I==CStyle ? i1 : i1-1);
	TMVAssert(jx<size());
	TMVAssert(i1x<=i2 && i2<=size());
	TMVAssert( jx<=i1x || i2<=jx+1 );
	if ((U==Upper && i2<=jx+1) || (U==Lower && jx<=i1x))
	  return VectorView<T,I>(ptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj FIRSTLAST); 
	else 
	  return VectorView<T,I>(ptr()+i1x*stepj()+jx*stepi(),
	      i2-i1x,stepj(),ConjOf(T,NonConj) FIRSTLAST); 
      }

      inline VectorView<T,I> diag()
      {
	return VectorView<T,I>(ptr(),size(),stepi()+stepj(),NonConj FIRSTLAST);
      }

      inline VectorView<T,I> diag(int i) 
      {
	ConjItType newct = 
	  ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
	i = abs(i);
	TMVAssert(i<=int(size())); 
	return VectorView<T,I>(ptr()+i*(U==Upper?stepj():stepi()),size()-i,
	    stepi()+stepj(),newct FIRSTLAST); 
      }

      inline VectorView<T,I> diag(int i, size_t j1, size_t j2) 
      {
	if (I == FortranStyle) { TMVAssert(j1>0); }
	const size_t j1x = (I==CStyle ? j1 : j1-1); 
	ConjItType newct = 
	  ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
	i = abs(i);
	TMVAssert(j1x <= j2);
	TMVAssert(j2 <= size()-abs(i));
	const int ds = stepi()+stepj();
	return VectorView<T,I>(ptr()+i*(U==Upper?stepj():stepi())+j1x*ds,
	    j2-j1x,ds,newct FIRSTLAST);
      }

      //
      // Modifying Functions
      //

      inline HermMatrix<T,U,S,I>& Zero() 
      { return SetAllTo(0); }

      inline HermMatrix<T,U,S,I>& SetAllTo(T x) 
      { 
	TMVAssert(IMAG(x) == RealType(T)(0));
	UpperTri().SetAllTo(x); 
	return *this; 
      }

      inline HermMatrix<T,U,S,I>& Clip(RealType(T) thresh) 
      { UpperTri().Clip(thresh); return *this; }

      inline HermMatrix<T,U,S,I>& ConjugateSelf() 
      { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

      inline HermMatrix<T,U,S,I>& TransposeSelf() 
      { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

      inline HermMatrix<T,U,S,I>& SetToIdentity(T x=T(1)) 
      { 
	TMVAssert(IMAG(x) == RealType(T)(0));
	Zero(); 
	diag().SetAllTo(x); 
	return *this; 
      }

      inline HermMatrix<T,U,S,I>& SwapRowsCols(size_t i1, size_t i2)
      { View().SwapRowsCols(i1,i2); return *this; }

      inline HermMatrix<T,U,S,I>& PermuteRowsCols(
	  const size_t* p, size_t i1, size_t i2)
      { View().PermuteRowsCols(p,i1,i2); return *this; }

      inline HermMatrix<T,U,S,I>& ReversePermuteRowsCols(
	  const size_t* p, size_t i1, size_t i2)
      { View().ReversePermuteRowsCols(p,i1,i2); return *this; }

      inline HermMatrix<T,U,S,I>& PermuteRowsCols(const size_t* p)
      { View().PermuteRowsCols(p); return *this; }

      inline HermMatrix<T,U,S,I>& ReversePermuteRowsCols(const size_t* p)
      { View().ReversePermuteRowsCols(p); return *this; }


      //
      // SubMatrix
      //

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	TMVAssert(i2<=j1x || j2<=i1x);
	if ((U==Upper && i2<=j1x) || (U==Lower && j2<=i1x))
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x,j2-j1x,stepi(),stepj(),S,NonConj);
	else
	  return ConstMatrixView<T,I>(cptr()+i1x*stepj()+j1x*stepi(),
	      i2-i1x,j2-j1x,stepj(),stepi(),TransOf(S),ConjOf(T,NonConj));
      }

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  istep == 1 ? ColMajor : NoMajor;
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	const int j2x = (I==CStyle ? j2 : j2-1+jstep); 
	TMVAssert(i2x<=j1x || j2x<=i1x);
	if ((U==Upper && i2x<=j1x) || (U==Lower && j2x<=i1x))
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	      newstor,NonConj);
	else
	  return ConstMatrixView<T,I>(cptr()+i1x*stepj()+j1x*stepi(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepj(),jstep*stepi(),
	      TransOf(newstor),ConjOf(T,NonConj));
      }

      inline ConstVectorView<T,I> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t n) const
      {
	TMVAssert(GenSymMatrix<T>::OKSubVector(i,j,istep,jstep,n));
	if ((U==Upper && i<=j) || (U==Lower && j<=i))
	  return ConstVectorView<T,I>(cptr()+i*stepi()+j*stepj(),n,
	      istep*stepi()+jstep*stepj(),NonConj);
	else
	  return ConstVectorView<T,I>(cptr()+i*stepj()+j*stepi(),n,
	      istep*stepj()+jstep*stepi(),ConjOf(T,NonConj));
      }

      inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	return ConstSymMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	    i2-i1x,stepi(),stepj(),Herm,U,S,NonConj);
      }

      inline ConstSymMatrixView<T,I> SubSymMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	return ConstSymMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	    (i2x-i1x)/istep,istep*stepi(),istep*stepj(),Herm,U,
	    istep==1 ? S : NoMajor, NonConj);
      }

      inline ConstUpperTriMatrixView<T,I> UpperTri(
	  DiagType dt = NonUnitDiag) const
      {
	return U==Upper ? 
	  ConstUpperTriMatrixView<T,I>(cptr(),size(),stepi(),stepj(),
	      dt,S,NonConj) :
	  ConstUpperTriMatrixView<T,I>(cptr(),size(),stepj(),stepi(),
	      dt,TransOf(S),ConjOf(T,NonConj));
      }

      inline ConstLowerTriMatrixView<T,I> LowerTri(
	  DiagType dt = NonUnitDiag) const
      {
	return U==Lower ? 
	  ConstLowerTriMatrixView<T,I>(cptr(),size(),stepi(),stepj(),
	      dt,S,NonConj) :
	  ConstLowerTriMatrixView<T,I>(cptr(),size(),stepj(),stepi(),
	      dt,TransOf(S),ConjOf(T,NonConj));
      }

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	TMVAssert(i2<=j1x || j2<=i1x);
	if ((U==Upper && i2<=j1x) || (U==Lower && j2<=i1x))
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x,j2-j1x,stepi(),stepj(),S,NonConj FIRSTLAST);
	else
	  return MatrixView<T,I>(ptr()+i1x*stepj()+j1x*stepi(),
	      i2-i1x,j2-j1x,stepj(),stepi(),TransOf(S),ConjOf(T,NonConj) 
	      FIRSTLAST);
      }

      inline MatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep)
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  istep == 1 ? ColMajor : NoMajor;
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int j1x = (I==CStyle ? j1 : j1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	const int j2x = (I==CStyle ? j2 : j2-1+jstep); 
	TMVAssert(i2x<=j1x || j2x<=i1x);
	if ((U==Upper && i2x<=j1x) || (U==Lower && j2x<=i1x))
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	      newstor,NonConj FIRSTLAST);
	else
	  return MatrixView<T,I>(ptr()+i1x*stepj()+j1x*stepi(),
	      (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepj(),jstep*stepi(),
	      TransOf(newstor),ConjOf(T,NonConj) FIRSTLAST);
      }

      inline VectorView<T,I> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t n)
      {
	TMVAssert(GenSymMatrix<T>::OKSubVector(i,j,istep,jstep,n));
	if ((U==Upper && i<=j) || (U==Lower && j<=i))
	  return VectorView<T,I>(ptr()+i*stepi()+j*stepj(),n,
	      istep*stepi()+jstep*stepj(),NonConj FIRSTLAST);
	else
	  return VectorView<T,I>(ptr()+i*stepj()+j*stepi(),n,
	      istep*stepj()+jstep*stepi(),ConjOf(T,NonConj) FIRSTLAST);
      }

      inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2)
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,1));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	return SymMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),
	    i2-i1x,stepi(),stepj(),Herm,U,S,NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) 
      {
	TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
	const int i1x = (I==CStyle ? i1 : i1-1); 
	const int i2x = (I==CStyle ? i2 : i2-1+istep); 
	return SymMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),
	    (i2x-i1x)/istep,istep*stepi(),istep*stepj(),Herm,U,
	    istep==1 ? S : NoMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> UpperTri(DiagType dt = NonUnitDiag)
      {
	return U==Upper ? 
	  UpperTriMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	      dt,S,NonConj FIRSTLAST) :
	  UpperTriMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	      dt,TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> LowerTri(DiagType dt = NonUnitDiag)
      {
	return U==Lower ? 
	  LowerTriMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	      dt,S,NonConj FIRSTLAST) :
	  LowerTriMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	      dt,TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
      }

      inline ConstSymMatrixView<T,I> View() const
      { 
	return ConstSymMatrixView<T,I>(cptr(),size(),stepi(),stepj(),
	    Herm,U,S,NonConj);
      }

      inline ConstSymMatrixView<T,I> Transpose() const
      {
	return ConstSymMatrixView<T,I>(cptr(),size(),stepj(),stepi(),
	    Herm,UTransOf(U),TransOf(S),NonConj);
      }

      inline ConstSymMatrixView<T,I> Conjugate() const
      { 
	return ConstSymMatrixView<T,I>(cptr(),size(),stepi(),stepj(),
	    Herm,U,S,ConjOf(T,NonConj));
      }

      inline ConstSymMatrixView<T,I> Adjoint() const
      {
	return ConstSymMatrixView<T,I>(cptr(),size(),stepj(),stepi(),
	    Herm,UTransOf(U),TransOf(S),ConjOf(T,NonConj));
      }

      inline ConstSymMatrixView<RealType(T),I> Real() const
      {
	return ConstSymMatrixView<RealType(T),I>(
	    reinterpret_cast<const RealType(T)*>(cptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    Herm,U,IsReal(T())?S:NoMajor,NonConj);
      }

      inline ConstSymMatrixView<RealType(T),I> Imag() const
      // The imaginary part of a Hermitian matrix is anti-symmetric
      // so this is illegal.
      { 
	TMVAssert(FALSE);
	return ConstSymMatrixView<RealType(T),I>(0,0,0,0,Herm,U,S,NonConj);
      }

      inline SymMatrixView<T,I> View() 
      { 
	return SymMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	    Herm,U,S,NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> Transpose() 
      {
	return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	    Herm,UTransOf(U),TransOf(S),NonConj FIRSTLAST);
      }

      inline SymMatrixView<T,I> Conjugate() 
      { 
	return SymMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
	    Herm,U,S,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline SymMatrixView<T,I> Adjoint() 
      {
	return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
	    Herm,UTransOf(U),TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
      }

      inline SymMatrixView<RealType(T),I> Real()
      {
	return SymMatrixView<RealType(T),I>(
	    reinterpret_cast<RealType(T)*>(ptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    Herm,U,IsReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline SymMatrixView<RealType(T),I> Imag()
      // The imaginary part of a Hermitian matrix is anti-symmetric
      // so this is illegal.
      { 
	TMVAssert(FALSE);
	return SymMatrixView<RealType(T),I>(0,0,0,0,Herm,U,S,NonConj
	    FIRSTLAST1(0,0) );
      }

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm.get(); }
      inline T* ptr() { return itsm.get(); }
      inline int stepi() const { return S==RowMajor ? itss : 1; }
      inline int stepj() const { return S==RowMajor ? 1 : itss; }
      inline SymType sym() const { return Herm; }
      inline UpLoType uplo() const { return U; }
      inline StorageType stor() const { return S; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isrm() const { return S==RowMajor; }
      inline bool iscm() const { return S==ColMajor; }
      inline bool isconj() const { return false; }
      inline bool isherm() const { return true; }
      inline bool issym() const { return IsReal(T()); }

    protected :

      const size_t itslen;
      auto_array<T> itsm;
      const size_t itss;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
#endif

      inline RefType(T) ref(size_t i, size_t j)
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
	  return REF(ptr() + i*stepi() + j*stepj(),NonConj);
	else 
	  return REF(ptr() + j*stepi() + i*stepj(),Conj);
      }

      inline T cref(size_t i, size_t j) const 
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
	  return *(cptr() + i*stepi() + j*stepj());
	else 
	  return CONJ(*(cptr() + j*stepi() + i*stepj()));
      }

  }; // HermMatrix

//---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   SymMatrixViewOf(m,uplo)
  //   HermMatrixViewOf(m,uplo)
  //   SymMatrixViewOf(mptr,size,uplo,stor)
  //   HermMatrixViewOf(mptr,size,uplo,stor)
  //

  template <class T> inline ConstSymMatrixView<T> SymMatrixViewOf(
      const GenMatrix<T>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymMatrixView<T>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
	Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T,I> SymMatrixViewOf(
	const ConstMatrixView<T,I>& m, UpLoType uplo)
    { 
      TMVAssert(m.colsize()==m.rowsize());
      return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
	  Sym,uplo,m.stor(),m.ct()); 
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> SymMatrixViewOf(
	const Matrix<T,S,I>& m, UpLoType uplo)
    {
      TMVAssert(m.colsize()==m.rowsize());
      return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
	  Sym,uplo,m.stor(),m.ct()); 
    }

  template <class T, IndexStyle I> inline SymMatrixView<T,I> SymMatrixViewOf(
      const MatrixView<T,I>& m, UpLoType uplo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
    return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
	Sym,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
  }

  template <class T, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> SymMatrixViewOf(
	Matrix<T,S,I>& m, UpLoType uplo)
    {
      TMVAssert(m.colsize()==m.rowsize());
      return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
	  Sym,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
    }

  template <class T> inline ConstSymMatrixView<T> HermMatrixViewOf(
      const GenMatrix<T>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
    return ConstSymMatrixView<T>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
	Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T> HermMatrixViewOf(
	const ConstMatrixView<T,I>& m, UpLoType uplo)
    { 
      TMVAssert(m.colsize()==m.rowsize());
      TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
      return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
	  Herm,uplo,m.stor(),m.ct()); 
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> HermMatrixViewOf(
	const Matrix<T,S,I>& m, UpLoType uplo)
    {
      TMVAssert(m.colsize()==m.rowsize());
      TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
      return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
	  Herm,uplo,m.stor(),m.ct()); 
    }

  template <class T, IndexStyle I> inline SymMatrixView<T> HermMatrixViewOf(
      const MatrixView<T,I>& m, UpLoType uplo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
    return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
	Herm,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
  }

  template <class T, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> HermMatrixViewOf(
	Matrix<T,S,I>& m, UpLoType uplo)
    {
      TMVAssert(m.colsize()==m.rowsize());
      TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
      return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
	  Herm,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
    }

  template <class T> inline ConstSymMatrixView<T> SymMatrixViewOf(
      const T* m, size_t size, UpLoType uplo, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstSymMatrixView<T>(m,size,size,1,Sym,uplo,RowMajor,NonConj);
    else
      return ConstSymMatrixView<T>(m,size,1,size,Sym,uplo,ColMajor,NonConj);
  }

  template <class T> inline ConstSymMatrixView<T> HermMatrixViewOf(
      const T* m, size_t size, UpLoType uplo, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    ConstSymMatrixView<T> ret = (stor == RowMajor ?
	ConstSymMatrixView<T>(m,size,size,1,Herm,uplo,RowMajor,NonConj) :
	ConstSymMatrixView<T>(m,size,1,size,Herm,uplo,ColMajor,NonConj) );
    TMVAssert(IsReal(T()) || ret.diag().Imag().NormInf() == RealType(T)(0));
    return ret;
  }

  template <class T> inline SymMatrixView<T> SymMatrixViewOf(
      T* m, size_t size, UpLoType uplo, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return SymMatrixView<T>(m,size,size,1,Sym,uplo,RowMajor,NonConj
	  FIRSTLAST1(m,m+size*size));
    else
      return SymMatrixView<T>(m,size,1,size,Sym,uplo,ColMajor,NonConj
	  FIRSTLAST1(m,m+size*size));
  }

  template <class T> inline SymMatrixView<T> HermMatrixViewOf(
      T* m, size_t size, UpLoType uplo, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    SymMatrixView<T> ret = (stor == RowMajor ?
	SymMatrixView<T>(m,size,size,1,Herm,uplo,RowMajor,NonConj
	  FIRSTLAST1(m,m+size*size)) :
	SymMatrixView<T>(m,size,1,size,Herm,uplo,ColMajor,NonConj
	  FIRSTLAST1(m,m+size*size)) );
    TMVAssert(IsReal(T()) || ret.diag().Imag().NormInf() == RealType(T)(0));
    return ret;
  }

  //
  // Swap Matrices
  //

  template <class T> inline void Swap(
      const SymMatrixView<T>& m1, const SymMatrixView<T>& m2)
  { Swap(m1.UpperTri(),m2.UpperTri()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const GenSymMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const GenSymMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const GenSymMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormSq(const GenSymMatrix<T>& m)
  { return m.NormSq(); }

  template <class T> inline RealType(T) NormF(const GenSymMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm1(const GenSymMatrix<T>& m)
  { return m.Norm1(); }

  template <class T> inline RealType(T) Norm2(const GenSymMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const GenSymMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> inline RealType(T) MaxAbsElement(const GenSymMatrix<T>& m)
  { return m.MaxAbsElement(); }

  template <class T> inline ConstSymMatrixView<T> Transpose(
      const GenSymMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> inline ConstSymMatrixView<T,I> Transpose(
      const ConstSymMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Transpose(const SymMatrix<T,U,S,I>& m)
    { return m.Transpose(); }

  template <class T, IndexStyle I> inline SymMatrixView<T,I> Transpose(
      const SymMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> Transpose(SymMatrix<T,U,S,I>& m)
    { return m.Transpose(); }

  template <class T> inline ConstSymMatrixView<T> Conjugate(
      const GenSymMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> inline ConstSymMatrixView<T,I> Conjugate(
      const ConstSymMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Conjugate(const SymMatrix<T,U,S,I>& m)
    { return m.Conjugate(); }

  template <class T, IndexStyle I> inline SymMatrixView<T,I> Conjugate(
      const SymMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> Conjugate(SymMatrix<T,U,S,I>& m)
    { return m.Conjugate(); }

  template <class T> inline ConstSymMatrixView<T> Adjoint(
      const GenSymMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline ConstSymMatrixView<T,I> Adjoint(
      const ConstSymMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Adjoint(const SymMatrix<T,U,S,I>& m)
    { return m.Adjoint(); }

  template <class T, IndexStyle I> inline SymMatrixView<T,I> Adjoint(
      const SymMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> Adjoint(SymMatrix<T,U,S,I>& m)
    { return m.Adjoint(); }

  template <class T> inline QuotXS<T,T> Inverse(const GenSymMatrix<T>& m)
  { return m.Inverse(); }

  //
  // SymMatrix ==, != SymMatrix
  //

  template <class T1, class T2> inline bool operator==(
      const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
  { return m1.UpperTri() == m2.UpperTri(); }
  template <class T1, class T2> inline bool operator!=(
      const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    istream& operator>>(istream& is, auto_ptr<SymMatrix<T,U,S,I> >& m);

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    istream& operator>>(istream& is, auto_ptr<HermMatrix<T,U,S,I> >& m);

  template <class T> istream& operator>>(
      istream& is, const SymMatrixView<T>& m);

  template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline istream& operator>>(istream& is, SymMatrix<T,U,S,I>& m)
    { return is>>m.View(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline istream& operator>>(istream& is, HermMatrix<T,U,S,I>& m)
    { return is>>m.View(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline string Type(const SymMatrix<T,U,S,I>& )
    { 
      return string("SymMatrix<")+Type(T())+","+Text(U)+","+Text(S)+
	","+Text(I)+">"; 
    }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline string Type(const HermMatrix<T,U,S,I>& )
    {
      return string("HermMatrix<")+Type(T())+","+Text(U)+","+Text(S)+
	","+Text(I)+">"; 
    }

  template <class T> inline string Type(const GenSymMatrix<T>& m)
  { 
    return string("GenSymMatrix<")+Type(T())+","+Text(m.sym())+","
      +Text(m.uplo())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }

  template <class T, IndexStyle I> inline string Type(
      const ConstSymMatrixView<T,I>& m)
  {
    return string("ConstSymMatrixView<")+Type(T())+","+Text(m.sym())+
      ","+Text(m.uplo())+ ","+Text(m.stor())+ ","+Text(I)+
      ","+Text(m.ct())+ ">"; 
  }

  template <class T, IndexStyle I> inline string Type(
      const SymMatrixView<T,I>& m)
  {
    return string("SymMatrixView<")+Type(T())+","+Text(m.sym())+
      ","+Text(m.uplo())+ ","+Text(m.stor())+ ","+Text(I)+
      ","+Text(m.ct())+ ">"; 
  }

} // namespace tmv

#endif
