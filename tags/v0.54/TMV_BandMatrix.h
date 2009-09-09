//---------------------------------------------------------------------------
//
// This file defines the TMV BandMatrix class.
//
// Constructors:
//
//    BandMatrices have 1 template parameter: T = the type of data.
//    The default is double if it is omitted.
//
//    A BandMatrix is only non-zero in a small number of diagonals around
//    the main diagonal.  Specifically, we store nhi super-diagonals 
//    (above the main diagonal) and nlo sub-diagonals (below the main).
//    
//    All constructors have the storage pattern as the last parameter.
//    RowMajor and ColMajor are as for a normal Matrix.
//    But there is a new possibility for BandMatrices: DiagMajor.
//    In this case the storage is in order of each diagonal starting 
//    with the lowest one, proceding along that whole diagonal, and then 
//    the next diagonal up, and so on.
//
//    Also, for each storage possibility, we store some extra elements 
//    in order to have the rows, columns, and diagonals all have constant
//    steps.  For example, a 6x6 ColMajor BandMatrix with nlo=2,nhi=3 has 
//    the following elements, numbered in order that they are stored.
//
//    [ 1   6  11  16         ]
//    [ 2   7  12  17  22     ]
//    [ 3   8  13  18  23  28 ]
//    [     9  14  19  24  29 ]
//    [        15  20  25  30 ]
//    [            21  26  31 ]
//
//    We do not ever use elements stored at the memory locations 4-5, 10, 27.
//    We need to skip those in order to get the rows and diagonals to have 
//    constant steps in the memory.  If we had started the second column
//    with memory location 4, and the third at 8, then the first row would 
//    be 1 4 8 13 which cannot be referenced by a VectorView.
//
//    Likewise, the same BandMatrix in RowMajor and DiagMajor storage
//    is as follows:
//
//    [  1   2   3   4         ]  [ 11  17  23  29         ]
//    [  6   7   8   9  10     ]  [  6  12  18  24  30     ]
//    [ 11  12  13  14  15  16 ]  [  1   7  13  19  25  31 ]
//    [     17  18  19  20  21 ]  [      2   8  14  20  26 ]
//    [         23  24  25  26 ]  [          3   9  15  21 ]
//    [             29  30  31 ]  [              4  10  16 ]
//
//    For square BandMatrices, the wasted storage is only 
//    (nlo-1)*nlo/2 + (nhi-1)*nhi/2 memory locations,
//    which, if nlo and nhi are small compared to the size N,
//    is negligible compared to the total memory allocated of 
//    (N-1)*(nlo+nhi+1)+1.
//    
//    Also, we don't actually require that the BandMatrix be square,
//    although that is the usual case.  Hopefully, the extension
//    of the above formats to non-square cases is obvious. 
//    The memory required for non-square BandMatrices is not quite as
//    simple a formula, and it depends on the StorageType.
//    The required memory (in units of sizeof(T)) can be obtained
//    by the function:
//
//    size_t BandStorageLength(stor, colsize, rowsize, nlo, nhi);
//
//
//    BandMatrix<T,stor>(colsize, rowsize, nlo, nhi)
//        Makes a BandMatrix with column size = row size = size 
//        with nhi non-zero superdiagonals and nlo non-zero subdiagonals
//        with _uninitialized_ values
//
//    BandMatrix<T,stor>(colsize, rowsize, nlo, nhi, x)
//        The same as above, but all values are initialized to x
//
//    BandMatrix<T,stor>(const Matrix<T>& m, nlo, nhi)
//        Makes a BandMatrix which copies the corresponding elements of m.
//
//    BandMatrix<T,stor>(colsize, rowsize, nlo, nhi, const T* m)
//    BandMatrix<T,stor>(colsize, rowsize, nlo, nhi, const valarray<T>& m)
//    BandMatrix<T,stor>(colsize, rowsize, nlo, nhi, const vector<T>& m)
//        Makes a BandMatrix which copies the copies the elements in m
//        according to the storage scheme (specified as stor) as 
//        described above.
//
// Special Constructors
//
//    UpperBiDiagMatrix(const Vector& v1, const Vector& v2)
//        Returns a (DiagMajor) BandMatrix with nlo=0, nhi=1, 
//        v1 on the main diagonal, and v2 on the superdiagonal
//
//    LowerBiDiagMatrix(const Vector& v1, const Vector& v2)
//        Returns a (DiagMajor) BandMatrix with nlo=1, nhi=0, 
//        v1 on the subdiagonal, and v2 on the main diagonal
//
//    TriDiagMatrix(const Vector& v1, const Vector& v2, const Vector& v3)
//        Returns a (DiagMajor) BandMatrix with nlo=1, nhi=1, 
//        v1 on the subdiagonal, v2 on the main diagonal, and 
//        v3 on the superdiagonal
//
//    ConstBandMatrixView BandMatrixViewOf(const Matrix<T>& m, nlo, nhi)
//    ConstBandMatrixView BandMatrixViewOf(const BandMatrix<T>& m, nlo, nhi)
//        Makes a constant BandMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//        For the second version, nlo,nhi must be <= the corresponding 
//        values in m.
//
//    BandMatrixView BandMatrixViewOf(Matrix<T>& m, nlo, nhi)
//    BandMatrixView BandMatrixViewOf(BandMatrix<T>& m, nlo, nhi)
//        Makes a modifiable BandMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the orginial Matrix.
//
//    ConstBandMatrixView BandMatrixViewOf(const T* m,
//            colsize, rowsize, nlo, nhi, stor)
//    BandMatrixView BandMatrixViewOf(T* m, colsize, rowsize, nlo, nhi, stor)
//        Makes a BandMatrixView of the elements in m using the actual
//        element m for the storage.  This is essentially the same as the 
//        constructor with (const T* m), except that the data isn't duplicated.
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//        Return the dimensions of the BandMatrix
//
//    int nlo() const
//    int nhi() const
//        Return the band dimensions
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the BandMatrix
//
//    VectorView row(size_t i, size_t j1, size_t j2)
//        Return a subset of the ith row of the BandMatrix
//        The range (i,j1)..(i,j2-1) must be entirely within the band.
//
//    VectorView col(size_t j, size_t i1, size_t i2)
//        Return a subset of the jth column of the BandMatrix
//        The range (i1,j)..(i2-1,j) must be entirely within the band.
//
//    VectorView diag()
//        Return the main diagonal of the BandMatrix
//
//    VectorView diag(int i)
//        Return the super- or sub-diagonal i
//        If i >= 0 return the super diagonal starting at (0,i)
//        If i <= 0 return the sub diagonal starting at (|i|,0)
//
//    VectorView diag(int i, int j1, int j2)
//        Return a subset of the (super/sub)-diagonal i. 
//        If i >= 0 the range is from (j1,i+j1)..(j2,i+j2)
//        If i <= 0 the range is from (|i|+j1,j1)..(|i|+j2,j2)
//
//
// Modifying Functions
//
//    BandMatrix& Zero()
//    BandMatrix<T>& TransposeSelf() 
//        Must be square, and have nhi=nlo for this function
//    BandMatrix& ConjugateSelf()
//    BandMatrix& SetToIdentity(x = 1)
//    void Swap(BandMatrix& m1, BandMatrix& m2)
//        Must be the same size and have the same band structure (nlo,nhi)
//
//
// Views of a BandMatrix:
//
//    SubBandMatrix(int i1, int i2, int j1, int j2, int nlo, int nhi,
//            int istep=1, int jstep=1)
//        Returns a BandMatrixView with (i1,j1) in the upper left corner,
//        (i2,j2) in the lower right corder, and with nlo and nhi 
//        sub- and super-diagonals respectively.
//        All members of the new SubBandMatrix must be within the 
//        original band.
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//    SubVector(int i1, int i2, int istep, int jstep, int size)
//        Just like the regular Matrix version, but all entries in the 
//        submatrix must be within the band.
//
//    Rows(int i1, int i2)
//        A bit different from the Matrix version, since the rowsize will
//        not be the full rowsize of the source BandMatrix.
//        Instead, this returns (as a BandMatrix) the nontrivial portion
//        of the BandMatrix in rows i1..i2 (not including i2)
//    Cols(int j1, int j2)
//        As with Rows, this returns the nontrivial portion of the 
//        BandMatrix in columns j1..j2 (not including j2)
//    Diags(int i1, int i2)
//        Returns a thinner BandMatrix including the diagonals from i1..i2
//        (not including i2)
//
//    m.View()
//    m.Transpose() or Transpose(m)
//    m.Conjugate() or Conjugate(m)
//    m.Adjoint() or Adjoint(m)
//
//
// Functions of BandMatrices:
//        (These are all both member functions and functions of a BandMatrix,
//         so Norm(m) and m.Norm() for example are equivalent.)
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
//    m.Inverse(minv)
//    m.InverseATA(invata)
//
//    m.NewTranspose()
//    m.NewConjugate()
//    m.NewAdjoint()
//    m.NewInverse()
//    m.NewView()
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
//          colsize rowsize nlo nhi
//          m(0,0) m(0,1) m(0,2) ... m(0,nhi) 
//          m(1,0) m(1,1) ... m(1,nhi+1)
//          ...
//          m(nlo,0) m(nlo,1) ... m(nlo,nhi+nlo)
//          ...
//          m(size-nhi,size-nlo-nhi) ... m(size-nhi,size)
//          ...
//          m(size-1,size-nlo-1) ... m(size-1,size) 
//          m(size,size-nlo) ... m(size,size)
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the BandMatrix to be read, you can
//        use this form where mptr is an auto_ptr to an undefined BandMatrix.
//
//
// Division Control Functions:
//
//    m.DivideUsing(dt) 
//    where dt is LU, QR, SV, SVS, SVU, SVV
// 
//    LUD(), QRD(), SVD() return the corresponding Divider classes.
//
//


#ifndef TMV_BandMatrix_H
#define TMV_BandMatrix_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenBandMatrix;
  template <class T, IndexStyle I=CStyle> class ConstBandMatrixView;
  template <class T, IndexStyle I=CStyle> class BandMatrixView;
  template <class T, StorageType S=RowMajor, IndexStyle I=CStyle> class BandMatrix;
  template <class T> class BandMatrixComposite;

  template <class T> class BandLUDiv;
  template <class T> class BandSVDiv;
  template <class T> class BandQRDiv;
  template <class T, class Tm> class QuotXB;

  template <class T> inline StorageType BaseStorOf(
      const GenBandMatrix<T>& m)
  {
    return (m.stor()==RowMajor || m.stor()==ColMajor || m.stor()==DiagMajor) ? 
      m.stor() : DiagMajor;
  }

  template <class T> class BandMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<BandMatrix<T> > m;
      char exp,got;
      size_t cs,rs;
      int lo,hi;
      bool is, iseof, isbad;

      inline BandMatrixReadError(
	  size_t _i, size_t _j, const GenBandMatrix<T>& _m,
	  istream& _is) throw() :
	ReadError("BandMatrix"),
	i(_i), j(_j), m(new BandMatrix<T>(_m)), exp(0), got(0),
	cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(istream& _is) throw() :
	ReadError("BandMatrix"),
	i(0), j(0), m(0), exp(0), got(0),
	cs(0), rs(0), lo(0), hi(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(
	  size_t _i, size_t _j, const GenBandMatrix<T>& _m,
	  istream& _is, char _e, char _g) throw() :
	ReadError("BandMatrix"),
	i(_i), j(_j), m(new BandMatrix<T>(_m)), exp(_e), got(_g),
	cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(istream& _is, char _e, char _g) throw() :
	ReadError("BandMatrix"),
	i(0), j(0), m(0), exp(_e), got(_g),
	cs(0), rs(0), lo(0), hi(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(const GenBandMatrix<T>& _m,
	  istream& _is, size_t _cs, size_t _rs, int _lo, int _hi) throw() :
	ReadError("BandMatrix"),
	i(0), j(0), m(new BandMatrix<T>(_m)), exp(0), got(0),
	cs(_cs), rs(_rs), lo(_lo), hi(_hi),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline BandMatrixReadError(const BandMatrixReadError<T>& rhs) :
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
	cs(rhs.cs), rs(rhs.rs), lo(rhs.lo), hi(rhs.hi),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~BandMatrixReadError() throw() {}

      virtual void Write(ostream& os) const throw();
  };
  

  size_t BandStorageLength(StorageType s, size_t cs, size_t rs,
      int lo, int hi);

  template <class T1, class T2> void Copy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2);

  template <class T> class GenBandMatrix : 
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline GenBandMatrix(StorageType s, ConjItType c, size_t ls) :
	itsstor(s), itsct(c), linsize(ls) {}
      inline GenBandMatrix(const GenBandMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs), itsstor(rhs.itsstor), itsct(rhs.itsct),
        linsize(rhs.linsize) {}
      inline ~GenBandMatrix() {}

      //
      // Access Functions
      //

      inline T operator()(size_t i,size_t j) const 
      { 
	TMVAssert(i< colsize());
	TMVAssert(j< rowsize());
	if (okij(i,j)) return cref(i,j); 
	else return T(0);
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2 && j2<=rowsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),
	    j2-j1,stepj(),itsct);
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2 && i2<=colsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),
	    i2-i1,stepi(),itsct);
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),
	    min(colsize(),rowsize()),diagstep(),itsct);
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i*stepj(),
	      diagsize,diagstep(),itsct);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return ConstVectorView<T>(cptr()-i*stepi(),
	      diagsize,diagstep(),itsct);
	}
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	TMVAssert(j1 <= j2);
	if (i >= 0) {
	  TMVAssert(j2<=min(rowsize()-i,colsize()));
	  return ConstVectorView<T>(cptr()+i*stepj()+j1*diagstep(),
	      j2-j1,diagstep(),itsct);
	} else {
	  TMVAssert(j2<=min(colsize()+i,rowsize()));
	  return ConstVectorView<T>(cptr()-i*stepi()+j1*diagstep(),
	      j2-j1,diagstep(),itsct);
	}
      }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      inline bool SameStorageAs(const GenBandMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      template <class T2> inline bool SameAs(const GenBandMatrix<T2>& ) const
      { return false; }

      inline bool SameAs(const GenBandMatrix<T>& m2) const
      { 
	return (this==&m2 || (cptr()==m2.cptr() && 
	    colsize()==m2.colsize() && rowsize()==m2.rowsize() && 
	    stepi()==m2.stepi() && stepj()==m2.stepj() &&
	    nhi()==m2.nhi() && nlo()==m2.nlo()));
      }

      inline void CopyToMatrix(const MatrixView<RealType(T)>& m2) const
      {
	TMVAssert(IsReal(T())); 
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	if (Real().SameStorageAs(m2)) {
	  UpperTriMatrixViewOf(m2.Cols(nhi()+1,rowsize())).Zero();
	  LowerTriMatrixViewOf(m2.Rows(nlo()+1,colsize())).Zero();
	} else {
	  m2.Zero();
	}
	BandMatrixViewOf(m2,nlo(),nhi()) = *this;
      }

      inline void CopyToMatrix(const MatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	if (Real().SameStorageAs(m2.Real())) {
	  UpperTriMatrixViewOf(m2.Cols(nhi()+1,rowsize())).Zero();
	  LowerTriMatrixViewOf(m2.Rows(nlo()+1,colsize())).Zero();
	} else {
	  m2.Zero();
	}
	BandMatrixViewOf(m2,nlo(),nhi()) = *this;
      }

      //
      // SubBandMatrix
      //

      bool OKSubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const;

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), itsstor, itsct);
      }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const 
      {
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
	    newstor,itsct);
      }

      bool OKSubVector(
	  size_t i, size_t j, int istep, int jstep, size_t size) const;

      inline ConstVectorView<T> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),itsct);
      }

      bool OKSubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const;

      inline ConstBandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(),itsstor,itsct);
      }

      inline ConstBandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
	    newstepi, newstepj, newstepi+newstepj, newstor, itsct);
      }

      inline ConstBandMatrixView<T> Rows(int i1, int i2) const
      {
	TMVAssert(i1 >= 0);
	TMVAssert(i1 < i2);
	TMVAssert(i2 <= int(colsize()));

	const int j1 = i1 > nlo() ? i1-nlo() : 0;
	const int j2 = min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1 < nlo() ? min(nlo(),i2-1) - i1 : 0;
	const int newnhi = min(nlo()+nhi()-newnlo,j2-j1-1);
	const size_t newlin = (ls() && isrm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin);
      }

      inline ConstBandMatrixView<T> Cols(int j1, int j2) const
      {
	TMVAssert(j1 >= 0);
	TMVAssert(j1 < j2);
	TMVAssert(j2 <= int(rowsize()));

	const int i1 = j1 > nhi() ? j1-nhi() : 0;
	const int i2 = min(j2 + nlo(),int(colsize()));
	const int newnhi = j1 < nhi() ? min(nhi(),j2-1) - j1 : 0;
	const int newnlo = min(nlo()+nhi()-newnhi,i2-i1-1);
	const size_t newlin = (ls() && iscm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin);
      }

      inline ConstBandMatrixView<T> Diags(int k1, int k2) const
      {
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 < 0 ? -k2 : 0;
	const int i2 = min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = min(rowsize(),colsize()+k2);
	const int newnlo = k2 < 0 ? k2-k1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k1 > 0 ? k2-k1-1 : k2 > 0 ? k2-1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct());
      }

      inline ConstBandMatrixView<T> UpperBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    min(colsize(),rowsize()),min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct());
      }

      inline ConstBandMatrixView<T> LowerBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    min(colsize(),rowsize()+nlo()),min(colsize(),rowsize()),
	    nlo(),0,stepi(),stepj(),diagstep(),stor(),ct());
      }
	
      inline ConstBandMatrixView<RealType(T)> Real() const
      {
	return ConstBandMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),
	    colsize(),rowsize(),nlo(),nhi(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? diagstep() : 2*diagstep(),
	    IsReal(T()) ? itsstor : NoMajor, NonConj);
      }

      inline ConstBandMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return ConstBandMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    colsize(),rowsize(),nlo(),nhi(),
	    2*stepi(),2*stepj(),2*diagstep(),NoMajor, NonConj);
      }

      //
      // Views
      //

      inline ConstBandMatrixView<T> View() const
      { 
	return ConstBandMatrixView<T>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),itsstor,itsct,
	    isdm()?0:linsize);
      }

      inline ConstBandMatrixView<T> Transpose() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(itsstor),itsct,
	    isdm()?0:linsize);
      }

      inline ConstBandMatrixView<T> Conjugate() const
      { 
	return ConstBandMatrixView<T>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),itsstor,ConjOf(T,itsct),
	    isdm()?0:linsize);
      }

      inline ConstBandMatrixView<T> Adjoint() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(itsstor),ConjOf(T,itsct),
	    isdm()?0:linsize);
      }

      inline BandMatrixView<T> NonConst() const
      {
	return BandMatrixView<T>(const_cast<T*>(cptr()),colsize(),rowsize(),
	    nlo(),nhi(),stepi(),stepj(),diagstep(),stor(),itsct,
	    isdm()?0:linsize
	    FIRSTLAST1(isdm() ? diag(-nlo()).begin().GetP() : cptr(),
	      isdm() ? diag(nhi()).end().GetP() :
	      ((diag().end()-1).GetP()+1)));
      }

      inline bool CanLinearize() const
      { 
	if (linsize == 1 && (rowsize() != 1 || colsize() != 1)) {
	  TMVAssert(isrm() || iscm());
	  if (isrm()) linsize = rowsize() + (colsize()-1)*stepi();
	  else linsize = colsize() + (rowsize()-1)*stepj();
	}
	//TMVAssert(linsize == 0 || linsize == 
	//    BandStorageLength(stor(),colsize(),rowsize(),nlo(),nhi()));
	return linsize > 0; 
      }

      // This needs to be virtual, since DiagMajor BandMatrices need to
      // start at itsm1, not cptr() (=itsm).  The BandMatrix version
      // of this gets that right.
      // Also, because of the virtual, we now need two names (hence the
      // Const), since the BandMatrixView::LinearView() function has a 
      // different return type, so it is not "covariant" with this
      // funtion, which it needs to be for virtual overriding.
      virtual inline ConstVectorView<T> ConstLinearView() const
      {
	TMVAssert(!isdm());
	TMVAssert(linsize != 1 || (rowsize() == 1 && colsize() == 1));
	// (To assure that next assert has no effect.)

	TMVAssert(CanLinearize());
	return ConstVectorView<T>(cptr(),linsize,1,itsct);
      }

      //
      // Functions of Matrix
      //

      inline T Trace() const
      { return diag().SumElements(); }

      // Default Norm = NormF()
      inline RealType(T) Norm() const
      { return NormF(); }

      // Frobenius norm = sqrt(sum_ij |a_ij|^2 )
      RealType(T) NormF() const;

      // NormF()^2
      RealType(T) NormSq() const;

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const;

      // 2-Norm defined in BaseMatrix.h
      
      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return Transpose().Norm1(); }

      // max_i,j (|a_ij|)
      RealType(T) MaxAbsElement() const;

      inline QuotXB<T,T> Inverse() const
      { return QuotXB<T,T>(T(1),*this); }
      
      using BaseMatrix<T>::Inverse;

      using BaseMatrix<T>::InverseATA;

      template <StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T,S,I>& minv) const
      { BaseMatrix<T>::Inverse(minv.View()); }

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& minv) const
      { BaseMatrix<T>::Inverse(minv.View()); }

      inline auto_ptr<BaseMatrix<T> > NewTranspose() const 
      {
	auto_ptr<BaseMatrix<T> > a(new ConstBandMatrixView<T>(Transpose())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewConjugate() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstBandMatrixView<T>(Conjugate())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewAdjoint() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstBandMatrixView<T>(Adjoint())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewView() const
      { 
	auto_ptr<BaseMatrix<T> > a(new ConstBandMatrixView<T>(View())); 
	return a;
      }

      inline auto_ptr<BaseMatrix<T> > NewCopy() const
      { 
	auto_ptr<BaseMatrix<T> > a;
	if (isrm()) a.reset(new BandMatrix<T,RowMajor>(*this)); 
	else if (iscm()) a.reset(new BandMatrix<T,ColMajor>(*this)); 
	else a.reset(new BandMatrix<T,DiagMajor>(*this)); 
	return a;
      }

      //
      // Division Control
      //

      using BaseMatrix<T>::GetDiv;

      inline void DivideUsing(DivType dt) const
      {
	TMVAssert(dt == LU || dt == QR || SV_Type(dt));
	BaseMatrix<T>::DivideUsing(dt);
      }

      inline const BandLUDiv<T>& LUD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const BandLUDiv<T>*>(GetDiv()));
	return *dynamic_cast<const BandLUDiv<T>*>(GetDiv());
      }

      inline const BandQRDiv<T>& QRD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const BandQRDiv<T>*>(GetDiv()));
	return *dynamic_cast<const BandQRDiv<T>*>(GetDiv());
      }

      inline const BandSVDiv<T>& SVD() const
      {
	TMVAssert(GetDiv());
	TMVAssert(dynamic_cast<const BandSVDiv<T>*>(GetDiv()));
	return *dynamic_cast<const BandSVDiv<T>*>(GetDiv());
      }


      //
      // I/O
      //

      void WriteCompact(ostream& fout) const;
      void Write(ostream& fout) const;
      void WriteCompact(ostream& fout, RealType(T) thresh) const;
      void Write(ostream& fout, RealType(T) thresh) const;

      using BaseMatrix<T>::colsize;
      using BaseMatrix<T>::rowsize;
      using BaseMatrix<T>::IsSquare;
      virtual int nlo() const =0;
      virtual int nhi() const =0;
      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual int diagstep() const = 0;
      virtual inline StorageType stor() const { return itsstor; }
      virtual inline ConjItType ct() const { return itsct; }
      virtual inline bool isrm() const { return itsstor == RowMajor; }
      virtual inline bool iscm() const { return itsstor == ColMajor; }
      virtual inline bool isdm() const { return itsstor == DiagMajor; }
      inline size_t ls() const { return linsize; }
      inline bool isconj() const
      {
	TMVAssert(IsComplex(T()) || itsct==NonConj);
	return IsComplex(T()) && itsct==Conj;
      }

    protected :

      inline bool okij(size_t i, size_t j) const
      { return (j+nlo() >= i && i+nhi() >= j); }

      virtual T cref(size_t i, size_t j) const;

      void NewDivider() const;

    private :

      const StorageType itsstor;
      const ConjItType itsct;
      mutable size_t linsize;

      inline void operator=(const GenBandMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenBandMatrix

  template <class T, IndexStyle I> class ConstBandMatrixView : 
    public GenBandMatrix<T>
  {
    public :

      inline ConstBandMatrixView(const ConstBandMatrixView<T,I>& rhs) :
	GenBandMatrix<T>(rhs), itsm(rhs.itsm),
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd) {}

      inline ConstBandMatrixView(const GenBandMatrix<T>& rhs) :
	GenBandMatrix<T>(rhs), itsm(rhs.cptr()),
	itscs(rhs.colsize()), itsrs(rhs.rowsize()), 
	itsnlo(rhs.nlo()), itsnhi(rhs.nhi()),
	itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()) {}

      inline ConstBandMatrixView(
	  const T* _m, size_t _cs, size_t _rs, 
	  int _lo, int _hi, int _si, int _sj, int _sd, 
	  StorageType instor, ConjItType inct, size_t ls=0) : 
	GenBandMatrix<T>(instor,inct,ls), itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd)
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ? _si==1 :
	    instor==DiagMajor ? _sd==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
	//TMVAssert(ls==0 || ls==1 || 
	//    ls==BandStorageLength(instor,_cs,_rs,_lo,_hi));
      }

      inline ~ConstBandMatrixView() {}

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline int diagstep() const { return itssd; }

    protected :

      const T*const itsm;
      const size_t itscs;
      const size_t itsrs;
      const int itsnlo;
      const int itsnhi;
      const int itssi;
      const int itssj;
      const int itssd;

    private :

      inline void operator=(const ConstBandMatrixView<T>&) 
      { TMVAssert(FALSE); }

  }; // ConstBandMatrixView

  template <class T> class ConstBandMatrixView<T,FortranStyle> : 
    public ConstBandMatrixView<T,CStyle>
  {
    public :

      inline ConstBandMatrixView(
	  const ConstBandMatrixView<T,FortranStyle>& rhs) :
	ConstBandMatrixView<T,CStyle>(rhs) {}

      inline ConstBandMatrixView(const GenBandMatrix<T>& rhs) :
	ConstBandMatrixView<T,CStyle>(rhs) {}

      inline ConstBandMatrixView(const T* _m, size_t _cs, size_t _rs, 
	  int _lo, int _hi, int _si, int _sj, int _sd, 
	  StorageType instor, ConjItType inct, size_t ls=0) : 
	ConstBandMatrixView<T,CStyle>(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,
	    instor,inct,ls) {}

      inline ~ConstBandMatrixView() {}

      //
      // Access Functions
      //

      inline T operator()(size_t i,size_t j) const 
      { 
	TMVAssert(i>0 && i<= colsize());
	TMVAssert(j>0 && j<= rowsize());
	if (okij(i-1,j-1)) return cref(i-1,j-1); 
	else return T(0);
      }

      inline ConstVectorView<T,FortranStyle> row(
	  size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i>0 && i<=colsize());
	TMVAssert(j1>0 && j1<=j2 && j2<=rowsize());
	TMVAssert(okij(i-1,j1-1));
	TMVAssert(okij(i-1,j2-1));
	return GenBandMatrix<T>::row(i-1,j1-1,j2);
      }

      inline ConstVectorView<T,FortranStyle> col(
	  size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j>0 && j<=rowsize());
	TMVAssert(i1>0 && i1<=i2 && i2<=colsize());
	TMVAssert(okij(i1-1,j-1));
	TMVAssert(okij(i2-1,j-1));
	return GenBandMatrix<T>::col(j-1,i1-1,i2);
      }

      inline ConstVectorView<T,FortranStyle> diag() const
      { return GenBandMatrix<T>::diag(); }

      inline ConstVectorView<T,FortranStyle> diag(int i) const
      { return GenBandMatrix<T>::diag(i); }

      inline ConstVectorView<T,FortranStyle> diag(
	  int i, size_t j1, size_t j2) const
      {
	TMVAssert(j1>0);
	return GenBandMatrix<T>::diag(i,j1-1,j2);
      }


      //
      // SubBandMatrix
      //

      bool OKSubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const;

      inline ConstMatrixView<T,FortranStyle> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return GenBandMatrix<T>::SubMatrix(i1-1,i2,j1-1,j2);
      }

      inline ConstMatrixView<T,FortranStyle> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const 
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return GenBandMatrix<T>::SubMatrix(i1-1,i2-1+istep,j1-1,j2-1+jstep,
	    istep,jstep);
      }

      bool OKSubVector(
	  size_t i, size_t j, int istep, int jstep, size_t size) const;

      inline ConstVectorView<T,FortranStyle> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return GenBandMatrix<T>::SubVector(i-1,j-1,istep,jstep,size);
      }

      bool OKSubBandMatrix(int i1, int i2, int j1, int j2,
	  int newnlo, int newnhi, int istep, int jstep) const;

      inline ConstBandMatrixView<T,FortranStyle> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	return GenBandMatrix<T>::SubBandMatrix(i1-1,i2,j1-1,j2,newnlo,newnhi);
      }

      inline ConstBandMatrixView<T,FortranStyle> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
	return GenBandMatrix<T>::SubBandMatrix(i1-1,i2-1+istep,j1-1,j2-1+jstep,
	    newnlo,newnhi,istep,jstep);
      }

      inline ConstBandMatrixView<T,FortranStyle> Rows(int i1, int i2) const
      {
	TMVAssert(i1 > 0 && i1 <= i2);
	TMVAssert(i2 <= int(colsize()));
	return GenBandMatrix<T>::Rows(i1-1,i2);
      }

      inline ConstBandMatrixView<T,FortranStyle> Cols(int j1, int j2) const
      {
	TMVAssert(j1 > 0 && j1 <= j2);
	TMVAssert(j2 <= int(rowsize()));
	return GenBandMatrix<T>::Cols(j1-1,j2);
      }

      inline ConstBandMatrixView<T,FortranStyle> Diags(int k1, int k2) const
      {
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);
	return GenBandMatrix<T>::Diags(k1,k2);
      }
	
      inline ConstBandMatrixView<T,FortranStyle> UpperBand() const
      { return GenBandMatrix<T>::UpperBand(); }

      inline ConstBandMatrixView<T,FortranStyle> LowerBand() const
      { return GenBandMatrix<T>::LowerBand(); }
	
      inline ConstBandMatrixView<RealType(T),FortranStyle> Real() const
      { return GenBandMatrix<T>::Real(); }

      inline ConstBandMatrixView<RealType(T),FortranStyle> Imag() const
      { return GenBandMatrix<T>::Imag(); }

      //
      // Views
      //

      inline ConstBandMatrixView<T,FortranStyle> View() const
      { return GenBandMatrix<T>::View(); }

      inline ConstBandMatrixView<T,FortranStyle> Transpose() const
      { return GenBandMatrix<T>::Transpose(); }

      inline ConstBandMatrixView<T,FortranStyle> Conjugate() const
      { return GenBandMatrix<T>::Conjugate(); }

      inline ConstBandMatrixView<T,FortranStyle> Adjoint() const
      { return GenBandMatrix<T>::Adjoint(); }

      inline BandMatrixView<T,FortranStyle> NonConst() const
      { return GenBandMatrix<T>::NonConst(); }

      using ConstBandMatrixView<T,CStyle>::colsize;
      using ConstBandMatrixView<T,CStyle>::rowsize;
      using ConstBandMatrixView<T,CStyle>::nlo;
      using ConstBandMatrixView<T,CStyle>::nhi;

    protected :

      using GenBandMatrix<T>::okij;
      using GenBandMatrix<T>::cref;

    private :

      inline void operator=(const ConstBandMatrixView<T,FortranStyle>&) 
      { TMVAssert(FALSE); }

  }; // FortranStyle ConstBandMatrixView

  template <class T, IndexStyle I> class BandMatrixView : 
    public GenBandMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline BandMatrixView(const BandMatrixView<T,I>& rhs) : 
	GenBandMatrix<T>(rhs), itsm(rhs.itsm), 
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd)
	  DEFFIRSTLAST(rhs.first,rhs.last) {}

      inline BandMatrixView(
	  T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType instor, ConjItType inct,
	  size_t ls PARAMFIRSTLAST(T) ) :
	GenBandMatrix<T>(instor,inct,ls), itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd)
	  DEFFIRSTLAST(_first,_last)
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : 
	    instor==ColMajor ? _si==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
	//TMVAssert(ls==0 || ls==1 || 
	//    ls==BandStorageLength(instor,_cs,_rs,_lo,_hi));
      }

      inline BandMatrixView(
	  T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType instor, ConjItType inct
	  PARAMFIRSTLAST(T) ) :
	GenBandMatrix<T>(instor,inct,0), itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd)
	  DEFFIRSTLAST(_first,_last)
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : 
	    instor==ColMajor ? _si==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      inline ~BandMatrixView() {}

      //
      // Op=
      //

      inline const BandMatrixView<T,I>& operator=(
	  const BandMatrixView<T,I>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	if (!SameAs(m2)) {
	  Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	  for(int i=-nlo();i<-m2.nlo();++i) diag(i).Zero();
	  for(int i=m2.nhi()+1;i<=nhi();++i) diag(i).Zero();
	}
	return *this; 
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenBandMatrix<T>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	if (!SameAs(m2)) {
	  Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	  for(int i=-nlo();i<-m2.nlo();++i) diag(i).Zero();
	  for(int i=m2.nhi()+1;i<=nhi();++i) diag(i).Zero();
	}
	return *this; 
      }

      template <class T2> inline const BandMatrixView<T,I>& operator=(
	    const GenBandMatrix<T2>& m2) const
      {
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	for(int i=-nlo();i<-m2.nlo();++i) diag(i).Zero();
	for(int i=m2.nhi()+1;i<=nhi();++i) diag(i).Zero();
	return *this; 
      }

      inline const BandMatrixView<T,I>& operator=(T x) const 
      {
	TMVAssert(IsSquare());
	return SetToIdentity(x); 
      }

      inline const BandMatrixView<T,I>& operator=(
	  const BandMatrixComposite<T>& mcomp) const
      { 
	TMVAssert(colsize() == mcomp.colsize());
	TMVAssert(rowsize() == mcomp.rowsize());
	TMVAssert(nlo() >= mcomp.nlo());
	TMVAssert(nhi() >= mcomp.nhi());
	mcomp.AssignTo(SubBandMatrix(0,colsize(),0,rowsize(),
	      mcomp.nlo(),mcomp.nhi()));
	for(int i=-nlo();i<-mcomp.nlo();++i) diag(i).Zero();
	for(int i=mcomp.nhi()+1;i<=nhi();++i) diag(i).Zero();
	return *this;
      }

      //
      // Access
      //
 
      inline RefType(T) operator()(size_t i,size_t j) const 
      { 
	TMVAssert(i<colsize());
	TMVAssert(j<rowsize());
	TMVAssert(okij(i,j));
	return ref(i,j); 
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1 <= j2 && j2 <= rowsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	const size_t start = i*stepi() + j1*stepj();
	const size_t rowlen = j2 - j1;
	return VectorView<T>(ptr()+start,rowlen,stepj(),ct() FIRSTLAST );
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1 <= i2 && i2 <= colsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	const size_t start = i1*stepi() + j*stepj();
	const size_t collen = i2 - i1;
	return VectorView<T>(ptr()+start,collen,stepi(),ct() FIRSTLAST );
      }

      inline VectorView<T> diag() const
      {
	return VectorView<T>(ptr(),min(colsize(),rowsize()),diagstep(),ct() 
	    FIRSTLAST);
      }

      inline VectorView<T> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diaglen = min(rowsize()-i,colsize());
	  return VectorView<T>(ptr()+i*stepj(),diaglen,diagstep(),ct() 
	      FIRSTLAST );
	} else {
	  const size_t diaglen = min(colsize()+i,rowsize());
	  return VectorView<T>(ptr()-i*stepi(),diaglen,diagstep(),ct() 
	      FIRSTLAST );
	}
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	TMVAssert(j1 <= j2);
	if (i >= 0) {
	  TMVAssert(j2<=min(rowsize()-i,colsize()));
	  return VectorView<T>(ptr()+i*stepj()+j1*diagstep(),
	      j2-j1, diagstep(),ct() FIRSTLAST );
	} else {
	  TMVAssert(j2<=min(colsize()+i,rowsize()));
	  return VectorView<T>(ptr()-i*stepi()+j1*diagstep(),
	      j2-j1, diagstep(),ct() FIRSTLAST );
	}
      }

      //
      // Modifying Functions
      //

      inline const BandMatrixView<T,I>& Zero() const 
      { return SetAllTo(T(0)); }

      const BandMatrixView<T,I>& Clip(RealType(T) thresh) const;

      const BandMatrixView<T,I>& SetAllTo(T x) const;

      inline const BandMatrixView<T,I>& TransposeSelf() const
      { 
	TMVAssert(IsSquare());
	TMVAssert(nlo() == nhi());
	for(int i=1;i<nhi();++i) Swap(diag(-i),diag(i));
	return *this;
      }

      const BandMatrixView<T,I>& ConjugateSelf() const;

      inline const BandMatrixView<T,I>& SetToIdentity(T x=T(1)) const 
      {
	TMVAssert(IsSquare());
	Zero(); diag().SetAllTo(x); return *this; 
      }

      //
      // SubBandMatrix
      //

      inline MatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(GenBandMatrix<T>::OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), stor(), ct() FIRSTLAST);
      }

      inline MatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const 
      {
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(GenBandMatrix<T>::OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, ct() FIRSTLAST);
      }

      inline VectorView<T> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t size) const
      {
	TMVAssert(GenBandMatrix<T>::OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(), ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
	      1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
	      istep,jstep));
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
	    stepi()*istep, stepj()*jstep, newstepi+newstepj, newstor, 
	    ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> Rows(int i1, int i2) const
      {
	TMVAssert(i1 >= 0);
	TMVAssert(i1 < i2);
	TMVAssert(i2 <= int(colsize()));

	const int j1 = i1 > nlo() ? i1-nlo() : 0;
	const int j2 = min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1 < nlo() ? min(nlo(),i2-1) - i1 : 0;
	const int newnhi = min(nlo()+nhi()-newnlo,j2-j1-1);
	const size_t newlin = (ls() && isrm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
	      1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin FIRSTLAST);
      }

      inline BandMatrixView<T,I> Cols(int j1, int j2) const
      {
	TMVAssert(j1 >= 0);
	TMVAssert(j1 < j2);
	TMVAssert(j2 <= int(rowsize()));

	const int i1 = j1 > nhi() ? j1-nhi() : 0;
	const int i2 = min(j2 + nlo(),int(colsize()));
	const int newnhi = j1 < nhi() ? min(nhi(),j2-1) - j1 : 0;
	const int newnlo = min(nlo()+nhi()-newnhi,i2-i1-1);
	const size_t newlin = (ls() && iscm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin FIRSTLAST);
      }

      inline BandMatrixView<T,I> Diags(int k1, int k2) const
      {
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 < 0 ? -k2 : 0;
	const int i2 = min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = min(rowsize(),colsize()+k2);
	const int newnlo = k2 < 0 ? k2-k1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k1 > 0 ? k2-k1-1 : k2 > 0 ? k2-1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct() FIRSTLAST);
      }
	
      inline BandMatrixView<T,I> UpperBand() const
      {
	return BandMatrixView<T,I>(ptr(),
	    min(colsize(),rowsize()),min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> LowerBand() const
      {
	return BandMatrixView<T,I>(ptr(),
	    min(colsize(),rowsize()+nlo()),min(colsize(),rowsize()),
	    nlo(),0,stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      }
	
      inline BandMatrixView<RealType(T),I> Real() const
      {
	return BandMatrixView<RealType(T),I>(
	    reinterpret_cast<RealType(T)*>(ptr()),
	    colsize(),rowsize(),nlo(),nhi(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? diagstep() : 2*diagstep(),
	    IsReal(T()) ? stor() : NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline BandMatrixView<RealType(T),I> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return BandMatrixView<RealType(T),I>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    colsize(),rowsize(),nlo(),nhi(),
	    2*stepi(),2*stepj(),2*diagstep(),NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      inline BandMatrixView<T,I> View() const
      { return *this; }

      inline BandMatrixView<T,I> Transpose() const
      { 
	return BandMatrixView<T,I>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(stor()),ct(),
	    isdm()?0:ls() FIRSTLAST);
      }

      inline BandMatrixView<T,I> Conjugate() const
      { 
	return BandMatrixView<T,I>(ptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),stor(),ConjOf(T,ct()),
	    isdm()?0:ls() FIRSTLAST);
      }

      inline BandMatrixView<T,I> Adjoint() const
      { 
	return BandMatrixView<T,I>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(stor()),ConjOf(T,ct()),
	    isdm()?0:ls()
	    FIRSTLAST);
      }

      inline VectorView<T> LinearView() const
      {
	TMVAssert(!isdm());
	TMVAssert(ls() != 1 || (rowsize() == 1 && colsize() == 1));
	  // To assure that next assert has no effect.

	TMVAssert(GenBandMatrix<T>::CanLinearize());
	return VectorView<T>(ptr(),ls(),1,ct() FIRSTLAST );
      }


      //
      // I/O
      //

      void Read(istream& fin) const;

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      using BaseMatrix<T>::IsSquare;
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline int diagstep() const { return itssd; }
      using GenBandMatrix<T>::stor;
      using GenBandMatrix<T>::isrm;
      using GenBandMatrix<T>::iscm;
      using GenBandMatrix<T>::isdm;
      using GenBandMatrix<T>::ct;
      using GenBandMatrix<T>::ls;
      using GenBandMatrix<T>::isconj;

      // This makes it easier for some things to compile.
      // The ones that work are overridden.
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

    protected:

      T*const itsm;
      const size_t itscs;
      const size_t itsrs;
      const int itsnlo;
      const int itsnhi;
      const int itssi;
      const int itssj;
      const int itssd;

#ifdef TMVFLDEBUG
    public :
      const T*const first;
      const T*const last;
    protected :
#endif

      using GenBandMatrix<T>::okij;
      RefType(T) ref(size_t i, size_t j) const;

  }; // BandMatrixView

  template <class T> class BandMatrixView<T,FortranStyle> : 
    public BandMatrixView<T,CStyle>
  {

    public:

      //
      // Constructors
      //

      inline BandMatrixView(const BandMatrixView<T,FortranStyle>& rhs) : 
	BandMatrixView<T,CStyle>(rhs) {}

      inline BandMatrixView(const BandMatrixView<T,CStyle>& rhs) : 
	BandMatrixView<T,CStyle>(rhs) {}

      inline BandMatrixView(
	  T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType instor, ConjItType inct,
	  size_t ls PARAMFIRSTLAST(T) ) :
	BandMatrixView<T,CStyle>(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,instor,inct,ls
	    FIRSTLAST1(_first,_last) ) {}

      inline BandMatrixView(
	  T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType instor, ConjItType inct
	  PARAMFIRSTLAST(T) ) :
	BandMatrixView<T,CStyle>(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,instor,inct
	    FIRSTLAST1(_first,_last) ) {}

      inline ~BandMatrixView() {}

      //
      // Op=
      //

      inline const BandMatrixView<T,FortranStyle>& operator=(
	  const BandMatrixView<T,FortranStyle>& m2) const
      {
	TMVAssert(colsize()==m2.colsize());
	TMVAssert(rowsize()==m2.rowsize());
	TMVAssert(nhi()==m2.nhi());
	TMVAssert(nlo()==m2.nlo());
	BandMatrixView<T,CStyle>::operator=(m2); 
	return *this; 
      }

      inline const BandMatrixView<T,FortranStyle>& operator=(
	  const GenBandMatrix<T>& m2) const
      { 
	TMVAssert(colsize()==m2.colsize());
	TMVAssert(rowsize()==m2.rowsize());
	TMVAssert(nhi()==m2.nhi());
	TMVAssert(nlo()==m2.nlo());
	BandMatrixView<T,CStyle>::operator=(m2); 
	return *this; 
      }

      template <class T2> 
	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenBandMatrix<T2>& m2) const
	{ 
	  TMVAssert(colsize()==m2.colsize());
	  TMVAssert(rowsize()==m2.rowsize());
	  TMVAssert(nhi()==m2.nhi());
	  TMVAssert(nlo()==m2.nlo());
	  BandMatrixView<T,CStyle>::operator=(m2); 
	  return *this; 
	}

      inline const BandMatrixView<T,FortranStyle>& operator=(T x) const 
      { 
	BandMatrixView<T,CStyle>::operator=(x); 
	return *this; 
      }

      inline const BandMatrixView<T,FortranStyle>& operator=(
	  const BandMatrixComposite<T>& mcomp) const
      { 
	TMVAssert(colsize()==mcomp.colsize());
	TMVAssert(rowsize()==mcomp.rowsize());
	TMVAssert(nhi()==mcomp.nhi());
	TMVAssert(nlo()==mcomp.nlo());
	BandMatrixView<T,CStyle>::operator=(mcomp); 
	return *this; 
      }

      //
      // Access
      //
 
      inline RefType(T) operator()(size_t i,size_t j) const 
      { 
	TMVAssert(i>0 && i<=colsize());
	TMVAssert(j>0 && j<=rowsize());
	TMVAssert(okij(i-1,j-1));
	return ref(i-1,j-1); 
      }

      inline VectorView<T,FortranStyle> row(
	  size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i>0 && i<=colsize());
	TMVAssert(j1 > 0 && j1 <= j2 && j2 <= rowsize());
	TMVAssert(okij(i-1,j1-1));
	TMVAssert(okij(i-1,j2-1));
	return BandMatrixView<T,CStyle>::row(i-1,j1-1,j2);
      }

      inline VectorView<T,FortranStyle> col(
	  size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j>0 && j<=rowsize());
	TMVAssert(i1 > 0 && i1 <= i2 && i2 <= colsize());
	TMVAssert(okij(i1-1,j-1));
	TMVAssert(okij(i2-1,j-1));
	return BandMatrixView<T,CStyle>::col(j-1,i1-1,i2);
      }

      inline VectorView<T,FortranStyle> diag() const
      { return BandMatrixView<T,CStyle>::diag(); }

      inline VectorView<T,FortranStyle> diag(int i) const
      { return BandMatrixView<T,CStyle>::diag(i); }

      inline VectorView<T,FortranStyle> diag(
	  int i, size_t j1, size_t j2) const
      {
	TMVAssert(j1>0);
	return BandMatrixView<T,CStyle>::diag(i,j1-1,j2); 
      }


      //
      // Modifying Functions
      //

      inline const BandMatrixView<T,FortranStyle>& Zero() const 
      { BandMatrixView<T,CStyle>::Zero(); return *this; }

      inline const BandMatrixView<T,FortranStyle>& Clip(
	  RealType(T) thresh) const
      { BandMatrixView<T,CStyle>::Clip(thresh); return *this; }

      inline const BandMatrixView<T,FortranStyle>& SetAllTo(T x) const
      { BandMatrixView<T,CStyle>::SetAllTo(x); return *this; }

      inline const BandMatrixView<T,FortranStyle>& TransposeSelf() const
      { BandMatrixView<T,CStyle>::TransposeSelf(); return *this; }

      inline const BandMatrixView<T,FortranStyle>& ConjugateSelf() const
      { BandMatrixView<T,CStyle>::ConjugateSelf(); return *this; }

      inline const BandMatrixView<T,FortranStyle>& SetToIdentity(
	  T x=T(1)) const 
      { BandMatrixView<T,CStyle>::SetToIdentity(x); return *this; }

      //
      // SubBandMatrix
      //

      inline bool OKSubVector(
	  int i, int j, int istep, int jstep, int size) const
      {
	return ConstBandMatrixView<T,FortranStyle>(*this).OKSubVector(
	    i,j,istep,jstep,size);
      }

      inline bool OKSubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	return ConstBandMatrixView<T,FortranStyle>(*this).OKSubMatrix(
	    i1,i2,j1,j2,istep,jstep);
      }

      inline bool OKSubBandMatrix(int i1, int i2, int j1, int j2,
	  int newnlo, int newnhi, int istep, int jstep) const
      {
	return ConstBandMatrixView<T,FortranStyle>(*this).OKSubBandMatrix(
	    i1,i2,j1,j2,newnlo,newnhi,istep,jstep);
      }

      inline MatrixView<T,FortranStyle> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return BandMatrixView<T,CStyle>::SubMatrix(i1-1,i2,j1-1,j2);
      }

      inline MatrixView<T,FortranStyle> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const 
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return BandMatrixView<T,CStyle>::SubMatrix(i1-1,i2-1+istep,
	    j1-1,j2-1+jstep,istep,jstep);
      }

      inline VectorView<T,FortranStyle> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return BandMatrixView<T,CStyle>::SubVector(i-1,j-1,istep,jstep,size);
      }

      inline BandMatrixView<T,FortranStyle> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,CStyle>::SubBandMatrix(i1-1,i2,j1-1,j2,
	    newnlo,newnhi);
      }

      inline BandMatrixView<T,FortranStyle> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
	return BandMatrixView<T,CStyle>::SubBandMatrix(i1-1,i2-1+istep,
	    j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
      }

      inline BandMatrixView<T,FortranStyle> Rows(int i1, int i2) const
      {
	TMVAssert(i1 > 0 && i1 <= i2);
	return BandMatrixView<T,CStyle>::Rows(i1-1,i2);
      }

      inline BandMatrixView<T,FortranStyle> Cols(int j1, int j2) const
      {
	TMVAssert(j1 > 0 && j1 <= j2);
	return BandMatrixView<T,CStyle>::Cols(j1-1,j2);
      }

      inline BandMatrixView<T,FortranStyle> Diags(int k1, int k2) const
      { return BandMatrixView<T,CStyle>::Diags(k1,k2); }
	
      inline BandMatrixView<T,FortranStyle> UpperBand() const
      { return BandMatrixView<T,CStyle>::UpperBand(); }

      inline BandMatrixView<T,FortranStyle> LowerBand() const
      { return BandMatrixView<T,CStyle>::LowerBand(); }
	
      inline BandMatrixView<RealType(T),FortranStyle> Real() const
      { return BandMatrixView<T,CStyle>::Real(); }

      inline BandMatrixView<RealType(T),FortranStyle> Imag() const
      { return BandMatrixView<T,CStyle>::Imag(); }

      inline BandMatrixView<T,FortranStyle> View() const
      { return BandMatrixView<T,CStyle>::View(); }

      inline BandMatrixView<T,FortranStyle> Transpose() const
      { return BandMatrixView<T,CStyle>::Transpose(); }

      inline BandMatrixView<T,FortranStyle> Conjugate() const
      { return BandMatrixView<T,CStyle>::Conjugate(); }

      inline BandMatrixView<T,FortranStyle> Adjoint() const
      { return BandMatrixView<T,CStyle>::Adjoint(); }

      inline VectorView<T,FortranStyle> LinearView() const
      { return BandMatrixView<T,CStyle>::LinearView(); }

      using BandMatrixView<T,CStyle>::colsize;
      using BandMatrixView<T,CStyle>::rowsize;
      using BandMatrixView<T,CStyle>::nlo;
      using BandMatrixView<T,CStyle>::nhi;

      template <class T2> inline void operator=(const BaseMatrix<T2>&) const 
      { TMVAssert(FALSE); }

    protected :

      using GenBandMatrix<T>::okij;
      using BandMatrixView<T,CStyle>::ref;

  }; // FortranStyle BandMatrixView

  template <class T, StorageType S, IndexStyle I> class BandMatrix : 
    public GenBandMatrix<T>
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(cs,rs,lo,hi) \
      GenBandMatrix<T>(S,NonConj,BandStorageLength(S,cs,rs,lo,hi)), \
      itsm1(new T[ls()]), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
      itssi(S==RowMajor ? lo+hi : S==ColMajor ? 1 : \
	  rs>= cs ? -int(cs)+1 : -int(rs) ), \
      itssj(S==RowMajor ? 1 : S==ColMajor ? lo+hi : -itssi+1), \
      itsds(S==RowMajor ? itssi+1 : S==ColMajor ? itssj+1 : 1), \
      itsm(S==DiagMajor ? itsm1.get() - lo*itssi : itsm1.get()) \
      DEFFIRSTLAST(itsm1.get(),itsm1.get()+ls())

      inline BandMatrix(size_t cs, size_t rs, int lo, int hi) :
	  NEW_SIZE(cs,rs,lo,hi) 
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(cs));
	TMVAssert(hi < int(rs));
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      inline BandMatrix(size_t cs, size_t rs, int lo, int hi, T x) :
	  NEW_SIZE(cs,rs,lo,hi) 
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(cs));
	TMVAssert(hi < int(rs));
	SetAllTo(x);
      }

      inline BandMatrix(
	  size_t cs, size_t rs, int lo, int hi, const valarray<T>& vv) :
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(cs));
	TMVAssert(hi < int(rs));
	TMVAssert(vv.size() == ls());
	T* vi=itsm1.get();
	for(size_t i=0;i<vv.size();++i,++vi) *vi=vv[i];
      }

      inline BandMatrix(size_t cs, size_t rs, int lo, int hi, const T* vv) :
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(cs));
	TMVAssert(hi < int(rs));
	memmove(itsm1.get(),vv,ls()*sizeof(T));
      }

      inline BandMatrix(
	  size_t cs, size_t rs, int lo, int hi, const vector<T>& vv) : 
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(cs));
	TMVAssert(hi < int(rs));
	TMVAssert(vv.size() == ls());
	T* vi = itsm1.get();
        typename vector<T>::const_iterator vvi = vv.begin();
	for(size_t i=vv.size();i>0;--i,++vi,++vvi) *vi=*vvi;
      }

      inline BandMatrix(const BandMatrix<T,S,I>& rhs) :
	GenBandMatrix<T>(S,NonConj,rhs.ls()), itsm1(new T[ls()]),
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itsds(rhs.itsds),
	itsm(S==DiagMajor ? itsm1.get()-itsnlo*itssi : itsm1.get())
	  DEFFIRSTLAST(itsm1.get(),itsm1.get()+ls())
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	memmove(itsm1.get(),rhs.itsm1.get(),ls()*sizeof(T));
      }

      template <IndexStyle I2> inline BandMatrix(
	  const BandMatrix<T,S,I2>& rhs) :
	GenBandMatrix<T>(S,NonConj,rhs.ls()), itsm1(new T[ls()]),
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itsds(rhs.itsds),
	itsm(S==DiagMajor ? itsm1.get()-itsnlo*itssi : itsm1.get())
	  DEFFIRSTLAST(itsm1.get(),itsm1.get()+ls())
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	memmove(itsm1.get(),rhs.start_mem(),ls()*sizeof(T));
      }

      inline BandMatrix(const GenBandMatrix<T>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),rhs.nlo(),rhs.nhi())
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	if (S == rhs.stor() && rhs.CanLinearize()) {
	  TMVAssert(stepi() == rhs.stepi());
	  TMVAssert(stepj() == rhs.stepj());
	  TMVAssert(ls() == rhs.ConstLinearView().size());
	  LinearView() = rhs.ConstLinearView();
	} else {
	  Copy(rhs,View()); 
	}
      }

      template <class T2> inline BandMatrix(const GenBandMatrix<T2>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),rhs.nlo(),rhs.nhi())
      { 
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	if (S == rhs.stor() && rhs.CanLinearize()) {
	  TMVAssert(stepi() == rhs.stepi());
	  TMVAssert(stepj() == rhs.stepj());
	  TMVAssert(ls() == rhs.ConstLinearView().size());
	  LinearView() = rhs.ConstLinearView();
	} else {
	  Copy(rhs,View()); 
	}
      }

      template <class T2> inline BandMatrix(
	  const GenBandMatrix<T2>& rhs, int lo, int hi) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),lo,hi)
      { 
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo <= rhs.nlo());
	TMVAssert(hi <= rhs.nhi());
	if (S == rhs.stor() && rhs.CanLinearize() && 
	    lo==rhs.nlo() && hi==rhs.nhi()) {
	  TMVAssert(stepi() == rhs.stepi());
	  TMVAssert(stepj() == rhs.stepj());
	  TMVAssert(ls() == rhs.ConstLinearView().size());
	  LinearView() = rhs.ConstLinearView();
	} else {
	  Copy(BandMatrixViewOf(rhs,lo,hi),View()); 
	  if (lo > rhs.nlo()) Diags(-nlo(),-rhs.nlo()).Zero();
	  if (hi > rhs.nhi()) Diags(rhs.nhi()+1,nhi()+1).Zero();
	}
      }

      inline BandMatrix(const GenMatrix<T>& rhs, int lo, int hi) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),lo,hi)
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rhs.colsize()));
	TMVAssert(hi < int(rhs.rowsize()));
	Copy(BandMatrixViewOf(rhs,lo,hi),View());
      }

      inline BandMatrix(const BandMatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.colsize(),mcomp.rowsize(),mcomp.nlo(),mcomp.nhi())
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	mcomp.AssignTo(View()); 
      }

#undef NEW_SIZE

      inline ~BandMatrix() { }


      //
      // Op=
      //

      inline BandMatrix<T,S,I>& operator=(const BandMatrix<T,S,I>& m2)
      {
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	if (&m2 != this) {
	  if (nlo() > m2.nlo()) Diags(-m2.nlo(),-nlo()).Zero();
	  if (nhi() > m2.nhi()) Diags(nhi()+1,m2.nhi()+1).Zero();
	  if (S == DiagMajor)
	    if (nlo() > m2.nlo())
	      memmove(itsm+m2.nlo()*stepi(),m2.itsm1.get(),m2.ls()*sizeof(T));
	    else 
	      memmove(itsm1.get(),m2.itsm1.get(),m2.ls()*sizeof(T));
	  else if (nlo() == m2.nlo() && nhi() == m2.nhi())
	    memmove(itsm1.get(),m2.itsm1.get(),m2.ls()*sizeof(T));
	  else
	    Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(const GenBandMatrix<T>& m2)
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	if (S == m2.stor() && m2.CanLinearize() && 
	    nlo()==m2.nlo() && nhi()==m2.nhi()) {
	  TMVAssert(stepi() == m2.stepi());
	  TMVAssert(stepj() == m2.stepj());
	  TMVAssert(ls() == m2.ConstLinearView().size());
	  LinearView() = m2.ConstLinearView();
	} else {
	  Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	  if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
	  if (nhi() > m2.nhi()) Diags(m2.nhi()+1,nhi()+1).Zero();
	}
	return *this;
      }

      template <class T2> inline BandMatrix<T,S,I>& operator=(
	  const GenBandMatrix<T2>& m2)
      { 
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	if (S == m2.stor() && m2.CanLinearize() && 
	    nlo()==m2.nlo() && nhi()==m2.nhi()) {
	  TMVAssert(stepi() == m2.stepi());
	  TMVAssert(stepj() == m2.stepj());
	  TMVAssert(ls() == m2.ConstLinearView().size());
	  LinearView() = m2.ConstLinearView();
	} else {
	  Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	  if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
	  if (nhi() > m2.nhi()) Diags(m2.nhi()+1,nhi()+1).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(T x) 
      {
	TMVAssert(IsSquare());
	return SetToIdentity(x); 
      }

      inline BandMatrix<T,S,I>& operator=(const BandMatrixComposite<T>& mcomp)
      {
	TMVAssert(colsize() == mcomp.colsize());
	TMVAssert(rowsize() == mcomp.rowsize());
	TMVAssert(nlo() >= mcomp.nlo());
	TMVAssert(nhi() >= mcomp.nhi());
	mcomp.AssignTo(SubBandMatrix(0,colsize(),0,rowsize(),
	      mcomp.nlo(),mcomp.nhi()));
	if (nlo() > mcomp.nlo()) Diags(-nlo(),-mcomp.nlo()).Zero();
	if (nhi() > mcomp.nhi()) Diags(mcomp.nhi()+1,nhi()+1).Zero();
	return *this;
      }

      //
      // Access
      //

      inline T operator()(size_t i,size_t j) const
      { 
	if (I==CStyle) {
	  TMVAssert(i < colsize());
	  TMVAssert(j < rowsize());
	  if (okij(i,j)) return cref(i,j); 
	  else return T(0);
	} else {
	  TMVAssert(i>0 && i <= colsize());
	  TMVAssert(j>0 && j <= rowsize());
	  if (okij(i-1,j-1)) return cref(i-1,j-1); 
	  else return T(0);
	}
      }

      inline T& operator()(size_t i,size_t j) 
      { 
	if (I==CStyle) {
	  TMVAssert(i < colsize());
	  TMVAssert(j < rowsize());
	  TMVAssert(okij(i,j));
	  return ref(i,j); 
	} else {
	  TMVAssert(i>0 && i <= colsize());
	  TMVAssert(j>0 && j <= rowsize());
	  TMVAssert(okij(i-1,j-1));
	  return ref(i-1,j-1);
	}
      }

      inline ConstVectorView<T,I> row(size_t i, size_t j1, size_t j2) const
      { 
	if (I == FortranStyle) { TMVAssert(i>0 && j1>0 && j1<=j2); }
	const size_t ix = (I == CStyle ? i : i-1);
	const size_t j1x = (I == CStyle ? j1 : j1-1);
	TMVAssert(ix<colsize());
	TMVAssert(j1x<=j2 && j2<=rowsize());
	TMVAssert(j1x==j2 || okij(ix,j1x));
	TMVAssert(j1x==j2 || okij(ix,j2-1));
	return ConstVectorView<T,I>(cptr()+ix*stepi()+j1x*stepj(),
	    j2-j1x,stepj(),NonConj);
      }

      inline ConstVectorView<T,I> col(size_t j, size_t i1, size_t i2) const
      {
	if (I == FortranStyle) { TMVAssert(j>0 && i1>0 && i1<=i2); }
	const size_t jx = (I == CStyle ? j : j-1);
	const size_t i1x = (I == CStyle ? i1 : i1-1);
	TMVAssert(jx<rowsize());
	TMVAssert(i1x<=i2 && i2<=colsize());
	TMVAssert(i1x==i2 || okij(i1x,jx));
	TMVAssert(i1x==i2 || okij(i2-1,jx));
	return ConstVectorView<T,I>(cptr()+i1x*stepi()+jx*stepj(),
	    i2-i1x,stepi(),NonConj);
      }

      inline ConstVectorView<T,I> diag() const
      {
	return ConstVectorView<T,I>(cptr(),
	    min(colsize(),rowsize()),diagstep(),NonConj);
      }

      inline ConstVectorView<T,I> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T,I>(cptr()+i*stepj(),
	      diagsize,diagstep(),NonConj);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return ConstVectorView<T,I>(cptr()-i*stepi(),
	      diagsize,diagstep(),NonConj);
	}
      }

      inline ConstVectorView<T,I> diag(int i, size_t j1, size_t j2) const
      {
	if (I==FortranStyle) { TMVAssert(j1>0); }
	const size_t j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  TMVAssert(j2<=min(rowsize()-i,colsize()));
	  return ConstVectorView<T,I>(cptr()+i*stepj()+j1x*diagstep(),
	      j2-j1x,diagstep(),NonConj);
	} else {
	  TMVAssert(j2<=min(colsize()+i,rowsize()));
	  return ConstVectorView<T,I>(cptr()-i*stepi()+j1x*diagstep(),
	      j2-j1x,diagstep(),NonConj);
	}
      }

      inline VectorView<T,I> row(size_t i, size_t j1, size_t j2)
      { 
	if (I == FortranStyle) { TMVAssert(i>0 && j1>0 && j1<=j2); }
	const size_t ix = (I == CStyle ? i : i-1);
	const size_t j1x = (I == CStyle ? j1 : j1-1);
	TMVAssert(ix<colsize());
	TMVAssert(j1x<=j2 && j2<=rowsize());
	TMVAssert(j1x==j2 || okij(ix,j1x));
	TMVAssert(j1x==j2 || okij(ix,j2-1));
	return VectorView<T,I>(ptr()+ix*stepi()+j1x*stepj(),
	    j2-j1x,stepj(),NonConj FIRSTLAST);
      }

      inline VectorView<T,I> col(size_t j, size_t i1, size_t i2)
      {
	if (I == FortranStyle) { TMVAssert(j>0 && i1>0 && i1<=i2); }
	const size_t jx = (I == CStyle ? j : j-1);
	const size_t i1x = (I == CStyle ? i1 : i1-1);
	TMVAssert(jx<rowsize());
	TMVAssert(i1x<=i2 && i2<=colsize());
	TMVAssert(i1x==i2 || okij(i1x,jx));
	TMVAssert(i1x==i2 || okij(i2-1,jx));
	return VectorView<T,I>(ptr()+i1x*stepi()+jx*stepj(),
	    i2-i1x,stepi(),NonConj FIRSTLAST);
      }

      inline VectorView<T,I> diag()
      {
	return VectorView<T,I>(ptr(),
	    min(colsize(),rowsize()),diagstep(),NonConj FIRSTLAST);
      }

      inline VectorView<T,I> diag(int i)
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return VectorView<T,I>(ptr()+i*stepj(),
	      diagsize,diagstep(),NonConj FIRSTLAST);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return VectorView<T,I>(ptr()-i*stepi(),
	      diagsize,diagstep(),NonConj FIRSTLAST);
	}
      }

      inline VectorView<T,I> diag(int i, size_t j1, size_t j2)
      {
	if (I==FortranStyle) { TMVAssert(j1>0); }
	const size_t j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  TMVAssert(j2<=min(rowsize()-i,colsize()));
	  return VectorView<T,I>(ptr()+i*stepj()+j1x*diagstep(),
	      j2-j1x,diagstep(),NonConj FIRSTLAST);
	} else {
	  TMVAssert(j2<=min(colsize()+i,rowsize()));
	  return VectorView<T,I>(ptr()-i*stepi()+j1x*diagstep(),
	      j2-j1x,diagstep(),NonConj FIRSTLAST);
	}
      }

      //
      // Modifying Functions
      //

      inline BandMatrix<T,S,I>& Zero() { return SetAllTo(0); }

      inline BandMatrix<T,S,I>& SetAllTo(T x) 
      { LinearView().SetAllTo(x); return *this; }

      inline BandMatrix<T,S,I>& Clip(RealType(T) thresh) 
      { LinearView().Clip(thresh); return *this; }

      inline BandMatrix<T,S,I>& TransposeSelf() 
      { 
	TMVAssert(IsSquare());
	TMVAssert(nlo()==nhi());
	View().TransposeSelf(); 
	return *this; 
      }

      inline BandMatrix<T,S,I>& ConjugateSelf() 
      { LinearView().ConjugateSelf(); return *this; }

      inline BandMatrix<T,S,I>& SetToIdentity(T x=T(1)) 
      { 
	TMVAssert(IsSquare());
	Zero(); diag().SetAllTo(x); 
	return *this; 
      }

      //
      // SubBandMatrix
      //

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(GenBandMatrix<T>::OKSubMatrix(i1x,i2,j1x,j2,1,1));
	return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, stepi(), stepj(), S, NonConj);
      }

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	TMVAssert(GenBandMatrix<T>::OKSubMatrix(i1x,i2x,j1x,j2x,istep,jstep));
	return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, NonConj);
      }

      inline ConstVectorView<T,I> SubVector(
	  size_t i, size_t j, int istep, int jstep, size_t size) const
      {
	const int ix = (I==CStyle ? i : i-1);
	const int jx = (I==CStyle ? j : j-1);
	TMVAssert(GenBandMatrix<T>::OKSubVector(ix,jx,istep,jstep,size));
	return ConstVectorView<T,I>(cptr()+ix*stepi()+jx*stepj(),size,
	    istep*stepi()+jstep*stepj(), NonConj);
      }

      inline ConstBandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1x,i2,j1x,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), S, NonConj);
      }

      inline ConstBandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : 
	  istep == 1 && jstep == 1 ? DiagMajor : NoMajor;
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1x,i2x,j1x,j2x,newnlo,newnhi,istep,jstep));
	return ConstBandMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    (i2x-i1x)/istep, (j2x-j1x)/jstep, newnlo, newnhi, 
	    newstepi, newstepj, newstepi+newstepj, newstor, NonConj);
      }

      inline ConstBandMatrixView<T,I> Rows(int i1, int i2) const
      {
	if (I==FortranStyle) { TMVAssert(i1>0); TMVAssert(i1<=i2); }
	else { TMVAssert(i1 >= 0); TMVAssert(i1 < i2); }
	TMVAssert(i2 <= int(colsize()));

	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1 = i1x > nlo() ? i1x-nlo() : 0;
	const int j2 = min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1x < nlo() ? min(nlo(),i2-1) - i1x : 0;
	const int newnhi = min(nlo()+nhi()-newnlo,j2-j1-1);
	const size_t newlin = (ls() && isrm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1x,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T,I>(cptr()+i1x*stepi()+j1*stepj(),
	    i2-i1x, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin);
      }

      inline ConstBandMatrixView<T,I> Cols(int j1, int j2) const
      {
	if (I==FortranStyle) { TMVAssert(j1>0); TMVAssert(j1 <= j2); }
	else { TMVAssert(j1 >= 0); TMVAssert(j1 < j2); }
	TMVAssert(j2 <= int(rowsize()));

	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i1 = j1x > nhi() ? j1x-nhi() : 0;
	const int i2 = min(j2 + nlo(),int(colsize()));
	const int newnhi = j1x < nhi() ? min(nhi(),j2-1) - j1x : 0;
	const int newnlo = min(nlo()+nhi()-newnhi,i2-i1-1);
	const size_t newlin = (ls() && iscm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1x,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T,I>(cptr()+i1*stepi()+j1x*stepj(),
	    i2-i1, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin);
      }

      inline ConstBandMatrixView<T,I> Diags(int k1, int k2) const
      {
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 < 0 ? -k2 : 0;
	const int i2 = min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = min(rowsize(),colsize()+k2);
	const int newnlo = k2 < 0 ? k2-k1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k1 > 0 ? k2-k1-1 : k2 > 0 ? k2-1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T,I>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct());
      }
	
      inline ConstBandMatrixView<T,I> UpperBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    min(colsize(),rowsize()),min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct());
      }

      inline ConstBandMatrixView<T,I> LowerBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    min(colsize(),rowsize()+nlo()),min(colsize(),rowsize()),
	    nlo(),0,stepi(),stepj(),diagstep(),stor(),ct());
      }
	
      inline ConstBandMatrixView<RealType(T),I> Real() const
      {
	return ConstBandMatrixView<RealType(T),I>(
	    reinterpret_cast<const RealType(T)*>(cptr()),
	    colsize(),rowsize(),nlo(),nhi(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? diagstep() : 2*diagstep(),
	    IsReal(T()) ? S : NoMajor, NonConj);
      }

      inline ConstBandMatrixView<RealType(T),I> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return ConstBandMatrixView<RealType(T),I>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    colsize(),rowsize(),nlo(),nhi(),
	    2*stepi(),2*stepj(),2*diagstep(),NoMajor,NonConj);
      }

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
      {
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(GenBandMatrix<T>::OKSubMatrix(i1x,i2,j1x,j2,1,1));
	return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, stepi(), stepj(), S, NonConj FIRSTLAST);
      }

      inline MatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep)
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	TMVAssert(GenBandMatrix<T>::OKSubMatrix(i1x,i2x,j1x,j2x,istep,jstep));
	return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, NonConj FIRSTLAST);
      }

      inline VectorView<T,I> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size)
      {
	const int ix = (I==CStyle ? i : i-1);
	const int jx = (I==CStyle ? j : j-1);
	TMVAssert(GenBandMatrix<T>::OKSubVector(ix,jx,istep,jstep,size));
	return VectorView<T,I>(ptr()+ix*stepi()+jx*stepj(),size,
	    istep*stepi()+jstep*stepj(), NonConj FIRSTLAST);
      }

      inline BandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi)
      {
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1x,i2,j1x,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), S, NonConj FIRSTLAST);
      }

      inline BandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep)
      {
	StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : 
	  istep == 1 && jstep == 1 ? DiagMajor : NoMajor;
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1x,i2x,j1x,j2x,newnlo,newnhi,istep,jstep));
	return BandMatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    (i2x-i1x)/istep, (j2x-j1x)/jstep, newnlo, newnhi, 
	    newstepi, newstepj, newstepi+newstepj, newstor, 
	    NonConj FIRSTLAST);
      }

      inline BandMatrixView<T,I> Rows(int i1, int i2)
      {
	if (I==FortranStyle) { TMVAssert(i1>0); TMVAssert(i1<=i2); }
	else { TMVAssert(i1 >= 0); TMVAssert(i1 < i2); }
	TMVAssert(i2 <= int(colsize()));

	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1 = i1x > nlo() ? i1x-nlo() : 0;
	const int j2 = min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1x < nlo() ? min(nlo(),i2-1) - i1x : 0;
	const int newnhi = min(nlo()+nhi()-newnlo,j2-j1-1);
	const size_t newlin = (ls() && isrm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1x,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1x*stepi()+j1*stepj(),
	    i2-i1x, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin FIRSTLAST);
      }

      inline BandMatrixView<T,I> Cols(int j1, int j2)
      {
	if (I==FortranStyle) { TMVAssert(j1>0); TMVAssert(j1 <= j2); }
	else { TMVAssert(j1 >= 0); TMVAssert(j1 < j2); }
	TMVAssert(j2 <= int(rowsize()));

	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i1 = j1x > nhi() ? j1x-nhi() : 0;
	const int i2 = min(j2 + nlo(),int(colsize()));
	const int newnhi = j1x < nhi() ? min(nhi(),j2-1) - j1x : 0;
	const int newnlo = min(nlo()+nhi()-newnhi,i2-i1-1);
	const size_t newlin = (ls() && iscm()) ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1x,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1x*stepj(),
	    i2-i1, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin FIRSTLAST);
      }

      inline BandMatrixView<T,I> Diags(int k1, int k2)
      {
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 < 0 ? -k2 : 0;
	const int i2 = min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = min(rowsize(),colsize()+k2);
	const int newnlo = k2 < 0 ? k2-k1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k1 > 0 ? k2-k1-1 : k2 > 0 ? k2-1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct()  FIRSTLAST);
      }
	
      inline BandMatrixView<T,I> UpperBand() 
      {
	return BandMatrixView<T,I>(ptr(),
	    min(colsize(),rowsize()),min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> LowerBand()
      {
	return BandMatrixView<T,I>(ptr(),
	    min(colsize(),rowsize()+nlo()),min(colsize(),rowsize()),
	    nlo(),0,stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      }
	
      inline BandMatrixView<RealType(T),I> Real()
      {
	return BandMatrixView<RealType(T),I>(
	    reinterpret_cast<RealType(T)*>(ptr()),
	    colsize(),rowsize(),nlo(),nhi(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? diagstep() : 2*diagstep(),
	    IsReal(T()) ? S : NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline BandMatrixView<RealType(T),I> Imag()
      {
	TMVAssert(IsComplex(T()));
	return BandMatrixView<RealType(T),I>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    colsize(),rowsize(),nlo(),nhi(),
	    2*stepi(),2*stepj(),2*diagstep(),NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // Views
      //

      // For these, we don't keep the linearview option for DiagMajor
      // BandMatrices, since the start of the vector needs to be at 
      // itsm1 in that case, rather than cptr()==itsm.  And anyway, 
      // the speed benefit of linearizing isn't nearly as good in this
      // case as it is for Row or ColMajor BandMatrices.
      inline ConstBandMatrixView<T,I> View() const
      { 
	return ConstBandMatrixView<T,I>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,NonConj,isdm()?0:ls());
      }

      inline ConstBandMatrixView<T,I> Transpose() const
      { 
	return ConstBandMatrixView<T,I>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),NonConj,isdm()?0:ls());
      }

      inline ConstBandMatrixView<T,I> Conjugate() const
      { 
	return ConstBandMatrixView<T,I>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,ConjOf(T,NonConj),isdm()?0:ls());
      }

      inline ConstBandMatrixView<T,I> Adjoint() const
      { 
	return ConstBandMatrixView<T,I>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj),
	    isdm()?0:ls());
      }

      inline ConstVectorView<T> ConstLinearView() const
      {
	return ConstVectorView<T>(itsm1.get(),ls(),1,NonConj);
      }

      inline BandMatrixView<T,I> View() 
      { 
	return BandMatrixView<T,I>(ptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,NonConj,isdm()?0:ls() FIRSTLAST);
      }

      inline BandMatrixView<T,I> Transpose()
      { 
	return BandMatrixView<T,I>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),NonConj,isdm()?0:ls() 
	    FIRSTLAST);
      }

      inline BandMatrixView<T,I> Conjugate()
      { 
	return BandMatrixView<T,I>(ptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,ConjOf(T,NonConj),isdm()?0:ls() 
	    FIRSTLAST);
      }

      inline BandMatrixView<T,I> Adjoint()
      { 
	return BandMatrixView<T,I>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj),
	    isdm()?0:ls() FIRSTLAST);
      }

      inline VectorView<T,I> LinearView()
      {
	return VectorView<T,I>(itsm1.get(),ls(),1,NonConj FIRSTLAST );
      }

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      using BaseMatrix<T>::IsSquare;
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline size_t mem_used() const { return ls(); }
      inline const T* start_mem() const { return itsm1.get(); }
      inline const T* cptr() const { return itsm; }
      inline T* ptr() { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline int diagstep() const { return itsds; }
      inline StorageType stor() const { return S; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isrm() const { return S==RowMajor; }
      inline bool iscm() const { return S==ColMajor; }
      inline bool isdm() const { return S==DiagMajor; }
      using GenBandMatrix<T>::ls;
      inline bool isconj() const { return false; }

    protected :

      auto_array<T> const itsm1;
      const size_t itscs;
      const size_t itsrs;
      const int itsnlo;
      const int itsnhi;
      const size_t itssi;
      const size_t itssj;
      const size_t itsds;
      T*const itsm;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      inline bool okij(size_t i, size_t j) const
      { return (j+nlo() >= i && i+nhi() >= j); }

      inline T cref(size_t i, size_t j) const
      {
	TMVAssert(i<colsize() && j<rowsize());
	return okij(i,j) ? *(cptr() + int(i)*stepi() + int(j)*stepj()) : T(0);
      }

      inline T& ref(size_t i, size_t j)
      {
	TMVAssert(GenBandMatrix<T>::okij(i,j));
	T* mi = ptr() + int(i)*stepi() + int(j)*stepj();
#ifdef TMVFLDEBUG
	TMVAssert(mi >= first);
	TMVAssert(mi < last);
#endif
	return *mi;
      }

  }; // BandMatrix

//---------------------------------------------------------------------------

  //
  // Special Constructors:
  //   UpperBiDiagMatrix(v1,v2)
  //        Returns a square BandMatrix with nlo=0, nhi=1, 
  //   LowerBiDiagMatrix(v1,v2)
  //        Returns a BandMatrix with nlo=1, nhi=0, 
  //   TriDiagMatrix(v1,v2,v3)
  //        Returns a BandMatrix with nhi=1, nlo=1, 
  //
  //   In all three cases, the vectors are given from bottom to top. 

  template <class T> BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2);

  template <class T> BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2);

  template <class T> BandMatrix<T,DiagMajor> TriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2,
      const GenVector<T>& v3);

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
	const GenMatrix<T>& m, int nlo, int nhi)
  { 
    return ConstBandMatrixView<T>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
      m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0); 
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstMatrixView<T,I>& m, int nlo, int nhi)
    { 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0); 
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const Matrix<T,S,I>& m, int nlo, int nhi)
    { 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0); 
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const MatrixView<T,I>& m, int nlo, int nhi)
  {  
    return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	Matrix<T,S,I>& m, int nlo, int nhi)
    {  
      return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0
	  FIRSTLAST1(m.first,m.last) ); 
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenBandMatrix<T>& m, int nlo, int nhi)
  { 
    TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
    return ConstBandMatrixView<T>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
      m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0); 
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstBandMatrixView<T,I>& m, int nlo, int nhi)
    { 
      TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0); 
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const BandMatrix<T,S,I>& m, int nlo, int nhi)
    { 
      TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0); 
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const BandMatrixView<T,I>& m, int nlo, int nhi)
  { 
    TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
    return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	BandMatrix<T,S,I>& m, int nlo, int nhi)
  { 
    TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
    return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct(),0
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenUpperTriMatrix<T>& m, int nhi)
  { 
    TMVAssert(!m.isunit());
    return ConstBandMatrixView<T>(m.cptr(),m.size(),m.size(),0,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	m.stor(),m.ct());
  }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> BandMatrixViewOf(
      const ConstUpperTriMatrixView<T,I>& m, int nhi)
  {
    TMVAssert(!m.isunit());
    return ConstBandMatrixView<T,I>(m.ptr(),m.size(),m.size(),0,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	m.stor(),m.ct());
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const UpperTriMatrix<T,D,S,I>& m, int nhi)
    {
      TMVAssert(D==NonUnitDiag);
      return ConstBandMatrixView<T,I>(m.ptr(),m.size(),m.size(),0,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	  m.stor(),m.ct());
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const UpperTriMatrixView<T,I>& m, int nhi)
  {
    TMVAssert(!m.isunit());
    return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),0,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	m.stor(),m.ct());
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	UpperTriMatrix<T,D,S,I>& m, int nhi)
    {
      TMVAssert(D==NonUnitDiag);
      return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),0,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	  m.stor(),m.ct());
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenUpperTriMatrix<T>& m)
  { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> BandMatrixViewOf(
      const ConstUpperTriMatrixView<T,I>& m)
  { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(const UpperTriMatrix<T,D,S,I>& m)
    { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const UpperTriMatrixView<T,I>& m)
  { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(UpperTriMatrix<T,D,S,I>& m)
    { return BandMatrixViewOf(m,m.size()-1); }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenLowerTriMatrix<T>& m, int nlo)
  { 
    TMVAssert(!m.isunit());
    return ConstBandMatrixView<T>(m.cptr(),m.size(),m.size(),nlo,0,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	m.stor(),m.ct());
  }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> BandMatrixViewOf(
      const ConstLowerTriMatrixView<T,I>& m, int nlo)
  {
    TMVAssert(!m.isunit());
    return ConstBandMatrixView<T,I>(m.ptr(),m.size(),m.size(),nlo,0,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	m.stor(),m.ct());
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const LowerTriMatrix<T,D,S,I>& m, int nlo)
    {
      TMVAssert(D==NonUnitDiag);
      return ConstBandMatrixView<T,I>(m.ptr(),m.size(),m.size(),nlo,0,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	  m.stor(),m.ct());
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const LowerTriMatrixView<T,I>& m, int nlo)
  {
    TMVAssert(!m.isunit());
    return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),nlo,0,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	m.stor(),m.ct());
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	LowerTriMatrix<T,D,S,I>& m, int nlo)
    {
      TMVAssert(D==NonUnitDiag);
      return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),nlo,0,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),
	  m.stor(),m.ct());
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenLowerTriMatrix<T>& m)
  { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> BandMatrixViewOf(
      const ConstLowerTriMatrixView<T,I>& m)
  { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(const LowerTriMatrix<T,D,S,I>& m)
    { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const LowerTriMatrixView<T,I>& m)
  { return BandMatrixViewOf(m,m.size()-1); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(LowerTriMatrix<T,D,S,I>& m)
    { return BandMatrixViewOf(m,m.size()-1); }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const T* m, size_t cs, size_t rs, int nlo, int nhi,
      StorageType stor)
  {
    TMVAssert(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    if (stor == DiagMajor) {
      int stepi = rs >= cs ? -int(cs)+1 : -int(rs);
      int stepj = rs >= cs ? int(cs) : int(rs)+1;
      return ConstBandMatrixView<T>(m-nlo*stepi,cs,rs,nlo,nhi,stepi,stepj,1,
	  DiagMajor,NonConj);
    } else {
      int lohi = nlo+nhi;
      if (stor == RowMajor)
	return ConstBandMatrixView<T>(m,cs,rs,nlo,nhi,lohi,1,lohi+1,
	  RowMajor,NonConj);
      else if (stor == ColMajor)
	return ConstBandMatrixView<T>(m,cs,rs,nlo,nhi,1,lohi,lohi+1,
	  ColMajor,NonConj);
    }
  }

  template <class T> inline BandMatrixView<T> BandMatrixViewOf(
      T* m, size_t cs, size_t rs, int nlo, int nhi,
      StorageType stor)
  {
    TMVAssert(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    if (stor == DiagMajor) {
      int stepi = rs >= cs ? -int(cs)+1 : -int(rs);
      int stepj = rs >= cs ? int(cs) : int(rs)+1;
      return BandMatrixView<T>(m-nlo*stepi,cs,rs,nlo,nhi,stepi,stepj,1,
	  DiagMajor,NonConj
	  FIRSTLAST1(m,m+BandStorageLength(DiagMajor,cs,rs,nlo,nhi)));
    } else {
      int lohi = nlo+nhi;
      if (stor == RowMajor)
	return BandMatrixView<T>(m,cs,rs,nlo,nhi,lohi,1,lohi+1,RowMajor,NonConj
	    FIRSTLAST1(m,m+BandStorageLength(RowMajor,cs,rs,nlo,nhi)));
      else
	return BandMatrixView<T>(m,cs,rs,nlo,nhi,1,lohi,lohi+1,ColMajor,NonConj
	    FIRSTLAST1(m,m+BandStorageLength(ColMajor,cs,rs,nlo,nhi)));
    }
  }
  
  //
  // Copy
  //

  template <class T1, class T2> inline void DoCopy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T2()));
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.nlo() == m1.nlo());
    TMVAssert(m2.nhi() == m1.nhi());
    TMVAssert(m2.colsize() > 0);
    TMVAssert(m2.rowsize() > 0);
    TMVAssert(m2.ct()==NonConj);
    TMVAssert(!m2.SameAs(m1));

    if (m1.stor() == m2.stor() && m1.CanLinearize() && m2.CanLinearize()) {
      TMVAssert(m1.ConstLinearView().size() == m2.LinearView().size());
      TMVAssert(m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj());
      m2.LinearView() = m1.ConstLinearView();
    }
    else {
      const int lo = m2.nlo();
      const int hi = m2.nhi();
      for(int k = -lo; k <= hi; ++k) m2.diag(k) = m1.diag(k);
    }
  }

  template <class T> inline void DoCopy(
      const GenBandMatrix<complex<T> >& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T1, class T2> inline void Copy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T2()));
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.nlo() == m1.nlo());
    TMVAssert(m2.nhi() == m1.nhi());
    if (m2.colsize() > 0 && m2.rowsize() > 0 && !m2.SameAs(m1)) {
      if (m2.isconj()) DoCopy(m1.Conjugate(),m2.Conjugate());
      else DoCopy(m1,m2);
    }
  }


  //
  // Swap Matrices
  //

  template <class T> inline void Swap(
      const BandMatrixView<T>& m1, const BandMatrixView<T>& m2);

  template <class T, StorageType S, IndexStyle I> inline void Swap(
      const BandMatrixView<T>& m1, BandMatrix<T,S,I>& m2)
  { Swap(m1,m2.View()); }

  template <class T, StorageType S, IndexStyle I> inline void Swap(
      BandMatrix<T,S,I>& m1, const BandMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, StorageType S1, IndexStyle I1, StorageType S2, IndexStyle I2> 
    inline void Swap(BandMatrix<T,S1,I1>& m1, BandMatrix<T,S2,I2>& m2)
    { Swap(m1.View(),m2.View()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const GenBandMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const GenBandMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const GenBandMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormSq(const GenBandMatrix<T>& m)
  { return m.NormSq(); }

  template <class T> inline RealType(T) NormF(const GenBandMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm1(const GenBandMatrix<T>& m)
  { return m.Norm1(); }

  template <class T> inline RealType(T) Norm2(const GenBandMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const GenBandMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> inline RealType(T) MaxAbsElement(
      const GenBandMatrix<T>& m)
  { return m.MaxAbsElement(); }

  template <class T> inline ConstBandMatrixView<T> Transpose(
      const GenBandMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> Transpose(
      const ConstBandMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Transpose(const BandMatrix<T,S,I>& m)
    { return m.Transpose(); }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> Transpose(
      const BandMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> Transpose(BandMatrix<T,S,I>& m)
    { return m.Transpose(); }

  template <class T> inline ConstBandMatrixView<T> Conjugate(
      const GenBandMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> Conjugate(
      const ConstBandMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Conjugate(const BandMatrix<T,S,I>& m)
    { return m.Conjugate(); }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> Conjugate(
      const BandMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> Conjugate(BandMatrix<T,S,I>& m)
    { return m.Conjugate(); }

  template <class T> inline ConstBandMatrixView<T> Adjoint(
      const GenBandMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline ConstBandMatrixView<T,I> Adjoint(
      const ConstBandMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Adjoint(const BandMatrix<T,S,I>& m)
    { return m.Adjoint(); }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> Adjoint(
      const BandMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> Adjoint(BandMatrix<T,S,I>& m)
    { return m.Adjoint(); }

  template <class T> inline QuotXB<T,T> Inverse(const GenBandMatrix<T>& m)
  { return m.Inverse(); }

  //
  // BandMatrix ==, != BandMatrix
  //

  template <class T1, class T2> bool operator==(
      const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2);

  template <class T1, class T2> inline bool operator!=(
      const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, StorageType S, IndexStyle I> istream& operator>>(
      istream& fin, auto_ptr<BandMatrix<T,S,I> >& m);

  template <class T> istream& operator>>(istream& fin,
      const BandMatrixView<T>& m);

  template <class T, StorageType S, IndexStyle I> inline istream& operator>>(
      istream& fin, BandMatrix<T,S,I>& m)
  { return fin >> m.View(); }

  template <class T, StorageType S, IndexStyle I> inline string Type(
      const BandMatrix<T,S,I>& )
  { return string("BandMatrix<")+Type(T())+","+Text(I)+","+Text(S)+">"; }

  template <class T> inline string Type(const GenBandMatrix<T>& m)
  {
    return string("GenBandMatrix<")+Type(T())+","+Text(m.ct())+
      ","+Text(m.stor())+">"; 
  }
  template <class T, IndexStyle I> inline string Type(
      const ConstBandMatrixView<T,I>& m)
  { 
    return string("ConstBandMatrixView<")+Type(T())+","+Text(I)+
      ","+Text(m.ct())+","+Text(m.stor())+">"; 
  }
  template <class T, IndexStyle I> inline string Type(
      const BandMatrixView<T,I>& m)
  {
    return string("BandMatrixView<")+Type(T())+","+Text(I)+","+Text(m.ct())+
      ","+Text(m.stor())+">"; 
  }

} // namespace tmv

#endif
