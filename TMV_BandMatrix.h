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
//       size_t BandStorageLength(colsize, rowsize, nlo, nhi, stor);
//
//
//    BandMatrix<T>(colsize, rowsize, nlo, nhi, stor)
//        Makes a BandMatrix with column size = row size = size 
//        with nhi non-zero superdiagonals and nlo non-zero subdiagonals
//        with _uninitialized_ values
//
//    BandMatrix<T>(const Matrix<T>& m, nlo, nhi, stor)
//        Makes a BandMatrix which copies the corresponding elements of m.
//
//    BandMatrix<T>(colsize, rowsize, nlo, nhi, const T* m,  stor)
//    BandMatrix<T>(colsize, rowsize, nlo, nhi, const valarray<T>& m, stor)
//    BandMatrix<T>(colsize, rowsize, nlo, nhi, const vector<T>& m, stor)
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
//    Inverse(m)
//    InverseATA(m)
//        For each of these inverses, when calling as a method
//        (eg. m.Inverse()), you can give the matrix to store the inverse
//        as a parameter, rather than taking it as a returned parameter.
//        ie. m.Inverse(minv) is equivalent to minv = m.Inverse, although
//        it may be more efficient.
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
//          size nhi nlo
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
//        use this form where mptr is a pointer to an undefined BandMatrix.
//        Then *mptr will be created with the right size using new,
//        so you should subsequently delete it.
//
//
// Division Control Functions:
//
//    m.DivideUsing(dt) 
//    where dt is LU, QR, SV, SVF, SVS, SVU, SVV
//    (however SVF and SV are currently identical)
// 
//    LUD(), QRD() and SVD() (or equivalently SVFD()) return the 
//      corresponding Divider classes.
//
//


#ifndef TMV_BandMatrix_H
#define TMV_BandMatrix_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenBandMatrix;
  template <class T> class ConstBandMatrixView;
  template <class T> class BandMatrixView;
  template <class T, StorageType S=RowMajor> class BandMatrix;
  template <class T> class BandMatrixComposite;

  template <class T> class BandLUDiv;
  template <class T> class BandSVDiv;
  template <class T> class BandQRDiv;

  template <class T> inline StorageType BaseStorOf(
      const GenBandMatrix<T>& m)
  {
    return (m.stor()==RowMajor || m.stor()==ColMajor || m.stor()==DiagMajor) ? 
      m.stor() : DiagMajor;
  }

  template <class T1, class T2> void Copy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2);

  template <class T> class GenBandMatrix : 
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      GenBandMatrix(StorageType s, ConjItType c) : itsstor(s), itsct(c) {}
      GenBandMatrix(const GenBandMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs), itsstor(rhs.itsstor), itsct(rhs.itsct) {}
      ~GenBandMatrix() {}

      //
      // Access Functions
      //

      using BaseMatrix<T>::colsize;
      using BaseMatrix<T>::rowsize;
      using BaseMatrix<T>::IsSquare;
      virtual int nlo() const =0;
      virtual int nhi() const =0;

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	    itsct);
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	    itsct);
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),min(colsize(),rowsize()),diagstep(),
	    itsct);
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i*stepj(),diagsize,diagstep(),itsct);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return ConstVectorView<T>(cptr()-i*stepi(),diagsize,diagstep(),itsct);
	}
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
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
	  const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      inline bool SameStorageAs(const GenBandMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      template <class T2> inline bool SameAs(const GenBandMatrix<T2>& m2) const
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
	m2.Zero();
	Copy(*this,BandMatrixViewOf(m2,nlo(),nhi()));
      }

      inline void CopyToMatrix(const MatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	m2.Zero();
	Copy(*this,BandMatrixViewOf(m2,nlo(),nhi()));
      }

      //
      // SubBandMatrix
      //

      bool OKSubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const;

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), itsstor, itsct);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const 
      {
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
	    newstor,itsct);
      }

      bool OKSubVector(size_t i, size_t j, int istep, int jstep,
	  size_t size) const;

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),itsct);
      }

      bool OKSubBandMatrix(int i1, int i2, int j1, int j2,
	  int newnlo, int newnhi, int istep, int jstep) const;

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
	const int j2 = min(i1 + nhi() + 1,int(rowsize()));
	const int newnlo = i1 < nlo() ? nlo() - i1 : 0;
	const int newnhi = nlo()+nhi()-newnlo;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }

      inline ConstBandMatrixView<T> Cols(int j1, int j2) const
      {
	TMVAssert(j1 >= 0);
	TMVAssert(j1 < j2);
	TMVAssert(j2 <= int(rowsize()));

	const int i1 = j1 > nhi() ? j1-nhi() : 0;
	const int i2 = min(j1 + nlo() + 1,int(colsize()));
	const int newnhi = j1 < nhi() ? nhi() - j1 : 0;
	const int newnlo = nlo()+nhi()-newnhi;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
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
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
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
	    stepi(),stepj(),diagstep(),itsstor,itsct);
      }

      inline ConstBandMatrixView<T> Transpose() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(itsstor),itsct);
      }

      inline ConstBandMatrixView<T> Conjugate() const
      { 
	return ConstBandMatrixView<T>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),itsstor,ConjOf(T,itsct));
      }

      inline ConstBandMatrixView<T> Adjoint() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(itsstor),ConjOf(T,itsct));
      }

      inline BandMatrixView<T> NonConst() const
      {
	return BandMatrixView<T>(const_cast<T*>(cptr()),colsize(),rowsize(),
	    nlo(),nhi(),stepi(),stepj(),diagstep(),stor(),itsct
	    FIRSTLAST1(isdm() ? diag(-nlo()).begin().GetP() : cptr(),
	      isdm() ? diag(nhi()).end().GetP() :
	      ((diag().end()-1).GetP()+1)));
      }

      //
      // Functions of BandMatrix
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

      using BaseMatrix<T>::Inverse;
      using BaseMatrix<T>::InverseATA;

      inline Matrix<T,ColMajor> Inverse() const
      { 
	Matrix<T,ColMajor> minv(rowsize(),colsize());
	Inverse(minv.View());
	return minv;
      }

      template <StorageType S> inline void Inverse(
	  const Matrix<T,S>& minv) const
      { Inverse(minv.View()); }

      inline Matrix<T,ColMajor> InverseATA() const
      { 
	size_t n = min(rowsize(),colsize());
	Matrix<T,ColMajor> minv(n,n);
	InverseATA(minv.View());
	return minv;
      }

      template <StorageType S> inline void InverseATA(
	  const Matrix<T,S>& minv) const
      { InverseATA(minv.View()); }

      inline BaseMatrix<T>* NewTranspose() const 
      { return new ConstBandMatrixView<T>(Transpose()); }

      inline BaseMatrix<T>* NewConjugate() const
      { return new ConstBandMatrixView<T>(Conjugate()); }

      inline BaseMatrix<T>* NewAdjoint() const
      { return new ConstBandMatrixView<T>(Adjoint()); }

      inline BaseMatrix<T>* NewView() const
      { return new ConstBandMatrixView<T>(View()); }

      inline BaseMatrix<T>* NewCopy() const
      { 
	if (isrm()) return new BandMatrix<T,RowMajor>(*this); 
	else if (iscm()) return new BandMatrix<T,ColMajor>(*this); 
	else return new BandMatrix<T,DiagMajor>(*this); 
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

      inline const BandSVDiv<T>& SVFD() const
      { return SVD(); }


      //
      // I/O
      //

      void WriteCompact(ostream& fout) const;
      void Write(ostream& fout) const;

      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual int diagstep() const = 0;
      virtual inline StorageType stor() const { return itsstor; }
      virtual inline bool isrm() const { return itsstor == RowMajor; }
      virtual inline bool iscm() const { return itsstor == ColMajor; }
      virtual inline bool isdm() const { return itsstor == DiagMajor; }
      virtual inline ConjItType ct() const { return itsct; }
      inline bool isconj() const
      {
	TMVAssert(IsComplex(T()) || itsct==NonConj);
	return IsComplex(T()) && itsct==Conj;
      }

    protected :

      inline bool okij(size_t i, size_t j) const
      { return (j+nlo() >= i && i+nhi() >= j); }

      T cref(size_t i, size_t j) const;

      void NewDivider() const;

    private :

      StorageType itsstor;
      ConjItType itsct;

      void operator=(const GenBandMatrix<T>&) { TMVAssert(false); }

  }; // GenBandMatrix

  template <class T> class ConstBandMatrixView : 
    public GenBandMatrix<T>
  {
    public :

      ConstBandMatrixView(const ConstBandMatrixView<T>& rhs) :
	GenBandMatrix<T>(rhs), itsm(rhs.itsm),
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd) {}

      ConstBandMatrixView(const GenBandMatrix<T>& rhs) :
	GenBandMatrix<T>(rhs), itsm(rhs.cptr()),
	itscs(rhs.colsize()), itsrs(rhs.rowsize()), 
	itsnlo(rhs.nlo()), itsnhi(rhs.nhi()),
	itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()) {}

      ConstBandMatrixView(const GenBandMatrix<T>& rhs, int _lo, int _hi) :
	GenBandMatrix<T>(rhs.stor(),rhs.ct()), itsm(rhs.cptr()),
	itscs(rhs.colsize()), itsrs(rhs.rowsize()), itsnlo(_lo), itsnhi(_hi),
	itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep())
      {
	TMVAssert(_lo <= rhs.nlo() && _hi <= rhs.nhi()); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      ConstBandMatrixView(const T* _m, size_t _cs, size_t _rs, 
	  int _lo, int _hi, int _si, int _sj, int _sd, 
	  StorageType instor, ConjItType inct) : 
	GenBandMatrix<T>(instor,inct), itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd)
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ? _si==1 :
	    instor==DiagMajor ? _sd==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      ConstBandMatrixView(const GenMatrix<T>& rhs, int _lo, int _hi) :
	GenBandMatrix<T>(rhs.stor(),rhs.ct()), itsm(rhs.cptr()), 
	itscs(rhs.colsize()), itsrs(rhs.rowsize()), itsnlo(_lo), itsnhi(_hi), 
	itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(itssi+itssj)
      { 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      ~ConstBandMatrixView() {}

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
      { TMVAssert(false); }

  }; // ConstBandMatrixView

  template <class T> class BandMatrixView : 
    public GenBandMatrix<T>
  {

    public:

      //
      // Constructors
      //

      BandMatrixView(const BandMatrixView<T>& rhs) : 
	GenBandMatrix<T>(rhs), itsm(rhs.itsm), 
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd)
	  DEFFIRSTLAST(rhs.first,rhs.last) {}

      BandMatrixView(const BandMatrixView<T>& rhs, int _lo, int _hi) :
	GenBandMatrix<T>(rhs.stor(),rhs.ct()), itsm(rhs.ptr()),
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(_lo), itsnhi(_hi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd)
	  DEFFIRSTLAST(rhs.first,rhs.last)
      { 
	TMVAssert(_lo <= rhs.nlo() && _hi <= rhs.nhi()); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      BandMatrixView(T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType instor, ConjItType inct
	  PARAMFIRSTLAST(T) ) :
	GenBandMatrix<T>(instor,inct), itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd)
	  DEFFIRSTLAST(_first,_last)
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : 
	    instor==ColMajor ? _si==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      BandMatrixView(const MatrixView<T>& rhs, int _lo, int _hi) :
	GenBandMatrix<T>(rhs.stor(),rhs.ct()), itsm(rhs.ptr()),
	itscs(rhs.colsize()), itsrs(rhs.rowsize()),
	itsnlo(_lo), itsnhi(_hi),
	itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(itssi+itssj)
	  DEFFIRSTLAST(rhs.first,rhs.last)
      { 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
      }

      ~BandMatrixView() {}

      //
      // Op=
      //

      inline const BandMatrixView<T>& operator=(
	  const BandMatrixView<T>& m2) const
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

      inline const BandMatrixView<T>& operator=(
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

      template <class T2> inline const BandMatrixView<T>& operator=(
	    const GenBandMatrix<T2>& m2) const
      {
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	for(int i=-nlo();i<-m2.nlo();++i) diag(i).Zero();
	for(int i=m2.nhi()+1;i<=nhi();++i) diag(i).Zero();
	return *this; 
      }

      inline const BandMatrixView<T>& operator=(T x) const 
      {
	TMVAssert(IsSquare());
	return SetToIdentity(x); 
      }

      inline const BandMatrixView<T>& operator=(
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
 
      inline VectorView<T> row(size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	TMVAssert(j1 <= j2);
	const size_t start = i*stepi() + j1*stepj();
	const size_t rowlen = j2 - j1;
	return VectorView<T>(ptr()+start,rowlen,stepj(),ct() FIRSTLAST );
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	TMVAssert(i1 <= i2);
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

      inline RefType(T) operator()(size_t i,size_t j) const 
      { return ref(i,j); }


      //
      // Modifying Functions
      //

      inline const BandMatrixView<T>& Zero() const 
      { return SetAllTo(T(0)); }

      const BandMatrixView<T>& Clip(RealType(T) thresh) const;

      const BandMatrixView<T>& SetAllTo(T x) const;

      inline const BandMatrixView<T>& TransposeSelf() const
      { 
	TMVAssert(IsSquare());
	TMVAssert(nlo() == nhi());
	for(int i=1;i<nhi();++i) Swap(diag(-i),diag(i));
	return *this;
      }

      const BandMatrixView<T>& ConjugateSelf() const;

      inline const BandMatrixView<T>& SetToIdentity(T x=T(1)) const 
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
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), stor(), ct() FIRSTLAST);
      }

      inline MatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const 
      {
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, ct() FIRSTLAST);
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(), ct() FIRSTLAST);
      }

      inline BandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(this->OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct() FIRSTLAST);
      }

      inline BandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(this->OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	return BandMatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
	    stepi()*istep, stepj()*jstep, newstepi+newstepj, newstor, 
	    ct() FIRSTLAST);
      }

      inline BandMatrixView<T> Rows(int i1, int i2) const
      {
	TMVAssert(i1 >= 0);
	TMVAssert(i1 < i2);
	TMVAssert(i2 <= int(colsize()));

	const int j1 = i1 > nlo() ? i1-nlo() : 0;
	const int j2 = min(i1 + nhi() + 1,int(rowsize()));
	const int newnlo = i1 < nlo() ? nlo() - i1 : 0;
	const int newnhi = nlo()+nhi()-newnlo;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }

      inline BandMatrixView<T> Cols(int j1, int j2) const
      {
	TMVAssert(j1 >= 0);
	TMVAssert(j1 < j2);
	TMVAssert(j2 <= int(rowsize()));

	const int i1 = j1 > nhi() ? j1-nhi() : 0;
	const int i2 = min(j1 + nlo() + 1,int(colsize()));
	const int newnhi = j1 < nhi() ? nhi() - j1 : 0;
	const int newnlo = nlo()+nhi()-newnhi;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }

      inline BandMatrixView<T> Diags(int k1, int k2) const
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
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }
	
      inline BandMatrixView<RealType(T)> Real() const
      {
	return BandMatrixView<RealType(T)>(
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

      inline BandMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return BandMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    colsize(),rowsize(),nlo(),nhi(),
	    2*stepi(),2*stepj(),2*diagstep(),NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      inline BandMatrixView<T> View() const
      { return *this; }

      inline BandMatrixView<T> Transpose() const
      { 
	return BandMatrixView<T>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(stor()),ct() FIRSTLAST);
      }

      inline BandMatrixView<T> Conjugate() const
      { 
	return BandMatrixView<T>(ptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),stor(),ConjOf(T,ct()) FIRSTLAST);
      }

      inline BandMatrixView<T> Adjoint() const
      { 
	return BandMatrixView<T>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(stor()),ConjOf(T,ct())
	    FIRSTLAST);
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
      using GenBandMatrix<T>::isconj;

      // This makes it easier for some things to compile.
      // The ones that work are overridden.
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

  template <StorageType S> size_t BandStorageLength(
      size_t cs, size_t rs, int lo, int hi);

  template <class T, StorageType S> class BandMatrix : 
    public GenBandMatrix<T>
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(cs,rs,lo,hi) \
      GenBandMatrix<T>(S,NonConj), \
      itslen(BandStorageLength<S>(cs,rs,lo,hi)), \
      itsm1(new T[itslen]), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
      itssi(S==RowMajor ? lo+hi : S==ColMajor ? 1 : \
	  rs>= cs ? -int(cs)+1 : -int(rs) ), \
      itssj(S==RowMajor ? 1 : S==ColMajor ? lo+hi : -itssi+1), \
      itsds(S==RowMajor ? itssi+1 : S==ColMajor ? itssj+1 : 1), \
      itsm(S==DiagMajor ? itsm1 - lo*itssi : itsm1) \
      DEFFIRSTLAST(itsm1,itsm1+itslen)

      BandMatrix(size_t cs, size_t rs, int lo, int hi) :
	  NEW_SIZE(cs,rs,lo,hi) 
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rs));
	TMVAssert(hi < int(cs));
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      BandMatrix(size_t cs, size_t rs, int lo, int hi, T x) :
	  NEW_SIZE(cs,rs,lo,hi) 
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rs));
	TMVAssert(hi < int(cs));
	SetAllTo(x);
      }

      BandMatrix(size_t cs, size_t rs, int lo, int hi, const valarray<T>& vv) :
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rs));
	TMVAssert(hi < int(cs));
	TMVAssert(vv.size() == itslen);
	T* vi=itsm1;
	for(size_t i=0;i<vv.size();++i,++vi) *vi=vv[i];
      }

      BandMatrix(size_t cs, size_t rs, int lo, int hi, const T* vv) :
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rs));
	TMVAssert(hi < int(cs));
	memmove(itsm1,vv,itslen*sizeof(T));
      }

      BandMatrix(size_t cs, size_t rs, int lo, int hi, const vector<T>& vv) : 
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rs));
	TMVAssert(hi < int(cs));
	TMVAssert(vv.size() == itslen);
	std::copy(vv.begin(),vv.end(),itsm1);
      }

      BandMatrix(const GenMatrix<T>& rhs, int lo, int hi) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),lo,hi)
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(rhs.rowsize()));
	TMVAssert(hi < int(rhs.colsize()));
	Copy(BandMatrixViewOf(rhs,lo,hi),View());
      }

      BandMatrix(const BandMatrix<T,S>& rhs) :
	GenBandMatrix<T>(S,NonConj), itslen(rhs.itslen), itsm1(new T[itslen]),
	itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itsds(rhs.itsds),
	itsm(S==DiagMajor ? itsm1-itsnlo*itssi : itsm1)
	  DEFFIRSTLAST(itsm1,itsm1+itslen)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	memmove(itsm1,rhs.itsm1,itslen*sizeof(T));
      }

      template <class T2> BandMatrix(const GenBandMatrix<T2>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),rhs.nlo(),rhs.nhi())
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	Copy(rhs,View()); 
      }

      template <class T2> BandMatrix(const GenBandMatrix<T2>& rhs, 
	  int lo, int hi) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),lo,hi)
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo <= rhs.nlo());
	TMVAssert(hi <= rhs.nhi());
	Copy(BandMatrixViewOf(rhs,lo,hi),View()); 
      }

      BandMatrix(const BandMatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.colsize(),mcomp.rowsize(),mcomp.nlo(),mcomp.nhi())
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	mcomp.AssignTo(View()); 
      }

#undef NEW_SIZE

      ~BandMatrix() { delete[] itsm1; }


      //
      // Op=
      //

      inline BandMatrix<T,S>& operator=(const BandMatrix<T,S>& m2)
      {
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	if (&m2 != this) {
	  if (S == DiagMajor && nlo() > m2.nlo())
	    memmove(itsm+m2.nlo()*stepi(),m2.itsm1,m2.itslen*sizeof(T));
	  else 
	    memmove(itsm1,m2.itsm1,m2.itslen*sizeof(T));
	  if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
	  if (nhi() > m2.nhi()) Diags(m2.nhi()+1,nhi()+1).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S>& operator=(const GenBandMatrix<T>& m2)
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
	if (nhi() > m2.nhi()) Diags(m2.nhi()+1,nhi()+1).Zero();
	return *this;
      }

      template <class T2> inline BandMatrix<T,S>& operator=(
	  const GenBandMatrix<T2>& m2)
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	Copy(m2,SubBandMatrix(0,colsize(),0,rowsize(),m2.nlo(),m2.nhi()));
	if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
	if (nhi() > m2.nhi()) Diags(m2.nhi()+1,nhi()+1).Zero();
	return *this;
      }

      inline BandMatrix<T,S>& operator=(T x) 
      {
	TMVAssert(IsSquare());
	return SetToIdentity(x); 
      }

      inline BandMatrix<T,S>& operator=(const BandMatrixComposite<T>& mcomp)
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

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	    NonConj);
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	    NonConj);
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),min(colsize(),rowsize()),diagstep(),
	    NonConj);
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i*stepj(),diagsize,diagstep(),
	      NonConj);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return ConstVectorView<T>(cptr()-i*stepi(),diagsize,diagstep(),
	      NonConj);
	}
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  TMVAssert(j2<=min(rowsize()-i,colsize()));
	  return ConstVectorView<T>(cptr()+i*stepj()+j1*diagstep(),
	      j2-j1,diagstep(),NonConj);
	} else {
	  TMVAssert(j2<=min(colsize()+i,rowsize()));
	  return ConstVectorView<T>(cptr()-i*stepi()+j1*diagstep(),
	      j2-j1,diagstep(),NonConj);
	}
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2)
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	return VectorView<T>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	    NonConj FIRSTLAST);
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2)
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	return VectorView<T>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	    NonConj FIRSTLAST);
      }

      inline VectorView<T> diag()
      {
	return VectorView<T>(ptr(),min(colsize(),rowsize()),diagstep(),NonConj
	    FIRSTLAST);
      }

      inline VectorView<T> diag(int i)
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return VectorView<T>(ptr()+i*stepj(),diagsize,diagstep(),
	      NonConj FIRSTLAST);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return VectorView<T>(ptr()-i*stepi(),diagsize,diagstep(),
	      NonConj FIRSTLAST);
	}
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2)
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  TMVAssert(j2<=min(rowsize()-i,colsize()));
	  return VectorView<T>(ptr()+i*stepj()+j1*diagstep(),
	      j2-j1,diagstep(),NonConj FIRSTLAST);
	} else {
	  TMVAssert(j2<=min(colsize()+i,rowsize()));
	  return VectorView<T>(ptr()-i*stepi()+j1*diagstep(),
	      j2-j1,diagstep(),NonConj FIRSTLAST);
	}
      }

      inline T& operator()(size_t i,size_t j) 
      { return ref(i,j); }

      //
      // Modifying Functions
      //

      inline BandMatrix<T,S>& Zero() { return SetAllTo(0); }
      inline BandMatrix<T,S>& SetAllTo(T x) 
      { std::fill(itsm1,itsm1+itslen,x); return *this; }
      inline BandMatrix<T,S>& Clip(RealType(T) thresh) 
      { VectorViewOf(ptr(),itslen).Clip(thresh); return *this; }
      inline BandMatrix<T,S>& TransposeSelf() 
      { 
	TMVAssert(IsSquare());
	TMVAssert(nlo()==nhi());
	View().TransposeSelf(); 
	return *this; 
      }
      inline BandMatrix<T,S>& ConjugateSelf() 
      { VectorViewOf(ptr(),itslen).ConjugateSelf(); return *this; }
      inline BandMatrix<T,S>& SetToIdentity(T x=T(1)) 
      { 
	TMVAssert(IsSquare());
	Zero(); diag().SetAllTo(x); 
	return *this; 
      }

      //
      // SubBandMatrix
      //

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), S, NonConj);
      }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, NonConj);
      }

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(), NonConj);
      }

      inline ConstBandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(this->OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), S, NonConj);
      }

      inline ConstBandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(this->OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
	StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : 
	  istep == 1 && jstep == 1 ? DiagMajor : NoMajor;
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
	    newstepi, newstepj, newstepi+newstepj, newstor, NonConj);
      }

      inline ConstBandMatrixView<T> Rows(int i1, int i2) const
      {
	TMVAssert(i1 >= 0);
	TMVAssert(i1 < i2);
	TMVAssert(i2 <= int(colsize()));

	const int j1 = i1 > nlo() ? i1-nlo() : 0;
	const int j2 = min(i1 + nhi() + 1,int(rowsize()));
	const int newnlo = i1 < nlo() ? nlo() - i1 : 0;
	const int newnhi = nlo()+nhi()-newnlo;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }

      inline ConstBandMatrixView<T> Cols(int j1, int j2) const
      {
	TMVAssert(j1 >= 0);
	TMVAssert(j1 < j2);
	TMVAssert(j2 <= int(rowsize()));

	const int i1 = j1 > nhi() ? j1-nhi() : 0;
	const int i2 = min(j1 + nlo() + 1, int(colsize()));
	const int newnhi = j1 < nhi() ? nhi() - j1 : 0;
	const int newnlo = nlo()+nhi()-newnhi;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
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
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }
	
      inline ConstBandMatrixView<RealType(T)> Real() const
      {
	return ConstBandMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),
	    colsize(),rowsize(),nlo(),nhi(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? diagstep() : 2*diagstep(),
	    IsReal(T()) ? S : NoMajor, NonConj);
      }

      inline ConstBandMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return ConstBandMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    colsize(),rowsize(),nlo(),nhi(),
	    2*stepi(),2*stepj(),2*diagstep(),NoMajor,NonConj);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2)
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), S, NonConj FIRSTLAST);
      }

      inline MatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep)
      {
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, NonConj FIRSTLAST);
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size)
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(), NonConj FIRSTLAST);
      }

      inline BandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi)
      {
	TMVAssert(this->OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), S, NonConj FIRSTLAST);
      }

      inline BandMatrixView<T> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep)
      {
	TMVAssert(this->OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
	StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : 
	  istep == 1 && jstep == 1 ? DiagMajor : NoMajor;
	const int newstepi = stepi()*istep;
	const int newstepj = stepj()*jstep;
	return BandMatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
	    newstepi, newstepj, newstepi+newstepj, newstor, 
	    NonConj FIRSTLAST);
      }

      inline BandMatrixView<T> Rows(int i1, int i2)
      {
	TMVAssert(i1 >= 0);
	TMVAssert(i1 < i2);
	TMVAssert(i2 <= int(colsize()));

	const int j1 = i1 > nlo() ? i1-nlo() : 0;
	const int j2 = min(i1 + nhi() + 1,int(rowsize()));
	const int newnlo = i1 < nlo() ? nlo() - i1 : 0;
	const int newnhi = nlo()+nhi()-newnlo;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }

      inline BandMatrixView<T> Cols(int j1, int j2)
      {
	TMVAssert(j1 >= 0);
	TMVAssert(j1 < j2);
	TMVAssert(j2 <= int(rowsize()));

	const int i1 = j1 > nhi() ? j1-nhi() : 0;
	const int i2 = min(j1 + nlo() + 1,int(colsize()));
	const int newnhi = j1 < nhi() ? nhi() - j1 : 0;
	const int newnlo = nlo()+nhi()-newnhi;
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }

      inline BandMatrixView<T> Diags(int k1, int k2)
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
	return SubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
      }
	
      inline BandMatrixView<RealType(T)> Real()
      {
	return BandMatrixView<RealType(T)>(
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

      inline BandMatrixView<RealType(T)> Imag()
      {
	TMVAssert(IsComplex(T()));
	return BandMatrixView<RealType(T)>(
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

      inline ConstBandMatrixView<T> View() const
      { 
	return ConstBandMatrixView<T>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,NonConj);
      }

      inline ConstBandMatrixView<T> Transpose() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),NonConj);
      }

      inline ConstBandMatrixView<T> Conjugate() const
      { 
	return ConstBandMatrixView<T>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,ConjOf(T,NonConj));
      }

      inline ConstBandMatrixView<T> Adjoint() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj));
      }

      inline BandMatrixView<T> View() 
      { 
	return BandMatrixView<T>(ptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      }

      inline BandMatrixView<T> Transpose()
      { 
	return BandMatrixView<T>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),NonConj FIRSTLAST);
      }

      inline BandMatrixView<T> Conjugate()
      { 
	return BandMatrixView<T>(ptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),S,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline BandMatrixView<T> Adjoint()
      { 
	return BandMatrixView<T>(ptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj)
	    FIRSTLAST);
      }

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      using BaseMatrix<T>::IsSquare;
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline size_t mem_used() const { return itslen; }
      inline const T* start_mem() const { return itsm1; }
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
      inline bool isconj() const { return false; }

    protected :

      size_t itslen;
      T*const itsm1;
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
	TMVAssert(this->okij(i,j));
	T* mi = ptr() + int(i)*stepi() + int(j)*stepj();
#ifdef TMVFLDEBUG
	TMVAssert(mi >= first);
	TMVAssert(mi < last);
#endif
	return *mi;
      }

  }; // BandMatrix<RowMajor>

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
  { return ConstBandMatrixView<T>(m,nlo,nhi); }

  template <class T> inline BandMatrixView<T> BandMatrixViewOf(
      const MatrixView<T>& m, int nlo, int nhi)
  { return BandMatrixView<T>(m,nlo,nhi); }

  template <class T, StorageType S> inline BandMatrixView<T> BandMatrixViewOf(
      Matrix<T,S>& m, int nlo, int nhi)
  { return BandMatrixView<T>(m.View(),nlo,nhi); }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenBandMatrix<T>& m, int nlo, int nhi)
  { return ConstBandMatrixView<T>(m,nlo,nhi); }

  template <class T> inline BandMatrixView<T> BandMatrixViewOf(
      const BandMatrixView<T>& m, int nlo, int nhi)
  { return BandMatrixView<T>(m,nlo,nhi); }

  template <class T, StorageType S> inline BandMatrixView<T> BandMatrixViewOf(
      BandMatrix<T,S>& m, int nlo, int nhi)
  { return BandMatrixView<T>(m.View(),nlo,nhi); }

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
	  FIRSTLAST1(m,m+BandStorageLength<DiagMajor>(cs,rs,nlo,nhi)));
    } else {
      int lohi = nlo+nhi;
      if (stor == RowMajor)
	return BandMatrixView<T>(m,cs,rs,nlo,nhi,lohi,1,lohi+1,RowMajor,NonConj
	    FIRSTLAST1(m,m+BandStorageLength<RowMajor>(cs,rs,nlo,nhi)));
      else
	return BandMatrixView<T>(m,cs,rs,nlo,nhi,1,lohi,lohi+1,ColMajor,NonConj
	    FIRSTLAST1(m,m+BandStorageLength<ColMajor>(cs,rs,nlo,nhi)));
    }
  }
  
  //
  // Copy
  //

  template <class T1, class T2> inline void DoCopy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.nlo() == m1.nlo());
    TMVAssert(m2.nhi() == m1.nhi());
    TMVAssert(m2.ct()==NonConj);
    TMVAssert(!(m1.isrm() && m2.isrm()));

    if (m1.iscm() && m2.iscm()) {
      size_t i1=0;
      size_t i2=m1.nlo()+1;
      size_t len=m1.nlo()+1;
      size_t k=m1.nhi();
      const T1* m1ptr=m1.cptr();
      T2* m2ptr=m2.ptr();
      for(size_t j=0;j<m1.rowsize();++j) {
	if (m1.isconj())
	  DoCopy(CVIt<T1,Unit,Conj>(m1ptr,1),
	      VIt<T2,Unit,NonConj>(m2ptr,1 FIRSTLAST1(m2.first,m2.last)),len);
	else
	  DoCopy(CVIt<T1,Unit,NonConj>(m1ptr,1),
	      VIt<T2,Unit,NonConj>(m2ptr,1 FIRSTLAST1(m2.first,m2.last)),len);
	if (k>0) { --k; ++len; m1ptr+=m1.stepj(); m2ptr+=m2.stepj(); } 
	else { ++i1; m1ptr+=m1.diagstep(); m2ptr+=m2.diagstep(); }
	if (i2<m1.colsize()) ++i2;
	else { --len; if (i1==m1.colsize()) break; }
      }
    } else {
      for(int i=-m1.nlo();i<=m1.nhi();++i) {
	m2.diag(i) = m1.diag(i);
      }
    }
  }

  template <class T1, class T2> inline void CallCopy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    if (m2.isrm() || (m1.isrm() && !m2.iscm()))
      DoCopy(m1.Transpose(),m2.Transpose());
    else DoCopy(m1,m2);
  }
  template <class T> inline void CallCopy(
      const GenBandMatrix<T>& m1, const BandMatrixView<complex<T> >& m2)
  {
    m2.Imag().Zero();
    CallCopy(m1,m2.Real());
  }
  template <class T> inline void CallCopy(
      const GenBandMatrix<complex<T> >& m1, const BandMatrixView<T>& m2)
  { TMVAssert(false); }

  template <class T1, class T2> inline void Copy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.nlo() == m1.nlo());
    TMVAssert(m2.nhi() == m1.nhi());
    if (m2.isconj()) Copy(m1.Conjugate(),m2.Conjugate());
    else CallCopy(m1,m2);
  }


  //
  // Swap Matrices
  //

  template <class T> inline void Swap(
      const BandMatrixView<T>& m1, const BandMatrixView<T>& m2);

  template <class T, StorageType S> inline void Swap(
      const BandMatrixView<T>& m1, BandMatrix<T,S>& m2)
  { Swap(m1,m2.View()); }

  template <class T, StorageType S> inline void Swap(
      BandMatrix<T,S>& m1, const BandMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, StorageType S1, StorageType S2> inline void Swap(
      BandMatrix<T,S1>& m1, BandMatrix<T,S2>& m2)
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

  template <class T> inline RealType(T) MaxAbsElement(const GenBandMatrix<T>& m)
  { return m.MaxAbsElement(); }

  template <class T> inline Matrix<T,ColMajor> Inverse(
      const GenBandMatrix<T>& m)
  { return m.Inverse(); }

  template <class T> inline Matrix<T,ColMajor> InverseATA(
      const GenBandMatrix<T>& m)
  { return m.InverseATA(); }

  template <class T> inline ConstBandMatrixView<T> Transpose(
      const GenBandMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline BandMatrixView<T> Transpose(
      const BandMatrixView<T>& m)
  { return m.Transpose(); }

  template <class T, StorageType S> inline BandMatrixView<T> Transpose(
      BandMatrix<T,S>& m)
  { return m.Transpose(); }

  template <class T> inline ConstBandMatrixView<T> Conjugate(
      const GenBandMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T> inline BandMatrixView<T> Conjugate(
      const BandMatrixView<T>& m)
  { return m.Conjugate(); }

  template <class T, StorageType S> inline BandMatrixView<T> Conjugate(
      BandMatrix<T,S>& m)
  { return m.Conjugate(); }

  template <class T> inline ConstBandMatrixView<T> Adjoint(
      const GenBandMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T> inline BandMatrixView<T> Adjoint(
      const BandMatrixView<T>& m)
  { return m.Adjoint(); }

  template <class T, StorageType S> inline BandMatrixView<T> Adjoint(
      BandMatrix<T,S>& m)
  { return m.Adjoint(); }

  //
  // BandMatrix ==, != BandMatrix
  //

  template <class T> bool operator==(
      const GenBandMatrix<T>& m1, const GenBandMatrix<T>& m2);

  template <class T> inline bool operator!=(
      const GenBandMatrix<T>& m1, const GenBandMatrix<T>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, StorageType S> istream& operator>>(istream& fin,
      BandMatrix<T,S>* m);

  template <class T> istream& operator>>(istream& fin,
      const BandMatrixView<T>& m);

  template <class T, StorageType S> istream& operator>>(istream& fin, 
      BandMatrix<T,S>& m)
  { return fin >> m.View(); }

  template <class T, StorageType S> inline std::string Type(
      const BandMatrix<T,S>& m)
  { return std::string("BandMatrix<")+Type(T())+","+Text(S)+">"; }

  template <class T> inline std::string Type(const GenBandMatrix<T>& m)
  {
    return std::string("GenBandMatrix<")+Type(T())+","+Text(m.ct())+
      ","+Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(const ConstBandMatrixView<T>& m)
  { 
    return std::string("ConstBandMatrixView<")+Type(T())+","+Text(m.ct())+
      ","+Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(const BandMatrixView<T>& m)
  {
    return std::string("BandMatrixView<")+Type(T())+","+Text(m.ct())+
      ","+Text(m.stor())+">"; 
  }

}; // namespace tmv

#endif
