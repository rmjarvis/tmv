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

#include "TMV_BaseBandMatrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include <vector>

namespace tmv {

  template <class T> class GenBandMatrix : 
    virtual public AssignableToBandMatrix<T>,
    virtual public AssignableToDiagMatrix<T>,
    virtual public AssignableToUpperTriMatrix<T>,
    virtual public AssignableToLowerTriMatrix<T>,
    public BaseMatrix<T>,
    private DivHelper<T>
  {

    public:

      //
      // Constructors
      //

      inline GenBandMatrix() {}
      inline GenBandMatrix(const GenBandMatrix<T>&) {}
      virtual inline ~GenBandMatrix() {}

      //
      // Access Functions
      //

      using AssignableToMatrix<T>::colsize;
      using AssignableToMatrix<T>::rowsize;
      using AssignableToBandMatrix<T>::nlo;
      using AssignableToBandMatrix<T>::nhi;
      // For Tri, Diag compatibility:
      inline size_t size() const 
      { TMVAssert(colsize() == rowsize()); return colsize(); }
      inline DiagType dt() const { return NonUnitDiag; }

      inline T operator()(size_t i,size_t j) const 
      { 
	TMVAssert(i< colsize());
	TMVAssert(j< rowsize());
	return okij(i,j) ? cref(i,j) : T(0);
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2 && j2<=rowsize());
	TMVAssert(j1==j2 || okij(i,j1));
	TMVAssert(j1==j2 || okij(i,j2-1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),
	    j2-j1,stepj(),ct());
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2 && i2<=colsize());
	TMVAssert(i1==i2 || okij(i1,j));
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),
	    i2-i1,stepi(),ct());
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),
	    std::min(colsize(),rowsize()),diagstep(),ct());
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	if (i >= 0) {
	  const size_t diagsize = std::min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i*stepj(),
	      diagsize,diagstep(),ct());
	} else {
	  const size_t diagsize = std::min(colsize()+i,rowsize());
	  return ConstVectorView<T>(cptr()-i*stepi(),
	      diagsize,diagstep(),ct());
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
	  TMVAssert(j2<=std::min(rowsize()-i,colsize()));
	  return ConstVectorView<T>(cptr()+i*stepj()+j1*diagstep(),
	      j2-j1,diagstep(),ct());
	} else {
	  TMVAssert(j2<=std::min(colsize()+i,rowsize()));
	  return ConstVectorView<T>(cptr()-i*stepi()+j1*diagstep(),
	      j2-j1,diagstep(),ct());
	}
      }

      template <class T2> inline bool SameAs(const GenBandMatrix<T2>& ) const
      { return false; }

      inline bool SameAs(const GenBandMatrix<T>& m2) const
      { 
	return (this==&m2 || (cptr()==m2.cptr() &&
	      colsize()==m2.colsize() && rowsize()==m2.rowsize() &&
	      stepi()==m2.stepi() && stepj()==m2.stepj() &&
	      nhi()==m2.nhi() && nlo()==m2.nlo() && isconj() == m2.isconj()));
      }

      inline void AssignToM(const MatrixView<RealType(T)>& m2) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	if (SameStorage(*this,m2)) {
	  UpperTriMatrixViewOf(m2.Cols(nhi()+1,rowsize())).Zero();
	  LowerTriMatrixViewOf(m2.Rows(nhi()+1,colsize())).Zero();
	} else {
	  m2.Zero();
	}
	AssignToB(BandMatrixViewOf(m2,nlo(),nhi()));
      }

      inline void AssignToM(const MatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	if (SameStorage(*this,m2)) {
	  UpperTriMatrixViewOf(m2.Cols(nhi()+1,rowsize())).Zero();
	  LowerTriMatrixViewOf(m2.Rows(nhi()+1,colsize())).Zero();
	} else {
	  m2.Zero();
	}
	AssignToB(BandMatrixViewOf(m2,nlo(),nhi()));
      }

      inline void AssignToB(const BandMatrixView<RealType(T)>& m2) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	TMVAssert(m2.nlo() >= nlo());
	TMVAssert(m2.nhi() >= nhi());
	if (!SameAs(m2)) Copy(*this,m2); 
      }

      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m2) const
      { 
	TMVAssert(m2.colsize() == colsize());
	TMVAssert(m2.rowsize() == rowsize());
	TMVAssert(m2.nlo() >= nlo());
	TMVAssert(m2.nhi() >= nhi());
	if (!SameAs(m2)) Copy(*this,m2); 
      }

      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m2) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m2.size() == colsize());
	TMVAssert(m2.size() == rowsize());
	TMVAssert(nlo() == 0); 
	AssignToB(BandMatrixViewOf(m2)); 
      }

      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.size() == colsize());
	TMVAssert(m2.size() == rowsize());
	TMVAssert(nlo() == 0); 
	AssignToB(BandMatrixViewOf(m2)); 
      }

      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m2) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m2.size() == colsize());
	TMVAssert(m2.size() == rowsize());
	TMVAssert(nhi() == 0); 
	AssignToB(BandMatrixViewOf(m2)); 
      }

      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.size() == colsize());
	TMVAssert(m2.size() == rowsize());
	TMVAssert(nhi() == 0); 
	AssignToB(BandMatrixViewOf(m2)); 
      }

      inline void AssignToD(const DiagMatrixView<RealType(T)>& m2) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m2.size() == colsize());
	TMVAssert(m2.size() == rowsize());
	TMVAssert(nhi() == 0 && nlo() == 0); 
	AssignToB(BandMatrixViewOf(m2)); 
      }

      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m2) const
      { 
	TMVAssert(m2.size() == colsize());
	TMVAssert(m2.size() == rowsize());
	TMVAssert(nhi() == 0 && nlo() == 0); 
	AssignToB(BandMatrixViewOf(m2)); 
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
	    i2-i1, j2-j1, stepi(), stepj(), stor(), ct());
      }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const 
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
	    newstor,ct());
      }

      bool OKSubVector(
	  int i, int j, int istep, int jstep, size_t size) const;

      inline ConstVectorView<T> SubVector(
	  int i, int j, int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),ct());
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
	    diagstep(),stor(),ct());
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
	    newstepi, newstepj, newstepi+newstepj, newstor, ct());
      }

      inline ConstBandMatrixView<T> Rows(int i1, int i2) const
      {
	TMVAssert(i1 >= 0);
	TMVAssert(i1 < i2);
	TMVAssert(i2 <= int(colsize()));

	const int j1 = i1 > nlo() ? i1-nlo() : 0;
	const int j2 = std::min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1 < nlo() ? std::min(nlo(),i2-1) - i1 : 0;
	const int newnhi = j1==j2 ? 0 : std::min(nlo()+nhi()-newnlo,j2-j1-1);
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
	const int i2 = std::min(j2 + nlo(),int(colsize()));
	const int newnhi = j1 < nhi() ? std::min(nhi(),j2-1) - j1 : 0;
	const int newnlo = i1==i2 ? 0 : std::min(nlo()+nhi()-newnhi,i2-i1-1);
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	const size_t newlin = (ls() && iscm()) ? 1 : 0;
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin);
      }

      inline ConstBandMatrixView<T> Diags(int k1, int k2) const
      {
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 <= 0 ? -k2+1 : 0;
	const int i2 = std::min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = std::min(rowsize(),colsize()+k2-1);
	const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct());
      }

      inline ConstBandMatrixView<T> UpperBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    std::min(colsize(),rowsize()),std::min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct());
      }

      inline ConstBandMatrixView<T> LowerBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    std::min(colsize(),rowsize()+nlo()),std::min(colsize(),rowsize()),
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
	    IsReal(T()) ? stor() : NoMajor, NonConj);
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
	    stepi(),stepj(),diagstep(),stor(),ct(),
	    isdm()?0:ls());
      }

      inline ConstBandMatrixView<T> Transpose() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(stor()),ct(),
	    isdm()?0:ls());
      }

      inline ConstBandMatrixView<T> Conjugate() const
      { 
	return ConstBandMatrixView<T>(cptr(),colsize(),rowsize(),nlo(),nhi(),
	    stepi(),stepj(),diagstep(),stor(),ConjOf(T,ct()),
	    isdm()?0:ls());
      }

      inline ConstBandMatrixView<T> Adjoint() const
      { 
	return ConstBandMatrixView<T>(cptr(),rowsize(),colsize(),nhi(),nlo(),
	    stepj(),stepi(),diagstep(),TransOf(stor()),ConjOf(T,ct()),
	    isdm()?0:ls());
      }

      inline BandMatrixView<T> NonConst() const
      {
	return BandMatrixView<T>(const_cast<T*>(cptr()),colsize(),rowsize(),
	    nlo(),nhi(),stepi(),stepj(),diagstep(),stor(),ct(),
	    isdm()?0:ls()
	    FIRSTLAST1(isdm() ? diag(-nlo()).begin().GetP() : cptr(),
	      isdm() ? diag(nhi()).end().GetP() :
	      ((diag().end()-1).GetP()+1)));
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
	TMVAssert(ls() != 1 || (rowsize() == 1 && colsize() == 1));
	// (To assure that next assert has no effect.)
	TMVAssert(CanLinearize());

	return ConstVectorView<T>(cptr(),ls(),1,ct());
      }

      //
      // Functions of Matrix
      //

      inline T Det() const
      { return DivHelper<T>::Det(); }

      inline T Trace() const
      { return diag().SumElements(); }

      inline RealType(T) Norm() const
      { return NormF(); }

      RealType(T) NormF() const;

      RealType(T) NormSq() const;

      RealType(T) Norm1() const;

      RealType(T) Norm2() const
      { return DivHelper<T>::Norm2(); }

      inline RealType(T) NormInf() const
      { return Transpose().Norm1(); }

      RealType(T) MaxAbsElement() const;

      QuotXB<T,T> QInverse() const;
      inline QuotXB<T,T> Inverse() const
      { return QInverse(); }

      using DivHelper<T>::Inverse;
      using DivHelper<T>::InverseATA;

      inline void Inverse(const MatrixView<T>& minv) const
      { DivHelper<T>::Inverse(minv); }

      inline void InverseATA(const MatrixView<T>& minv) const
      { DivHelper<T>::InverseATA(minv); }

      inline bool Singular() const
      { return DivHelper<T>::Singular(); }

      inline RealType(T) Condition() const
      { return DivHelper<T>::Condition(); }

      auto_ptr<BaseMatrix<T> > NewCopy() const;
      auto_ptr<BaseMatrix<T> > NewView() const;
      auto_ptr<BaseMatrix<T> > NewTranspose() const;
      auto_ptr<BaseMatrix<T> > NewConjugate() const;
      auto_ptr<BaseMatrix<T> > NewAdjoint() const;
      auto_ptr<BaseMatrix<T> > NewInverse() const;

      //
      // Division Control
      //

      using DivHelper<T>::DivideInPlace;
      using DivHelper<T>::SaveDiv;
      using DivHelper<T>::SetDiv;
      using DivHelper<T>::UnSetDiv;
      using DivHelper<T>::ReSetDiv;
      using DivHelper<T>::GetDiv;
      using DivHelper<T>::CheckDecomp;

      inline void DivideUsing(DivType dt) const
      {
	TMVAssert(dt == LU || dt == QR || SV_Type(dt));
	DivHelper<T>::DivideUsing(dt);
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

      using DivHelper<T>::LDivEq;
      using DivHelper<T>::RDivEq;
      using DivHelper<T>::LDiv;
      using DivHelper<T>::RDiv;

      //
      // I/O
      //

      void WriteCompact(std::ostream& fout) const;
      void Write(std::ostream& fout) const;
      void WriteCompact(std::ostream& fout, RealType(T) thresh) const;
      void Write(std::ostream& fout, RealType(T) thresh) const;

      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual int diagstep() const = 0;
      virtual size_t ls() const  = 0;
      virtual inline bool isrm() const { return stor() == RowMajor; }
      virtual inline bool iscm() const { return stor() == ColMajor; }
      virtual inline bool isdm() const { return stor() == DiagMajor; }
      inline bool isconj() const
      {
	TMVAssert(IsComplex(T()) || ct()==NonConj);
	return IsComplex(T()) && ct()==Conj;
      }
      virtual StorageType stor() const = 0;
      virtual ConjType ct() const = 0;

      virtual bool CanLinearize() const = 0;

    protected :

      inline bool okij(size_t i, size_t j) const
      { return (j+nlo() >= i && i+nhi() >= j); }

      virtual T cref(size_t i, size_t j) const;

      void NewDivider() const;
      inline const BaseMatrix<T>& GetMatrix() const { return *this; }

    private :

      inline void operator=(const GenBandMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenBandMatrix

  template <class T, IndexStyle I> class ConstBandMatrixView : 
    public GenBandMatrix<T>
  {
    public :

      inline ConstBandMatrixView(const ConstBandMatrixView<T,I>& rhs) :
	itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
	itsstor(rhs.itsstor), itsct(rhs.itsct), linsize(rhs.linsize) 
      { TMVAssert(!(isdm() && linsize != 0)); }

      inline ConstBandMatrixView(const GenBandMatrix<T>& rhs) :
	itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
	itsnlo(rhs.nlo()), itsnhi(rhs.nhi()),
	itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
	itsstor(rhs.stor()), itsct(rhs.ct()), 
	linsize(rhs.isdm()?0:rhs.ls()) 
      { TMVAssert(!(isdm() && linsize != 0)); }

      inline ConstBandMatrixView(
	  const T* _m, size_t _cs, size_t _rs, 
	  int _lo, int _hi, int _si, int _sj, int _sd, 
	  StorageType _stor, ConjType _ct, size_t _ls=0) : 
	itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
	itsstor(_stor), itsct(_ct), linsize(_ls)
      { 
	TMVAssert(itsstor==RowMajor ? itssj==1 : itsstor==ColMajor ? itssi==1 :
	    itsstor==DiagMajor ? itssd==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
	TMVAssert(nhi() >= 0);
	TMVAssert(nlo() >= 0);
	TMVAssert(!(isdm() && linsize!=0));
#ifdef XTEST
	TMVAssert(linsize==0 || linsize==1 || 
	    linsize==BandStorageLength(itsstor,itscs,itsrs,itsnlo,itsnhi));
#endif
      }

      virtual inline ~ConstBandMatrixView() {}

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline int diagstep() const { return itssd; }
      inline size_t ls() const { return linsize; }
      inline StorageType stor() const { return itsstor; }
      inline ConjType ct() const { return itsct; }
      using GenBandMatrix<T>::isdm;

      bool CanLinearize() const;

    protected :

      const T*const itsm;
      const size_t itscs;
      const size_t itsrs;
      const int itsnlo;
      const int itsnhi;
      const int itssi;
      const int itssj;
      const int itssd;

      const StorageType itsstor;
      const ConjType itsct;
      mutable size_t linsize;

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

	inline ConstBandMatrixView(
	    const GenBandMatrix<T>& rhs) :
	  ConstBandMatrixView<T,CStyle>(rhs) {}

	inline ConstBandMatrixView(const T* _m, size_t _cs, size_t _rs, 
	    int _lo, int _hi, int _si, int _sj, int _sd, 
	    StorageType instor, ConjType inct, size_t ls=0) : 
	  ConstBandMatrixView<T,CStyle>(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,
	      instor,inct,ls) {}

	virtual inline ~ConstBandMatrixView() {}

	//
	// Access Functions
	//

	inline T operator()(size_t i,size_t j) const 
	{ 
	  TMVAssert(i>0 && i<= colsize());
	  TMVAssert(j>0 && j<= rowsize());
	  return okij(i-1,j-1) ? cref(i-1,j-1) : T(0);
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
	    int i, int j, int istep, int jstep, size_t size) const;

	inline ConstVectorView<T,FortranStyle> SubVector(
	    int i, int j, int istep, int jstep, size_t size) const
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
	  TMVAssert(k1 <= k2);
	  TMVAssert(k2 <= nhi());
	  return GenBandMatrix<T>::Diags(k1,k2+1);
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
	itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
	itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
	itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
	itsstor(rhs.stor()), itsct(rhs.ct()), linsize(rhs.ls())
	  DEFFIRSTLAST(rhs.first,rhs.last) 
      { 
	TMVAssert(itsstor==RowMajor ? itssj==1 : itsstor==ColMajor ? itssi==1 :
	    itsstor==DiagMajor ? itssd==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
	TMVAssert(nhi() >= 0);
	TMVAssert(nlo() >= 0);
	TMVAssert(!(isdm() && linsize != 0));
	TMVAssert(itsstor==DiagMajor ? linsize == 0 : true);
#ifdef XTEST
	TMVAssert(linsize==0 || linsize==1 || 
	    linsize==BandStorageLength(itsstor,itscs,itsrs,itsnlo,itsnhi));
#endif
      }

      inline BandMatrixView(
	  T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType _stor, ConjType _ct,
	  size_t _ls PARAMFIRSTLAST(T) ) :
	itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
	itsstor(_stor), itsct(_ct), linsize(_ls)
	  DEFFIRSTLAST(_first,_last)
      { 
	TMVAssert(itsstor==RowMajor ? itssj==1 : itsstor==ColMajor ? itssi==1 :
	    itsstor==DiagMajor ? itssd==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
	TMVAssert(nhi() >= 0);
	TMVAssert(nlo() >= 0);
	TMVAssert(!(isdm() && linsize != 0));
	TMVAssert(itsstor==DiagMajor ? linsize == 0 : true);
#ifdef XTEST
	TMVAssert(linsize==0 || linsize==1 || 
	    linsize==BandStorageLength(itsstor,itscs,itsrs,itsnlo,itsnhi));
#endif
      }

      inline BandMatrixView(
	  T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	  int _si, int _sj, int _sd, StorageType _stor, ConjType _ct
	  PARAMFIRSTLAST(T) ) :
	itsm(_m), itscs(_cs), itsrs(_rs),
	itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
	itsstor(_stor), itsct(_ct), linsize(0)
	  DEFFIRSTLAST(_first,_last)
      { 
	TMVAssert(itsstor==RowMajor ? itssj==1 : itsstor==ColMajor ? itssi==1 :
	    itsstor==DiagMajor ? itssd==1 : true); 
	TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
	TMVAssert(colsize() == 0 || nlo() < int(colsize()));
	TMVAssert(nhi() >= 0);
	TMVAssert(nlo() >= 0);
	TMVAssert(!(isdm() && linsize != 0));
	TMVAssert(itsstor==DiagMajor ? linsize == 0 : true);
	TMVAssert(linsize == 0);
      }

      virtual inline ~BandMatrixView() {}

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
	m2.AssignToB(*this);
	return *this; 
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenBandMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(*this);
	return *this; 
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenBandMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(*this);
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
	if (!SameAs(m2)) Copy(m2,View());
	return *this; 
      }

      inline const BandMatrixView<T,I>& operator=(T x) const 
      {
	TMVAssert(colsize() == rowsize());
	return SetToIdentity(x); 
      }

      inline const BandMatrixView<T,I>& operator=(
	  const AssignableToBandMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(*this);
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const AssignableToBandMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(*this);
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenDiagMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	m2.AssignToD(DiagMatrixViewOf(diag()));
	if (nhi() > 0) Diags(1,nhi()+1).Zero();
	if (nlo() > 0) Diags(-nlo(),0).Zero();
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenDiagMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	m2.AssignToD(DiagMatrixViewOf(diag()));
	if (nhi() > 0) Diags(1,nhi()+1).Zero();
	if (nlo() > 0) Diags(-nlo(),0).Zero();
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenUpperTriMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nhi() == int(rowsize())-1);
	m2.AssignToU(UpperTriMatrixView<T>(ptr(),colsize(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:stor(),ct() FIRSTLAST));
	if (nlo() > 0) Diags(-nlo(),0).Zero();
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenUpperTriMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nhi() == int(rowsize())-1);
	m2.AssignToU(UpperTriMatrixView<T>(ptr(),colsize(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:stor(),ct() FIRSTLAST));
	if (nlo() > 0) Diags(-nlo(),0).Zero();
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenLowerTriMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() == int(rowsize())-1);
	m2.AssignToL(LowerTriMatrixView<T>(ptr(),colsize(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:stor(),ct() FIRSTLAST));
	if (nhi() > 0) Diags(1,nhi()+1).Zero();
	return *this;
      }

      inline const BandMatrixView<T,I>& operator=(
	  const GenLowerTriMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() == int(rowsize())-1);
	m2.AssignToL(LowerTriMatrixView<T>(ptr(),colsize(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:stor(),ct() FIRSTLAST));
	if (nhi() > 0) Diags(1,nhi()+1).Zero();
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
	return VectorView<T>(ptr(),std::min(colsize(),rowsize()),
	    diagstep(),ct() FIRSTLAST);
      }

      inline VectorView<T> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diaglen = std::min(rowsize()-i,colsize());
	  return VectorView<T>(ptr()+i*stepj(),diaglen,diagstep(),ct() 
	      FIRSTLAST );
	} else {
	  const size_t diaglen = std::min(colsize()+i,rowsize());
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
	  TMVAssert(j2<=std::min(rowsize()-i,colsize()));
	  return VectorView<T>(ptr()+i*stepj()+j1*diagstep(),
	      j2-j1, diagstep(),ct() FIRSTLAST );
	} else {
	  TMVAssert(j2<=std::min(colsize()+i,rowsize()));
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

      void DoTransposeSelf() const;
      inline const BandMatrixView<T,I>& TransposeSelf() const
      { 
	TMVAssert(colsize() == rowsize());
	TMVAssert(nlo() == nhi());
	DoTransposeSelf();
	return *this;
      }

      const BandMatrixView<T,I>& ConjugateSelf() const;

      inline const BandMatrixView<T,I>& SetToIdentity(T x=T(1)) const 
      {
	TMVAssert(colsize() == rowsize());
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
	  int i, int j, int istep, int jstep, size_t size) const
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
	const int j2 = std::min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1 < nlo() ? std::min(nlo(),i2-1) - i1 : 0;
	const int newnhi = j1==j2 ? 0 : std::min(nlo()+nhi()-newnlo,j2-j1-1);
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
	const int i2 = std::min(j2 + nlo(),int(colsize()));
	const int newnhi = j1 < nhi() ? std::min(nhi(),j2-1) - j1 : 0;
	const int newnlo = i1==i2 ? 0 : std::min(nlo()+nhi()-newnhi,i2-i1-1);
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

	const int i1 = k2 <= 0 ? -k2+1 : 0;
	const int i2 = std::min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = std::min(rowsize(),colsize()+k2-1);
	const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> UpperBand() const
      {
	return BandMatrixView<T,I>(ptr(),
	    std::min(colsize(),rowsize()),std::min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> LowerBand() const
      {
	return BandMatrixView<T,I>(ptr(),
	    std::min(colsize(),rowsize()+nlo()),std::min(colsize(),rowsize()),
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
	// (To assure that next assert has no effect.)
	TMVAssert(CanLinearize());

	return VectorView<T>(ptr(),ls(),1,ct() FIRSTLAST );
      }


      //
      // I/O
      //

      void Read(std::istream& fin) const;

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      inline int nlo() const { return itsnlo; }
      inline int nhi() const { return itsnhi; }
      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline int diagstep() const { return itssd; }
      using GenBandMatrix<T>::isrm;
      using GenBandMatrix<T>::iscm;
      using GenBandMatrix<T>::isdm;
      using GenBandMatrix<T>::isconj;
      inline size_t ls() const { return linsize; }
      inline StorageType stor() const { return itsstor; }
      inline ConjType ct() const { return itsct; }

      bool CanLinearize() const;

#if 0
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
#endif

    protected:

      T*const itsm;
      const size_t itscs;
      const size_t itsrs;
      const int itsnlo;
      const int itsnhi;
      const int itssi;
      const int itssj;
      const int itssd;

      const StorageType itsstor;
      const ConjType itsct;
      mutable size_t linsize;

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
	    int _si, int _sj, int _sd, StorageType instor, ConjType inct,
	    size_t ls PARAMFIRSTLAST(T) ) :
	  BandMatrixView<T,CStyle>(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,instor,inct,ls
	      FIRSTLAST1(_first,_last) ) {}

	inline BandMatrixView(
	    T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
	    int _si, int _sj, int _sd, StorageType instor, ConjType inct
	    PARAMFIRSTLAST(T) ) :
	  BandMatrixView<T,CStyle>(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,instor,inct
	      FIRSTLAST1(_first,_last) ) {}

	virtual inline ~BandMatrixView() {}

	//
	// Op=
	//

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const BandMatrixView<T,FortranStyle>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenBandMatrix<RealType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenBandMatrix<ComplexType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	template <class T2> 
	  inline const BandMatrixView<T,FortranStyle>& operator=(
	      const GenBandMatrix<T2>& m2) const
	  { BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(T x) const 
	{ BandMatrixView<T,CStyle>::operator=(x); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const AssignableToBandMatrix<RealType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const AssignableToBandMatrix<ComplexType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenDiagMatrix<RealType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenDiagMatrix<ComplexType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenUpperTriMatrix<RealType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenUpperTriMatrix<ComplexType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenLowerTriMatrix<RealType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const BandMatrixView<T,FortranStyle>& operator=(
	    const GenLowerTriMatrix<ComplexType(T)>& m2) const
	{ BandMatrixView<T,CStyle>::operator=(m2); return *this; }

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

	inline VectorView<T,FortranStyle> SubVector(int i, int j,
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
	{ return BandMatrixView<T,CStyle>::Diags(k1,k2+1); }

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
      linsize(BandStorageLength(S,cs,rs,lo,hi)), \
      itsm1(new T[linsize]), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
      itssi(S==RowMajor ? lo+hi : S==ColMajor ? 1 : \
	  rs>= cs ? -int(cs)+1 : -int(rs) ), \
      itssj(S==RowMajor ? 1 : S==ColMajor ? lo+hi : -itssi+1), \
      itsds(S==RowMajor ? itssi+1 : S==ColMajor ? itssj+1 : 1), \
      itsm(S==DiagMajor ? itsm1.get() - lo*itssi : itsm1.get()) \
      DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize)

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

      inline BandMatrix(size_t cs, size_t rs, int lo, int hi,
	  const std::vector<T>& vv) :
	NEW_SIZE(cs,rs,lo,hi)
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	TMVAssert(lo >= 0);
	TMVAssert(hi >= 0);
	TMVAssert(lo < int(cs));
	TMVAssert(hi < int(rs));
	TMVAssert(vv.size() == ls());
	T* vi = itsm1.get();
	typename std::vector<T>::const_iterator vvi = vv.begin();
	for(size_t i=vv.size();i>0;--i,++vi,++vvi) *vi=*vvi;
      }

      inline BandMatrix(const BandMatrix<T,S,I>& rhs) :
	linsize(rhs.ls()), itsm1(new T[linsize]),
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
	linsize(rhs.ls()), itsm1(new T[linsize]),
	itscs(rhs.colsize()), itsrs(rhs.rowsize()),
	itsnlo(rhs.nlo()), itsnhi(rhs.nhi()),
	itssi(rhs.stepi()), itssj(rhs.stepj()), itsds(rhs.diagstep()),
	itsm(S==DiagMajor ? itsm1.get()-itsnlo*itssi : itsm1.get())
	  DEFFIRSTLAST(itsm1.get(),itsm1.get()+ls())
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	memmove(itsm1.get(),rhs.start_mem(),ls()*sizeof(T));
      }

      inline BandMatrix(const GenBandMatrix<RealType(T)>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),rhs.nlo(),rhs.nhi())
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	rhs.AssignToB(View());
      }

      inline BandMatrix(const GenBandMatrix<ComplexType(T)>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),rhs.nlo(),rhs.nhi())
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	rhs.AssignToB(View());
      }

      template <class T2> inline BandMatrix(const GenBandMatrix<T2>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize(),rhs.nlo(),rhs.nhi())
      { 
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	Copy(rhs,View()); 
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
	Copy(BandMatrixViewOf(rhs,lo,hi),View()); 
	if (I==CStyle) {
	  if (lo > rhs.nlo()) Diags(-lo,-rhs.nlo()).Zero();
	  if (hi > rhs.nhi()) Diags(rhs.nhi()+1,hi+1).Zero();
	} else {
	  if (lo > rhs.nlo()) Diags(-lo,-rhs.nlo()-1).Zero();
	  if (hi > rhs.nhi()) Diags(rhs.nhi()+1,hi).Zero();
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

      inline BandMatrix(const AssignableToBandMatrix<RealType(T)>& m2) :
	NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
      {
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToB(View());
      }

      inline BandMatrix(const AssignableToBandMatrix<ComplexType(T)>& m2) :
	NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToB(View());
      }

      inline BandMatrix(const GenDiagMatrix<RealType(T)>& m2) :
	NEW_SIZE(m2.size(),m2.size(),0,0)
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToD(DiagMatrixViewOf(diag()));
      }

      inline BandMatrix(const GenDiagMatrix<ComplexType(T)>& m2) :
	NEW_SIZE(m2.size(),m2.size(),0,0)
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToD(DiagMatrixViewOf(diag()));
      }

      inline BandMatrix(const GenUpperTriMatrix<RealType(T)>& m2) :
	NEW_SIZE(m2.size(),m2.size(),0,m2.size()-1)
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToU(UpperTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
      }

      inline BandMatrix(const GenUpperTriMatrix<ComplexType(T)>& m2) :
	NEW_SIZE(m2.size(),m2.size(),0,m2.size()-1)
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToU(UpperTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
      }

      inline BandMatrix(const GenLowerTriMatrix<RealType(T)>& m2) :
	NEW_SIZE(m2.size(),m2.size(),m2.size()-1,0)
      { 
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToL(LowerTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
      }

      inline BandMatrix(const GenLowerTriMatrix<ComplexType(T)>& m2) :
	NEW_SIZE(m2.size(),m2.size(),m2.size()-1,0)
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
	m2.AssignToL(LowerTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
      }

#undef NEW_SIZE

      virtual inline ~BandMatrix() {}

      //
      // Op=
      //

      inline BandMatrix<T,S,I>& operator=(const BandMatrix<T,S,I>& m2)
      {
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignTo(View());
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenBandMatrix<RealType(T)>& m2)
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(View());
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenBandMatrix<ComplexType(T)>& m2)
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(View());
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
	Copy(m2,View());
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(T x) 
      {
	TMVAssert(colsize() == rowsize());
	return SetToIdentity(x); 
      }

      inline BandMatrix<T,S,I>& operator=(
	  const AssignableToBandMatrix<RealType(T)>& m2)
      {
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(View());
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const AssignableToBandMatrix<ComplexType(T)>& m2)
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() >= m2.nlo());
	TMVAssert(nhi() >= m2.nhi());
	m2.AssignToB(View());
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenDiagMatrix<RealType(T)>& m2)
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	m2.AssignToD(DiagMatrixViewOf(diag()));
	if (I==CStyle) {
	  if (nhi() > 0) Diags(1,nhi()+1).Zero();
	  if (nlo() > 0) Diags(-nlo(),0).Zero();
	} else {
	  if (nhi() > 0) Diags(1,nhi()).Zero();
	  if (nlo() > 0) Diags(-nlo(),-1).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenDiagMatrix<ComplexType(T)>& m2) 
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	m2.AssignToD(DiagMatrixViewOf(diag()));
	if (I==CStyle) {
	  if (nhi() > 0) Diags(1,nhi()+1).Zero();
	  if (nlo() > 0) Diags(-nlo(),0).Zero();
	} else {
	  if (nhi() > 0) Diags(1,nhi()).Zero();
	  if (nlo() > 0) Diags(-nlo(),-1).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenUpperTriMatrix<RealType(T)>& m2) 
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nhi() == rowsize()-1);
	m2.AssignToU(UpperTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
	if (I==CStyle) {
	  if (nlo() > 0) Diags(-nlo(),0).Zero();
	} else {
	  if (nlo() > 0) Diags(-nlo(),-1).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenUpperTriMatrix<ComplexType(T)>& m2)
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nhi() == rowsize()-1);
	m2.AssignToU(UpperTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
	if (I==CStyle) {
	  if (nlo() > 0) Diags(-nlo(),0).Zero();
	} else {
	  if (nlo() > 0) Diags(-nlo(),-1).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenLowerTriMatrix<RealType(T)>& m2) 
      { 
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() == rowsize()-1);
	m2.AssignToL(LowerTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
	if (I==CStyle) {
	  if (nhi() > 0) Diags(1,nhi()+1).Zero();
	} else {
	  if (nhi() > 0) Diags(1,nhi()).Zero();
	}
	return *this;
      }

      inline BandMatrix<T,S,I>& operator=(
	  const GenLowerTriMatrix<ComplexType(T)>& m2) 
      { 
	TMVAssert(IsComplex(T()));
	TMVAssert(colsize() == m2.colsize());
	TMVAssert(rowsize() == m2.rowsize());
	TMVAssert(nlo() == rowsize()-1);
	m2.AssignToL(LowerTriMatrixView<T>(ptr(),m2.size(),stepi(),stepj(),
	      NonUnitDiag,isdm()?NoMajor:S,NonConj FIRSTLAST));
	if (I==CStyle) {
	  if (nhi() > 0) Diags(1,nhi()+1).Zero();
	} else {
	  if (nhi() > 0) Diags(1,nhi()).Zero();
	}
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
	  return okij(i,j) ? cref(i,j) : T(0); 
	} else {
	  TMVAssert(i>0 && i <= colsize());
	  TMVAssert(j>0 && j <= rowsize());
	  return okij(i-1,j-1) ? cref(i-1,j-1) : T(0); 
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
	    std::min(colsize(),rowsize()),diagstep(),NonConj);
      }

      inline ConstVectorView<T,I> diag(int i) const
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = std::min(colsize(),rowsize()-i);
	  return ConstVectorView<T,I>(cptr()+i*stepj(),
	      diagsize,diagstep(),NonConj);
	} else {
	  const size_t diagsize = std::min(colsize()+i,rowsize());
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
	  TMVAssert(j2<=std::min(rowsize()-i,colsize()));
	  return ConstVectorView<T,I>(cptr()+i*stepj()+j1x*diagstep(),
	      j2-j1x,diagstep(),NonConj);
	} else {
	  TMVAssert(j2<=std::min(colsize()+i,rowsize()));
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
	    std::min(colsize(),rowsize()),diagstep(),NonConj FIRSTLAST);
      }

      inline VectorView<T,I> diag(int i)
      {
	TMVAssert(i<=nhi());
	TMVAssert(-i<=nlo());
	if (i >= 0) {
	  const size_t diagsize = std::min(colsize(),rowsize()-i);
	  return VectorView<T,I>(ptr()+i*stepj(),
	      diagsize,diagstep(),NonConj FIRSTLAST);
	} else {
	  const size_t diagsize = std::min(colsize()+i,rowsize());
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
	  TMVAssert(j2<=std::min(rowsize()-i,colsize()));
	  return VectorView<T,I>(ptr()+i*stepj()+j1x*diagstep(),
	      j2-j1x,diagstep(),NonConj FIRSTLAST);
	} else {
	  TMVAssert(j2<=std::min(colsize()+i,rowsize()));
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
	TMVAssert(colsize() == rowsize());
	TMVAssert(nlo()==nhi());
	View().TransposeSelf(); 
	return *this; 
      }

      inline BandMatrix<T,S,I>& ConjugateSelf() 
      { LinearView().ConjugateSelf(); return *this; }

      inline BandMatrix<T,S,I>& SetToIdentity(T x=T(1)) 
      { 
	TMVAssert(colsize() == rowsize());
	Zero(); diag().SetAllTo(x); 
	return *this; 
      }

      //
      // SubBandMatrix
      //

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, stepi(), stepj(), S, NonConj);
      }

      inline ConstMatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, NonConj);
      }

      inline ConstVectorView<T,I> SubVector(
	  int i, int j, int istep, int jstep, size_t size) const
      {
	TMVAssert(View().OKSubVector(i,j,istep,jstep,size));
	const int ix = (I==CStyle ? i : i-1);
	const int jx = (I==CStyle ? j : j-1);
	return ConstVectorView<T,I>(cptr()+ix*stepi()+jx*stepj(),size,
	    istep*stepi()+jstep*stepj(), NonConj);
      }

      inline ConstBandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
      {
	TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	return ConstBandMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), S, NonConj);
      }

      inline ConstBandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep) const
      {
	TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,
	      newnlo,newnhi,istep,jstep));
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
	const int j2 = std::min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1x < nlo() ? std::min(nlo(),i2-1) - i1x : 0;
	const int newnhi = j1==j2 ? 0 : std::min(nlo()+nhi()-newnlo,j2-j1-1);
	const size_t newlin = isrm() ? 1 : 0;
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
	const int i2 = std::min(j2 + nlo(),int(colsize()));
	const int newnhi = j1x < nhi() ? std::min(nhi(),j2-1) - j1x : 0;
	const int newnlo = i1==i2 ? 0 : std::min(nlo()+nhi()-newnhi,i2-i1-1);
	const size_t newlin = iscm() ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1x,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T,I>(cptr()+i1*stepi()+j1x*stepj(),
	    i2-i1, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin);
      }

      inline ConstBandMatrixView<T,I> Diags(int k1, int k2) const
      {
	if (I==FortranStyle) ++k2;
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 <= 0 ? -k2+1 : 0;
	const int i2 = std::min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = std::min(rowsize(),colsize()+k2-1);
	const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return ConstBandMatrixView<T,I>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct());
      }

      inline ConstBandMatrixView<T,I> UpperBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    std::min(colsize(),rowsize()),std::min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct());
      }

      inline ConstBandMatrixView<T,I> LowerBand() const
      {
	return ConstBandMatrixView<T>(cptr(),
	    std::min(colsize(),rowsize()+nlo()),std::min(colsize(),rowsize()),
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
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, stepi(), stepj(), S, NonConj FIRSTLAST);
      }

      inline MatrixView<T,I> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep)
      {
	TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	const StorageType newstor = S==RowMajor ?
	  jstep == 1 ? RowMajor : NoMajor :
	  S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	const int i2x = (I==CStyle ? i2 : i2-1+istep);
	const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    (i2x-i1x)/istep,(j2x-j1x)/jstep,istep*stepi(),jstep*stepj(),
	    newstor, NonConj FIRSTLAST);
      }

      inline VectorView<T,I> SubVector(int i, int j,
	  int istep, int jstep, size_t size)
      {
	TMVAssert(View().OKSubVector(i,j,istep,jstep,size));
	const int ix = (I==CStyle ? i : i-1);
	const int jx = (I==CStyle ? j : j-1);
	return VectorView<T,I>(ptr()+ix*stepi()+jx*stepj(),size,
	    istep*stepi()+jstep*stepj(), NonConj FIRSTLAST);
      }

      inline BandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi)
      {
	TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
	const int i1x = (I==CStyle ? i1 : i1-1);
	const int j1x = (I==CStyle ? j1 : j1-1);
	return BandMatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	    i2-i1x, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), S, NonConj FIRSTLAST);
      }

      inline BandMatrixView<T,I> SubBandMatrix(
	  int i1, int i2, int j1, int j2, int newnlo, int newnhi,
	  int istep, int jstep)
      {
	TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,
	      newnlo,newnhi,istep,jstep));
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
	const int j2 = std::min(i2 + nhi(),int(rowsize()));
	const int newnlo = i1x < nlo() ? std::min(nlo(),i2-1) - i1x : 0;
	const int newnhi = j1==j2 ? 0 : std::min(nlo()+nhi()-newnlo,j2-j1-1);
	const size_t newlin = isrm() ? 1 : 0;
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
	const int i2 = std::min(j2 + nlo(),int(colsize()));
	const int newnhi = j1x < nhi() ? std::min(nhi(),j2-1) - j1x : 0;
	const int newnlo = i1==i2 ? 0 : std::min(nlo()+nhi()-newnhi,i2-i1-1);
	const size_t newlin = iscm() ? 1 : 0;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1x,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1x*stepj(),
	    i2-i1, j2-j1x, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct(), newlin FIRSTLAST);
      }

      inline BandMatrixView<T,I> Diags(int k1, int k2)
      {
	if (I==FortranStyle) ++k2;
	TMVAssert(k1 >= -nlo());
	TMVAssert(k1 < k2);
	TMVAssert(k2 <= nhi()+1);

	const int i1 = k2 <= 0 ? -k2+1 : 0;
	const int i2 = std::min(rowsize()-k1,colsize());
	const int j1 = k1 <= 0 ? 0 : k1;
	const int j2 = std::min(rowsize(),colsize()+k2-1);
	const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
	const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
	TMVAssert(GenBandMatrix<T>::OKSubBandMatrix(
	      i1,i2,j1,j2,newnlo,newnhi,1,1));
	return BandMatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
	    diagstep(), stor(), ct()  FIRSTLAST);
      }

      inline BandMatrixView<T,I> UpperBand() 
      {
	return BandMatrixView<T,I>(ptr(),
	    std::min(colsize(),rowsize()),std::min(colsize()+nhi(),rowsize()),
	    0,nhi(),stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      }

      inline BandMatrixView<T,I> LowerBand()
      {
	return BandMatrixView<T,I>(ptr(),
	    std::min(colsize(),rowsize()+nlo()),std::min(colsize(),rowsize()),
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
      inline ConjType ct() const { return NonConj; }
      inline bool isrm() const { return S==RowMajor; }
      inline bool iscm() const { return S==ColMajor; }
      inline bool isdm() const { return S==DiagMajor; }
      inline bool isconj() const { return false; }
      inline size_t ls() const { return linsize; }

      inline bool CanLinearize() const
      { return true; }

    protected :

      const size_t linsize;
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
	TMVAssert(GenBandMatrix<T>::okij(i,j));
	return *(cptr() + int(i)*stepi() + int(j)*stepj());
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
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstMatrixView<T,I>& m, int nlo, int nhi)
    { 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const Matrix<T,S,I>& m, int nlo, int nhi)
    { 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const MatrixView<T,I>& m, int nlo, int nhi)
  {  
    return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	Matrix<T,S,I>& m, int nlo, int nhi)
    {  
      return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
	  FIRSTLAST1(m.first,m.last) ); 
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenBandMatrix<T>& m, int nlo, int nhi)
  { 
    TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
    return ConstBandMatrixView<T>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstBandMatrixView<T,I>& m, int nlo, int nhi)
    { 
      TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const BandMatrix<T,S,I>& m, int nlo, int nhi)
    { 
      TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
      return ConstBandMatrixView<T,I>(m.cptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const BandMatrixView<T,I>& m, int nlo, int nhi)
  { 
    TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
    return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	BandMatrix<T,S,I>& m, int nlo, int nhi)
    { 
      TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
      return BandMatrixView<T,I>(m.ptr(),m.colsize(),m.rowsize(),nlo,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
	  FIRSTLAST1(m.first,m.last) ); 
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenUpperTriMatrix<T>& m, int nhi=-1)
  { 
    if (nhi < 0) nhi = m.size()-1;
    TMVAssert(!m.isunit());
    return ConstBandMatrixView<T>(m.cptr(),m.size(),m.size(),0,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstUpperTriMatrixView<T,I>& m, int nhi=-1)
    {
      if (nhi < 0) nhi = m.size()-1;
      TMVAssert(!m.isunit());
      return ConstBandMatrixView<T,I>(m.cptr(),m.size(),m.size(),0,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
    }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const UpperTriMatrix<T,D,S,I>& m, int nhi=-1)
    {
      if (nhi < 0) nhi = m.size()-1;
      TMVAssert(D==NonUnitDiag);
      return ConstBandMatrixView<T,I>(m.cptr(),m.size(),m.size(),0,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const UpperTriMatrixView<T,I>& m, int nhi=-1)
  {
    if (nhi < 0) nhi = m.size()-1;
    TMVAssert(!m.isunit());
    return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),0,nhi,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct() 
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	UpperTriMatrix<T,D,S,I>& m, int nhi=-1)
    {
      if (nhi < 0) nhi = m.size()-1;
      TMVAssert(D==NonUnitDiag);
      return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),0,nhi,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct() 
	  FIRSTLAST1(m.first,m.last) ); 
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenLowerTriMatrix<T>& m, int nlo=-1)
  { 
    if (nlo < 0) nlo = m.size() - 1;
    TMVAssert(!m.isunit());
    return ConstBandMatrixView<T>(m.cptr(),m.size(),m.size(),nlo,0,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstLowerTriMatrixView<T,I>& m, int nlo=-1)
    {
      if (nlo < 0) nlo = m.size() - 1;
      TMVAssert(!m.isunit());
      return ConstBandMatrixView<T,I>(m.cptr(),m.size(),m.size(),nlo,0,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
    }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const LowerTriMatrix<T,D,S,I>& m, int nlo=-1)
    {
      if (nlo < 0) nlo = m.size() - 1;
      TMVAssert(D==NonUnitDiag);
      return ConstBandMatrixView<T,I>(m.cptr(),m.size(),m.size(),nlo,0,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const LowerTriMatrixView<T,I>& m, int nlo=-1)
  {
    if (nlo < 0) nlo = m.size() - 1;
    TMVAssert(!m.isunit());
    return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),nlo,0,
	m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
	FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
	LowerTriMatrix<T,D,S,I>& m, int nlo=-1)
    {
      if (nlo < 0) nlo = m.size() - 1;
      TMVAssert(D==NonUnitDiag);
      return BandMatrixView<T,I>(m.ptr(),m.size(),m.size(),nlo,0,
	  m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct() 
	  FIRSTLAST1(m.first,m.last) ); 
    }

  template <class T> inline ConstBandMatrixView<T> BandMatrixViewOf(
      const GenDiagMatrix<T>& m)
  { 
    return ConstBandMatrixView<T>(m.diag().cptr(),m.size(),m.size(),0,0,
	m.diag().step()-1,1,m.diag().step(),
	m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct());
  }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const ConstDiagMatrixView<T,I>& m)
    {
      return ConstBandMatrixView<T,I>(m.diag().cptr(),m.size(),m.size(),0,0,
	m.diag().step()-1,1,m.diag().step(),
	m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct());
    }

  template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
	const DiagMatrix<T,I>& m)
    {
      return ConstBandMatrixView<T,I>(m.diag().cptr(),m.size(),m.size(),0,0,
	m.diag().step()-1,1,m.diag().step(),
	m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct());
    }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      const DiagMatrixView<T,I>& m)
  {
    return BandMatrixView<T,I>(m.diag().ptr(),m.size(),m.size(),0,0,
	m.diag().step()-1,1,m.diag().step(),
	m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct()
	FIRSTLAST1(m.diag().first,m.diag().last) ); 
  }

  template <class T, IndexStyle I> inline BandMatrixView<T,I> BandMatrixViewOf(
      DiagMatrix<T,I>& m)
  {
    return BandMatrixView<T,I>(m.diag().ptr(),m.size(),m.size(),0,0,
	m.diag().step()-1,1,m.diag().step(),
	m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct()
	FIRSTLAST1(m.diag().first,m.diag().last) ); 
  }

  template <class T> ConstBandMatrixView<T> BandMatrixViewOf(
      const T* m, size_t cs, size_t rs, int nlo, int nhi, StorageType stor);

  template <class T> BandMatrixView<T> BandMatrixViewOf(
      T* m, size_t cs, size_t rs, int nlo, int nhi, StorageType stor);

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
      const GenBandMatrix<std::complex<T> >& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T1, class T2> inline void DoCopy1(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T2()));
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.nlo() == m1.nlo());
    TMVAssert(m2.nhi() == m1.nhi());
    if (m2.colsize() > 0 && m2.rowsize() > 0) {
      if (SameStorage(m1,m2)) {
	if (m2.SameAs(m1)) {} // Do Nothing
	else if (m2.nlo() == m2.nhi() &&  m2.Transpose().SameAs(m1)) 
	  m2.TransposeSelf();
	else if (m1.isconj() != m2.isconj() && m2.Conjugate().SameAs(m1)) 
	  m2.ConjugateSelf();
	else if (m1.isrm()) 
	  DoCopy1(BandMatrix<T1,RowMajor>(m1),m2);
	else if (m1.iscm()) 
	  DoCopy1(BandMatrix<T1,ColMajor>(m1),m2);
	else 
	  DoCopy1(BandMatrix<T1,DiagMajor>(m1),m2);
      }
      else if (m2.isconj()) DoCopy(m1.Conjugate(),m2.Conjugate());
      else DoCopy(m1,m2);
    }
  }

  template <class T1, class T2> inline void Copy(
      const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T2()));
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.nlo() >= m1.nlo());
    TMVAssert(m2.nhi() >= m1.nhi());

    DoCopy1(m1,m2.SubBandMatrix(0,m2.colsize(),0,m2.rowsize(),
	  m1.nlo(),m1.nhi()));
    if (m2.nhi() > m1.nhi())
      m2.Diags(m1.nhi()+1,m2.nhi()+1).Zero();
    if (m2.nlo() > m1.nlo())
      m2.Diags(-m2.nlo(),-m1.nlo()).Zero();
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

  template <class T, StorageType S, IndexStyle I> std::istream& operator>>(
      std::istream& fin, auto_ptr<BandMatrix<T,S,I> >& m);

  template <class T> std::istream& operator>>(std::istream& fin,
      const BandMatrixView<T>& m);

  template <class T, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(std::istream& fin, BandMatrix<T,S,I>& m)
    { return fin >> m.View(); }

} // namespace tmv

#endif
