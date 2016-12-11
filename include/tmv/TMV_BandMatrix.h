///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// This file defines the TMV BandMatrix class.
//
// A BandMatrix is only non-zero in a (typically small) number of
// diagonals around the main diagonal.
// Specifically, we store nhi super-diagonals (above the main diagonal)
// and nlo sub-diagonals (below the main).
//
// As with the regular Matrix class, there are two template arguments.
// The first is simply the type of the data.
// The second, which is optional, specifies the known attributes of
// the matrix.  The valid attributes are:
// - ColMajor or RowMajor or DiagMajor
// - CStyle or FortranStyle
// - WithDivider or NoDivider
//
// The attributes are treated as a bit field, so you | them together to
// get the complete value, just like with a regular Matrix.
// The default values are ColMajor, CStyle and WithDivider, so you
// only need to specify changes to that.
//
// The RowMajor and ColMajor storage options are the same as for a normal
// Matrix.  But there is a new possibility for BandMatrices -- DiagMajor.
// In this case the storage is in order of each diagonal starting
// with the lowest one, proceding along that whole diagonal, and then
// the next diagonal up, and so on.
//
// Also, for each storage possibility, we store some extra elements
// in order to have the rows, columns, and diagonals all have constant
// steps.  For example, a 6x6 ColMajor BandMatrix with nlo=2,nhi=3 has
// the following elements, numbered in order that they are stored.
//
// [ 1   6  11  16         ]
// [ 2   7  12  17  22     ]
// [ 3   8  13  18  23  28 ]
// [     9  14  19  24  29 ]
// [        15  20  25  30 ]
// [            21  26  31 ]
//
// We do not ever use elements stored at the memory locations 4-5, 10, 27.
// We need to skip those in order to get the rows and diagonals to have
// constant steps in the memory.  If we had started the second column
// with memory location 4, and the third at 8, then the first row would
// be 1 4 8 13 which cannot be referenced by a VectorView.
//
// Likewise, the same BandMatrix in RowMajor and DiagMajor storage
// is as follows:
//
// [  1   2   3   4         ]  [ 11  17  23  29         ]
// [  6   7   8   9  10     ]  [  6  12  18  24  30     ]
// [ 11  12  13  14  15  16 ]  [  1   7  13  19  25  31 ]
// [     17  18  19  20  21 ]  [      2   8  14  20  26 ]
// [         23  24  25  26 ]  [          3   9  15  21 ]
// [             29  30  31 ]  [              4  10  16 ]
//
// For square BandMatrices, the wasted storage is only
// (nlo-1)*nlo/2 + (nhi-1)*nhi/2 memory locations,
// which, if nlo and nhi are small compared to the size N,
// is negligible compared to the total memory allocated of
// (N-1)*(nlo+nhi+1)+1.
//
// Also, we don't actually require that the BandMatrix be square,
// although that is the usual case.  Hopefully, the extension
// of the above formats to non-square cases is obvious.
// The memory required for non-square BandMatrices is not quite as
// simple a formula, and it depends on the StorageType.
// The required memory (in units of sizeof(T)) can be obtained
// by the function:
//
// int BandStorageLength(stor, colsize, rowsize, nlo, nhi);
//
// Constructors:
//
//    BandMatrix<T,A>(colsize, rowsize, nlo, nhi)
//        Makes a BandMatrix with column size = row size = size
//        with nhi non-zero superdiagonals and nlo non-zero subdiagonals
//        with _uninitialized_ values
//
//    BandMatrix<T,A>(colsize, rowsize, nlo, nhi, x)
//        The same as above, but all values are initialized to x
//
//    BandMatrix<T,A>(const Matrix<T>& m, nlo, nhi)
//        Makes a BandMatrix which copies the corresponding elements of m.
//
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
//    ConstBandMatrixView BandMatrixViewOf(const DiagMatrix<T>& m)
//    ConstBandMatrixView BandMatrixViewOf(const TriMatrix<T>& m, nlo=n-1)
//        Makes a constant BandMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//        For the second version, nlo,nhi must be <= the corresponding
//        values in m.
//        For the fourth version, nlo is optional.
//
//    BandMatrixView BandMatrixViewOf(Matrix<T>& m, nlo, nhi)
//    BandMatrixView BandMatrixViewOf(BandMatrix<T>& m, nlo, nhi)
//    BandMatrixView BandMatrixViewOf(DiagMatrix<T>& m)
//    BandMatrixView BandMatrixViewOf(TriMatrix<T>& m, nhi=n-1)
//        Makes a modifiable BandMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the orginial Matrix.
//
//    ConstBandMatrixView BandMatrixViewOf(const T* m,
//            colsize, rowsize, nlo, nhi, stor)
//    BandMatrixView BandMatrixViewOf(T* m,
//            colsize, rowsize, nlo, nhi, stor)
//    ConstBandMatrixView BandMatrixViewOf(const T* m,
//            colsize, rowsize, nlo, nhi, stepi, stepj)
//    BandMatrixView BandMatrixViewOf(T* m,
//            colsize, rowsize, nlo, nhi, stepi, stepj)
//        Makes a BandMatrixView of the elements in m using the actual
//        element m for the storage.  This is essentially the same as the
//        constructor with (const T* m), except that the data isn't duplicated.
//
// Access Functions
//
//    int colsize() const
//    int rowsize() const
//        Return the dimensions of the BandMatrix
//
//    int nlo() const
//    int nhi() const
//        Return the band dimensions
//
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the BandMatrix
//
//    VectorView row(int i, int j1, int j2)
//        Return a subset of the ith row of the BandMatrix
//        The range (i,j1)..(i,j2-1) must be entirely within the band.
//
//    VectorView col(int j, int i1, int i2)
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
//    BandMatrix& setZero()
//    BandMatrix& setAllTo(T x)
//    BandMatrix& addToAll(T x)
//    BandMatrix<T>& transposeSelf()
//        Must be square, and have nhi=nlo for this function
//    BandMatrix& conjugateSelf()
//    BandMatrix& setToIdentity(T x = 1)
//    void Swap(BandMatrix& m1, BandMatrix& m2)
//        Must be the same size and have the same band structure (nlo,nhi)
//
//
// Views of a BandMatrix:
//
//    subBandMatrix(int i1, int i2, int j1, int j2, int nlo, int nhi,
//            int istep=1, int jstep=1)
//        Returns a BandMatrixView with (i1,j1) in the upper left corner,
//        (i2,j2) in the lower right corder, and with nlo and nhi
//        sub- and super-diagonals respectively.
//        All members of the new subBandMatrix must be within the
//        original band.
//
//    subBandMatrix(int i1, int i2, int j1, int j2)
//        This is normally equivalenet to
//        b.subBandMatrix(i1,i2,j1,j2,b.nlo(),b.nhi())
//        However, it lowers the new nlo, nhi as necessary for small
//        i2-i1 or j2-j1 if the full band-width doesn't fit into the
//        new matrix size.
//
//    subMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//    subVector(int i1, int i2, int istep, int jstep, int size)
//        Just like the regular Matrix version, but all entries in the
//        submatrix must be within the band.
//
//    rowRange(int i1, int i2)
//        A bit different from the Matrix version, since the rowsize will
//        not be the full rowsize of the source BandMatrix.
//        Instead, this returns (as a BandMatrix) the nontrivial portion
//        of the BandMatrix in rows i1..i2 (not including i2)
//    colRange(int j1, int j2)
//        As with rowRange, this returns the nontrivial portion of the
//        BandMatrix in columns j1..j2 (not including j2)
//    diagRange(int i1, int i2)
//        Returns a thinner BandMatrix including the diagonals from i1..i2
//        (not including i2)
//
//    upperBand()
//    lowerBand()
//        Returns a BandMatrix of only the upper or lower portion of the
//        matrix.
//
//    upperBandOff()
//    lowerBandOff()
//        Returns a BandMatrix of only the strictly upper or lower portion
//        of the matrix.  (ie. not including the diagonal.)
//        Conceptually equivalent to upperBand().offDiag().
//        Except that a BandMatrix doesn't have an offDiag() function.
//
//    m.view()
//    m.transpose() or Transpose(m)
//    m.conjugate() or Conjugate(m)
//    m.adjoint() or Adjoint(m)
//
//
// Functions of BandMatrices:
//        (These are all both member functions and functions of a BandMatrix,
//         so Norm(m) and m.norm() for example are equivalent.)
//
//    m.det() or Det(m)
//    m.logDet() or m.logDet(T* sign) or LogDet(m)
//    m.trace() or Trace(m)
//    m.norm() or m.normF() or Norm(m) or NormF(m)
//    m.sumElements() or SumElements(m)
//    m.sumAbsElements() or SumAbsElements(m)
//    m.sumAbs2Elements() or SumAbs2Elements(m)
//    m.normSq() or NormSq(m)
//    m.norm1() or Norm1(m)
//    m.norm2() or Norm2(m)
//    m.normInf() or NormInf(m)
//    m.maxAbsElement() or MaxAbsElements(m)
//    m.maxAbs2Element() or MaxAbs2Elements(m)
//
//    m.inverse() or Inverse(m)
//    m.makeInverse(minv)
//    m.makeInverseATA(invata)
//
//
// I/O:
//
//    os << m
//        Writes m to ostream os in the usual Matrix format
//
//    os << CompactIO() << m
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
//    is >> CompactIO() >> m
//        Reads m from istream is in either format
//
//
// Division Control Functions:
//
//    m.divideUsing(dt)
//    where dt is LU, QR, or SV
//
//    m.lud(), m.qrd(), m.svd() return the corresponding Divider classes.
//
//


#ifndef TMV_BandMatrix_H
#define TMV_BandMatrix_H

#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Array.h"

#ifdef TMV_USE_VALGRIND
#include <vector>
#endif

namespace tmv {

    template <typename T>
    class GenBandMatrix :
        virtual public AssignableToBandMatrix<T>,
        virtual public AssignableToDiagMatrix<T>,
        virtual public AssignableToUpperTriMatrix<T>,
        virtual public AssignableToLowerTriMatrix<T>,
        public BaseMatrix<T>,
        public DivHelper<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenBandMatrix<T> type;
        typedef BandMatrix<T> copy_type;
        typedef ConstBandMatrixView<T> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef ConstMatrixView<T> const_rec_type;
        typedef ConstBandMatrixView<RT> const_realpart_type;
        typedef BandMatrixView<T> nonconst_type;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef CDMIt<type> const_diagmajor_iterator;

        //
        // Constructors
        //

        inline GenBandMatrix() {}
        inline GenBandMatrix(const type&) {}
        virtual inline ~GenBandMatrix() {}

        //
        // Access Functions
        //

        using AssignableToMatrix<T>::colsize;
        using AssignableToMatrix<T>::rowsize;
        using AssignableToBandMatrix<T>::nlo;
        using AssignableToBandMatrix<T>::nhi;
        // For Tri, Diag compatibility:
        inline ptrdiff_t size() const
        { TMVAssert(colsize() == rowsize()); return colsize(); }
        inline DiagType dt() const { return NonUnitDiag; }

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j>=0 && j<rowsize());
            return cref(i,j);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return const_vec_type(
                cptr()+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),j2-j1,stepj(),ct());
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return const_vec_type(
                cptr()+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),i2-i1,stepi(),ct());
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                cptr(),TMV_MIN(colsize(),rowsize()),diagstep(),ct());
        }

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    cptr()+i*ptrdiff_t(stepj()),diagsize,diagstep(),ct());
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    cptr()-i*ptrdiff_t(stepi()),diagsize,diagstep(),ct());
            }
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            TMVAssert(j1>=0 && j1-j2<=0);
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(rowsize()-i,colsize()));
                return const_vec_type(
                    cptr()+i*ptrdiff_t(stepj())+j1*diagstep(),j2-j1,diagstep(),ct());
            } else {
                TMVAssert(j2<=TMV_MIN(colsize()+i,rowsize()));
                return const_vec_type(
                    cptr()-i*ptrdiff_t(stepi())+j1*diagstep(),j2-j1,diagstep(),ct());
            }
        }

        template <typename T2>
        inline bool isSameAs(const GenBandMatrix<T2>& ) const
        { return false; }

        inline bool isSameAs(const type& m2) const
        {
            return ( this==&m2 ||
                     ( cptr()==m2.cptr() &&
                       colsize()==m2.colsize() && rowsize()==m2.rowsize() &&
                       stepi()==m2.stepi() && stepj()==m2.stepj() &&
                       nhi()==m2.nhi() && nlo()==m2.nlo() &&
                       isconj() == m2.isconj() ) );
        }

        inline void assignToM(MatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            assignToB(BandMatrixView<RT>(m2,nlo(),nhi()));
            if (rowsize() > nhi()+1)
                BandMatrixView<RT>(
                    m2.colRange(nhi()+1,rowsize()),0,rowsize()-nhi()-2).setZero();
            if (colsize() > nlo()+1)
                BandMatrixView<RT>(
                    m2.rowRange(nlo()+1,colsize()),colsize()-nlo()-2,0).setZero();
        }

        inline void assignToM(MatrixView<CT> m2) const
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            assignToB(BandMatrixView<CT>(m2,nlo(),nhi()));
            if (rowsize() > nhi()+1)
                BandMatrixView<CT>(
                    m2.colRange(nhi()+1,rowsize()),0,rowsize()-nhi()-2).setZero();
            if (colsize() > nlo()+1)
                BandMatrixView<CT>(
                    m2.rowRange(nlo()+1,colsize()),colsize()-nlo()-2,0).setZero();
        }

        inline void assignToB(BandMatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(m2.nhi() >= nhi());
            if (!isSameAs(m2)) Copy(*this,m2);
        }

        inline void assignToB(BandMatrixView<CT> m2) const
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(m2.nhi() >= nhi());
            if (!isSameAs(m2)) Copy(*this,m2);
        }

        inline void assignToU(UpperTriMatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nlo() == 0);
            assignToB(BandMatrixViewOf(m2));
        }

        inline void assignToU(UpperTriMatrixView<CT> m2) const
        {
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nlo() == 0);
            assignToB(BandMatrixViewOf(m2));
        }

        inline void assignToL(LowerTriMatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0);
            assignToB(BandMatrixViewOf(m2));
        }

        inline void assignToL(LowerTriMatrixView<CT> m2) const
        {
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0);
            assignToB(BandMatrixViewOf(m2));
        }

        inline void assignToD(DiagMatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0 && nlo() == 0);
            assignToB(BandMatrixViewOf(m2));
        }

        inline void assignToD(DiagMatrixView<CT> m2) const
        {
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0 && nlo() == 0);
            assignToB(BandMatrixViewOf(m2));
        }

        //
        // subBandMatrix
        //

        bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return const_rec_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, stepi(), stepj(), ct());
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return const_rec_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                ct());
        }

        bool hasSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const;

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return const_vec_type(
                cptr()+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),size,
                istep*stepi()+jstep*stepj(),ct());
        }

        bool hasSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(),ct());
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            const ptrdiff_t newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const ptrdiff_t newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            const ptrdiff_t newstepi = stepi()*istep;
            const ptrdiff_t newstepj = stepj()*jstep;
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi,
                newstepi, newstepj, newstepi+newstepj, ct());
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());

            const ptrdiff_t j1 = i1 > nlo() ? i1-nlo() : 0;
            const ptrdiff_t j2 = TMV_MIN(i2 + nhi(),rowsize());
            const ptrdiff_t newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const ptrdiff_t newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const ptrdiff_t newlin = (ls() && isrm()) ? -1 : 0;
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());

            const ptrdiff_t i1 = j1 > nhi() ? j1-nhi() : 0;
            const ptrdiff_t i2 = TMV_MIN(j2 + nlo(),colsize());
            const ptrdiff_t newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const ptrdiff_t newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            const ptrdiff_t newlin = (ls() && iscm()) ? -1 : 0;
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin);
        }

        inline const_view_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);

            const ptrdiff_t i1 = k2 <= 0 ? -k2+1 : 0;
            const ptrdiff_t i2 = TMV_MIN(rowsize()-k1,colsize());
            const ptrdiff_t j1 = k1 <= 0 ? 0 : k1;
            const ptrdiff_t j2 = TMV_MIN(rowsize(),colsize()+k2-1);
            const ptrdiff_t newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const ptrdiff_t newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct());
        }

        inline const_view_type upperBand() const
        {
            return const_view_type(
                cptr(),TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),ct());
        }

        inline const_view_type lowerBand() const
        {
            return const_view_type(
                cptr(),TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),ct());
        }

        inline const_view_type upperBandOff() const
        {
            return const_view_type(
                cptr()+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),ct());
        }

        inline const_view_type lowerBandOff() const
        {
            return const_view_type(
                cptr()+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),ct());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NonConj);
        }

        //
        // Views
        //

        inline const_view_type view() const
        {
            return const_view_type(
                cptr(),colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),ct(),
                isdm()?0:ls());
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                cptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),ct(),
                isdm()?0:ls());
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                cptr(),colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),TMV_ConjOf(T,ct()),
                isdm()?0:ls());
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                cptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),
                TMV_ConjOf(T,ct()),isdm()?0:ls());
        }

        // This needs to be virtual, since DiagMajor BandMatrices need to
        // start at itsm1, not cptr() (=itsm).  The BandMatrix version
        // of this gets that right.
        // Also, because of the virtual, we now need two names (hence the
        // Const), since the BandMatrixView::linearView() function has a
        // different return type, so it is not "covariant" with this
        // funtion, which it needs to be for virtual overriding.
        virtual inline const_vec_type constLinearView() const
        {
            TMVAssert(iscm() || isrm());
            TMVAssert(ls() != -1);
            // (To assure that next assert has no effect.)
            TMVAssert(canLinearize());
#ifdef TMV_USE_VALGRIND
            std::vector<bool> ok(ls(),false);
            for (ptrdiff_t i=0;i<colsize();++i) {
                for (ptrdiff_t j=0;j<rowsize();++j) {
                    if (okij(i,j)) {
                        ptrdiff_t k = i*stepi() + j*stepj();
                        ok[k] = true;
                    }
                }
            }
            for (ptrdiff_t k=0;k<ls();++k) if (!ok[k]) {
                const_cast<T*>(cptr())[k] = T(-777);
            }
#endif
            return const_vec_type(cptr(),ls(),1,ct());
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<T*>(cptr()),colsize(),rowsize(),
                nlo(),nhi(),stepi(),stepj(),diagstep(),ct(),
                isdm()?0:ls()
                TMV_FIRSTLAST1(isdm() ? diag(-nlo()).begin().getP() : cptr(),
                               isdm() ? diag(nhi()).end().getP() :
                               ((diag().end()-1).getP()+1)));
        }

        //
        // Functions of Matrix
        //

        T det() const;

        RT logDet(T* sign=0) const;

        bool isSingular() const;

        inline T trace() const
        { return diag().sumElements(); }

        T sumElements() const;

        RT sumAbsElements() const;

        RT sumAbs2Elements() const;

        inline RT norm() const
        { return normF(); }

        RT normF() const;

        RT normSq(const RT scale = RT(1)) const;

        RT norm1() const;


        inline RT normInf() const
        { return transpose().norm1(); }

        RT maxAbsElement() const;
        RT maxAbs2Element() const;

        RT doNorm2() const;
        RT norm2() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::norm2();
            else return doNorm2();
        }

        RT doCondition() const;
        inline RT condition() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::condition();
            else return doCondition();
        }


        QuotXB<T,T> QInverse() const;
        inline QuotXB<T,T> inverse() const
        { return QInverse(); }

        //
        // Division Control
        //

        void setDiv() const;

        inline void divideUsing(DivType dt) const
        {
            TMVAssert(dt == LU || dt == QR || dt == SV);
            DivHelper<T>::divideUsing(dt);
        }

        inline const BandLUDiv<T>& lud() const
        {
            divideUsing(LU);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsLUDiv());
            return static_cast<const BandLUDiv<T>&>(*this->getDiv());
        }

        inline const BandQRDiv<T>& qrd() const
        {
            divideUsing(QR);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsQRDiv());
            return static_cast<const BandQRDiv<T>&>(*this->getDiv());
        }

        inline const BandSVDiv<T>& svd() const
        {
            divideUsing(SV);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsSVDiv());
            return static_cast<const BandSVDiv<T>&>(*this->getDiv());
        }


        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        virtual const T* cptr() const = 0;
        virtual ptrdiff_t stepi() const = 0;
        virtual ptrdiff_t stepj() const = 0;
        virtual ptrdiff_t diagstep() const = 0;
        virtual ptrdiff_t ls() const  = 0;
        virtual inline bool isrm() const { return stepj() == 1; }
        virtual inline bool iscm() const { return stepi() == 1; }
        virtual inline bool isdm() const { return diagstep() == 1; }
        inline bool isconj() const
        {
            TMVAssert(isComplex(T()) || ct()==NonConj);
            return isComplex(T()) && ct()==Conj;
        }
        virtual ConjType ct() const = 0;

        virtual bool canLinearize() const = 0;
        virtual T cref(ptrdiff_t i, ptrdiff_t j) const;

        inline ptrdiff_t rowstart(ptrdiff_t i) const
        { return TMV_MAX(ptrdiff_t(0),i-nlo()); }
        inline ptrdiff_t rowend(ptrdiff_t i) const
        { return TMV_MIN(rowsize(),i+nhi()+1); }
        inline ptrdiff_t colstart(ptrdiff_t j) const
        { return TMV_MAX(ptrdiff_t(0),j-nhi()); }
        inline ptrdiff_t colend(ptrdiff_t j) const
        { return TMV_MIN(colsize(),j+nlo()+1); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        {
            return const_rowmajor_iterator(
                this,TMV_MIN(colsize(),rowsize()+nlo()),
                rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        {
            return const_colmajor_iterator(
                this,colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }

        inline const_diagmajor_iterator diagmajor_begin() const
        { return const_diagmajor_iterator(this,nlo(),0); }
        inline const_diagmajor_iterator diagmajor_end() const
        { return const_diagmajor_iterator(this,0,nhi()+1); }


    protected :

        inline bool okij(ptrdiff_t i, ptrdiff_t j) const
        { return (j+nlo() >= i && i+nhi() >= j); }

        inline const BaseMatrix<T>& getMatrix() const { return *this; }

    private :

        type& operator=(const type&);

        bool divIsLUDiv() const;
        bool divIsQRDiv() const;
        bool divIsSVDiv() const;

    }; // GenBandMatrix

    template <typename T, int A>
    class ConstBandMatrixView : public GenBandMatrix<T>
    {
    public :

        typedef GenBandMatrix<T> base;
        typedef ConstBandMatrixView<T,A> type;

        inline ConstBandMatrixView(const type& rhs) :
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
            itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
            itsct(rhs.itsct), linsize(rhs.linsize)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(!(isdm() && linsize != 0));
        }

        inline ConstBandMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(rhs.nlo()), itsnhi(rhs.nhi()),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
            itsct(rhs.ct()),
            linsize(rhs.isdm()?0:rhs.ls())
            {
                TMVAssert(Attrib<A>::viewok);
                TMVAssert(!(isdm() && linsize != 0));
            }

        inline ConstBandMatrixView(
            const T* _m, ptrdiff_t _cs, ptrdiff_t _rs,
            ptrdiff_t _lo, ptrdiff_t _hi, ptrdiff_t _si, ptrdiff_t _sj, ptrdiff_t _sd,
            ConjType _ct, ptrdiff_t _ls=0) :
            itsm(_m), itscs(_cs), itsrs(_rs),
            itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
            itsct(_ct), linsize(_ls)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
            TMVAssert(linsize==0 || linsize==-1 ||
                      ( (itssi == 1 || itssj == 1) &&
                        linsize==BandStorageLength(
                            (itssi == 1 ? ColMajor : RowMajor),
                            itscs,itsrs,itsnlo,itsnhi)));
        }

        // These two work slightly differently than the BandMatrixViewOf
        // commands when the rhs matrix is not square.
        // These two constructors copy the size of rhs viewing only
        // the relevant rows.
        // In contrast, BandMatrixViewOf shrinks colsize or rowsize down to
        // only the rows and columns which include the bands.
        //   e.g. if rhs is 10 x 8, then:
        //   BandMatrixView(rhs,0,2) will have cs = rs = 8
        //   BandMatrixViewOf(rhs,0,2) will have cs = 10, rs = 8
        inline ConstBandMatrixView(const base& rhs, ptrdiff_t lo, ptrdiff_t hi) :
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
            itsct(rhs.ct()), linsize(0)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        inline ConstBandMatrixView(const GenMatrix<T>& rhs, ptrdiff_t lo, ptrdiff_t hi) :
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(itssi+itssj),
            itsct(rhs.ct()), linsize(0)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        virtual inline ~ConstBandMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        inline ptrdiff_t colsize() const { return itscs; }
        inline ptrdiff_t rowsize() const { return itsrs; }
        inline ptrdiff_t nlo() const { return itsnlo; }
        inline ptrdiff_t nhi() const { return itsnhi; }
        inline const T* cptr() const { return itsm; }
        inline ptrdiff_t stepi() const { return itssi; }
        inline ptrdiff_t stepj() const { return itssj; }
        inline ptrdiff_t diagstep() const { return itssd; }
        inline ptrdiff_t ls() const { return linsize; }
        inline ConjType ct() const { return itsct; }
        using base::isdm;

        bool canLinearize() const;

    protected :

        const T*const itsm;
        const ptrdiff_t itscs;
        const ptrdiff_t itsrs;
        const ptrdiff_t itsnlo;
        const ptrdiff_t itsnhi;
        const ptrdiff_t itssi;
        const ptrdiff_t itssj;
        const ptrdiff_t itssd;

        const ConjType itsct;
        mutable ptrdiff_t linsize;

    private :

        type& operator=(const type&);

    }; // ConstBandMatrixView

    template <typename T>
    class ConstBandMatrixView<T,FortranStyle> :
        public ConstBandMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef GenBandMatrix<T> base;
        typedef ConstBandMatrixView<T,FortranStyle> type;
        typedef ConstBandMatrixView<T,CStyle> c_type;
        typedef ConstBandMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstBandMatrixView<RT,FortranStyle> const_realpart_type;
        typedef BandMatrixView<T,FortranStyle> nonconst_type;

        inline ConstBandMatrixView(const type& rhs) : c_type(rhs) {}

        inline ConstBandMatrixView(const base& rhs) : c_type(rhs) {}

        inline ConstBandMatrixView(
            const T* _m, ptrdiff_t _cs, ptrdiff_t _rs,
            ptrdiff_t _lo, ptrdiff_t _hi, ptrdiff_t _si, ptrdiff_t _sj, ptrdiff_t _sd,
            ConjType _ct, ptrdiff_t ls=0) :
            c_type(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,_ct,ls) {}

        inline ConstBandMatrixView(const base& rhs, ptrdiff_t lo, ptrdiff_t hi) :
            c_type(rhs,lo,hi) {}

        inline ConstBandMatrixView(const GenMatrix<T>& rhs, ptrdiff_t lo, ptrdiff_t hi) :
            c_type(rhs,lo,hi) {}

        virtual inline ~ConstBandMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            TMVAssert(i>0 && i<=c_type::colsize());
            TMVAssert(j>0 && j<=c_type::rowsize());
            return c_type::cref(i-1,j-1);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>0 && i<=c_type::colsize());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=c_type::rowsize());
            TMVAssert(c_type::okij(i-1,j1-1));
            TMVAssert(c_type::okij(i-1,j2-1));
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>0 && j<=c_type::rowsize());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=c_type::colsize());
            TMVAssert(c_type::okij(i1-1,j-1));
            TMVAssert(c_type::okij(i2-1,j-1));
            return base::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag() const
        { return base::diag(); }

        inline const_vec_type diag(ptrdiff_t i) const
        { return base::diag(i); }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>0);
            return base::diag(i,j1-1,j2);
        }


        //
        // subBandMatrix
        //

        bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::subMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::subMatrix(i1-1,i2-1+istep,j1-1,j2-1+jstep,
                                   istep,jstep);
        }

        bool hasSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const;

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return base::subVector(i-1,j-1,istep,jstep,size);
        }

        bool hasSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2,
            ptrdiff_t newnlo, ptrdiff_t newnhi, ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return base::subBandMatrix(i1-1,i2,j1-1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            const ptrdiff_t newnlo = TMV_MIN(c_type::nlo()+j1-i1,i2-i1);
            const ptrdiff_t newnhi = TMV_MIN(c_type::nhi()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            return base::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=c_type::colsize());
            return base::rowRange(i1-1,i2);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>0 && j1-j2<=0 && j2<=c_type::rowsize());
            return base::colRange(j1-1,j2);
        }

        inline const_view_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        {
            TMVAssert(k1>=-c_type::nlo() && k1<=k2 && k2<=c_type::nhi());
            return base::diagRange(k1,k2+1);
        }

        inline const_view_type upperBand() const
        { return base::upperBand(); }

        inline const_view_type lowerBand() const
        { return base::lowerBand(); }

        inline const_view_type upperBandOff() const
        { return base::upperBandOff(); }

        inline const_view_type lowerBandOff() const
        { return base::lowerBandOff(); }

        inline const_realpart_type realPart() const
        { return base::realPart(); }

        inline const_realpart_type imagPart() const
        { return base::imagPart(); }

        //
        // Views
        //

        inline const_view_type view() const
        { return base::view(); }

        inline const_view_type transpose() const
        { return base::transpose(); }

        inline const_view_type conjugate() const
        { return base::conjugate(); }

        inline const_view_type adjoint() const
        { return base::adjoint(); }

        inline nonconst_type nonConst() const
        { return base::nonConst(); }

    private :

        type& operator=(const type&);

    }; // FortranStyle ConstBandMatrixView

    template <typename T, int A>
    class BandMatrixView : public GenBandMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenBandMatrix<T> base;
        typedef BandMatrixView<T,A> type;
        typedef BandMatrixView<T,A> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,A> vec_type;
        typedef MatrixView<T,A> rec_type;
        typedef BandMatrixView<RT,A> realpart_type;
        typedef ConstBandMatrixView<T,A> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,A> const_vec_type;
        typedef ConstMatrixView<T,A> const_rec_type;
        typedef ConstBandMatrixView<RT,A> const_realpart_type;
        typedef typename RefHelper<T>::reference reference;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef DMIt<type> diagmajor_iterator;
        typedef RMIt<const type> const_rowmajor_iterator;
        typedef CMIt<const type> const_colmajor_iterator;
        typedef DMIt<const type> const_diagmajor_iterator;

        //
        // Constructors
        //

        inline BandMatrixView(const type& rhs) :
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
            itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
            itsct(rhs.ct()), linsize(rhs.ls())
            TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        inline BandMatrixView(
            T* _m, ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _lo, ptrdiff_t _hi,
            ptrdiff_t _si, ptrdiff_t _sj, ptrdiff_t _sd, ConjType _ct,
            ptrdiff_t _ls=0 TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itscs(_cs), itsrs(_rs),
            itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
            itsct(_ct), linsize(_ls)
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
            TMVAssert(linsize==0 || linsize==-1 ||
                      ( (itssi == 1 || itssj == 1) &&
                        linsize==BandStorageLength(
                            (itssi == 1 ? ColMajor : RowMajor),
                            itscs,itsrs,itsnlo,itsnhi)));
        }

        inline BandMatrixView(BandMatrixView<T> rhs, ptrdiff_t lo, ptrdiff_t hi) :
            itsm(rhs.ptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
            itsct(rhs.ct()), linsize(0)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        inline BandMatrixView(MatrixView<T> rhs, ptrdiff_t lo, ptrdiff_t hi) :
            itsm(rhs.ptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(itssi+itssj),
            itsct(rhs.ct()), linsize(0)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(rowsize() == 0 || nhi() < rowsize());
            TMVAssert(colsize() == 0 || nlo() < colsize());
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        virtual inline ~BandMatrixView()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            const_cast<T*&>(itsm) = 0;
#endif
        }

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        inline type& operator=(const GenBandMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        inline type& operator=(const GenBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const GenBandMatrix<T2>& m2)
        {
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            if (!this->isSameAs(m2)) Copy(m2,view());
            return *this;
        }

        inline type& operator=(const T& x)
        {
            TMVAssert(colsize() == rowsize());
            return setToIdentity(x);
        }

        inline type& operator=(const AssignableToBandMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        inline type& operator=(const AssignableToBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignToD(DiagMatrixViewOf(this->diag()));
            if (this->nhi() > 0) upperBandOff().setZero();
            if (this->nlo() > 0) lowerBandOff().setZero();
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignToD(DiagMatrixViewOf(this->diag()));
            if (this->nhi() > 0) upperBandOff().setZero();
            if (this->nlo() > 0) lowerBandOff().setZero();
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == rowsize()-1);
            m2.assignToU(
                UpperTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
            if (nlo() > 0) diagRange(-nlo(),0).setZero();
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == rowsize()-1);
            m2.assignToU(
                UpperTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
            if (nlo() > 0) diagRange(-nlo(),0).setZero();
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == rowsize()-1);
            m2.assignToL(
                LowerTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
            if (nhi() > 0) diagRange(1,nhi()+1).setZero();
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == rowsize()-1);
            m2.assignToL(
                LowerTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
            if (nhi() > 0) diagRange(1,nhi()+1).setZero();
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        {
            const ptrdiff_t n = BandNumElements(colsize(),rowsize(),nlo(),nhi());
            return MyListAssigner(rowmajor_begin(),n,x);
        }

        //
        // Access
        //

        inline reference operator()(ptrdiff_t i, ptrdiff_t j)
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(okij(i,j));
            return ref(i,j);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return vec_type(ptr()+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                            j2-j1,stepj(),ct() TMV_FIRSTLAST);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return vec_type(ptr()+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),
                            i2-i1,stepi(),ct() TMV_FIRSTLAST );
        }

        inline vec_type diag()
        {
            return vec_type(ptr(),TMV_MIN(colsize(),rowsize()),
                            diagstep(),ct() TMV_FIRSTLAST);
        }

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            if (i >= 0) {
                const ptrdiff_t diaglen = TMV_MIN(rowsize()-i,colsize());
                return vec_type(ptr()+i*ptrdiff_t(stepj()),diaglen,diagstep(),ct()
                                TMV_FIRSTLAST );
            } else {
                const ptrdiff_t diaglen = TMV_MIN(colsize()+i,rowsize());
                return vec_type(ptr()-i*ptrdiff_t(stepi()),diaglen,diagstep(),ct()
                                TMV_FIRSTLAST );
            }
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            TMVAssert(j1>=0 && j1-j2<=0);
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(rowsize()-i,colsize()));
                return vec_type(ptr()+i*ptrdiff_t(stepj())+j1*diagstep(),
                                j2-j1, diagstep(),ct() TMV_FIRSTLAST );
            } else {
                TMVAssert(j2<=TMV_MIN(colsize()+i,rowsize()));
                return vec_type(ptr()-i*ptrdiff_t(stepi())+j1*diagstep(),
                                j2-j1, diagstep(),ct() TMV_FIRSTLAST );
            }
        }

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        { return base::operator()(i,j); }
        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::row(i,j1,j2); }
        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        { return base::col(j,i1,i2); }
        inline const_vec_type diag() const
        { return base::diag(); }
        inline const_vec_type diag(ptrdiff_t i) const
        { return base::diag(i); }
        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::diag(i,j1,j2); }

        //
        // Modifying Functions
        //

        type& setZero();

        type& setAllTo(const T& x);

        type& addToAll(const T& x);

        type& clip(RT thresh);

        void doTransposeSelf();
        inline type& transposeSelf()
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(nlo() == nhi());
            doTransposeSelf();
            return *this;
        }

        type& conjugateSelf();

        inline type& setToIdentity(const T& x=T(1))
        {
            TMVAssert(colsize() == rowsize());
            setZero(); diag().setAllTo(x); return *this;
        }


        //
        // subBandMatrix
        //

        inline rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            return rec_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, stepi(), stepj(), ct() TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return rec_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                ct() TMV_FIRSTLAST);
        }

        inline vec_type subVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size)
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,size));
            return vec_type(
                ptr()+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),size,
                istep*stepi()+jstep*stepj(), ct() TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi)
        {
            TMVAssert(base::hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct() TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            const ptrdiff_t newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const ptrdiff_t newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            const ptrdiff_t newstepi = stepi()*istep;
            const ptrdiff_t newstepj = stepj()*jstep;
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi,
                stepi()*istep, stepj()*jstep, newstepi+newstepj, ct()
                TMV_FIRSTLAST);
        }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());

            const ptrdiff_t j1 = i1 > nlo() ? i1-nlo() : 0;
            const ptrdiff_t j2 = TMV_MIN(i2 + nhi(),rowsize());
            const ptrdiff_t newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const ptrdiff_t newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const ptrdiff_t newlin = (ls() && isrm()) ? -1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());

            const ptrdiff_t i1 = j1 > nhi() ? j1-nhi() : 0;
            const ptrdiff_t i2 = TMV_MIN(j2 + nlo(),colsize());
            const ptrdiff_t newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const ptrdiff_t newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            const ptrdiff_t newlin = (ls() && iscm()) ? -1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type diagRange(ptrdiff_t k1, ptrdiff_t k2)
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);

            const ptrdiff_t i1 = k2 <= 0 ? -k2+1 : 0;
            const ptrdiff_t i2 = TMV_MIN(rowsize()-k1,colsize());
            const ptrdiff_t j1 = k1 <= 0 ? 0 : k1;
            const ptrdiff_t j2 = TMV_MIN(rowsize(),colsize()+k2-1);
            const ptrdiff_t newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const ptrdiff_t newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct() TMV_FIRSTLAST);
        }

        inline view_type upperBand()
        {
            return view_type(
                ptr(),TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBand()
        {
            return view_type(
                ptr(),TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline view_type upperBandOff()
        {
            return view_type(
                ptr()+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBandOff()
        {
            return view_type(
                ptr()+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(), NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        {
            TMVAssert(isComplex(T()));
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline view_type view()
        { return *this; }

        inline view_type transpose()
        {
            return view_type(
                ptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),ct(),
                isdm()?0:ls() TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return view_type(
                ptr(),colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),TMV_ConjOf(T,ct()),
                isdm()?0:ls() TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return view_type(
                ptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),
                TMV_ConjOf(T,ct()),isdm()?0:ls() TMV_FIRSTLAST);
        }

        inline vec_type linearView()
        {
            TMVAssert(iscm() || isrm());
            TMVAssert(ls() != -1);
            // (To assure that next assert has no effect.)
            TMVAssert(canLinearize());

#ifdef TMV_USE_VALGRIND
            // Valgrind will complain about using a BandMatrix linearView
            // since there are some values that are not initialized.
            // These are ok to be not initialized, since they are also
            // never used for anything, so this ifdef will assign -777
            // to each of these values to make sure valgrind doesn't give
            // a conditional jump error when they are accessed.
            // I didn't try to make this at all efficient, since it's not
            // ever used in production versions of the code.
            // e.g. there should be a guard like if (isSetupForValgrind) {...}
            // but I haven't bothered to do that.
            std::vector<bool> ok(ls(),false);
            for (ptrdiff_t i=0;i<colsize();++i) {
                for (ptrdiff_t j=0;j<rowsize();++j) {
                    if (okij(i,j)) {
                        ptrdiff_t k = i*stepi() + j*stepj();
                        ok[k] = true;
                    }
                }
            }
            for (ptrdiff_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return vec_type(ptr(),ls(),1,ct() TMV_FIRSTLAST );
        }


        inline const_rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::subMatrix(i1,i2,j1,j2); }
        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::subMatrix(i1,i2,j1,j2,istep,jstep); }
        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        { return base::subVector(i,j,istep,jstep,size); }
        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        { return base::subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::subBandMatrix(i1,i2,j1,j2); }
        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::rowRange(i1,i2); }
        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        { return base::colRange(j1,j2); }
        inline const_view_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        { return base::diagRange(k1,k2); }
        inline const_view_type upperBand() const
        { return base::upperBand(); }
        inline const_view_type lowerBand() const
        { return base::lowerBand(); }
        inline const_view_type upperBandOff() const
        { return base::upperBandOff(); }
        inline const_view_type lowerBandOff() const
        { return base::lowerBandOff(); }
        inline const_realpart_type realPart() const
        { return base::realPart(); }
        inline const_realpart_type imagPart() const
        { return base::imagPart(); }
        inline const_view_type view() const
        { return base::view(); }
        inline const_view_type transpose() const
        { return base::transpose(); }
        inline const_view_type conjugate() const
        { return base::conjugate(); }
        inline const_view_type adjoint() const
        { return base::adjoint(); }
        inline const_vec_type linearView() const
        { return base::linearView(); }


        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        inline ptrdiff_t colsize() const { return itscs; }
        inline ptrdiff_t rowsize() const { return itsrs; }
        inline ptrdiff_t nlo() const { return itsnlo; }
        inline ptrdiff_t nhi() const { return itsnhi; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }
        inline ptrdiff_t stepi() const { return itssi; }
        inline ptrdiff_t stepj() const { return itssj; }
        inline ptrdiff_t diagstep() const { return itssd; }
        using base::isrm;
        using base::iscm;
        using base::isdm;
        using base::isconj;
        inline ptrdiff_t ls() const { return linsize; }
        inline ConjType ct() const { return itsct; }

        bool canLinearize() const;

        reference ref(ptrdiff_t i, ptrdiff_t j);

        inline rowmajor_iterator rowmajor_begin()
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end()
        {
            return rowmajor_iterator(
                this,TMV_MIN(colsize(),rowsize()+nlo()),
                this->rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }

        inline colmajor_iterator colmajor_begin()
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end()
        {
            return colmajor_iterator(
                this,this->colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }

        inline diagmajor_iterator diagmajor_begin()
        { return diagmajor_iterator(this,nlo(),0); }
        inline diagmajor_iterator diagmajor_end()
        { return diagmajor_iterator(this,0,nhi()+1); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        {
            return const_rowmajor_iterator(
                this,TMV_MIN(colsize(),rowsize()+nlo()),
                this->rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        {
            return const_colmajor_iterator(
                this,this->colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }

        inline const_diagmajor_iterator diagmajor_begin() const
        { return const_diagmajor_iterator(this,nlo(),0); }
        inline const_diagmajor_iterator diagmajor_end() const
        { return const_diagmajor_iterator(this,0,nhi()+1); }

    protected:

        T*const itsm;
        const ptrdiff_t itscs;
        const ptrdiff_t itsrs;
        const ptrdiff_t itsnlo;
        const ptrdiff_t itsnhi;
        const ptrdiff_t itssi;
        const ptrdiff_t itssj;
        const ptrdiff_t itssd;

        const ConjType itsct;
        mutable ptrdiff_t linsize;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
    protected :
#endif

        using base::okij;

    }; // BandMatrixView

    template <typename T>
    class BandMatrixView<T,FortranStyle> : public BandMatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenBandMatrix<T> base;
        typedef BandMatrixView<T,FortranStyle> type;
        typedef BandMatrixView<T,CStyle> c_type;
        typedef BandMatrixView<T,FortranStyle> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,FortranStyle> vec_type;
        typedef MatrixView<T,FortranStyle> rec_type;
        typedef BandMatrixView<RT,FortranStyle> realpart_type;
        typedef ConstBandMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstBandMatrixView<RT,FortranStyle> const_realpart_type;
        typedef typename RefHelper<T>::reference reference;

        //
        // Constructors
        //

        inline BandMatrixView(const type& rhs) : c_type(rhs) {}

        inline BandMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline BandMatrixView(
            T* _m, ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _lo, ptrdiff_t _hi,
            ptrdiff_t _si, ptrdiff_t _sj, ptrdiff_t _sd, ConjType _ct,
            ptrdiff_t _ls=0 TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,_ct,_ls
                   TMV_FIRSTLAST1(_first,_last) ) {}

        inline BandMatrixView(BandMatrixView<T> rhs, ptrdiff_t lo, ptrdiff_t hi) :
            c_type(rhs,lo,hi) {}

        inline BandMatrixView(MatrixView<T> rhs, ptrdiff_t lo, ptrdiff_t hi) :
            c_type(rhs,lo,hi) {}

        virtual inline ~BandMatrixView() {}

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const c_type& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenBandMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenBandMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        template <typename T2>
        inline type& operator=(const GenBandMatrix<T2>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const T& x)
        { c_type::operator=(x); return *this; }

        inline type& operator=(const AssignableToBandMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const AssignableToBandMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenDiagMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenUpperTriMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenUpperTriMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenLowerTriMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenLowerTriMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        typedef typename c_type::MyListAssigner MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return c_type::operator<<(x); }

        //
        // Access
        //

        inline reference operator()(ptrdiff_t i,ptrdiff_t j)
        {
            TMVAssert(i>0 && i<=c_type::colsize());
            TMVAssert(j>0 && j<=c_type::rowsize());
            TMVAssert(c_type::okij(i-1,j-1));
            return c_type::ref(i-1,j-1);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>0 && i<=c_type::colsize());
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= c_type::rowsize());
            TMVAssert(c_type::okij(i-1,j1-1));
            TMVAssert(c_type::okij(i-1,j2-1));
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(j>0 && j<=c_type::rowsize());
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= c_type::colsize());
            TMVAssert(c_type::okij(i1-1,j-1));
            TMVAssert(c_type::okij(i2-1,j-1));
            return c_type::col(j-1,i1-1,i2);
        }

        inline vec_type diag()
        { return c_type::diag(); }

        inline vec_type diag(ptrdiff_t i)
        { return c_type::diag(i); }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1>0);
            return c_type::diag(i,j1-1,j2);
        }

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            TMVAssert(i>0 && i<=c_type::colsize());
            TMVAssert(j>0 && j<=c_type::rowsize());
            TMVAssert(c_type::okij(i-1,j-1));
            return c_type::cref(i-1,j-1);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>0 && i<=c_type::colsize());
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= c_type::rowsize());
            TMVAssert(c_type::okij(i-1,j1-1));
            TMVAssert(c_type::okij(i-1,j2-1));
            return c_type::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>0 && j<=c_type::rowsize());
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= c_type::colsize());
            TMVAssert(c_type::okij(i1-1,j-1));
            TMVAssert(c_type::okij(i2-1,j-1));
            return c_type::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag() const
        { return c_type::diag(); }

        inline const_vec_type diag(ptrdiff_t i) const
        { return c_type::diag(i); }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>0);
            return c_type::diag(i,j1-1,j2);
        }


        //
        // Modifying Functions
        //

        inline type& setZero()
        { c_type::setZero(); return *this; }

        inline type& setAllTo(const T& x)
        { c_type::setAllTo(x); return *this; }

        inline type& addToAll(const T& x)
        { c_type::addToAll(x); return *this; }

        inline type& clip(RT thresh)
        { c_type::clip(thresh); return *this; }

        inline type& transposeSelf()
        { c_type::transposeSelf(); return *this; }

        inline type& conjugateSelf()
        { c_type::conjugateSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1))
        { c_type::setToIdentity(x); return *this; }

        //
        // subBandMatrix
        //

        inline bool hasSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        {
            return const_view_type(*this).hasSubVector(i,j,istep,jstep,size);
        }

        inline bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_view_type(*this).hasSubMatrix(
                i1,i2,j1,j2,istep,jstep);
        }

        inline bool hasSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2,
            ptrdiff_t newnlo, ptrdiff_t newnhi, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_view_type(*this).hasSubBandMatrix(
                i1,i2,j1,j2,newnlo,newnhi,istep,jstep);
        }

        inline rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::subMatrix(i1-1,i2,j1-1,j2);
        }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::subMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size)
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return c_type::subVector(i-1,j-1,istep,jstep,size);
        }

        inline view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi)
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return c_type::subBandMatrix(i1-1,i2,j1-1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            const ptrdiff_t newnlo = TMV_MIN(c_type::nlo()+j1-i1,i2-i1);
            const ptrdiff_t newnhi = TMV_MIN(c_type::nhi()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            return c_type::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
        }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=c_type::colsize());
            return c_type::rowRange(i1-1,i2);
        }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1>0 && j1-j2<=0 && j2<=c_type::rowsize());
            return c_type::colRange(j1-1,j2);
        }

        inline view_type diagRange(ptrdiff_t k1, ptrdiff_t k2)
        { return c_type::diagRange(k1,k2+1); }

        inline view_type upperBand()
        { return c_type::upperBand(); }

        inline view_type lowerBand()
        { return c_type::lowerBand(); }

        inline view_type upperBandOff()
        { return c_type::upperBandOff(); }

        inline view_type lowerBandOff()
        { return c_type::lowerBandOff(); }

        inline realpart_type realPart()
        { return c_type::realPart(); }

        inline realpart_type imagPart()
        { return c_type::imagPart(); }

        inline view_type view()
        { return c_type::view(); }

        inline view_type transpose()
        { return c_type::transpose(); }

        inline view_type conjugate()
        { return c_type::conjugate(); }

        inline view_type adjoint()
        { return c_type::adjoint(); }

        inline vec_type linearView()
        { return c_type::linearView(); }


        inline const_rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::subMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::subMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return c_type::subVector(i-1,j-1,istep,jstep,size);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return c_type::subBandMatrix(i1-1,i2,j1-1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            const ptrdiff_t newnlo = TMV_MIN(c_type::nlo()+j1-i1,i2-i1);
            const ptrdiff_t newnhi = TMV_MIN(c_type::nhi()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            return c_type::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=c_type::colsize());
            return c_type::rowRange(i1-1,i2);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>0 && j1-j2<=0 && j2<=c_type::rowsize());
            return c_type::colRange(j1-1,j2);
        }

        inline const_view_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        { return c_type::diagRange(k1,k2+1); }

        inline const_view_type upperBand() const
        { return c_type::upperBand(); }

        inline const_view_type lowerBand() const
        { return c_type::lowerBand(); }

        inline const_view_type upperBandOff() const
        { return c_type::upperBandOff(); }

        inline const_view_type lowerBandOff() const
        { return c_type::lowerBandOff(); }

        inline const_realpart_type realPart() const
        { return c_type::realPart(); }

        inline const_realpart_type imagPart() const
        { return c_type::imagPart(); }

        inline const_view_type view() const
        { return c_type::view(); }

        inline const_view_type transpose() const
        { return c_type::transpose(); }

        inline const_view_type conjugate() const
        { return c_type::conjugate(); }

        inline const_view_type adjoint() const
        { return c_type::adjoint(); }

        inline const_vec_type linearView() const
        { return c_type::linearView(); }

    }; // FortranStyle BandMatrixView

    template <typename T, int A>
    class BandMatrix : public GenBandMatrix<T>
    {
    public:

        enum { S = A & AllStorageType };
        enum { I = A & FortranStyle };
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenBandMatrix<T> base;
        typedef BandMatrix<T,A> type;
        typedef BandMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef BandMatrixView<RT,I> realpart_type;
        typedef ConstBandMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstBandMatrixView<RT,I> const_realpart_type;
        typedef ConstVectorView<T> const_c_vec_type;
        typedef T& reference;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef DMIt<type> diagmajor_iterator;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef CDMIt<type> const_diagmajor_iterator;

        //
        // Constructors
        //

#define NEW_SIZE(cs,rs,lo,hi) \
        linsize(BandStorageLength(static_cast<StorageType>(S),cs,rs,lo,hi)), \
        itsm1(linsize), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
        itssi(S==int(RowMajor) ? lo+hi : S==int(ColMajor) ? 1 : \
              rs>= cs ? 1-cs : -rs ), \
        itssj(S==int(RowMajor) ? 1 : S==int(ColMajor) ? lo+hi : -itssi+1), \
        itsds(S==int(RowMajor) ? itssi+1 : S==int(ColMajor) ? itssj+1 : 1), \
        itsm(S==int(DiagMajor) ? itsm1.get() - lo*ptrdiff_t(itssi) : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize)

#define NEW_SIZE2(ls,cs,rs,lo,hi) \
        linsize(ls), \
        itsm1(linsize), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
        itssi(S==int(RowMajor) ? lo+hi : S==int(ColMajor) ? 1 : \
              rs>= cs ? 1-cs : -rs ), \
        itssj(S==int(RowMajor) ? 1 : S==int(ColMajor) ? lo+hi : -itssi+1), \
        itsds(S==int(RowMajor) ? itssi+1 : S==int(ColMajor) ? itssj+1 : 1), \
        itsm(S==int(DiagMajor) ? itsm1.get() - lo*ptrdiff_t(itssi) : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize)

        inline BandMatrix() : NEW_SIZE(0,0,0,0)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
        }

        inline BandMatrix(ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi) :
            NEW_SIZE(cs,rs,lo,hi)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(cs >= 0 && rs >= 0);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < cs);
            TMVAssert(hi < rs);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline BandMatrix(ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi, const T& x) :
            NEW_SIZE(cs,rs,lo,hi)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(cs >= 0 && rs >= 0);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < cs);
            TMVAssert(hi < rs);
            setAllTo(x);
        }

        inline BandMatrix(const BandMatrix<T,A>& m2) :
            NEW_SIZE2(m2.ls(),m2.itscs,m2.itsrs,m2.itsnlo,m2.itsnhi)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            std::copy(m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
        }

        inline BandMatrix(const GenBandMatrix<RT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            m2.assignToB(view());
        }

        inline BandMatrix(const GenBandMatrix<CT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isComplex(T()));
            m2.assignToB(view());
        }

        template <typename T2>
        inline BandMatrix(const GenBandMatrix<T2>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isReal(T2()) || isComplex(T()));
            Copy(m2,view());
        }

        template <typename T2>
        inline BandMatrix(const GenBandMatrix<T2>& m2, ptrdiff_t lo, ptrdiff_t hi) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),lo,hi)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo <= m2.nlo());
            TMVAssert(hi <= m2.nhi());
            Copy(ConstBandMatrixView<T2>(m2,lo,hi),view());
            if (I==int(CStyle)) {
                if (lo > m2.nlo()) diagRange(-lo,-m2.nlo()).setZero();
                if (hi > m2.nhi()) diagRange(m2.nhi()+1,hi+1).setZero();
            } else {
                if (lo > m2.nlo()) diagRange(-lo,-m2.nlo()-1).setZero();
                if (hi > m2.nhi()) diagRange(m2.nhi()+1,hi).setZero();
            }
        }

        inline BandMatrix(const GenMatrix<T>& m2, ptrdiff_t lo, ptrdiff_t hi) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),lo,hi)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < m2.colsize());
            TMVAssert(hi < m2.rowsize());
            Copy(ConstBandMatrixView<T>(m2,lo,hi),view());
        }

        inline BandMatrix(const GenUpperTriMatrix<T>& m2, ptrdiff_t hi) :
            NEW_SIZE(m2.size(),m2.size(),0,hi)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(hi >= 0);
            TMVAssert(hi < m2.size());
            Copy(BandMatrixViewOf(m2,hi),view());
        }

        inline BandMatrix(const GenLowerTriMatrix<T>& m2, ptrdiff_t lo) :
            NEW_SIZE(m2.size(),m2.size(),lo,0)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(lo >= 0);
            TMVAssert(lo < m2.size());
            Copy(BandMatrixViewOf(m2,lo),view());
        }

        inline BandMatrix(const AssignableToBandMatrix<RT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            m2.assignToB(view());
        }

        inline BandMatrix(const AssignableToBandMatrix<CT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isComplex(T()));
            m2.assignToB(view());
        }

        inline BandMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,0)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        inline BandMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,0)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isComplex(T()));
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        inline BandMatrix(const GenUpperTriMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,m2.size()-1)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            setZero();
            m2.assignToU(
                UpperTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
        }

        inline BandMatrix(const GenUpperTriMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,m2.size()-1)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isComplex(T()));
            setZero();
            m2.assignToU(
                UpperTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
        }

        inline BandMatrix(const GenLowerTriMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),m2.size()-1,0)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            setZero();
            m2.assignToL(
                LowerTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
        }

        inline BandMatrix(const GenLowerTriMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),m2.size()-1,0)
        {
            TMVAssert(Attrib<A>::bandmatrixok);
            TMVAssert(isComplex(T()));
            setZero();
            m2.assignToL(
                LowerTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,ct() TMV_FIRSTLAST));
        }

#undef NEW_SIZE
#undef NEW_SIZE2

        virtual inline ~BandMatrix()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(999));
#endif
            itsm = 0;
        }


        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(view());
            return *this;
        }

        inline type& operator=(const GenBandMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(view());
            return *this;
        }

        inline type& operator=(const GenBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(view());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const GenBandMatrix<T2>& m2)
        {
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            Copy(m2,view());
            return *this;
        }

        inline type& operator=(const T& x)
        {
            TMVAssert(colsize() == rowsize());
            return setToIdentity(x);
        }

        inline type& operator=(const AssignableToBandMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(view());
            return *this;
        }

        inline type& operator=(const AssignableToBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(view());
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == rowsize()-1);
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == rowsize()-1);
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == rowsize()-1);
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == rowsize()-1);
            view() = m2;
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        {
            const ptrdiff_t n = BandNumElements(colsize(),rowsize(),nlo(),nhi());
            return MyListAssigner(rowmajor_begin(),n,x);
        }

        //
        // Access
        //

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            if (I==int(CStyle)) {
                TMVAssert(i>=0 && i<colsize());
                TMVAssert(j>=0 && j<rowsize());
                return cref(i,j);
            } else {
                TMVAssert(i>0 && i<=colsize());
                TMVAssert(j>0 && j<=rowsize());
                return cref(i-1,j-1);
            }
        }

        inline T& operator()(ptrdiff_t i,ptrdiff_t j)
        {
            if (I==int(CStyle)) {
                TMVAssert(i>=0 && i<colsize());
                TMVAssert(j>=0 && j<rowsize());
                TMVAssert(okij(i,j));
                return ref(i,j);
            } else {
                TMVAssert(i>0 && i<=colsize());
                TMVAssert(j>0 && j<=rowsize());
                TMVAssert(okij(i-1,j-1));
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I == int(FortranStyle)) {
                TMVAssert(i>0 && i<=colsize()); --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); --j1;
            } else {
                TMVAssert(i>=0 && i<colsize());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return const_vec_type(
                itsm+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                j2-j1,stepj(),NonConj);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I == int(FortranStyle)) {
                TMVAssert(j>0 && j<=rowsize()); --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1;
            } else {
                TMVAssert(j>=0 && j<rowsize());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return const_vec_type(
                itsm+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),
                i2-i1,stepi(),NonConj);
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                itsm,TMV_MIN(colsize(),rowsize()),diagstep(),NonConj);
        }

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    itsm+i*ptrdiff_t(stepj()),diagsize,diagstep(),NonConj);
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    itsm-i*ptrdiff_t(stepi()),diagsize,diagstep(),NonConj);
            }
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0);
            }
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(rowsize()-i,colsize()));
                return const_vec_type(
                    itsm+i*ptrdiff_t(stepj())+j1*ptrdiff_t(diagstep()),
                    j2-j1,diagstep(),NonConj);
            } else {
                TMVAssert(j2<=TMV_MIN(colsize()+i,rowsize()));
                return const_vec_type(
                    itsm-i*ptrdiff_t(stepi())+j1*ptrdiff_t(diagstep()),
                    j2-j1,diagstep(),NonConj);
            }
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(FortranStyle)) {
                TMVAssert(i>0 && i<=colsize()); --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); --j1;
            } else {
                TMVAssert(i>=0 && i<colsize());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return vec_type(
                itsm+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                j2-j1,stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(FortranStyle)) {
                TMVAssert(j>0 && j<=rowsize()); --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1;
            } else {
                TMVAssert(j>=0 && j<rowsize());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return vec_type(
                itsm+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),
                i2-i1,stepi(),NonConj TMV_FIRSTLAST );
        }

        inline vec_type diag()
        {
            return vec_type(
                itsm,TMV_MIN(colsize(),rowsize()),diagstep(),NonConj
                TMV_FIRSTLAST);
        }

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return vec_type(
                    itsm+i*ptrdiff_t(stepj()),diagsize,diagstep(),NonConj TMV_FIRSTLAST);
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return vec_type(
                    itsm-i*ptrdiff_t(stepi()),diagsize,diagstep(),NonConj TMV_FIRSTLAST);
            }
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-colsize() && i<=rowsize());
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0);
            }
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(rowsize()-i,colsize()));
                return vec_type(
                    itsm+i*ptrdiff_t(stepj())+j1*ptrdiff_t(diagstep()),j2-j1,diagstep(),NonConj
                    TMV_FIRSTLAST);
            } else {
                TMVAssert(j2<=TMV_MIN(colsize()+i,rowsize()));
                return vec_type(
                    itsm-i*ptrdiff_t(stepi())+j1*ptrdiff_t(diagstep()),j2-j1,diagstep(),NonConj
                    TMV_FIRSTLAST);
            }
        }

        //
        // Modifying Functions
        //

        inline type& setZero()
        { linearView().setZero(); return *this; }

        inline type& setAllTo(const T& x)
        { linearView().setAllTo(x); return *this; }

        inline type& addToAll(const T& x)
        { linearView().addToAll(x); return *this; }

        inline type& clip(RT thresh)
        { linearView().clip(thresh); return *this; }

        inline type& transposeSelf()
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(nlo()==nhi());
            view().transposeSelf();
            return *this;
        }

        inline type& conjugateSelf()
        { linearView().conjugateSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1))
        {
            TMVAssert(colsize() == rowsize());
            setZero(); diag().setAllTo(x);
            return *this;
        }

        //
        // subBandMatrix
        //

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I == int(FortranStyle)) { --i1; --j1; }
            return const_rec_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()), i2-i1, j2-j1,
                stepi(), stepj(), NonConj);
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I == int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return const_rec_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(),jstep*stepj(), NonConj);
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I == int(FortranStyle)) { --i; --j; }
            return const_vec_type(
                itsm+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()), size,
                istep*stepi()+jstep*stepj(), NonConj);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I == int(FortranStyle)) { --i1; --j1; }
            return const_view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()), i2-i1, j2-j1,
                newnlo, newnhi, stepi(), stepj(), diagstep(), NonConj);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            const ptrdiff_t newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==int(CStyle)?1:0));
            const ptrdiff_t newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-(I==int(CStyle)?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(view().hasSubBandMatrix(
                    i1,i2,j1,j2, newnlo,newnhi,istep,jstep));
            if (I == int(FortranStyle)) {
                --i1; --j1; i2+=istep-1; j2+=jstep-1;
            }
            const ptrdiff_t newstepi = stepi()*istep;
            const ptrdiff_t newstepj = stepj()*jstep;
            return const_view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi,
                newstepi, newstepj, newstepi+newstepj, NonConj);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1;
            } else {
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }

            const ptrdiff_t j1 = i1 > nlo() ? i1-nlo() : 0;
            const ptrdiff_t j2 = TMV_MIN(i2 + nhi(),rowsize());
            const ptrdiff_t newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const ptrdiff_t newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const ptrdiff_t newlin = isrm() ? -1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }

            const ptrdiff_t i1 = j1 > nhi() ? j1-nhi() : 0;
            const ptrdiff_t i2 = TMV_MIN(j2 + nlo(),colsize());
            const ptrdiff_t newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const ptrdiff_t newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            const ptrdiff_t newlin = iscm() ? -1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin);
        }

        inline const_view_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(k1>=-nlo() && k1<=k2 && k2<=nhi()); ++k2;
            } else {
                TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);
            }

            const ptrdiff_t i1 = k2 <= 0 ? -k2+1 : 0;
            const ptrdiff_t i2 = TMV_MIN(rowsize()-k1,colsize());
            const ptrdiff_t j1 = k1 <= 0 ? 0 : k1;
            const ptrdiff_t j2 = TMV_MIN(rowsize(),colsize()+k2-1);
            const ptrdiff_t newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const ptrdiff_t newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct());
        }

        inline const_view_type upperBand() const
        {
            return const_view_type(
                itsm,TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),ct());
        }

        inline const_view_type lowerBand() const
        {
            return const_view_type(
                itsm,TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),ct());
        }

        inline const_view_type upperBandOff() const
        {
            return const_view_type(
                itsm+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),ct());
        }

        inline const_view_type lowerBandOff() const
        {
            return const_view_type(
                itsm+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),ct());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep());
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm)+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NonConj);
        }

        inline rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I == int(FortranStyle)) { --i1; --j1; }
            return rec_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, stepi(), stepj(), NonConj TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I == int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return rec_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                NonConj TMV_FIRSTLAST);
        }

        inline vec_type subVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I == int(FortranStyle)) { --i; --j; }
            return vec_type(
                itsm+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),size,
                istep*stepi()+jstep*stepj(), NonConj TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I==int(FortranStyle)) { --i1; --j1; }
            return view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), NonConj TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            const ptrdiff_t newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==int(CStyle)?1:0));
            const ptrdiff_t newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-(I==int(CStyle)?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(view().hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            if (I==int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const ptrdiff_t newstepi = stepi()*istep;
            const ptrdiff_t newstepj = stepj()*jstep;
            return view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi,
                newstepi, newstepj, newstepi+newstepj, NonConj TMV_FIRSTLAST);
        }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1;
            } else {
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            const ptrdiff_t j1 = i1 > nlo() ? i1-nlo() : 0;
            const ptrdiff_t j2 = TMV_MIN(i2 + nhi(),rowsize());
            const ptrdiff_t newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const ptrdiff_t newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const ptrdiff_t newlin = isrm() ? -1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }

            const ptrdiff_t i1 = j1 > nhi() ? j1-nhi() : 0;
            const ptrdiff_t i2 = TMV_MIN(j2 + nlo(),colsize());
            const ptrdiff_t newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const ptrdiff_t newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            const ptrdiff_t newlin = iscm() ? -1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type diagRange(ptrdiff_t k1, ptrdiff_t k2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(k1>=-nlo() && k1<=k2 && k2<=nhi()); ++k2;
            } else {
                TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);
            }

            const ptrdiff_t i1 = k2 <= 0 ? -k2+1 : 0;
            const ptrdiff_t i2 = TMV_MIN(rowsize()-k1,colsize());
            const ptrdiff_t j1 = k1 <= 0 ? 0 : k1;
            const ptrdiff_t j2 = TMV_MIN(rowsize(),colsize()+k2-1);
            const ptrdiff_t newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const ptrdiff_t newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                itsm+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), ct()  TMV_FIRSTLAST);
        }

        inline view_type upperBand()
        {
            return view_type(
                itsm,TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBand()
        {
            return view_type(
                itsm,TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline view_type upperBandOff()
        {
            return view_type(
                itsm+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBandOff()
        {
            return view_type(
                itsm+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),ct() TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(), NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        {
            TMVAssert(isComplex(T()));
            return realpart_type(
                reinterpret_cast<RT*>(itsm)+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
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
        inline const_view_type view() const
        {
            return const_view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),NonConj,isdm()?0:linsize);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),NonConj,isdm()?0:linsize);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),TMV_ConjOf(T,NonConj),
                isdm()?0:linsize);
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),
                TMV_ConjOf(T,NonConj),isdm()?0:linsize);
        }

        inline const_c_vec_type constLinearView() const
        {
#ifdef TMV_USE_VALGRIND
            std::vector<bool> ok(linsize,false);
            for (ptrdiff_t i=0;i<colsize();++i) {
                for (ptrdiff_t j=0;j<rowsize();++j) {
                    if (okij(i,j)) {
                        ptrdiff_t k = i*stepi() + j*stepj() +
                            ptrdiff_t(cptr()-itsm1.get());
                        ok[k] = true;
                    }
                }
            }
            for (ptrdiff_t k=0;k<linsize;++k) if (!ok[k]) {
                const_cast<T*>(itsm1.get())[k] = T(-777);
            }
#endif
            return const_c_vec_type(itsm1.get(),linsize,1,NonConj);
        }

        inline view_type view()
        {
            return view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),NonConj,isdm()?0:linsize
                TMV_FIRSTLAST);
        }

        inline view_type transpose()
        {
            return view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),NonConj,isdm()?0:linsize
                TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),TMV_ConjOf(T,NonConj),
                isdm()?0:linsize TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),
                TMV_ConjOf(T,NonConj),isdm()?0:linsize TMV_FIRSTLAST);
        }

        inline vec_type linearView()
        {
#ifdef TMV_USE_VALGRIND
            std::vector<bool> ok(linsize,false);
            for (ptrdiff_t i=0;i<colsize();++i) {
                for (ptrdiff_t j=0;j<rowsize();++j) {
                    if (okij(i,j)) {
                        ptrdiff_t k = i*stepi() + j*stepj() +
                            ptrdiff_t(cptr()-itsm1.get());
                        ok[k] = true;
                    }
                }
            }
            for (ptrdiff_t k=0;k<linsize;++k) if (!ok[k]) {
                itsm1.get()[k] = T(-777);
            }
#endif
            return vec_type(itsm1.get(),linsize,1,NonConj TMV_FIRSTLAST );
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        inline ptrdiff_t colsize() const { return itscs; }
        inline ptrdiff_t rowsize() const { return itsrs; }
        inline ptrdiff_t nlo() const { return itsnlo; }
        inline ptrdiff_t nhi() const { return itsnhi; }
        inline ptrdiff_t mem_used() const { return linsize; }
        inline const T* start_mem() const { return itsm1.get(); }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }
        inline ptrdiff_t stepi() const { return itssi; }
        inline ptrdiff_t stepj() const { return itssj; }
        inline ptrdiff_t diagstep() const { return itsds; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==int(RowMajor); }
        inline bool iscm() const { return S==int(ColMajor); }
        inline bool isdm() const { return S==int(DiagMajor); }
        inline bool isconj() const { return false; }
        inline ptrdiff_t ls() const { return linsize; }

        inline bool canLinearize() const
        { return true; }

        inline T cref(ptrdiff_t i, ptrdiff_t j) const
        {
            return okij(i,j) ?
                itsm[i*ptrdiff_t(itssi) + j*ptrdiff_t(itssj)] : T(0);
        }

        inline T& ref(ptrdiff_t i, ptrdiff_t j)
        { return itsm[i*ptrdiff_t(itssi) + j*ptrdiff_t(itssj)]; }

        inline void resize(ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi)
        {
            TMVAssert(cs >= 0 && rs >= 0);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < cs);
            TMVAssert(hi < rs);
            linsize = BandStorageLength(
                static_cast<StorageType>(S),cs,rs,lo,hi);
            itsm1.resize(linsize);
            itscs = cs;
            itsrs = rs;
            itsnlo = lo;
            itsnhi = hi;
            itssi =
                S==int(RowMajor) ? lo+hi :
                S==int(ColMajor) ? 1 :
                rs>= cs ? 1-cs : -rs;
            itssj =
                S==int(RowMajor) ? 1 :
                S==int(ColMajor) ? lo+hi :
                -itssi+1;
            itsds =
                S==int(RowMajor) ? itssi+1 :
                S==int(ColMajor) ? itssj+1 :
                1;
            itsm = S==int(DiagMajor) ?
                itsm1.get()-lo*ptrdiff_t(itssi) : itsm1.get();
            DivHelper<T>::resetDivType();
#ifdef TMVFLDEBUG
            _first = itsm1.get();
            _last = _first + linsize;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline rowmajor_iterator rowmajor_begin()
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end()
        {
            return rowmajor_iterator(
                this,TMV_MIN(colsize(),rowsize()+nlo()),
                this->rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        {
            return const_rowmajor_iterator(
                this,TMV_MIN(colsize(),rowsize()+nlo()),
                this->rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }

        inline colmajor_iterator colmajor_begin()
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end()
        {
            return colmajor_iterator(
                this,this->colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        {
            return const_colmajor_iterator(
                this,this->colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }

        inline diagmajor_iterator diagmajor_begin()
        { return diagmajor_iterator(this,nlo(),0); }
        inline diagmajor_iterator diagmajor_end()
        { return diagmajor_iterator(this,0,nhi()+1); }

        inline const_diagmajor_iterator diagmajor_begin() const
        { return const_diagmajor_iterator(this,nlo(),0); }
        inline const_diagmajor_iterator diagmajor_end() const
        { return const_diagmajor_iterator(this,0,nhi()+1); }


    protected :

        ptrdiff_t linsize;
        AlignedArray<T> itsm1;
        ptrdiff_t itscs;
        ptrdiff_t itsrs;
        ptrdiff_t itsnlo;
        ptrdiff_t itsnhi;
        ptrdiff_t itssi;
        ptrdiff_t itssj;
        ptrdiff_t itsds;
        T* itsm;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected:
#endif

        inline bool okij(ptrdiff_t i, ptrdiff_t j) const
        { return (j+nlo() >= i && i+nhi() >= j); }

        friend void Swap(BandMatrix<T,A>& m1, BandMatrix<T,A>& m2)
        {
            TMVAssert(m1.colsize() == m2.colsize());
            TMVAssert(m1.rowsize() == m2.rowsize());
            TMVAssert(m1.nlo() == m2.nlo());
            TMVAssert(m1.nhi() == m2.nhi());
            m1.itsm1.swapWith(m2.itsm1);
            TMV_SWAP(m1.itsm,m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // BandMatrix

    //-------------------------------------------------------------------------

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

    template <typename T>
    BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2);

    template <typename T>
    BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2);

    template <typename T>
    BandMatrix<T,DiagMajor> TriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2,
        const GenVector<T>& v3);

    //
    // More special constructors:
    //   BandMatrixViewOf(Matrix m, nlo, nhi)
    //   BandMatrixViewOf(BandMatrix m, nlo, nhi)
    //   BandMatrixViewOf(DiagMatrix m)
    //   BandMatrixViewOf(UpperTriMatrix m)
    //   BandMatrixViewOf(UpperTriMatrix m, nhi)
    //   BandMatrixViewOf(LowerTriMatrix m)
    //   BandMatrixViewOf(LowerTriMatrix m, nlo)
    //   BandMatrixViewOf(T* m, cs, rs, nlo, nhi, stor)
    //   BandMatrixViewOf(T* m, cs, rs, nlo, nhi, stepi, stepj)
    //

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenMatrix<T>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        return ConstBandMatrixView<T>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> BandMatrixViewOf(
        const ConstMatrixView<T,A>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        return ConstBandMatrixView<T,A>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        const Matrix<T,A>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        return ConstBandMatrixView<T,A&FortranStyle>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline BandMatrixView<T,A> BandMatrixViewOf(
        MatrixView<T,A> m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        return BandMatrixView<T,A>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        Matrix<T,A>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        return BandMatrixView<T,A&FortranStyle>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenBandMatrix<T>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi());
        return ConstBandMatrixView<T>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> BandMatrixViewOf(
        const ConstBandMatrixView<T,A>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi());
        return ConstBandMatrixView<T,A>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        const BandMatrix<T,A>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi());
        return ConstBandMatrixView<T,A&FortranStyle>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline BandMatrixView<T,A> BandMatrixViewOf(
        BandMatrixView<T,A> m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi());
        return BandMatrixView<T,A>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        BandMatrix<T,A>& m, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi());
        return BandMatrixView<T,A&FortranStyle>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(const GenDiagMatrix<T>& m)
    {
        return ConstBandMatrixView<T>(
            m.diag().cptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),m.diag().ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> BandMatrixViewOf(
        const ConstDiagMatrixView<T,A>& m)
    {
        return ConstBandMatrixView<T,A>(
            m.diag().cptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),m.diag().ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        const DiagMatrix<T,A>& m)
    {
        return ConstBandMatrixView<T,A&FortranStyle>(
            m.diag().cptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),m.diag().ct());
    }

    template <typename T, int A>
    inline BandMatrixView<T,A> BandMatrixViewOf(DiagMatrixView<T,A> m)
    {
        return BandMatrixView<T,A>(
            m.diag().ptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),m.diag().ct()
            TMV_FIRSTLAST1(m.diag()._first,m.diag()._last) );
    }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        DiagMatrix<T,A>& m)
    {
        return BandMatrixView<T,A&FortranStyle>(
            m.diag().ptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),m.diag().ct()
            TMV_FIRSTLAST1(m.diag()._first,m.diag()._last) );
    }

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenUpperTriMatrix<T>& m, ptrdiff_t nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T>(
            m.cptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> BandMatrixViewOf(
        const ConstUpperTriMatrixView<T,A>& m, ptrdiff_t nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T,A>(
            m.cptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        const UpperTriMatrix<T,A>& m, ptrdiff_t nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T,A&FortranStyle>(
            m.cptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline BandMatrixView<T,A> BandMatrixViewOf(
        UpperTriMatrixView<T,A> m, ptrdiff_t nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return BandMatrixView<T,A>(
            m.ptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        UpperTriMatrix<T,A>& m, ptrdiff_t nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return BandMatrixView<T,A&FortranStyle>(
            m.ptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenLowerTriMatrix<T>& m, ptrdiff_t nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T>(
            m.cptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> BandMatrixViewOf(
        const ConstLowerTriMatrixView<T,A>& m, ptrdiff_t nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T,A>(
            m.cptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        const LowerTriMatrix<T,A>& m, ptrdiff_t nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T,A&FortranStyle>(
            m.cptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct());
    }

    template <typename T, int A>
    inline BandMatrixView<T,A> BandMatrixViewOf(
        LowerTriMatrixView<T,A> m, ptrdiff_t nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return BandMatrixView<T,A>(
            m.ptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> BandMatrixViewOf(
        LowerTriMatrix<T,A>& m, ptrdiff_t nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return BandMatrixView<T,A&FortranStyle>(
            m.ptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) );
    }

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t nlo, ptrdiff_t nhi, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo >= 0);
        TMVAssert(nhi >= 0);
        TMVAssert(nlo < cs);
        TMVAssert(nhi < rs);
        const ptrdiff_t stepi = (
            stor == RowMajor ? nlo+nhi :
            stor == ColMajor ? 1 :
            rs >= cs ? -cs+1 : -rs );
        const ptrdiff_t stepj = (
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo+nhi :
            rs >= cs ? cs : rs+1 );
        T* m0 = (stor == DiagMajor) ? m - nlo*ptrdiff_t(stepi) : m;
        return ConstBandMatrixView<T>(
            m0,cs,rs,nlo,nhi,stepi,stepj,NonConj);
    }

    template <typename T>
    inline BandMatrixView<T> BandMatrixViewOf(
        T* m, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t nlo, ptrdiff_t nhi, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo >= 0);
        TMVAssert(nhi >= 0);
        TMVAssert(nlo < cs);
        TMVAssert(nhi < rs);
        const ptrdiff_t stepi = (
            stor == RowMajor ? nlo+nhi :
            stor == ColMajor ? 1 :
            rs >= cs ? -cs+1 : -rs );
        const ptrdiff_t stepj = (
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo+nhi :
            rs >= cs ? cs : rs+1 );
        T* m0 = (stor == DiagMajor) ? m - nlo*ptrdiff_t(stepi) : m;
        return BandMatrixView<T>(
            m0,cs,rs,nlo,nhi,stepi,stepj,stepi+stepj,NonConj
            TMV_FIRSTLAST1(m,m+BandStorageLength(stor,cs,rs,nlo,nhi)));
    }

    template <typename T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t nlo, ptrdiff_t nhi, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo >= 0);
        TMVAssert(nhi >= 0);
        TMVAssert(nlo < cs);
        TMVAssert(nhi < rs);
        return BandMatrixView<T>(
            m,cs,rs,nlo,nhi,stepi,stepj,stepi+stepj,NonConj);
    }

    template <typename T>
    inline BandMatrixView<T> BandMatrixViewOf(
        T* m, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t nlo, ptrdiff_t nhi, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo >= 0);
        TMVAssert(nhi >= 0);
        TMVAssert(nlo < cs);
        TMVAssert(nhi < rs);
        return BandMatrixView<T>(
            m,cs,rs,nlo,nhi,stepi,stepj,stepi+stepj,NonConj
            TMV_FIRSTLAST1(
                m + (stepi < 0 ?
                     (stepj < 0 ? stepi*(cs-1) + stepj*(rs-1) : stepi*nlo) :
                     (stepj < 0 ? stepj*nhi : 0)),
                m + (stepi > 0 ?
                     (stepj > 0 ? stepi*(cs-1) + stepj*(rs-1) : stepi*nlo) :
                     (stepj > 0 ? stepj*nhi : 0)) + 1) );
    }


    //
    // Copy
    //

    template <typename T1, typename T2>
    inline void DoCopy(
        const GenBandMatrix<T1>& m1, BandMatrixView<T2> m2)
    {
        TMVAssert(isReal(T1()) || isComplex(T2()));
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.nlo() == m1.nlo());
        TMVAssert(m2.nhi() == m1.nhi());
        TMVAssert(m2.colsize() > 0);
        TMVAssert(m2.rowsize() > 0);
        TMVAssert(m2.ct()==NonConj);
        TMVAssert(!m2.isSameAs(m1));

        if (m1.canLinearize() && m2.canLinearize() &&
            m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) {
            TMVAssert(m1.constLinearView().size() == m2.linearView().size());
            m2.linearView() = m1.constLinearView();
        } else {
            const ptrdiff_t lo = m2.nlo();
            const ptrdiff_t hi = m2.nhi();
            for(ptrdiff_t k = -lo; k <= hi; ++k) m2.diag(k) = m1.diag(k);
        }
    }

    template <typename T>
    inline void DoCopy(
        const GenBandMatrix<std::complex<T> >& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T1, typename T2>
    inline void DoCopy1(
        const GenBandMatrix<T1>& m1, BandMatrixView<T2> m2)
    {
        TMVAssert(isReal(T1()) || isComplex(T2()));
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.nlo() == m1.nlo());
        TMVAssert(m2.nhi() == m1.nhi());
        if (m2.colsize() > 0 && m2.rowsize() > 0) {
            if (SameStorage(m1,m2)) {
                if (m2.isSameAs(m1)) {} // Do Nothing
                else if (m2.nlo() == m2.nhi() &&
                         m2.transpose().isSameAs(m1))
                    m2.transposeSelf();
                else if (m1.isconj() != m2.isconj() &&
                         m2.conjugate().isSameAs(m1))
                    m2.conjugateSelf();
                else if (m1.isrm())
                    DoCopy1(BandMatrix<T1,RowMajor>(m1),m2);
                else if (m1.iscm())
                    DoCopy1(BandMatrix<T1,ColMajor>(m1),m2);
                else
                    DoCopy1(BandMatrix<T1,DiagMajor>(m1),m2);
            } else if (m2.isconj()) {
                DoCopy(m1.conjugate(),m2.conjugate());
            } else {
                DoCopy(m1,m2);
            }
        }
    }

    template <typename T1, typename T2>
    inline void Copy(
        const GenBandMatrix<T1>& m1, BandMatrixView<T2> m2)
    {
        TMVAssert(isReal(T1()) || isComplex(T2()));
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.nlo() >= m1.nlo());
        TMVAssert(m2.nhi() >= m1.nhi());

        DoCopy1(m1,m2.subBandMatrix(
                0,m2.colsize(),0,m2.rowsize(),m1.nlo(),m1.nhi()));
        if (m2.nhi() > m1.nhi())
            m2.diagRange(m1.nhi()+1,m2.nhi()+1).setZero();
        if (m2.nlo() > m1.nlo())
            m2.diagRange(-m2.nlo(),-m1.nlo()).setZero();
    }


    //
    // Swap Matrices
    //

    template <typename T>
    inline void Swap(
        BandMatrixView<T> m1, BandMatrixView<T> m2);

    template <typename T, int A>
    inline void Swap(
        BandMatrixView<T> m1, BandMatrix<T,A>& m2)
    { Swap(m1,m2.view()); }

    template <typename T, int A>
    inline void Swap(
        BandMatrix<T,A>& m1, BandMatrixView<T> m2)
    { Swap(m1.view(),m2); }

    template <typename T, int A1, int A2>
    inline void Swap(BandMatrix<T,A1>& m1, BandMatrix<T,A2>& m2)
    { Swap(m1.view(),m2.view()); }


    //
    // Views:
    //

    template <typename T>
    inline ConstBandMatrixView<T> Transpose(const GenBandMatrix<T>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> Transpose(const ConstBandMatrixView<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> Transpose(
        const BandMatrix<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline BandMatrixView<T,A> Transpose(BandMatrixView<T,A> m)
    { return m.transpose(); }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> Transpose(BandMatrix<T,A>& m)
    { return m.transpose(); }

    template <typename T>
    inline ConstBandMatrixView<T> Conjugate(const GenBandMatrix<T>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> Conjugate(const ConstBandMatrixView<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> Conjugate(
        const BandMatrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline BandMatrixView<T,A> Conjugate(BandMatrixView<T,A> m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> Conjugate(BandMatrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T>
    inline ConstBandMatrixView<T> Adjoint(const GenBandMatrix<T>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A> Adjoint(const ConstBandMatrixView<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstBandMatrixView<T,A&FortranStyle> Adjoint(
        const BandMatrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline BandMatrixView<T,A> Adjoint(BandMatrixView<T,A> m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline BandMatrixView<T,A&FortranStyle> Adjoint(BandMatrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T>
    inline QuotXB<T,T> Inverse(const GenBandMatrix<T>& m)
    { return m.inverse(); }



    //
    // BandMatrix ==, != BandMatrix
    //

    template <typename T1, typename T2>
    bool operator==(
        const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2);

    template <typename T1, typename T2>
    inline bool operator!=(
        const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <typename T1, typename T2>
    bool operator==(
        const GenBandMatrix<T1>& m1, const GenMatrix<T2>& m2);

    template <typename T1, typename T2>
    inline bool operator==(
        const GenMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    { return m2 == m1; }

    template <typename T1, typename T2>
    inline bool operator!=(
        const GenBandMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <typename T1, typename T2>
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <typename T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, BandMatrixView<T> m)
    { m.read(reader); return reader.getis(); }

    template <typename T, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, BandMatrix<T,A>& m)
    { m.read(reader); return reader.getis(); }

    template <typename T>
    inline std::istream& operator>>(
        std::istream& is, BandMatrixView<T> m)
    { return is >> IOStyle() >> m; }

    template <typename T, int A>
    inline std::istream& operator>>(std::istream& is, BandMatrix<T,A>& m)
    { return is >> IOStyle() >> m; }


} // namespace tmv

#endif
