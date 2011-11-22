///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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
// along with this program in the file LICENSE.                              //
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
//    m.newView()
//    m.newTranspose()
//    m.newConjugate()
//    m.newAdjoint()
//    m.newInverse()
//    m.newCopy()
//
//
// I/O: 
//
//    os << m 
//        Writes m to ostream os in the usual Matrix format
//
//    m.writeCompact(os)
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

    template <class T>
    class GenBandMatrix : 
        virtual public AssignableToBandMatrix<T>,
        virtual public AssignableToDiagMatrix<T>,
        virtual public AssignableToUpperTriMatrix<T>,
        virtual public AssignableToLowerTriMatrix<T>,
        public BaseMatrix<T>,
        private DivHelper<T>
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
        inline size_t size() const 
        { TMVAssert(colsize() == rowsize()); return colsize(); }
        inline DiagType dt() const { return NonUnitDiag; }

        inline T operator()(int i,int j) const 
        { 
            TMVAssert(i>=0 && i<int(colsize()));
            TMVAssert(j>=0 && j<int(rowsize()));
            return okij(i,j) ? cref(i,j) : T(0);
        }

        inline const_vec_type row(int i, int j1, int j2) const
        { 
            TMVAssert(i>=0 && i<int(colsize()));
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize()));
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return const_vec_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct());
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<int(rowsize()));
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return const_vec_type(
                cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct());
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                cptr(),TMV_MIN(colsize(),rowsize()),diagstep(),ct());
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            if (i >= 0) {
                const size_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    cptr()+i*stepj(),diagsize,diagstep(),ct());
            } else {
                const size_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    cptr()-i*stepi(),diagsize,diagstep(),ct());
            }
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            TMVAssert(j1>=0 && j1-j2<=0);
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(int(rowsize())-i,int(colsize())));
                return const_vec_type(
                    cptr()+i*stepj()+j1*diagstep(),j2-j1,diagstep(),ct());
            } else {
                TMVAssert(j2<=TMV_MIN(int(colsize())+i,int(rowsize())));
                return const_vec_type(
                    cptr()-i*stepi()+j1*diagstep(),j2-j1,diagstep(),ct());
            }
        }

        template <class T2> 
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

        template <class T2>
        TMV_DEPRECATED(bool SameAs(const BaseMatrix<T2>& m2) const);
        TMV_DEPRECATED(bool SameAs(const type& m2) const)
        { return isSameAs(m2); }

        inline void assignToM(const MatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            assignToB(BandMatrixView<RT>(m2,nlo(),nhi()));
            if (int(rowsize()) > nhi()+1)
                BandMatrixView<RT>(
                    m2.colRange(nhi()+1,rowsize()),0,rowsize()-nhi()-2).setZero();
            if (int(colsize()) > nlo()+1)
                BandMatrixView<RT>(
                    m2.rowRange(nlo()+1,colsize()),colsize()-nlo()-2,0).setZero();
        }

        inline void assignToM(const MatrixView<CT>& m2) const
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            assignToB(BandMatrixView<CT>(m2,nlo(),nhi()));
            if (int(rowsize()) > nhi()+1)
                BandMatrixView<CT>(
                    m2.colRange(nhi()+1,rowsize()),0,rowsize()-nhi()-2).setZero();
            if (int(colsize()) > nlo()+1)
                BandMatrixView<CT>(
                    m2.rowRange(nlo()+1,colsize()),colsize()-nlo()-2,0).setZero();
        }

        inline void assignToB(const BandMatrixView<RT>& m2) const
        { 
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(m2.nhi() >= nhi());
            if (!isSameAs(m2)) Copy(*this,m2); 
        }

        inline void assignToB(const BandMatrixView<CT>& m2) const
        { 
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(m2.nhi() >= nhi());
            if (!isSameAs(m2)) Copy(*this,m2); 
        }

        inline void assignToU(const UpperTriMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nlo() == 0); 
            assignToB(BandMatrixViewOf(m2)); 
        }

        inline void assignToU(const UpperTriMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nlo() == 0); 
            assignToB(BandMatrixViewOf(m2)); 
        }

        inline void assignToL(const LowerTriMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0); 
            assignToB(BandMatrixViewOf(m2)); 
        }

        inline void assignToL(const LowerTriMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0); 
            assignToB(BandMatrixViewOf(m2)); 
        }

        inline void assignToD(const DiagMatrixView<RT>& m2) const
        { 
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == colsize());
            TMVAssert(m2.size() == rowsize());
            TMVAssert(nhi() == 0 && nlo() == 0); 
            assignToB(BandMatrixViewOf(m2)); 
        }

        inline void assignToD(const DiagMatrixView<CT>& m2) const
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
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return const_rec_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(), stepj(), stor(), ct());
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const 
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            const StorageType newstor = 
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return const_rec_type(
                cptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,ct());
        }

        bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const;

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return const_vec_type(
                cptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),ct());
        }

        bool hasSubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const;

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(),stor(),ct());
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            const StorageType newstor = 
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            const int newstepi = stepi()*istep;
            const int newstepj = stepj()*jstep;
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
                newstepi, newstepj, newstepi+newstepj, newstor, ct());
        }

        inline const_view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));

            const int j1 = i1 > nlo() ? i1-nlo() : 0;
            const int j2 = TMV_MIN(i2 + nhi(),int(rowsize()));
            const int newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const int newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const size_t newlin = (ls() && isrm()) ? 1 : 0;
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin);
        }

        inline const_view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize()));

            const int i1 = j1 > nhi() ? j1-nhi() : 0;
            const int i2 = TMV_MIN(j2 + nlo(),int(colsize()));
            const int newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const int newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            const size_t newlin = (ls() && iscm()) ? 1 : 0;
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin);
        }

        inline const_view_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);

            const int i1 = k2 <= 0 ? -k2+1 : 0;
            const int i2 = TMV_MIN(int(rowsize())-k1,int(colsize()));
            const int j1 = k1 <= 0 ? 0 : k1;
            const int j2 = TMV_MIN(int(rowsize()),int(colsize())+k2-1);
            const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct());
        }

        inline const_view_type upperBand() const
        {
            return const_view_type(
                cptr(),TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_view_type lowerBand() const
        {
            return const_view_type(
                cptr(),TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_view_type upperBandOff() const
        {
            return const_view_type(
                cptr()+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_view_type lowerBandOff() const
        {
            return const_view_type(
                cptr()+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                isReal(T()) ? stor() : NoMajor, NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NoMajor, NonConj);
        }

        TMV_DEPRECATED(const_rec_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(const_rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(const_vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(const_view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_DEPRECATED(const_view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_DEPRECATED(const_view_type Cols(int j1, int j2) const)
        { return colRange(j1,j2); }
        TMV_DEPRECATED(const_view_type Rows(int i1, int i2) const)
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(const_view_type Diags(int k1, int k2) const)
        { return diagRange(k1,k2); }
        TMV_DEPRECATED(const_view_type UpperBand(
                DiagType dt=NonUnitDiag) const)
        { return upperBand(dt); }
        TMV_DEPRECATED(const_view_type LowerBand(
                DiagType dt=NonUnitDiag) const)
        { return lowerBand(dt); }
        TMV_DEPRECATED(const_view_type UpperBandOff(
                DiagType dt=NonUnitDiag) const)
        { return upperBandOff(dt); }
        TMV_DEPRECATED(const_view_type LowerBandOff(
                DiagType dt=NonUnitDiag) const)
        { return lowerBandOff(dt); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }


        //
        // Views
        //

        inline const_view_type view() const
        { 
            return const_view_type(
                cptr(),colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),stor(),ct(),
                isdm()?0:ls());
        }

        inline const_view_type transpose() const
        { 
            return const_view_type(
                cptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(stor()),ct(),
                isdm()?0:ls());
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                cptr(),colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),stor(),TMV_ConjOf(T,ct()),
                isdm()?0:ls());
        }

        inline const_view_type adjoint() const
        { 
            return const_view_type(
                cptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
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
            TMVAssert(!isdm());
            TMVAssert(ls() != 1 || (rowsize() == 1 && colsize() == 1));
            // (To assure that next assert has no effect.)
            TMVAssert(canLinearize());
#ifdef TMV_USE_VALGRIND
            std::vector<bool> ok(ls(),false);
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) 
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = true;
                }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                const_cast<T*>(cptr())[k] = T(-777);
            }
#endif
            return const_vec_type(cptr(),ls(),1,ct());
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<T*>(cptr()),colsize(),rowsize(),
                nlo(),nhi(),stepi(),stepj(),diagstep(),stor(),ct(),
                isdm()?0:ls()
                TMV_FIRSTLAST1(isdm() ? diag(-nlo()).begin().getP() : cptr(),
                               isdm() ? diag(nhi()).end().getP() :
                               ((diag().end()-1).getP()+1)));
        }

        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(const_vec_type ConstLinearView() const)
        { return constLinearView(); }
        TMV_DEPRECATED(nonconst_type NonConst() const)
        { return nonConst(); }

        //
        // Functions of Matrix
        //

        T det() const;

        RT logDet(T* sign=0) const;

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

        RT doNorm2() const;
        RT norm2() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::norm2(); 
            TMV_Warning("Calling BandMatrix::norm2 without previously calling "
                        "divideUsing(SV)");
            return doNorm2();
        }

        inline RT normInf() const
        { return transpose().norm1(); }

        RT maxAbsElement() const;
        RT maxAbs2Element() const;

        QuotXB<T,T> QInverse() const;
        inline QuotXB<T,T> inverse() const
        { return QInverse(); }

        inline void makeInverse(const MatrixView<T>& minv) const
        { DivHelper<T>::makeInverse(minv); }

        template <class T1> 
        inline void makeInverse(const MatrixView<T1>& minv) const
        { DivHelper<T>::makeInverse(minv); }

        template <class T1, StorageType S, IndexStyle I> 
        inline void makeInverse(Matrix<T1,S,I>& minv) const
        { DivHelper<T>::makeInverse(minv); }

        inline void makeInverseATA(const MatrixView<T>& ata) const
        { DivHelper<T>::makeInverseATA(ata); }

        template <StorageType S, IndexStyle I> 
        inline void makeInverseATA(Matrix<T,S,I>& ata) const
        { DivHelper<T>::makeInverseATA(ata); }

        inline bool isSingular() const
        { return DivHelper<T>::isSingular(); }

        RT doCondition() const;
        inline RT condition() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::condition();
            TMV_Warning("Calling BandMatrix::condition without previously "
                        "calling divideUsing(SV)");
            return doCondition();
        }

        auto_ptr<BaseMatrix<T> > newCopy() const;
        auto_ptr<BaseMatrix<T> > newView() const;
        auto_ptr<BaseMatrix<T> > newTranspose() const;
        auto_ptr<BaseMatrix<T> > newConjugate() const;
        auto_ptr<BaseMatrix<T> > newAdjoint() const;
        auto_ptr<BaseMatrix<T> > newInverse() const;

        typedef QuotXB<T,T> MyQuotXB;
        TMV_DEPRECATED(MyQuotXB Inverse() const)
        { return inverse(); }
        template <class T1>
        TMV_DEPRECATED(void Inverse(const MatrixView<T1>& minv) const);
        template <class T1, StorageType S, IndexStyle I>
        TMV_DEPRECATED(void Inverse(Matrix<T1,S,I>& minv) const);
        template <StorageType S, IndexStyle I>
        TMV_DEPRECATED(void InverseATA(Matrix<T,S,I>& ata) const);

        //
        // Division Control
        //

        using DivHelper<T>::divideInPlace;
        using DivHelper<T>::saveDiv;
        using DivHelper<T>::setDiv;
        using DivHelper<T>::unsetDiv;
        using DivHelper<T>::resetDiv;
        using DivHelper<T>::divIsSet;
        using DivHelper<T>::checkDecomp;

        inline void divideUsing(DivType dt) const
        {
            TMVAssert(dt == LU || dt == QR || dt == SV);
            DivHelper<T>::divideUsing(dt);
        }

        inline const BandLUDiv<T>& lud() const
        {
            divideUsing(LU);
            setDiv();
            TMVAssert(getDiv());
            TMVAssert(divIsLUDiv());
            return static_cast<const BandLUDiv<T>&>(*getDiv());
        }

        inline const BandQRDiv<T>& qrd() const
        {
            divideUsing(QR);
            setDiv();
            TMVAssert(getDiv());
            TMVAssert(divIsQRDiv());
            return static_cast<const BandQRDiv<T>&>(*getDiv());
        }

        inline const BandSVDiv<T>& svd() const
        {
            divideUsing(SV);
            setDiv();
            TMVAssert(getDiv());
            TMVAssert(divIsSVDiv());
            return static_cast<const BandSVDiv<T>&>(*getDiv());
        }

        template <class T1> 
        inline void LDivEq(const VectorView<T1>& v) const 
        { DivHelper<T>::LDivEq(v); }
        template <class T1> 
        inline void LDivEq(const MatrixView<T1>& m) const 
        { DivHelper<T>::LDivEq(m); }
        template <class T1> 
        inline void RDivEq(const VectorView<T1>& v) const 
        { DivHelper<T>::RDivEq(v); }
        template <class T1> 
        inline void RDivEq(const MatrixView<T1>& m) const 
        { DivHelper<T>::RDivEq(m); }
        template <class T1, class T0> 
        inline void LDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { DivHelper<T>::LDiv(v1,v0); }
        template <class T1, class T0> 
        inline void LDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { DivHelper<T>::LDiv(m1,m0); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { DivHelper<T>::RDiv(v1,v0); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { DivHelper<T>::RDiv(m1,m0); }

        using DivHelper<T>::DivideInPlace;
        using DivHelper<T>::SaveDiv;
        using DivHelper<T>::SetDiv;
        using DivHelper<T>::UnSetDiv;
        using DivHelper<T>::ReSetDiv;
        using DivHelper<T>::DivIsSet;
        using DivHelper<T>::CheckDecomp;

        TMV_DEPRECATED(void DivideUsing(DivType dt) const)
        { divideUsing(dt); }
        TMV_DEPRECATED(const BandLUDiv<T>& LUD() const)
        { return lud(); }
        TMV_DEPRECATED(const BandQRDiv<T>& QRD() const)
        { return qrd(); }
        TMV_DEPRECATED(const BandSVDiv<T>& SVD() const)
        { return svd(); }


        //
        // I/O
        //

        void writeCompact(std::ostream& fout) const;
        void write(std::ostream& fout) const;
        void writeCompact(std::ostream& fout, RT thresh) const;
        void write(std::ostream& fout, RT thresh) const;

        TMV_DEPRECATED(void WriteCompact(std::ostream& fout) const)
        { writeCompact(fout); }
        TMV_DEPRECATED(void WriteCompact(std::ostream& fout, RT thresh) const)
        { writeCompact(fout,thresh); }

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
            TMVAssert(isComplex(T()) || ct()==NonConj);
            return isComplex(T()) && ct()==Conj;
        }
        virtual StorageType stor() const = 0;
        virtual ConjType ct() const = 0;

        virtual bool canLinearize() const = 0;
        TMV_DEPRECATED(bool CanLinearize() const)
        { return canLinearize(); }

        virtual T cref(int i, int j) const;

        inline int rowstart(int i) const { return TMV_MAX(0,i-nlo()); }
        inline int rowend(int i) const 
        { return TMV_MIN(int(rowsize()),i+nhi()+1); }
        inline int colstart(int j) const { return TMV_MAX(0,j-nhi()); }
        inline int colend(int j) const 
        { return TMV_MIN(int(colsize()),j+nlo()+1); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        {
            return const_rowmajor_iterator(
                this,rowstart(TMV_MIN(colsize(),rowsize()+nlo())),0); 
        }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        {
            return const_colmajor_iterator(
                this,colstart(TMV_MIN(rowsize(),colsize()+nhi())),0); 
        }

        inline const_diagmajor_iterator diagmajor_begin() const
        { return const_diagmajor_iterator(this,nlo(),0); }
        inline const_diagmajor_iterator diagmajor_end() const
        { return const_diagmajor_iterator(this,0,nhi()+1); }


    protected :

        inline bool okij(int i, int j) const
        { return (j+nlo() >= i && i+nhi() >= j); }

        using DivHelper<T>::getDiv;

        void newDivider() const;
        inline const BaseMatrix<T>& getMatrix() const { return *this; }

    private :

        type& operator=(const type&);

        bool divIsLUDiv() const;
        bool divIsQRDiv() const;
        bool divIsSVDiv() const;

    }; // GenBandMatrix

    template <class T>
    template <class T2>
    inline bool GenBandMatrix<T>::SameAs(const BaseMatrix<T2>& m2) const
    { return isSameAs(m2); }

    template <class T>
    template <class T1>
    inline void GenBandMatrix<T>::Inverse(const MatrixView<T1>& minv) const
    { makeInverse(minv); }

    template <class T>
    template <class T1, StorageType S, IndexStyle I>
    inline void GenBandMatrix<T>::Inverse(Matrix<T1,S,I>& minv) const
    { makeInverse(minv); }

    template <class T>
    template <StorageType S, IndexStyle I>
    inline void GenBandMatrix<T>::InverseATA(Matrix<T,S,I>& ata) const
    { makeInverseATA(ata); }



    template <class T, IndexStyle I> 
    class ConstBandMatrixView : public GenBandMatrix<T>
    {
    public :

        typedef GenBandMatrix<T> base;
        typedef ConstBandMatrixView<T,I> type;

        inline ConstBandMatrixView(const type& rhs) :
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
            itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
            itsstor(rhs.itsstor), itsct(rhs.itsct), linsize(rhs.linsize) 
        { TMVAssert(!(isdm() && linsize != 0)); }

        inline ConstBandMatrixView(const base& rhs) :
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
            TMVAssert(itsstor==RowMajor ? itssj==1 :
                      itsstor==ColMajor ? itssi==1 :
                      itsstor==DiagMajor ? itssd==1 : true); 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
            TMVAssert(!(isdm() && linsize!=0));
            TMVAssert(linsize==0 || linsize==1 || 
                      linsize==BandStorageLength(
                          itsstor,itscs,itsrs,itsnlo,itsnhi));
        }

        // These two work slightly differently than the BandMatrixViewOf
        // commands when the rhs matrix is not square.  
        // These two constructors copy the size of rhs viewing only 
        // the relevant rows.
        // In contrast, BandMatrixViewOf shrinks colsize of rowsize down to
        // only the rows and columns which include the bands.
        //   e.g. if rhs is 10 x 8, then:
        //   BandMatrixView(rhs,0,2) will have cs = rs = 8
        //   BandMatrixViewOf(rhs,0,2) will have cs = 10, rs = 8
        inline ConstBandMatrixView(const base& rhs, int lo, int hi) : 
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
            itsstor(rhs.stor()), itsct(rhs.ct()), linsize(0)
        { 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        inline ConstBandMatrixView(const GenMatrix<T>& rhs, int lo, int hi) : 
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(itssi+itssj),
            itsstor(rhs.stor()), itsct(rhs.ct()), linsize(0)
        { 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        virtual inline ~ConstBandMatrixView() 
        {
#ifdef TMVDEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

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
        using base::isdm;

        bool canLinearize() const;

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

        type& operator=(const type&);

    }; // ConstBandMatrixView

    template <class T> 
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
            const T* _m, size_t _cs, size_t _rs, 
            int _lo, int _hi, int _si, int _sj, int _sd, 
            StorageType instor, ConjType inct, size_t ls=0) : 
            c_type(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd, instor,inct,ls) {}

        inline ConstBandMatrixView(const base& rhs, int lo, int hi) : 
            c_type(rhs,lo,hi) {}

        inline ConstBandMatrixView(const GenMatrix<T>& rhs, int lo, int hi) : 
            c_type(rhs,lo,hi) {}

        virtual inline ~ConstBandMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(int i,int j) const 
        { 
            TMVAssert(i>0 && i<=int(colsize()));
            TMVAssert(j>0 && j<=int(rowsize()));
            return okij(i-1,j-1) ? cref(i-1,j-1) : T(0);
        }

        inline const_vec_type row(int i, int j1, int j2) const
        { 
            TMVAssert(i>0 && i<=int(colsize()));
            TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize()));
            TMVAssert(okij(i-1,j1-1));
            TMVAssert(okij(i-1,j2-1));
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=int(rowsize()));
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize()));
            TMVAssert(okij(i1-1,j-1));
            TMVAssert(okij(i2-1,j-1));
            return base::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag() const
        { return base::diag(); }

        inline const_vec_type diag(int i) const
        { return base::diag(i); }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(j1>0);
            return base::diag(i,j1-1,j2);
        }


        //
        // subBandMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::subMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const 
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::subMatrix(i1-1,i2-1+istep,j1-1,j2-1+jstep,
                                   istep,jstep);
        }

        bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const;

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return base::subVector(i-1,j-1,istep,jstep,size);
        }

        bool hasSubBandMatrix(
            int i1, int i2, int j1, int j2,
            int newnlo, int newnhi, int istep, int jstep) const;

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return base::subBandMatrix(i1-1,i2,j1-1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1);
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            return base::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
        }

        inline const_view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize()));
            return base::rowRange(i1-1,i2);
        }

        inline const_view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize()));
            return base::colRange(j1-1,j2);
        }

        inline const_view_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<=k2 && k2<=nhi());
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

        TMV_DEPRECATED(const_rec_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(const_rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(const_vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(const_view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_DEPRECATED(const_view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_DEPRECATED(const_view_type Cols(int j1, int j2) const)
        { return colRange(j1,j2); }
        TMV_DEPRECATED(const_view_type Rows(int i1, int i2) const)
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(const_view_type Diags(int k1, int k2) const)
        { return diagRange(k1,k2); }
        TMV_DEPRECATED(const_view_type UpperBand(
                DiagType dt=NonUnitDiag) const)
        { return upperBand(dt); }
        TMV_DEPRECATED(const_view_type LowerBand(
                DiagType dt=NonUnitDiag) const)
        { return lowerBand(dt); }
        TMV_DEPRECATED(const_view_type UpperBandOff(
                DiagType dt=NonUnitDiag) const)
        { return upperBandOff(dt); }
        TMV_DEPRECATED(const_view_type LowerBandOff(
                DiagType dt=NonUnitDiag) const)
        { return lowerBandOff(dt); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }


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

        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(nonconst_type NonConst() const)
        { return nonConst(); }

        using c_type::colsize;
        using c_type::rowsize;
        using c_type::nlo;
        using c_type::nhi;

        using base::cref;

    protected :

        using base::okij;

    private :

        type& operator=(const type&);

    }; // FortranStyle ConstBandMatrixView

    template <class T, IndexStyle I> 
    class BandMatrixView : public GenBandMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenBandMatrix<T> base;
        typedef BandMatrixView<T,I> type;
        typedef BandMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef BandMatrixView<RT,I> realpart_type;
        typedef TMV_RefType(T) reference;
        typedef RMIt<const type> rowmajor_iterator;
        typedef CMIt<const type> colmajor_iterator;
        typedef DMIt<const type> diagmajor_iterator;

        //
        // Constructors
        //

        inline BandMatrixView(const type& rhs) : 
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itsnlo(rhs.itsnlo), itsnhi(rhs.itsnhi),
            itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
            itsstor(rhs.stor()), itsct(rhs.ct()), linsize(rhs.ls())
            TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { 
            TMVAssert(itsstor==RowMajor ? itssj==1 :
                      itsstor==ColMajor ? itssi==1 :
                      itsstor==DiagMajor ? itssd==1 : true); 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
            TMVAssert(!(isdm() && linsize != 0));
            TMVAssert(itsstor==DiagMajor ? linsize == 0 : true);
            TMVAssert(linsize==0 || linsize==1 || 
                      linsize==BandStorageLength(
                          itsstor,itscs,itsrs,itsnlo,itsnhi));
        }

        inline BandMatrixView(
            T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
            int _si, int _sj, int _sd, StorageType _stor, ConjType _ct,
            size_t _ls TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itscs(_cs), itsrs(_rs),
            itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
            itsstor(_stor), itsct(_ct), linsize(_ls)
            TMV_DEFFIRSTLAST(_first,_last)
        { 
            TMVAssert(itsstor==RowMajor ? itssj==1 :
                      itsstor==ColMajor ? itssi==1 :
                      itsstor==DiagMajor ? itssd==1 : true); 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
            TMVAssert(!(isdm() && linsize != 0));
            TMVAssert(itsstor==DiagMajor ? linsize == 0 : true);
            TMVAssert(linsize==0 || linsize==1 || 
                      linsize==BandStorageLength(
                          itsstor,itscs,itsrs,itsnlo,itsnhi));
        }

        inline BandMatrixView(
            T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
            int _si, int _sj, int _sd, StorageType _stor, ConjType _ct
            TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itscs(_cs), itsrs(_rs),
            itsnlo(_lo), itsnhi(_hi), itssi(_si), itssj(_sj), itssd(_sd),
            itsstor(_stor), itsct(_ct), linsize(0)
            TMV_DEFFIRSTLAST(_first,_last)
        { 
            TMVAssert(itsstor==RowMajor ? itssj==1 :
                      itsstor==ColMajor ? itssi==1 :
                      itsstor==DiagMajor ? itssd==1 : true); 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
            TMVAssert(!(isdm() && linsize != 0));
            TMVAssert(itsstor==DiagMajor ? linsize == 0 : true);
            TMVAssert(linsize == 0);
        }

        inline BandMatrixView(const BandMatrixView<T>& rhs, int lo, int hi) : 
            itsm(rhs.ptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
            itsstor(rhs.stor()), itsct(rhs.ct()), linsize(0)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { 
            TMVAssert(itsstor==RowMajor ? itssj==1 :
                      itsstor==ColMajor ? itssi==1 :
                      true); 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        inline BandMatrixView(const MatrixView<T>& rhs, int lo, int hi) : 
            itsm(rhs.ptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itsnlo(lo), itsnhi(hi),
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(itssi+itssj),
            itsstor(rhs.stor()), itsct(rhs.ct()), linsize(0)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { 
            TMVAssert(itsstor==RowMajor ? itssj==1 :
                      itsstor==ColMajor ? itssi==1 :
                      true); 
            TMVAssert(rowsize() == 0 || nhi() < int(rowsize()));
            TMVAssert(colsize() == 0 || nlo() < int(colsize()));
            TMVAssert(nhi() >= 0);
            TMVAssert(nlo() >= 0);
        }

        virtual inline ~BandMatrixView() 
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMVDEBUG
            const_cast<T*&>(itsm) = 0;
#endif
        }

        //
        // Op=
        //

        inline const type& operator=(const type& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this; 
        }

        inline const type& operator=(const type& m2)
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this; 
        }

        inline const type& operator=(const GenBandMatrix<RT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this; 
        }

        inline const type& operator=(const GenBandMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenBandMatrix<T2>& m2) const
        {
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            if (!this->isSameAs(m2)) Copy(m2,view());
            return *this; 
        }

        inline const type& operator=(const T& x) const 
        {
            TMVAssert(colsize() == rowsize());
            return setToIdentity(x); 
        }

        inline const type& operator=(
            const AssignableToBandMatrix<RT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        inline const type& operator=(
            const AssignableToBandMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
            m2.assignToB(*this);
            return *this;
        }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignToD(DiagMatrixViewOf(this->diag()));
            if (this->nhi() > 0) upperBandOff().setZero();
            if (this->nlo() > 0) lowerBandOff().setZero();
            return *this;
        }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignToD(DiagMatrixViewOf(this->diag()));
            if (this->nhi() > 0) upperBandOff().setZero();
            if (this->nlo() > 0) lowerBandOff().setZero();
            return *this;
        }

        inline const type& operator=(const GenUpperTriMatrix<RT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == int(rowsize())-1);
            m2.assignToU(
                UpperTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
            if (nlo() > 0) diagRange(-nlo(),0).setZero();
            return *this;
        }

        inline const type& operator=(
            const GenUpperTriMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == int(rowsize())-1);
            m2.assignToU(
                UpperTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
            if (nlo() > 0) diagRange(-nlo(),0).setZero();
            return *this;
        }

        inline const type& operator=(
            const GenLowerTriMatrix<RT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == int(rowsize())-1);
            m2.assignToL(
                LowerTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
            if (nhi() > 0) diagRange(1,nhi()+1).setZero();
            return *this;
        }

        inline const type& operator=(
            const GenLowerTriMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == int(rowsize())-1);
            m2.assignToL(
                LowerTriMatrixView<T>(
                    ptr(),colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
            if (nhi() > 0) diagRange(1,nhi()+1).setZero();
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        {
            const int n = BandNumElements(colsize(),rowsize(),nlo(),nhi());
            return MyListAssigner(rowmajor_begin(),n,x); 
        }

        TMV_DEPRECATED(inline MyListAssigner operator=(ListInitClass) const)
        {
            const int n = BandNumElements(colsize(),rowsize(),nlo(),nhi());
            return MyListAssigner(rowmajor_begin(),n); 
        }


        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        { 
            TMVAssert(i>=0 && i<int(colsize()));
            TMVAssert(j>=0 && j<int(rowsize()));
            TMVAssert(okij(i,j));
            return ref(i,j); 
        }

        inline vec_type row(int i, int j1, int j2) const
        { 
            TMVAssert(i>=0 && i<int(colsize()));
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize()));
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return vec_type(ptr()+i*stepi()+j1*stepj(),
                            j2-j1,stepj(),ct() TMV_FIRSTLAST);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<int(rowsize()));
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return vec_type(ptr()+i1*stepi()+j*stepj(),
                            i2-i1,stepi(),ct() TMV_FIRSTLAST );
        }

        inline vec_type diag() const
        {
            return vec_type(ptr(),TMV_MIN(colsize(),rowsize()),
                            diagstep(),ct() TMV_FIRSTLAST);
        }

        inline vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            if (i >= 0) {
                const int diaglen = TMV_MIN(int(rowsize())-i,int(colsize()));
                return vec_type(ptr()+i*stepj(),diaglen,diagstep(),ct() 
                                TMV_FIRSTLAST );
            } else {
                const int diaglen = TMV_MIN(int(colsize())+i,int(rowsize()));
                return vec_type(ptr()-i*stepi(),diaglen,diagstep(),ct() 
                                TMV_FIRSTLAST );
            }
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            TMVAssert(j1>=0 && j1-j2<=0);
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(int(rowsize())-i,int(colsize())));
                return vec_type(ptr()+i*stepj()+j1*diagstep(),
                                j2-j1, diagstep(),ct() TMV_FIRSTLAST );
            } else {
                TMVAssert(j2<=TMV_MIN(int(colsize())+i,int(rowsize())));
                return vec_type(ptr()-i*stepi()+j1*diagstep(),
                                j2-j1, diagstep(),ct() TMV_FIRSTLAST );
            }
        }

        //
        // Modifying Functions
        //

        const type& setZero() const;

        const type& setAllTo(const T& x) const;

        const type& addToAll(const T& x) const;

        const type& clip(RT thresh) const;

        void doTransposeSelf() const;
        inline const type& transposeSelf() const
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(nlo() == nhi());
            doTransposeSelf();
            return *this;
        }

        const type& conjugateSelf() const;

        inline const type& setToIdentity(const T& x=T(1)) const 
        {
            TMVAssert(colsize() == rowsize());
            setZero(); diag().setAllTo(x); return *this; 
        }

        TMV_DEPRECATED(const type& Zero() const)
        { return setZero(); }
        TMV_DEPRECATED(const type& SetAllTo(const T& x) const)
        { return setAllTo(x); }
        TMV_DEPRECATED(const type& Clip(RT thresh) const)
        { return clip(thresh); }
        TMV_DEPRECATED(const type& TransposeSelf() const)
        { return transposeSelf(); }
        TMV_DEPRECATED(const type& ConjugateSelf() const)
        { return conjugateSelf(); }
        TMV_DEPRECATED(const type& SetToIdentity(const T& x=T(1)) const)
        { return setToIdentity(x); }


        //
        // subBandMatrix
        //

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            return rec_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(), stepj(), stor(), ct() TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const 
        {
            const StorageType newstor = 
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return rec_type(
                ptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                newstor, ct() TMV_FIRSTLAST);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,size));
            return vec_type(
                ptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(), ct() TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(base::hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
                                             1,1));
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct() TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(base::hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
                                             istep,jstep));
            const StorageType newstor = 
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            const int newstepi = stepi()*istep;
            const int newstepj = stepj()*jstep;
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
                stepi()*istep, stepj()*jstep, newstepi+newstepj, newstor, 
                             ct() TMV_FIRSTLAST);
        }

        inline view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));

            const int j1 = i1 > nlo() ? i1-nlo() : 0;
            const int j2 = TMV_MIN(i2 + nhi(),int(rowsize()));
            const int newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const int newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const size_t newlin = (ls() && isrm()) ? 1 : 0;
            TMVAssert(base::hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
                                             1,1));
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize()));

            const int i1 = j1 > nhi() ? j1-nhi() : 0;
            const int i2 = TMV_MIN(j2 + nlo(),int(colsize()));
            const int newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const int newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            const size_t newlin = (ls() && iscm()) ? 1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);

            const int i1 = k2 <= 0 ? -k2+1 : 0;
            const int i2 = TMV_MIN(int(rowsize())-k1,int(colsize()));
            const int j1 = k1 <= 0 ? 0 : k1;
            const int j2 = TMV_MIN(int(rowsize()),int(colsize())+k2-1);
            const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct() TMV_FIRSTLAST);
        }

        inline view_type upperBand() const
        {
            return view_type(
                ptr(),TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBand() const
        {
            return view_type(
                ptr(),TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type upperBandOff() const
        {
            return view_type(
                ptr()+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBandOff() const
        {
            return view_type(
                ptr()+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline realpart_type realPart() const
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                isReal(T()) ? stor() : NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline view_type view() const
        { return *this; }

        inline view_type transpose() const
        { 
            return view_type(
                ptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(stor()),ct(),
                isdm()?0:ls() TMV_FIRSTLAST);
        }

        inline view_type conjugate() const
        { 
            return view_type(
                ptr(),colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),stor(),TMV_ConjOf(T,ct()),
                isdm()?0:ls() TMV_FIRSTLAST);
        }

        inline view_type adjoint() const
        { 
            return view_type(
                ptr(),rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                TMV_ConjOf(T,ct()),isdm()?0:ls() TMV_FIRSTLAST);
        }

        inline vec_type linearView() const
        {
            TMVAssert(!isdm());
            TMVAssert(ls() != 1 || (rowsize() == 1 && colsize() == 1));
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
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) 
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = true;
                }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return vec_type(ptr(),ls(),1,ct() TMV_FIRSTLAST );
        }

        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_DEPRECATED(view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_DEPRECATED(view_type Cols(int j1, int j2) const)
        { return colRange(j1,j2); }
        TMV_DEPRECATED(view_type Rows(int i1, int i2) const)
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(view_type Diags(int k1, int k2) const)
        { return diagRange(k1,k2); }
        TMV_DEPRECATED(view_type UpperBand(
                DiagType dt=NonUnitDiag) const)
        { return upperBand(dt); }
        TMV_DEPRECATED(view_type LowerBand(
                DiagType dt=NonUnitDiag) const)
        { return lowerBand(dt); }
        TMV_DEPRECATED(view_type UpperBandOff(
                DiagType dt=NonUnitDiag) const)
        { return upperBandOff(dt); }
        TMV_DEPRECATED(view_type LowerBandOff(
                DiagType dt=NonUnitDiag) const)
        { return lowerBandOff(dt); }
        TMV_DEPRECATED(realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(view_type View() const)
        { return view(); }
        TMV_DEPRECATED(view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(vec_type LinearView() const)
        { return linearView(); }




        //
        // I/O
        //

        void read(std::istream& fin) const;

        inline size_t colsize() const { return itscs; }
        inline size_t rowsize() const { return itsrs; }
        inline int nlo() const { return itsnlo; }
        inline int nhi() const { return itsnhi; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline int diagstep() const { return itssd; }
        using base::isrm;
        using base::iscm;
        using base::isdm;
        using base::isconj;
        inline size_t ls() const { return linsize; }
        inline StorageType stor() const { return itsstor; }
        inline ConjType ct() const { return itsct; }

        bool canLinearize() const;

        reference ref(int i, int j) const;

        inline rowmajor_iterator rowmajor_begin() const
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end() const
        {
            return rowmajor_iterator(
                this,rowstart(TMV_MIN(colsize(),rowsize()+nlo())),0); 
        }

        inline colmajor_iterator colmajor_begin() const
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end() const
        {
            return colmajor_iterator(
                this,colstart(TMV_MIN(rowsize(),colsize()+nhi())),0); 
        }

        inline diagmajor_iterator diagmajor_begin() const
        { return diagmajor_iterator(this,nlo(),0); }
        inline diagmajor_iterator diagmajor_end() const
        { return diagmajor_iterator(this,0,nhi()+1); }

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
        const T* _first;
        const T* _last;
    protected :
#endif

        using base::okij;

    }; // BandMatrixView

    template <class T> 
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

        //
        // Constructors
        //

        inline BandMatrixView(const type& rhs) : c_type(rhs) {}

        inline BandMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline BandMatrixView(
            T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
            int _si, int _sj, int _sd, StorageType instor, ConjType inct,
            size_t ls TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,instor,inct,ls
                   TMV_FIRSTLAST1(_first,_last) ) {}

        inline BandMatrixView(
            T* _m, size_t _cs, size_t _rs, int _lo, int _hi,
            int _si, int _sj, int _sd, StorageType instor, ConjType inct
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_cs,_rs,_lo,_hi,_si,_sj,_sd,instor,inct
                   TMV_FIRSTLAST1(_first,_last) ) {}

        inline BandMatrixView(const BandMatrixView<T>& rhs, int lo, int hi) : 
            c_type(rhs,lo,hi) {}

        inline BandMatrixView(const MatrixView<T>& rhs, int lo, int hi) : 
            c_type(rhs,lo,hi) {}

        virtual inline ~BandMatrixView() {}

        //
        // Op=
        //

        inline const type& operator=(const type& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const type& m2)
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const c_type& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const c_type& m2)
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenBandMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenBandMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenBandMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(
            const AssignableToBandMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(
            const AssignableToBandMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenUpperTriMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenUpperTriMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenLowerTriMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenLowerTriMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        typedef typename c_type::MyListAssigner MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return c_type::operator<<(x); }

        TMV_DEPRECATED(inline MyListAssigner operator=(ListInitClass li) const)
        { return c_type::operator=(li); }

        //
        // Access
        //

        inline TMV_RefType(T) operator()(int i,int j) const 
        { 
            TMVAssert(i>0 && i<=int(colsize()));
            TMVAssert(j>0 && j<=int(rowsize()));
            TMVAssert(okij(i-1,j-1));
            return ref(i-1,j-1); 
        }

        inline vec_type row(int i, int j1, int j2) const
        { 
            TMVAssert(i>0 && i<=int(colsize()));
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= int(rowsize()));
            TMVAssert(okij(i-1,j1-1));
            TMVAssert(okij(i-1,j2-1));
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=int(rowsize()));
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= int(colsize()));
            TMVAssert(okij(i1-1,j-1));
            TMVAssert(okij(i2-1,j-1));
            return c_type::col(j-1,i1-1,i2);
        }

        inline vec_type diag() const
        { return c_type::diag(); }

        inline vec_type diag(int i) const
        { return c_type::diag(i); }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(j1>0);
            return c_type::diag(i,j1-1,j2); 
        }


        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { c_type::setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        { c_type::setAllTo(x); return *this; }

        inline const type& addToAll(const T& x) const
        { c_type::addToAll(x); return *this; }

        inline const type& clip(RT thresh) const
        { c_type::clip(thresh); return *this; }

        inline const type& transposeSelf() const
        { c_type::transposeSelf(); return *this; }

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const 
        { c_type::setToIdentity(x); return *this; }

        TMV_DEPRECATED(const type& Zero() const)
        { return setZero(); }
        TMV_DEPRECATED(const type& SetAllTo(const T& x) const)
        { return setAllTo(x); }
        TMV_DEPRECATED(const type& Clip(RT thresh) const)
        { return clip(thresh); }
        TMV_DEPRECATED(const type& TransposeSelf() const)
        { return transposeSelf(); }
        TMV_DEPRECATED(const type& ConjugateSelf() const)
        { return conjugateSelf(); }
        TMV_DEPRECATED(const type& SetToIdentity(const T& x=T(1)) const)
        { return setToIdentity(x); }

        //
        // subBandMatrix
        //

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            return const_view_type(*this).hasSubVector(i,j,istep,jstep,size);
        }

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            return const_view_type(*this).hasSubMatrix(
                i1,i2,j1,j2,istep,jstep);
        }

        inline bool hasSubBandMatrix(
            int i1, int i2, int j1, int j2,
            int newnlo, int newnhi, int istep, int jstep) const
        {
            return const_view_type(*this).hasSubBandMatrix(
                i1,i2,j1,j2,newnlo,newnhi,istep,jstep);
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::subMatrix(i1-1,i2,j1-1,j2);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const 
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::subMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return c_type::subVector(i-1,j-1,istep,jstep,size);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return c_type::subBandMatrix(i1-1,i2,j1-1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1);
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            return c_type::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
        }

        inline view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize()));
            return c_type::rowRange(i1-1,i2);
        }

        inline view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize()));
            return c_type::colRange(j1-1,j2);
        }

        inline view_type diagRange(int k1, int k2) const
        { return c_type::diagRange(k1,k2+1); }

        inline view_type upperBand() const
        { return c_type::upperBand(); }

        inline view_type lowerBand() const
        { return c_type::lowerBand(); }

        inline view_type upperBandOff() const
        { return c_type::upperBandOff(); }

        inline view_type lowerBandOff() const
        { return c_type::lowerBandOff(); }

        inline realpart_type realPart() const
        { return c_type::realPart(); }

        inline realpart_type imagPart() const
        { return c_type::imagPart(); }

        inline view_type view() const
        { return c_type::view(); }

        inline view_type transpose() const
        { return c_type::transpose(); }

        inline view_type conjugate() const
        { return c_type::conjugate(); }

        inline view_type adjoint() const
        { return c_type::adjoint(); }

        inline vec_type linearView() const
        { return c_type::linearView(); }

        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_DEPRECATED(view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_DEPRECATED(view_type Cols(int j1, int j2) const)
        { return colRange(j1,j2); }
        TMV_DEPRECATED(view_type Rows(int i1, int i2) const)
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(view_type Diags(int k1, int k2) const)
        { return diagRange(k1,k2); }
        TMV_DEPRECATED(view_type UpperBand(
                DiagType dt=NonUnitDiag) const)
        { return upperBand(dt); }
        TMV_DEPRECATED(view_type LowerBand(
                DiagType dt=NonUnitDiag) const)
        { return lowerBand(dt); }
        TMV_DEPRECATED(view_type UpperBandOff(
                DiagType dt=NonUnitDiag) const)
        { return upperBandOff(dt); }
        TMV_DEPRECATED(view_type LowerBandOff(
                DiagType dt=NonUnitDiag) const)
        { return lowerBandOff(dt); }
        TMV_DEPRECATED(realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(view_type View() const)
        { return view(); }
        TMV_DEPRECATED(view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(vec_type LinearView() const)
        { return linearView(); }



        using c_type::colsize;
        using c_type::rowsize;
        using c_type::nlo;
        using c_type::nhi;

    protected :

        using base::okij;
        using c_type::ref;

    }; // FortranStyle BandMatrixView

    template <class T, StorageType S, IndexStyle I> 
    class BandMatrix : public GenBandMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenBandMatrix<T> base;
        typedef BandMatrix<T,S,I> type;
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
        linsize(BandStorageLength(S,cs,rs,lo,hi)), \
        itsm1(linsize), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
        itssi(S==RowMajor ? lo+hi : S==ColMajor ? 1 : \
              rs>= cs ? 1-int(cs) : -int(rs) ), \
        itssj(S==RowMajor ? 1 : S==ColMajor ? lo+hi : -itssi+1), \
        itsds(S==RowMajor ? itssi+1 : S==ColMajor ? itssj+1 : 1), \
        itsm(S==DiagMajor ? itsm1.get() - lo*itssi : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize)

#define NEW_SIZE2(ls,cs,rs,lo,hi) \
        linsize(ls), \
        itsm1(linsize), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi), \
        itssi(S==RowMajor ? lo+hi : S==ColMajor ? 1 : \
              rs>= cs ? 1-int(cs) : -int(rs) ), \
        itssj(S==RowMajor ? 1 : S==ColMajor ? lo+hi : -itssi+1), \
        itsds(S==RowMajor ? itssi+1 : S==ColMajor ? itssj+1 : 1), \
        itsm(S==DiagMajor ? itsm1.get() - lo*itssi : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize)

        inline BandMatrix(size_t cs, size_t rs, int lo, int hi) :
            NEW_SIZE(cs,rs,lo,hi) 
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < int(cs));
            TMVAssert(hi < int(rs));
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        inline BandMatrix(size_t cs, size_t rs, int lo, int hi, const T& x) :
            NEW_SIZE(cs,rs,lo,hi)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < int(cs));
            TMVAssert(hi < int(rs));
            setAllTo(x);
        }

        TMV_DEPRECATED(inline BandMatrix(
                size_t cs, size_t rs, int lo, int hi, const T* vv)) :
            NEW_SIZE(cs,rs,lo,hi)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < int(cs));
            TMVAssert(hi < int(rs));
            std::copy(vv,vv+linsize,itsm1.get());
        }

        TMV_DEPRECATED(inline BandMatrix(
            size_t cs, size_t rs, int lo, int hi, const std::vector<T>& vv)) :
            NEW_SIZE(cs,rs,lo,hi)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < int(cs));
            TMVAssert(hi < int(rs));
            TMVAssert(vv.size() == linsize);
            std::copy(vv.begin(),vv.end(),itsm1.get());
        }

        inline BandMatrix(const BandMatrix<T,S,I>& m2) :
            NEW_SIZE2(m2.ls(),m2.itscs,m2.itsrs,m2.itsnlo,m2.itsnhi)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            std::copy(m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
        }

        template <IndexStyle I2> 
        inline BandMatrix(const BandMatrix<T,S,I2>& m2) :
            NEW_SIZE2(m2.ls(),m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            std::copy(m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
        }

        inline BandMatrix(const GenBandMatrix<RT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            m2.assignToB(view());
        }

        inline BandMatrix(const GenBandMatrix<CT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            m2.assignToB(view());
        }

        template <class T2> 
        inline BandMatrix(const GenBandMatrix<T2>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            Copy(m2,view()); 
        }

        template <class T2> 
        inline BandMatrix(const GenBandMatrix<T2>& m2, int lo, int hi) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),lo,hi)
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo <= m2.nlo());
            TMVAssert(hi <= m2.nhi());
            Copy(ConstBandMatrixView<T2>(m2,lo,hi),view()); 
            if (I==CStyle) {
                if (lo > m2.nlo()) diagRange(-lo,-m2.nlo()).setZero();
                if (hi > m2.nhi()) diagRange(m2.nhi()+1,hi+1).setZero();
            } else {
                if (lo > m2.nlo()) diagRange(-lo,-m2.nlo()-1).setZero();
                if (hi > m2.nhi()) diagRange(m2.nhi()+1,hi).setZero();
            }
        }

        inline BandMatrix(const GenMatrix<T>& m2, int lo, int hi) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),lo,hi)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(hi >= 0);
            TMVAssert(lo < int(m2.colsize()));
            TMVAssert(hi < int(m2.rowsize()));
            Copy(ConstBandMatrixView<T>(m2,lo,hi),view());
        }

        inline BandMatrix(const GenUpperTriMatrix<T>& m2, int hi) :
            NEW_SIZE(m2.size(),m2.size(),0,hi)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(hi >= 0);
            TMVAssert(hi < int(m2.size()));
            Copy(BandMatrixViewOf(m2,hi),view());
        }

        inline BandMatrix(const GenLowerTriMatrix<T>& m2, int lo) :
            NEW_SIZE(m2.size(),m2.size(),lo,0)
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            TMVAssert(lo >= 0);
            TMVAssert(lo < int(m2.size()));
            Copy(BandMatrixViewOf(m2,lo),view());
        }

        inline BandMatrix(const AssignableToBandMatrix<RT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            m2.assignToB(view());
        }

        inline BandMatrix(const AssignableToBandMatrix<CT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize(),m2.nlo(),m2.nhi())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            m2.assignToB(view());
        }

        inline BandMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,0)
        { 
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        inline BandMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,0)
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        inline BandMatrix(const GenUpperTriMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,m2.size()-1)
        { 
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            setZero();
            m2.assignToU(
                UpperTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
        }

        inline BandMatrix(const GenUpperTriMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),0,m2.size()-1)
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            setZero();
            m2.assignToU(
                UpperTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
        }

        inline BandMatrix(const GenLowerTriMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),m2.size()-1,0)
        { 
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            setZero();
            m2.assignToL(
                LowerTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
        }

        inline BandMatrix(const GenLowerTriMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.size(),m2.size()-1,0)
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(S == RowMajor || S == ColMajor || S == DiagMajor);
            setZero();
            m2.assignToL(
                LowerTriMatrixView<T>(
                    itsm,colsize(),stepi(),stepj(),
                    NonUnitDiag,isdm()?NoMajor:stor(),ct() TMV_FIRSTLAST));
        }

#undef NEW_SIZE
#undef NEW_SIZE2

        virtual inline ~BandMatrix() 
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMVDEBUG
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

        template <IndexStyle I2>
        inline type& operator=(const BandMatrix<T,S,I2>& m2)
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

        template <class T2> 
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
            TMVAssert(nhi() == int(rowsize())-1);
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<CT>& m2)
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nhi() == int(rowsize())-1);
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<RT>& m2) 
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == int(rowsize())-1);
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<CT>& m2) 
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() == int(rowsize())-1);
            view() = m2;
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x) 
        {
            const int n = BandNumElements(colsize(),rowsize(),nlo(),nhi());
            return MyListAssigner(rowmajor_begin(),n,x); 
        }

        TMV_DEPRECATED(inline MyListAssigner operator=(ListInitClass))
        {
            const int n = BandNumElements(colsize(),rowsize(),nlo(),nhi());
            return MyListAssigner(rowmajor_begin(),n); 
        }

        //
        // Access
        //

        inline T operator()(int i,int j) const
        { 
            if (I==CStyle) {
                TMVAssert(i>=0 && i<int(colsize()));
                TMVAssert(j>=0 && j<int(rowsize()));
                return okij(i,j) ? cref(i,j) : T(0); 
            } else {
                TMVAssert(i>0 && i<=int(colsize()));
                TMVAssert(j>0 && j<=int(rowsize()));
                return okij(i-1,j-1) ? cref(i-1,j-1) : T(0); 
            }
        }

        inline T& operator()(int i,int j) 
        { 
            if (I==CStyle) {
                TMVAssert(i>=0 && i<int(colsize()));
                TMVAssert(j>=0 && j<int(rowsize()));
                TMVAssert(okij(i,j));
                return ref(i,j); 
            } else {
                TMVAssert(i>0 && i<=int(colsize()));
                TMVAssert(j>0 && j<=int(rowsize()));
                TMVAssert(okij(i-1,j-1));
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(int i, int j1, int j2) const
        { 
            if (I == FortranStyle) {
                TMVAssert(i>0 && i<=int(colsize())); --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize())); --j1;
            } else {
                TMVAssert(i>=0 && i<int(colsize()));
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize()));
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return const_vec_type(
                itsm+i*stepi()+j1*stepj(),j2-j1,stepj(),NonConj);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I == FortranStyle) {
                TMVAssert(j>0 && j<=int(rowsize())); --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize())); --i1;
            } else {
                TMVAssert(j>=0 && j<int(rowsize()));
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return const_vec_type(
                itsm+i1*stepi()+j*stepj(),i2-i1,stepi(),NonConj);
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                itsm,TMV_MIN(colsize(),rowsize()),diagstep(),NonConj);
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            if (i >= 0) {
                const size_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    itsm+i*stepj(),diagsize,diagstep(),NonConj);
            } else {
                const size_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    itsm-i*stepi(),diagsize,diagstep(),NonConj);
            }
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            if (I==FortranStyle) { 
                TMVAssert(j1>0 && j1-j2<=0); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0);
            }
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(int(rowsize())-i,int(colsize())));
                return const_vec_type(
                    itsm+i*stepj()+j1*diagstep(),j2-j1,diagstep(),NonConj);
            } else {
                TMVAssert(j2<=TMV_MIN(int(colsize())+i,int(rowsize())));
                return const_vec_type(
                    itsm-i*stepi()+j1*diagstep(),j2-j1,diagstep(),NonConj);
            }
        }

        inline vec_type row(int i, int j1, int j2)
        { 
            if (I == FortranStyle) {
                TMVAssert(i>0 && i<=int(colsize())); --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize())); --j1;
            } else {
                TMVAssert(i>=0 && i<int(colsize()));
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize()));
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            return vec_type(
                itsm+i*stepi()+j1*stepj(),j2-j1,stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I == FortranStyle) {
                TMVAssert(j>0 && j<=int(rowsize())); --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize())); --i1;
            } else {
                TMVAssert(j>=0 && j<int(rowsize()));
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            return vec_type(
                itsm+i1*stepi()+j*stepj(),i2-i1,stepi(),NonConj TMV_FIRSTLAST );
        }

        inline vec_type diag()
        {
            return vec_type(
                itsm,TMV_MIN(colsize(),rowsize()),diagstep(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i)
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            if (i >= 0) {
                const size_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return vec_type(
                    itsm+i*stepj(),diagsize,diagstep(),NonConj TMV_FIRSTLAST);
            } else {
                const size_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return vec_type(
                    itsm-i*stepi(),diagsize,diagstep(),NonConj TMV_FIRSTLAST);
            }
        }

        inline vec_type diag(int i, int j1, int j2)
        {
            TMVAssert(i>=-nlo() && i<=nhi());
            TMVAssert(i>=-int(colsize()) && i<=int(rowsize()));
            if (I==FortranStyle) { 
                TMVAssert(j1>0 && j1-j2<=0); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0);
            }
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(int(rowsize())-i,int(colsize())));
                return vec_type(
                    itsm+i*stepj()+j1*diagstep(),j2-j1,diagstep(),NonConj 
                    TMV_FIRSTLAST);
            } else {
                TMVAssert(j2<=TMV_MIN(int(colsize())+i,int(rowsize())));
                return vec_type(
                    itsm-i*stepi()+j1*diagstep(),j2-j1,diagstep(),NonConj 
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

        TMV_DEPRECATED(type& Zero())
        { return setZero(); }
        TMV_DEPRECATED(type& SetAllTo(const T& x))
        { return setAllTo(x); }
        TMV_DEPRECATED(type& Clip(RT thresh))
        { return clip(thresh); }
        TMV_DEPRECATED(type& TransposeSelf())
        { return transposeSelf(); }
        TMV_DEPRECATED(type& ConjugateSelf())
        { return conjugateSelf(); }
        TMV_DEPRECATED(type& SetToIdentity(const T& x=T(1)))
        { return setToIdentity(x); }

        //
        // subBandMatrix
        //

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I == FortranStyle) { --i1; --j1; }
            return const_rec_type(
                itsm+i1*stepi()+j1*stepj(), i2-i1, j2-j1,
                stepi(), stepj(), S, NonConj);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
            if (I == FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return const_rec_type(
                itsm+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(),jstep*stepj(), newstor, NonConj);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I == FortranStyle) { --i; --j; }
            return const_vec_type(
                itsm+i*stepi()+j*stepj(), size,
                istep*stepi()+jstep*stepj(), NonConj);
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I == FortranStyle) { --i1; --j1; }
            return const_view_type(
                itsm+i1*stepi()+j1*stepj(), i2-i1, j2-j1,
                newnlo, newnhi, stepi(), stepj(), diagstep(), S, NonConj);
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==CStyle?1:0));
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-(I==CStyle?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(view().hasSubBandMatrix(
                    i1,i2,j1,j2, newnlo,newnhi,istep,jstep));
            StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                S==ColMajor ? istep == 1 ? ColMajor : NoMajor : 
                istep == 1 && jstep == 1 ? DiagMajor : NoMajor;
            if (I == FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const int newstepi = stepi()*istep;
            const int newstepj = stepj()*jstep;
            return const_view_type(
                itsm+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
                newstepi, newstepj, newstepi+newstepj, newstor, NonConj);
        }

        inline const_view_type rowRange(int i1, int i2) const
        {
            if (I==FortranStyle) { 
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize())); --i1; 
            } else {
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));
            }

            const int j1 = i1 > nlo() ? i1-nlo() : 0;
            const int j2 = TMV_MIN(i2 + nhi(),int(rowsize()));
            const int newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const int newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const size_t newlin = isrm() ? 1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin);
        }

        inline const_view_type colRange(int j1, int j2) const
        {
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize())); --j1;
            } else { 
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize())); 
            }

            const int i1 = j1 > nhi() ? j1-nhi() : 0;
            const int i2 = TMV_MIN(j2 + nlo(),int(colsize()));
            const int newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const int newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            const size_t newlin = iscm() ? 1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin);
        }

        inline const_view_type diagRange(int k1, int k2) const
        {
            if (I==FortranStyle) {
                TMVAssert(k1>=-nlo() && k1<=k2 && k2<=nhi()); ++k2;
            } else {
                TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);
            }

            const int i1 = k2 <= 0 ? -k2+1 : 0;
            const int i2 = TMV_MIN(int(rowsize())-k1,int(colsize()));
            const int j1 = k1 <= 0 ? 0 : k1;
            const int j2 = TMV_MIN(int(rowsize()),int(colsize())+k2-1);
            const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return const_view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct());
        }

        inline const_view_type upperBand() const
        {
            return const_view_type(
                itsm,TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_view_type lowerBand() const
        {
            return const_view_type(
                itsm,TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_view_type upperBandOff() const
        {
            return const_view_type(
                itsm+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_view_type lowerBandOff() const
        {
            return const_view_type(
                itsm+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                isReal(T()) ? S : NoMajor, NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm)+1,
                colsize(),rowsize(),nlo(),nhi(),
                2*stepi(),2*stepj(),2*diagstep(),NoMajor,NonConj);
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I == FortranStyle) { --i1; --j1; }
            return rec_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(), stepj(), S, NonConj TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                S==ColMajor ? istep == 1 ? ColMajor : NoMajor : NoMajor;
            if (I == FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return rec_type(
                itsm+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                newstor, NonConj TMV_FIRSTLAST);
        }

        inline vec_type subVector(int i, int j, int istep, int jstep, int size)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I == FortranStyle) { --i; --j; }
            return vec_type(
                itsm+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(), NonConj TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            return view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), S, NonConj TMV_FIRSTLAST);
        }

        inline view_type subBandMatrix(int i1, int i2, int j1, int j2)
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==CStyle?1:0));
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-(I==CStyle?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline view_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,
                                              newnlo,newnhi,istep,jstep));
            StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                S==ColMajor ? istep == 1 ? ColMajor : NoMajor : 
                istep == 1 && jstep == 1 ? DiagMajor : NoMajor;
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const int newstepi = stepi()*istep;
            const int newstepj = stepj()*jstep;
            return view_type(
                itsm+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, newnlo, newnhi, 
                newstepi, newstepj, newstepi+newstepj, newstor, 
                                       NonConj TMV_FIRSTLAST);
        }

        inline view_type rowRange(int i1, int i2)
        {
            if (I==FortranStyle) { 
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(colsize())); --i1; 
            } else {
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(colsize()));
            }
            const int j1 = i1 > nlo() ? i1-nlo() : 0;
            const int j2 = TMV_MIN(i2 + nhi(),int(rowsize()));
            const int newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const int newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            const size_t newlin = isrm() ? 1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type colRange(int j1, int j2)
        {
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(rowsize())); --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(rowsize())); 
            }

            const int i1 = j1 > nhi() ? j1-nhi() : 0;
            const int i2 = TMV_MIN(j2 + nlo(),int(colsize()));
            const int newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const int newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            const size_t newlin = iscm() ? 1 : 0;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct(), newlin TMV_FIRSTLAST);
        }

        inline view_type diagRange(int k1, int k2)
        {
            if (I==FortranStyle) {
                TMVAssert(k1>=-nlo() && k1<=k2 && k2<=nhi()); ++k2;
            } else {
                TMVAssert(k1>=-nlo() && k1<k2 && k2<=nhi()+1);
            }

            const int i1 = k2 <= 0 ? -k2+1 : 0;
            const int i2 = TMV_MIN(int(rowsize())-k1,int(colsize()));
            const int j1 = k1 <= 0 ? 0 : k1;
            const int j2 = TMV_MIN(int(rowsize()),int(colsize())+k2-1);
            const int newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const int newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            return view_type(
                itsm+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, newnlo, newnhi, stepi(), stepj(),
                diagstep(), stor(), ct()  TMV_FIRSTLAST);
        }

        inline view_type upperBand() 
        {
            return view_type(
                itsm,TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBand()
        {
            return view_type(
                itsm,TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type upperBandOff() 
        {
            return view_type(
                itsm+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type lowerBandOff()
        {
            return view_type(
                itsm+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm),
                colsize(),rowsize(),nlo(),nhi(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                isReal(T()) ? S : NoMajor, NonConj
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
                2*stepi(),2*stepj(),2*diagstep(),NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        TMV_DEPRECATED(const_rec_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(const_rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(const_vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(const_view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_DEPRECATED(const_view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const)
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_DEPRECATED(const_view_type Cols(int j1, int j2) const)
        { return colRange(j1,j2); }
        TMV_DEPRECATED(const_view_type Rows(int i1, int i2) const)
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(const_view_type Diags(int k1, int k2) const)
        { return diagRange(k1,k2); }
        TMV_DEPRECATED(const_view_type UpperBand(
                DiagType dt=NonUnitDiag) const)
        { return upperBand(dt); }
        TMV_DEPRECATED(const_view_type LowerBand(
                DiagType dt=NonUnitDiag) const)
        { return lowerBand(dt); }
        TMV_DEPRECATED(const_view_type UpperBandOff(
                DiagType dt=NonUnitDiag) const)
        { return upperBandOff(dt); }
        TMV_DEPRECATED(const_view_type LowerBandOff(
                DiagType dt=NonUnitDiag) const)
        { return lowerBandOff(dt); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }

        TMV_DEPRECATED(rec_type SubMatrix(int i1, int i2, int j1, int j2))
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep))
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s))
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi))
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_DEPRECATED(view_type SubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep))
        { return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_DEPRECATED(view_type Cols(int j1, int j2))
        { return colRange(j1,j2); }
        TMV_DEPRECATED(view_type Rows(int i1, int i2))
        { return rowRange(i1,i2); }
        TMV_DEPRECATED(view_type Diags(int k1, int k2))
        { return diagRange(k1,k2); }
        TMV_DEPRECATED(view_type UpperBand(DiagType dt=NonUnitDiag))
        { return upperBand(dt); }
        TMV_DEPRECATED(view_type LowerBand(DiagType dt=NonUnitDiag))
        { return lowerBand(dt); }
        TMV_DEPRECATED(view_type UpperBandOff(DiagType dt=NonUnitDiag))
        { return upperBandOff(dt); }
        TMV_DEPRECATED(view_type LowerBandOff(DiagType dt=NonUnitDiag))
        { return lowerBandOff(dt); }
        TMV_DEPRECATED(realpart_type Real())
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag())
        { return imagPart(); }


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
                stepi(),stepj(),diagstep(),S,NonConj,isdm()?0:linsize);
        }

        inline const_view_type transpose() const
        { 
            return const_view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(S),
                NonConj,isdm()?0:linsize);
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),S,TMV_ConjOf(T,NonConj),
                isdm()?0:linsize);
        }

        inline const_view_type adjoint() const
        { 
            return const_view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(S),
                TMV_ConjOf(T,NonConj),isdm()?0:linsize);
        }

        inline const_c_vec_type constLinearView() const
        {
#ifdef TMV_USE_VALGRIND
            std::vector<bool> ok(linsize,false);
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) 
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj() + int(cptr()-itsm1.get());
                    ok[k] = true;
                }
            for (size_t k=0;k<linsize;++k) if (!ok[k]) {
                const_cast<T*>(itsm1.get())[k] = T(-777);
            }
#endif
            return const_c_vec_type(itsm1.get(),linsize,1,NonConj); 
        }

        inline view_type view() 
        { 
            return view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),S,NonConj,isdm()?0:linsize 
                TMV_FIRSTLAST);
        }

        inline view_type transpose()
        { 
            return view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(S),
                NonConj,isdm()?0:linsize TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        { 
            return view_type(
                itsm,colsize(),rowsize(),nlo(),nhi(),
                stepi(),stepj(),diagstep(),S,TMV_ConjOf(T,NonConj),
                isdm()?0:linsize TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        { 
            return view_type(
                itsm,rowsize(),colsize(),nhi(),nlo(),
                stepj(),stepi(),diagstep(),TMV_TransOf(S),
                TMV_ConjOf(T,NonConj),isdm()?0:linsize TMV_FIRSTLAST);
        }

        inline vec_type linearView()
        {
#ifdef TMV_USE_VALGRIND
            std::vector<bool> ok(linsize,false);
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) 
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj() + int(cptr()-itsm1.get());
                    ok[k] = true;
                }
            for (size_t k=0;k<linsize;++k) if (!ok[k]) {
                itsm1.get()[k] = T(-777);
            }
#endif
            return vec_type(itsm1.get(),linsize,1,NonConj TMV_FIRSTLAST ); 
        }

        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(const_vec_type ConstLinearView() const)
        { return constLinearView(); }

        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Transpose())
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint())
        { return adjoint(); }
        TMV_DEPRECATED(vec_type LinearView())
        { return linearView(); }

        inline size_t colsize() const { return itscs; }
        inline size_t rowsize() const { return itsrs; }
        inline int nlo() const { return itsnlo; }
        inline int nhi() const { return itsnhi; }
        inline size_t mem_used() const { return linsize; }
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

        inline bool canLinearize() const
        { return true; }

        inline T cref(int i, int j) const
        { return itsm[i*itssi + j*itssj]; }

        inline T& ref(int i, int j)
        { return itsm[i*itssi + j*itssj]; }

        inline void resize(size_t cs, size_t rs, int lo, int hi)
        {
            linsize = BandStorageLength(S,cs,rs,lo,hi);
            itsm1.resize(linsize);
            itscs = cs;
            itsrs = rs;
            itsnlo = lo;
            itsnhi = hi;
            itssi = 
                S==RowMajor ? lo+hi : 
                S==ColMajor ? 1 : 
                rs>= cs ? 1-int(cs) : -int(rs);
            itssj = 
                S==RowMajor ? 1 :
                S==ColMajor ? lo+hi :
                -itssi+1;
            itsds = 
                S==RowMajor ? itssi+1 :
                S==ColMajor ? itssj+1 :
                1;
            itsm = S==DiagMajor ? itsm1.get()-lo*itssi : itsm1.get();
#ifdef TMVFLDEBUG
            _first = itsm.get();
            _last = _first + linsize;
#endif
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        inline rowmajor_iterator rowmajor_begin()
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end()
        {
            return rowmajor_iterator(
                this,rowstart(TMV_MIN(colsize(),rowsize()+nlo())),0); 
        }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        {
            return const_rowmajor_iterator(
                this,rowstart(TMV_MIN(colsize(),rowsize()+nlo())),0); 
        }

        inline colmajor_iterator colmajor_begin()
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end()
        {
            return colmajor_iterator(
                this,colstart(TMV_MIN(rowsize(),colsize()+nhi())),0); 
        }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        {
            return const_colmajor_iterator(
                this,colstart(TMV_MIN(rowsize(),colsize()+nhi())),0); 
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

        size_t linsize;
        AlignedArray<T> itsm1;
        size_t itscs;
        size_t itsrs;
        int itsnlo;
        int itsnhi;
        int itssi;
        int itssj;
        int itsds;
        T* itsm;

#ifdef TMVFLDEBUG
    public:
        T* _first;
        T* _last;
    protected:
#endif

        inline bool okij(int i, int j) const
        { return (j+nlo() >= i && i+nhi() >= j); }

        template <IndexStyle I2>
        friend void Swap(BandMatrix<T,S,I>& m1, BandMatrix<T,S,I2>& m2)
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

    template <class T> 
    BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2);

    template <class T> 
    BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2);

    template <class T> 
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
    //
    
    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenMatrix<T>& m, int nlo, int nhi)
    { 
        return ConstBandMatrixView<T>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const ConstMatrixView<T,I>& m, int nlo, int nhi)
    { 
        return ConstBandMatrixView<T,I>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const Matrix<T,S,I>& m, int nlo, int nhi)
    { 
        return ConstBandMatrixView<T,I>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        const MatrixView<T,I>& m, int nlo, int nhi)
    {  
        return BandMatrixView<T,I>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        Matrix<T,S,I>& m, int nlo, int nhi)
    {  
        return BandMatrixView<T,I>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenBandMatrix<T>& m, int nlo, int nhi)
    { 
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
        return ConstBandMatrixView<T>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const ConstBandMatrixView<T,I>& m, int nlo, int nhi)
    { 
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
        return ConstBandMatrixView<T,I>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const BandMatrix<T,S,I>& m, int nlo, int nhi)
    { 
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
        return ConstBandMatrixView<T,I>(
            m.cptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        const BandMatrixView<T,I>& m, int nlo, int nhi)
    { 
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
        return BandMatrixView<T,I>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        BandMatrix<T,S,I>& m, int nlo, int nhi)
    { 
        TMVAssert(nlo <= m.nlo() && nhi <= m.nhi()); 
        return BandMatrixView<T,I>(
            m.ptr(),
            TMV_MIN(m.colsize(),m.rowsize()+nlo),
            TMV_MIN(m.rowsize(),m.colsize()+nhi),nlo,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(const GenDiagMatrix<T>& m)
    { 
        return ConstBandMatrixView<T>(
            m.diag().cptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),
            m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct());
    }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const ConstDiagMatrixView<T,I>& m)
    {
        return ConstBandMatrixView<T,I>(
            m.diag().cptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),
            m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct());
    }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(const DiagMatrix<T,I>& m)
    {
        return ConstBandMatrixView<T,I>(
            m.diag().cptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),
            m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct());
    }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(const DiagMatrixView<T,I>& m)
    {
        return BandMatrixView<T,I>(
            m.diag().ptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),
            m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct()
            TMV_FIRSTLAST1(m.diag()._first,m.diag()._last) ); 
    }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(DiagMatrix<T,I>& m)
    {
        return BandMatrixView<T,I>(
            m.diag().ptr(),m.size(),m.size(),0,0,
            m.diag().step()-1,1,m.diag().step(),
            m.diag().step()==1?DiagMajor:RowMajor,m.diag().ct()
            TMV_FIRSTLAST1(m.diag()._first,m.diag()._last) ); 
    }

    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenUpperTriMatrix<T>& m, int nhi=-1)
    { 
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T>(
            m.cptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
    }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const ConstUpperTriMatrixView<T,I>& m, int nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T,I>(
            m.cptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
    }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const UpperTriMatrix<T,D,S,I>& m, int nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(D==NonUnitDiag);
        return ConstBandMatrixView<T,I>(
            m.cptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        const UpperTriMatrixView<T,I>& m, int nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(!m.isunit());
        return BandMatrixView<T,I>(
            m.ptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct() 
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        UpperTriMatrix<T,D,S,I>& m, int nhi=-1)
    {
        if (nhi < 0) nhi = m.size()-1;
        TMVAssert(D==NonUnitDiag);
        return BandMatrixView<T,I>(
            m.ptr(),m.size(),m.size(),0,nhi,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct() 
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const GenLowerTriMatrix<T>& m, int nlo=-1)
    { 
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T>(
            m.cptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
    }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const ConstLowerTriMatrixView<T,I>& m, int nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return ConstBandMatrixView<T,I>(
            m.cptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct());
    }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> BandMatrixViewOf(
        const LowerTriMatrix<T,D,S,I>& m, int nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(D==NonUnitDiag);
        return ConstBandMatrixView<T,I>(
            m.cptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        const LowerTriMatrixView<T,I>& m, int nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(!m.isunit());
        return BandMatrixView<T,I>(
            m.ptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> BandMatrixViewOf(
        LowerTriMatrix<T,D,S,I>& m, int nlo=-1)
    {
        if (nlo < 0) nlo = m.size() - 1;
        TMVAssert(D==NonUnitDiag);
        return BandMatrixView<T,I>(
            m.ptr(),m.size(),m.size(),nlo,0,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.stor(),m.ct() 
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, size_t cs, size_t rs, int nlo, int nhi, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo < int(cs));
        TMVAssert(nhi < int(rs));
        const int stepi = (
            stor == RowMajor ? nlo+nhi : 
            stor == ColMajor ? 1 : 
            rs >= cs ? -int(cs)+1 : -int(rs) );
        const int stepj = (
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo+nhi : 
            rs >= cs ? int(cs) : int(rs)+1 );
        return ConstBandMatrixView<T>(
            m,cs,rs,nlo,nhi,stepi,stepj,stor,NonConj);
    }

    template <class T> 
    BandMatrixView<T> BandMatrixViewOf(
        T* m, size_t cs, size_t rs, int nlo, int nhi, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo < int(cs));
        TMVAssert(nhi < int(rs));
        const int stepi = (
            stor == RowMajor ? nlo+nhi : 
            stor == ColMajor ? 1 : 
            rs >= cs ? -int(cs)+1 : -int(rs) );
        const int stepj = (
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo+nhi : 
            rs >= cs ? int(cs) : int(rs)+1 );
        return BandMatrixView<T>(
            m,cs,rs,nlo,nhi,stepi,stepj,stepi+stepj,stor,NonConj
            TMV_FIRSTLAST1(m,m+BandStorageLength(stor,cs,rs,nlo,nhi)));
    }

    template <class T> 
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, size_t cs, size_t rs, int nlo, int nhi,
        int stepi, int stepj)
    {
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo < int(cs));
        TMVAssert(nhi < int(rs));
        const StorageType stor = (
            stepi == 1 ? ColMajor : stepj == 1 ? RowMajor :
            stepi+stepj == 1 ? DiagMajor : NoMajor);
        return BandMatrixView<T>(
            m,cs,rs,nlo,nhi,stepi,stepj,stepi+stepj,stor,NonConj);
    }

    template <class T> 
    inline BandMatrixView<T> BandMatrixViewOf(
        T* m, size_t cs, size_t rs, int nlo, int nhi,
        int stepi, int stepj)
    {
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(nlo < int(cs));
        TMVAssert(nhi < int(rs));
        const StorageType stor = (
            stepi == 1 ? ColMajor : stepj == 1 ? RowMajor :
            stepi+stepj == 1 ? DiagMajor : NoMajor);
        return BandMatrixView<T>(
            m,cs,rs,nlo,nhi,stepi,stepj,stepi+stepj,stor,NonConj);
    }


    //
    // Copy
    //

    template <class T1, class T2> 
    inline void DoCopy(
        const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
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

        if (m1.stor() == m2.stor() && m1.canLinearize() && m2.canLinearize()) {
            TMVAssert(m1.constLinearView().size() == m2.linearView().size());
            TMVAssert(m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj());
            m2.linearView() = m1.constLinearView();
        }
        else {
            const int lo = m2.nlo();
            const int hi = m2.nhi();
            for(int k = -lo; k <= hi; ++k) m2.diag(k) = m1.diag(k);
        }
    }

    template <class T> 
    inline void DoCopy(
        const GenBandMatrix<std::complex<T> >& , const BandMatrixView<T>& )
    { TMVAssert(TMV_FALSE); }

    template <class T1, class T2> 
    inline void DoCopy1(
        const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
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
            }
            else if (m2.isconj()) DoCopy(m1.conjugate(),m2.conjugate());
            else DoCopy(m1,m2);
        }
    }

    template <class T1, class T2> 
    inline void Copy(
        const GenBandMatrix<T1>& m1, const BandMatrixView<T2>& m2)
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

    template <class T> 
    inline void Swap(
        const BandMatrixView<T>& m1, const BandMatrixView<T>& m2);

    template <class T, StorageType S, IndexStyle I> 
    inline void Swap(
        const BandMatrixView<T>& m1, BandMatrix<T,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, StorageType S, IndexStyle I> 
    inline void Swap(
        BandMatrix<T,S,I>& m1, const BandMatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, StorageType S1, IndexStyle I1, StorageType S2, IndexStyle I2> 
    inline void Swap(BandMatrix<T,S1,I1>& m1, BandMatrix<T,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }


    //
    // Views:
    //

    template <class T> 
    inline ConstBandMatrixView<T> Transpose(const GenBandMatrix<T>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Transpose(const ConstBandMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Transpose(const BandMatrix<T,S,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> Transpose(const BandMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> Transpose(BandMatrix<T,S,I>& m)
    { return m.transpose(); }

    template <class T> 
    inline ConstBandMatrixView<T> Conjugate(const GenBandMatrix<T>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Conjugate(const ConstBandMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Conjugate(const BandMatrix<T,S,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> Conjugate(const BandMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> Conjugate(BandMatrix<T,S,I>& m)
    { return m.conjugate(); }

    template <class T> 
    inline ConstBandMatrixView<T> Adjoint(const GenBandMatrix<T>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Adjoint(const ConstBandMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstBandMatrixView<T,I> Adjoint(const BandMatrix<T,S,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline BandMatrixView<T,I> Adjoint(const BandMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, StorageType S, IndexStyle I> 
    inline BandMatrixView<T,I> Adjoint(BandMatrix<T,S,I>& m)
    { return m.adjoint(); }

    template <class T> 
    inline QuotXB<T,T> Inverse(const GenBandMatrix<T>& m)
    { return m.inverse(); }



    //
    // BandMatrix ==, != BandMatrix
    //

    template <class T1, class T2> 
    bool operator==(
        const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2);

    template <class T1, class T2> 
    inline bool operator!=(
        const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    bool operator==(
        const GenBandMatrix<T1>& m1, const GenMatrix<T2>& m2);

    template <class T1, class T2> 
    inline bool operator==(
        const GenMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    { return m2 == m1; }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenBandMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& fin, auto_ptr<BandMatrix<T,S,I> >& m);

    template <class T> 
    std::istream& operator>>(
        std::istream& fin, const BandMatrixView<T>& m);

    template <class T, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(std::istream& fin, BandMatrix<T,S,I>& m)
    { return fin >> m.view(); }

} // namespace tmv

#endif
