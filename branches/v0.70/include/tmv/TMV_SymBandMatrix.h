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
// This file defines the TMV SymBandMatrix and HermBandMatrix classes.
//
// Constructors:
//
//    SymBandMatrix is used for symmetric band matrices, and 
//    HermBandMatrix is used for Hermitian band matrices.  
//    For real matrices, these are the same thing:
//        A = A.transpose().
//    But for complex, they are different:
//        A_sym = A_sym.transpose()
//        A_herm = A_herm.adjoint()
//
//    For these notes, I will always write SymBandMatrix, but (except where
//    otherwise indicated) everything applies the same for Sym and Herm.
//    Also, the Views keep track of sym/herm difference with a parameter,
//    so it is always a GenSymBandMatrix, ConstSymBandMatrixView, or 
//    SymBandMatrixView - never Herm in any of these.
//
//    Caveat: Complex Hermitian matrices are such that A = At, which 
//    implies that their diagonal elements are real.  Many routines
//    involving HermBandMatrixes assume the reality of the diagonal.
//    However, it is possible to assign a non-real value to a diagonal
//    element.  If the user does this, the results are undefined.
//
//    In addition to the type template parameter (T), SymBandMatrixes have two
//    additional template parameters:
//        UpLoType uplo = Upper || Lower
//        StorageType stor = RowMajor || ColMajor || DiagMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {Upper, RowMajor}
//
//        The first, uplo, refers to which triangular half stores the actual 
//        data.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    SymBandMatrix<T,uplo,stor>(int n, int nlo)
//        Makes a SymBandMatrix with column size = row size = n
//        and nlo (=nhi) off-diagonals with _uninitialized_ values.
//
//    SymBandMatrix<T,uplo,stor>(int n, int nlo, T x)
//        Makes a SymBandMatrix with values all initialized to x
//        For Hermitian matrixces, x must be real.
//
//    SymBandMatrix<T,uplo,stor>(const Matrix<T>& m)
//    SymBandMatrix<T,uplo,stor>(const SymMatrix<T>& m)
//    SymBandMatrix<T,uplo,stor>(const BandMatrix<T>& m)
//        Makes a SymBandMatrix which copies the corresponding elements of m.
//
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(
//            const Matrix<T>& m,  uplo, int nlo)
//    ConstSymBandMatrixView<T> HermBandMatrixViewOf(
//            const Matrix<T>& m, uplo, int nlo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(
//            const BandMatrix<T>& m, uplo)
//    ConstSymBandMatrixView<T> HermBandMatrixViewOf(
//            const BandMatrix<T>& m, uplo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(
//            const BandMatrix<T>& m, uplo, int nlo)
//    ConstSymBandMatrixView<T> HermBandMatrixViewOf(
//            const BandMatrix<T>& m, uplo, int nlo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(
//            const SymMatrix<T>& m, int nlo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(
//            const HermMatrix<T>& m, int nlo)
//        Makes a constant SymBandMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//
//    SymMatrixView<T> SymMatrixViewOf(Matrix<T>& m, uplo, int nlo)
//    SymMatrixView<T> HermMatrixViewOf(Matrix<T>& m, uplo, int nlo)
//    [ ... ]
//        Makes a modifiable SymMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the original Matrix.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(
//            const T* m, int size, int nlo, uplo, stor)
//    ConstSymMatrixView<T> HermMatrixViewOf(
//            const T* m, int size, int nlo, uplo, stor)
//    SymMatrixView<T> SymMatrixViewOf(
//            T* m, int size, int nlo, uplo, stor)
//    SymMatrixView<T> HermMatrixViewOf(
//            T* m, int size, int nlo, uplo, stor)
//        View the actual memory pointed to by m as a 
//        SymBandMatrix/HermBandMatrix with the given size and storage.
//
// Access Functions
//
//    int colsize() const
//    int rowsize() const
//    int size() const
//    int nlo() const
//    int nhi() const
//        Return the dimensions of the SymMatrix.  
//        (colsize()==rowsize()==size() and nlo()==nhi().)
//
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row(int i, int j1, int j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row
//        and be entirely in the upper or lower band.
//
//    Vector& col(int j, int i1, int i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column
//        and be entirely in the upper or lower band.
//
//    Vector& diag()
//        Return the main diagonal
//
//    Vector& diag(int i, int j1, int j2)
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
//    setZero()
//    setAllTo(T x) 
//    addToAll(T x) 
//        For HermBandMatrix, x must be real in both setAllTo and addToAll.
//    clip(RT thresh)
//    conjugateSelf()
//    transposeSelf()
//    setToIdentity(x = 1)
//    void Swap(SymBandMatrix& m1, SymBandMatrix& m2)
//        The SymBandMatrices must be the same size and Hermitianity.
//
// Views of a SymBandMatrix:
//
//    subMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//    subBandMatrix(int i1, int i2, int j1, int j2, int nlo, int nhi, 
//        int istep=1, int jstep=1)
//    subBandMatrix(int i1, int i2, int j1, int j2)
//    subVector(int i, int j, int istep, int jstep, int size)
//        Just like for a regular band matrix except that the 
//        sub-matrix/bandmatrix/vector must be completely contained within 
//        either the upper or lower band.
//
//    subSymBandMatrix(int i1, int i2, int nlo, int istep)
//    subSymMatrix(int i1, int i2, int istep)
//        Returns the SymMatrix or SymBandMatrix which runs from i1 to i2 
//        along the diagonal (not including i2) with an optional step, 
//        and includes the off diagonals in the same rows/cols.
//        The first version allows for having fewer off-diagonals than the 
//        original.
//        The parameter nlo may be omitted in which case it is taken to be the 
//        same number of off-diagonals as the original matrix.
//        The second version must have i2-i1 <= nlo+1 (ie. the SymMatrix is 
//        fully contained within the band structure).
//
//    symDiagRange(int nlo)
//        Returns a SymBandMatrix which has the same size as the original
//        but fewer off-diagonals.
//
//    diagRange(int k1, int k2)
//        Returns a BandMatrix of a range of diagonals of the SymBandMatrix.
//
//    lowerBand()
//    upperBand()
//        Returns a BandMatrix of either the upper or lower portion of
//        the SymBandMatrix.
//
//    lowerBandOff()
//    upperBandOff()
//        Returns a BandMatrix of either the strictly upper or lower 
//        portion of the SymBandMatrix.  ie. not including the main diagonal.
//
//    view(m)
//    transpose(m)
//    adjoint(m)
//    conjugate(m)
//
//
// Functions of Matrixs:
//
//    m.det() or Det(m)
//    m.logDet() or m.logDet(T* sign) or LogDet(m)
//    m.trace() or Trace(m)
//    m.sumElements() or SumElements(m)
//    m.sumAbsElements() or SumAbsElements(m)
//    m.sumAbs2Elements() or SumAbs2Elements(m)
//    m.norm() or m.normF() or Norm(m) or NormF(m)
//    m.normSq() or NormSq(m)
//    m.norm1() or Norm1(m)
//    m.norm2() or Norm2(m)
//    m.normInf() or NormInf(m)
//    m.maxAbsElement() or MaxAbsElements(m)
//    m.maxAbs2Element() or MaxAbs2Elements(m)
//
//    m.inverse() or Inverse(m)
//    m.makeInverse(minv) // Takes either a SymMatrix or Matrix argument
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
//          sB/hB size nlo
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(nlo,0) ... m(nlo,nlo) )
//          ( m(size-1,size-nlo-1) ... m(size-1,size-1) )
//
//    is >> m
//    is >> CompactIO() >> m
//        Reads m from istream is in either format.
//
//
// Division Control Functions:
//
//    m.divideUsing(dt)
//    where dt is LU, CH, or SV
//     
//    m.lud(), m.chd(), m.svd(), and m.symsvd() return the corresponding 
//        Divider classes.
//
//    As for SymMatrix, LU actually does an LDLT decomposition.
//        ie. L(Lower Triangle) * Diagonal * L.transpose()
//        For non-positive definite matrices, D may have 2x2 blocks along
//        the diagonal, which can then be further decomposed via normal LU
//        decomposition.
//        
//    Remember that Cholskey decomposition is only appropriate if you know 
//        that your HermMatrix is positive definite.  In this case, the 
//        HermmMatrix can be decomposed into L*L.adjoint().  
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


#ifndef TMV_SymBandMatrix_H
#define TMV_SymBandMatrix_H

#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_SymBandMatrixArithFunc.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Array.h"
#include <vector>

namespace tmv {

    template <class T>
    struct SymBandCopyHelper // real
    { typedef SymBandMatrix<T> type; };
    template <class T>
    struct SymBandCopyHelper<std::complex<T> > // complex
    { typedef BandMatrix<std::complex<T> > type; };

    template <class T> 
    class GenSymBandMatrix : 
        virtual public AssignableToSymBandMatrix<T>,
        virtual public AssignableToDiagMatrix<T>,
        public BaseMatrix<T>,
        private DivHelper<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenSymBandMatrix<T> type;
        typedef typename SymBandCopyHelper<T>::type copy_type;
        typedef ConstSymBandMatrixView<T> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef ConstMatrixView<T> const_rec_type;
        typedef ConstSymMatrixView<T> const_sym_type;
        typedef ConstBandMatrixView<T> const_band_type;
        typedef ConstSymBandMatrixView<RT> const_realpart_type;
        typedef SymBandMatrixView<T> nonconst_type;

        //
        // Constructors
        //

        inline GenSymBandMatrix() {}
        inline GenSymBandMatrix(const type& rhs) {}
        virtual inline ~GenSymBandMatrix() {}

        //
        // Access Functions
        //

        using AssignableToSymMatrix<T>::size;
        inline int colsize() const { return size(); }
        inline int rowsize() const { return size(); }
        using AssignableToSymMatrix<T>::sym;
        using AssignableToBandMatrix<T>::nlo;
        inline int nhi() const { return nlo(); }

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            return cref(i,j);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
                return const_vec_type(
                    cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
                    ct()); 
            else
                return const_vec_type(
                    cptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
                return const_vec_type(
                    cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
                    ct()); 
            else 
                return const_vec_type(
                    cptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_vec_type diag() const
        { return const_vec_type(cptr(),size(),diagstep(),ct()); }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            if (i>=0)
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()+i*stepj(),size()-i,
                        diagstep(),ct()); 
                else
                    return const_vec_type(
                        cptr()+i*stepi(),size()-i,
                        diagstep(),issym()?ct():TMV_ConjOf(T,ct()));
            else
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()-i*stepj(),size()+i,
                        diagstep(),issym()?ct():TMV_ConjOf(T,ct()));
                else
                    return const_vec_type(
                        cptr()-i*stepi(),size()+i,
                        diagstep(),ct()); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            if (i>=0)
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()+i*stepj()+j1*diagstep(),
                        j2-j1,diagstep(),ct()); 
                else
                    return const_vec_type(
                        cptr()+i*stepi()+j1*diagstep(),
                        j2-j1,diagstep(),issym()?ct():TMV_ConjOf(T,ct()));
            else
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()-i*stepj()+j1*diagstep(),
                        j2-j1,diagstep(),issym()?ct():TMV_ConjOf(T,ct()));
                else
                    return const_vec_type(
                        cptr()-i*stepi()+j1*diagstep(),
                        j2-j1,diagstep(),ct()); 
        }

        template <class T2> 
        inline bool isSameAs(const BaseMatrix<T2>& ) const
        { return false; }

        inline bool isSameAs(const type& m2) const
        {
            if (this == &m2) return true;
            else if (cptr()==m2.cptr() && size()==m2.size() && 
                     nlo()==m2.nlo() && sym() == m2.sym()) {
                if (uplo() == m2.uplo())
                    return (stepi() == m2.stepi() && stepj() == m2.stepj() && 
                            ct() == m2.ct());
                else
                    return (stepi() == m2.stepj() && stepj() == m2.stepi() && 
                            issym() == (ct()==m2.ct()));
            } else return false;
        }

        inline void assignToM(const MatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            assignToB(BandMatrixViewOf(m2,nlo(),nlo()));
            if (size() > nlo()+1) {
                m2.upperTri().offDiag(nlo()+1).setZero();
                m2.lowerTri().offDiag(nlo()+1).setZero();
            }
        }

        inline void assignToM(const MatrixView<CT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            assignToB(BandMatrixViewOf(m2,nlo(),nlo()));
            if (size() > nlo()+1) {
                m2.upperTri().offDiag(nlo()+1).setZero();
                m2.lowerTri().offDiag(nlo()+1).setZero();
            }
        }

        inline void assignToD(const DiagMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == size());
            TMVAssert(nlo() == 0);
            m2.diag() = diag();
        }

        inline void assignToD(const DiagMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(nlo() == 0);
            m2.diag() = diag();
            if (!issym()) m2.diag().imagPart().setZero();
        }

        inline void assignToB(const BandMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(m2.nhi() >= nlo());
            assignTosB(SymBandMatrixViewOf(m2,Upper,nlo()));
            if (nlo() > 0) 
                m2.diagRange(-nlo(),0) = m2.diagRange(1,nlo()+1).transpose();
            if (m2.nlo() > nlo()) m2.diagRange(-m2.nlo(),-nlo()).setZero();
            if (m2.nhi() > nlo()) m2.diagRange(nlo()+1,m2.nhi()+1).setZero();
        }

        inline void assignToB(const BandMatrixView<CT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(m2.nhi() >= nlo());
            if (issym()) {
                assignTosB(SymBandMatrixViewOf(m2,Upper,nlo()));
                if (nlo() > 0) 
                    m2.diagRange(-nlo(),0) =
                        m2.diagRange(1,nlo()+1).transpose();
            } else {
                m2.diag().imagPart().setZero();
                assignTosB(HermBandMatrixViewOf(m2,Upper,nlo()));
                if (nlo() > 0) 
                    m2.diagRange(-nlo(),0) = m2.diagRange(1,nlo()+1).adjoint();
            }
            if (m2.nlo() > nlo()) m2.diagRange(-m2.nlo(),-nlo()).setZero();
            if (m2.nhi() > nlo()) m2.diagRange(nlo()+1,m2.nhi()+1).setZero();
        }

        inline void assignToS(const SymMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == size());
            assignTosB(SymBandMatrixViewOf(m2,nlo()));
            if (size() > nlo()+1) m2.upperTri().offDiag(nlo()+1).setZero();
        }

        inline void assignToS(const SymMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(m2.sym() == sym());
            if (issym()) assignTosB(SymBandMatrixViewOf(m2,nlo()));
            else assignTosB(HermBandMatrixViewOf(m2,nlo()));
            if (size() > nlo()+1) m2.upperTri().offDiag(nlo()+1).setZero();
        }

        inline void assignTosB(const SymBandMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == size());
            TMVAssert(m2.nlo() >= nlo());

            if (!isSameAs(m2)) m2.upperBand() = upperBand();
            if (m2.nlo() > nlo()) m2.diagRange(-m2.nlo(),-nlo()).setZero();
        }

        inline void assignTosB(const SymBandMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(m2.nlo() >= nlo());
            TMVAssert(isReal(T()) || m2.sym() == sym());

            if (!isSameAs(m2)) m2.upperBand() = upperBand();
            if (m2.nlo() > nlo()) m2.diagRange(-m2.nlo(),-nlo()).setZero();
        }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            if ( (uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1) )
                return const_rec_type(
                    cptr()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),stor(),ct());
            else
                return const_rec_type(
                    cptr()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if ( (uplo()==Upper && i2-j1<=istep) || 
                 (uplo()==Lower && j2-i1<=jstep) )
                return const_rec_type(
                    cptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,ct());
            else 
                return const_rec_type(
                    cptr()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),issym()?ct():TMV_ConjOf(T,ct()));
        }

        bool hasSubBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const;

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if ((i1+newnlo-j1<=0 && uplo()==Upper) || 
                (j1+newnhi-i1<=0 && uplo()==Lower))
                return const_band_type(
                    cptr()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
                    stor(),ct());
            else
                return const_band_type(
                    cptr()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
                    TMV_TransOf(stor()),issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const int newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            const StorageType newstor=
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            if ((i1+newnlo*istep<=j1 && uplo()==Upper) || 
                (j1+newnhi*jstep<=i1 && uplo()==Lower)) {
                const int newstepi = stepi()*istep;
                const int newstepj = stepj()*jstep;
                return const_band_type(
                    cptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,newstor,ct());
            } else {
                const int newstepi = stepj()*istep;
                const int newstepj = stepi()*jstep;
                return const_band_type(
                    cptr()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,
                    TMV_TransOf(newstor),issym()?ct():TMV_ConjOf(T,ct()));
            }
        }

        bool hasSubVector(int i, int j, int istep, int jstep, int n) const;

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            if ((i-j<=0 && uplo()==Upper) || (j-i<=0 && uplo()==Lower))
                return const_vec_type(
                    cptr()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),ct());
            else
                return const_vec_type(
                    cptr()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),issym()?ct():TMV_ConjOf(T,ct()));
        }

        bool hasSubSymMatrix(int i1, int i2, int istep) const;

        inline const_sym_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return const_sym_type(
                cptr()+i1*diagstep(),i2-i1,
                stepi(),stepj(),sym(),uplo(),stor(),ct());
        }

        inline const_sym_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return const_sym_type(
                cptr()+i1*diagstep(),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),
                sym(),uplo(),istep==1 ? stor() : NoMajor,ct());
        }

        bool hasSubSymBandMatrix(int i1, int i2, int newnlo, int istep) const;

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo) const
        {
            TMVAssert(hasSubSymBandMatrix(i1,i2,newnlo,1));
            return const_view_type(
                cptr()+i1*diagstep(),i2-i1,newnlo,
                stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct());
        }

        inline const_view_type subSymBandMatrix(int i1, int i2) const
        { return subSymBandMatrix(i1,i2,TMV_MIN(nlo(),i2-i1-1)); }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep) const
        {
            TMVAssert(hasSubSymBandMatrix(i1,i2,newnlo,istep));
            return const_view_type(
                cptr()+i1*(diagstep()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),istep*diagstep(),
                sym(),uplo(),istep==1 ? stor() : NoMajor,ct());
        }

        inline const_band_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
            TMVAssert(k2 <= 1 || k1 >= 0);

            if (k1 >= 0) {
                const int newsize = size()-k1;
                const int newnhi = k2-k1-1;
                if (uplo() == Upper)
                    return const_band_type(
                        cptr()+k1*stepj(), newsize, newsize, 0, newnhi,
                        stepi(), stepj(), diagstep(), stor(), ct());
                else
                    return const_band_type(
                        cptr()+k1*stepi(), newsize, newsize, 0, newnhi,
                        stepj(), stepi(), diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()));
            } else {
                const int newsize = size()+k2-1;
                const int newnlo = k2-k1-1;
                if (uplo() == Lower)
                    return const_band_type(
                        cptr()-k2*stepi(), newsize, newsize, newnlo, 0,
                        stepi(), stepj(), diagstep(), stor(), ct());
                else
                    return const_band_type(
                        cptr()-k2*stepj(), newsize, newsize, newnlo, 0,
                        stepj(), stepi(), diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()));
            }
        }

        inline const_view_type symDiagRange(int newnlo) const
        {
            TMVAssert(newnlo>=0 && newnlo <= nlo());
            return const_view_type(
                cptr(),size(),newnlo,
                stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct());
        }

        inline const_band_type upperBand() const
        {
            if (uplo() == Upper)
                return const_band_type(
                    cptr(),size(),size(),0,nlo(),
                    stepi(),stepj(),diagstep(),stor(),ct());
            else
                return const_band_type(
                    cptr(),size(),size(),0,nlo(),
                    stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_band_type upperBandOff() const
        {
            TMVAssert(size()>0);
            TMVAssert(nlo()>0);
            if (uplo() == Upper)
                return const_band_type(
                    cptr()+stepj(),size()-1,size()-1,
                    0,nlo()-1,stepi(),stepj(),diagstep(),stor(),ct());
            else
                return const_band_type(
                    cptr()+stepi(),size()-1,size()-1,
                    0,nlo()-1,stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_band_type lowerBand() const
        {
            if (uplo() == Lower)
                return const_band_type(
                    cptr(),size(),size(),nlo(),0,
                    stepi(),stepj(),diagstep(),stor(),ct());
            else
                return const_band_type(
                    cptr(),size(),size(),nlo(),0,
                    stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_band_type lowerBandOff() const
        {
            TMVAssert(size()>0);
            TMVAssert(nlo()>0);
            if (uplo() == Lower)
                return const_band_type(
                    cptr()+stepi(),size()-1,size()-1,
                    nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct());
            else
                return const_band_type(
                    cptr()+stepj(),size()-1,size()-1,
                    nlo()-1,0,stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),size(),nlo(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                Sym, uplo(), isReal(T()) ? stor() : NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(issym());
            // The imaginary part of a Hermitian matrix is anti-symmetric
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                size(),nlo(),2*stepi(),2*stepj(),2*diagstep(),
                Sym,uplo(),NoMajor,NonConj);
        }

        // 
        // Views
        //

        inline const_view_type view() const
        {
            return const_view_type(
                cptr(),size(),nlo(),
                stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct());
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                cptr(),size(),nlo(),
                stepj(),stepi(),diagstep(),sym(),
                TMV_UTransOf(uplo()),TMV_TransOf(stor()),ct());
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                cptr(),size(),nlo(), stepi(),stepj(),diagstep(),
                sym(),uplo(),stor(),TMV_ConjOf(T,ct()));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                cptr(),size(),nlo(), stepj(),stepi(),diagstep(),sym(),
                TMV_UTransOf(uplo()),TMV_TransOf(stor()),TMV_ConjOf(T,ct()));
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<T*>(cptr()),size(),nlo(),
                stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct()
                TMV_FIRSTLAST1(
                    cptr(),row(colsize()-1,0,colsize()).end().getP()));
        }

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
        inline RT norm2() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::norm2(); 
            else return doNorm2();
        }

        inline RT normInf() const
        { return norm1(); }

        RT maxAbsElement() const;
        RT maxAbs2Element() const;

        RT doCondition() const;
        inline RT condition() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::condition();
            else return doCondition();
        }

        inline bool isSingular() const
        { return DivHelper<T>::isSingular(); }

        template <class T1> 
        void doMakeInverse(const SymMatrixView<T1>& sinv) const;

        template <class T1> 
        inline void makeInverse(const SymMatrixView<T1>& minv) const
        {
            TMVAssert(minv.size() == size());
            TMVAssert(isherm() == minv.isherm());
            TMVAssert(issym() == minv.issym());
            doMakeInverse(minv);
        }

        template <UpLoType U, StorageType S, IndexStyle I> 
        inline void makeInverse(SymMatrix<T,U,S,I>& sinv) const
        {
            TMVAssert(issym());
            makeInverse(sinv.view()); 
        }

        template <UpLoType U, StorageType S, IndexStyle I> 
        inline void makeInverse(HermMatrix<T,U,S,I>& sinv) const
        {
            TMVAssert(isherm());
            makeInverse(sinv.view()); 
        }

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

        QuotXsB<T,T> QInverse() const;
        inline QuotXsB<T,T> inverse() const
        { return QInverse(); }

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
            TMVAssert(dt == CH || dt == LU || dt == SV);
            TMVAssert(isherm() || dt != CH);
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

        inline const HermBandCHDiv<T>& chd() const
        {
            TMVAssert(isherm());
            divideUsing(CH);
            setDiv();
            TMVAssert(getDiv());
            TMVAssert(divIsCHDiv());
            return static_cast<const HermBandCHDiv<T>&>(*getDiv());
        }

        inline const HermBandSVDiv<T>& svd() const
        {
            TMVAssert(isherm());
            divideUsing(SV);
            setDiv();
            TMVAssert(getDiv());
            TMVAssert(divIsHermSVDiv());
            return static_cast<const HermBandSVDiv<T>&>(*getDiv());
        }

        inline const SymBandSVDiv<T>& symsvd() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(issym());
            divideUsing(SV);
            setDiv();
            TMVAssert(getDiv());
            TMVAssert(divIsSymSVDiv());
            return static_cast<const SymBandSVDiv<T>&>(*getDiv());
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

        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        virtual const T* cptr() const = 0;
        virtual int stepi() const = 0;
        virtual int stepj() const = 0;
        virtual int diagstep() const = 0;
        virtual UpLoType uplo() const = 0;
        virtual StorageType stor() const = 0;
        virtual ConjType ct() const = 0;
        virtual inline bool isrm() const { return stor() == RowMajor; }
        virtual inline bool iscm() const { return stor() == ColMajor; }
        virtual inline bool isdm() const { return stor() == DiagMajor; }
        inline bool isconj() const
        {
            TMVAssert(isComplex(T()) || ct()==NonConj);
            return isComplex(T()) && ct() == Conj; 
        }
        inline bool isherm() const { return isReal(T()) || sym() == Herm; }
        inline bool issym() const { return isReal(T()) || sym() == Sym; }
        inline bool isupper() const { return uplo() == Upper; }

        inline bool isHermOK() const
        {
            if (issym()) return true;
            else return diag().imagPart().normInf() == RT(0);
        }

        virtual T cref(int i, int j) const;

    protected :

        using DivHelper<T>::getDiv;

        inline bool okij(int i, int j) const
        { return (j+nlo() >= i && i+nlo() >= j); }

        void newDivider() const;
        inline const BaseMatrix<T>& getMatrix() const { return *this; }

    private :

        type& operator=(const type&);

        bool divIsLUDiv() const;
        bool divIsCHDiv() const;
        bool divIsHermSVDiv() const;
        bool divIsSymSVDiv() const;

    }; // GenSymBandMatrix

    template <class T, IndexStyle I> 
    class ConstSymBandMatrixView : public GenSymBandMatrix<T>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef ConstSymBandMatrixView<T,I> type;
        typedef GenSymBandMatrix<T> base;

        inline ConstSymBandMatrixView(const type& rhs) :
            itsm(rhs.itsm), itss(rhs.itss), itslo(rhs.itslo),
            itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
            itssym(rhs.itssym), itsuplo(rhs.itsuplo), itsstor(rhs.itsstor),
            itsct(rhs.itsct) {}

        inline ConstSymBandMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itss(rhs.size()), itslo(rhs.nlo()), 
            itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
            itssym(rhs.sym()), itsuplo(rhs.uplo()), itsstor(rhs.stor()),
            itsct(rhs.ct()) {}

        inline ConstSymBandMatrixView(
            const T* _m, int _s, int _lo, int _si, int _sj, int _sd,
            SymType _sym, UpLoType _uplo, StorageType _stor,
            ConjType _ct) : 
            itsm(_m), itss(_s), itslo(_lo), itssi(_si), itssj(_sj), itssd(_sd),
            itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
        {
            TMVAssert(_stor==RowMajor ? _sj==1 : _stor==ColMajor ?
                      _si==1 : _stor==DiagMajor ? _sd==1 : true); 
            TMVAssert(size()==0 || nlo() < size()); 
        }

        virtual inline ~ConstSymBandMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        inline int size() const { return itss; }
        inline int nlo() const { return itslo; }
        inline const T* cptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline int diagstep() const { return itssd; }
        inline SymType sym() const { return itssym; }
        inline UpLoType uplo() const { return itsuplo; }
        inline StorageType stor() const { return itsstor; }
        inline ConjType ct() const { return itsct; }

    protected :

        const T*const itsm;
        const int itss;
        const int itslo;
        const int itssi;
        const int itssj;
        const int itssd;

        const SymType itssym;
        const UpLoType itsuplo;
        const StorageType itsstor;
        const ConjType itsct;

    private :

        type& operator=(const type&);

    }; // ConstSymBandMatrixView

    template <class T> 
    class ConstSymBandMatrixView<T,FortranStyle> : 
        public ConstSymBandMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef ConstSymBandMatrixView<T,FortranStyle> type;
        typedef ConstSymBandMatrixView<T,CStyle> c_type;
        typedef GenSymBandMatrix<T> base;
        typedef ConstSymBandMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstSymMatrixView<T,FortranStyle> const_sym_type;
        typedef ConstBandMatrixView<T,FortranStyle> const_band_type;
        typedef ConstSymBandMatrixView<RT,FortranStyle> const_realpart_type;
        typedef SymBandMatrixView<T,FortranStyle> nonconst_type;

        inline ConstSymBandMatrixView(const type& rhs) : c_type(rhs) {}

        inline ConstSymBandMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline ConstSymBandMatrixView(const base& rhs) : c_type(rhs) {}

        inline ConstSymBandMatrixView(
            const T* _m, int _s, int _lo,
            int _si, int _sj, int _sd,
            SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct) : 
            c_type(_m,_s,_lo,_si,_sj,_sd, _sym,_uplo,_stor,_ct) {}

        virtual inline ~ConstSymBandMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j>0 && j<=this->size());
            return base::operator()(i-1,j-1);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size());
            TMVAssert(this->okij(i-1,j1-1));
            TMVAssert(this->okij(i-1,j2-1));
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=this->size());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
            TMVAssert(this->okij(i1-1,j-1));
            TMVAssert(this->okij(i2-1,j-1));
            return base::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag() const
        { return base::diag(); }

        inline const_vec_type diag(int i) const
        { return base::diag(i); }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-this->nlo() && i<=this->nlo());
            TMVAssert(i>=-this->size() && i<=this->size()); 
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size()-std::abs(i));
            return base::diag(i,j1-1,j2);
        }

        //
        // SubMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        bool hasSubBandMatrix(
            int i1, int i2, int j1, int j2, 
            int newnlo, int newnhi, int istep, int jstep) const;

        bool hasSubVector(int i, int j, int istep, int jstep, int n) const;

        bool hasSubSymMatrix(int i1, int i2, int istep) const;

        bool hasSubSymBandMatrix(int i1, int i2, int newnlo, int istep) const;

        inline const_rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::subMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::subMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            return base::subBandMatrix(
                i1-1,i2,j1-1,j2,newnlo,newnhi);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            return base::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(this->nlo()+j1-i1,i2-i1);
            const int newnhi = TMV_MIN(this->nlo()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return base::subVector(i-1,j-1,istep,jstep,n);
        }

        inline const_sym_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return base::subSymMatrix(i1-1,i2);
        }

        inline const_sym_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return base::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo) const
        {
            TMVAssert(hasSubSymBandMatrix(i1,i2,newnlo,1));
            return base::subSymBandMatrix(i1-1,i2,newnlo);
        }

        inline const_view_type subSymBandMatrix(int i1, int i2) const
        { return subSymBandMatrix(i1,i2,TMV_MIN(this->nlo(),i2-i1)); }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep) const
        {
            TMVAssert(hasSubSymBandMatrix(i1,i2,newnlo,istep));
            return base::subSymBandMatrix(i1-1,i2-1+istep,newnlo,istep);
        }

        inline const_band_type diagRange(int k1, int k2) const
        { return base::diagRange(k1,k2); }

        inline const_view_type symDiagRange(int newnlo) const
        { return base::symDiagRange(newnlo); }

        inline const_band_type upperBand() const
        { return base::upperBand(); }

        inline const_band_type upperBandOff() const
        { return base::upperBandOff(); }

        inline const_band_type lowerBand() const
        { return base::lowerBand(); }

        inline const_band_type lowerBandOff() const
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

        inline nonconst_type nonConst() const
        { return base::nonConst(); }

    private :

        type& operator=(const type&);

    }; // FortranStyle ConstSymBandMatrixView

    template <class T, IndexStyle I> 
    class SymBandMatrixView : public GenSymBandMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymBandMatrixView<T,I> type;
        typedef GenSymBandMatrix<T> base;
        typedef SymBandMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef SymMatrixView<T,I> sym_type;
        typedef BandMatrixView<T,I> band_type;
        typedef SymBandMatrixView<RT,I> realpart_type;
        typedef TMV_RefType(T) reference;

        //
        // Constructors
        //

        inline SymBandMatrixView(const type& rhs) : 
            itsm(rhs.itsm), itss(rhs.itss), itslo(rhs.itslo), 
            itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
            itssym(rhs.itssym), itsuplo(rhs.itsuplo), itsstor(rhs.itsstor),
            itsct(rhs.itsct) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        inline SymBandMatrixView(
            T* _m, int _s, int _lo, int _si, int _sj, int _sd,
            SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
            TMV_PARAMFIRSTLAST(T) 
        ) :
            itsm(_m), itss(_s), itslo(_lo), itssi(_si), itssj(_sj), itssd(_sd),
            itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(_stor==RowMajor ? _sj==1 : _stor==ColMajor ?
                      _si==1 : _stor==DiagMajor ? _sd==1 : true); 
            TMVAssert(size()==0 || nlo() < size()); 
        }

        virtual inline ~SymBandMatrixView()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            const_cast<T*&>(itsm) = 0;
#endif
        }

        //
        // Op=
        //

        inline const type& operator=(const type& m2) const
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T()) || m2.sym() == sym());
            TMVAssert(nlo() >= m2.nlo());
            if (!this->isSameAs(m2)) {
                upperBand() = m2.upperBand();
                if (nlo() > m2.nlo()) diagRange(-nlo(),-m2.nlo()).setZero();
            }
            return *this;
        }

        inline const type& operator=(const type& m2) 
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T()) || m2.sym() == sym());
            TMVAssert(nlo() >= m2.nlo());
            if (!this->isSameAs(m2)) {
                upperBand() = m2.upperBand();
                if (nlo() > m2.nlo()) diagRange(-nlo(),-m2.nlo()).setZero();
            }
            return *this;
        }

        inline const type& operator=(const GenSymBandMatrix<RT>& m2) const
        {
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            if (!this->isSameAs(m2)) m2.assignTosB(*this); 
            return *this; 
        }

        inline const type& operator=(const GenSymBandMatrix<CT>& m2) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.sym() == sym());
            TMVAssert(nlo() >= m2.nlo());
            if (!this->isSameAs(m2)) m2.assignTosB(*this); 
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenSymBandMatrix<T2>& m2) const
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(isReal(T2()) || m2.sym() == sym());
            TMVAssert(nlo() >= m2.nlo());
            upperBand() = m2.upperBand();
            if (nlo() > m2.nlo()) diagRange(-nlo(),-m2.nlo()).setZero();
            return *this; 
        }

        inline const type& operator=(const T& x) const 
        {
            TMVAssert(issym() || TMV_IMAG(x) == RT(0));
            return setToIdentity(x); 
        }

        inline const type& operator=(
            const AssignableToSymBandMatrix<RT>& m2) const
        {
            TMVAssert(size() == m2.size());
            m2.assignTosB(view());
            return *this;
        }

        inline const type& operator=(
            const AssignableToSymBandMatrix<CT>& m2) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.sym() == sym());
            m2.assignTosB(view());
            return *this;
        }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        {
            TMVAssert(size() == m2.size());
            m2.assignToD(DiagMatrixViewOf(diag()));
            if (nlo() > 0) upperBandOff().setZero();
            return *this;
        }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            m2.assignToD(DiagMatrixViewOf(diag()));
            if (nlo() > 0) upperBandOff().setZero();
            TMVAssert(this->isHermOK());
            return *this;
        }

        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            TMVAssert(okij(i,j));
            return ref(i,j); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
                return vec_type(
                    ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
                    ct() TMV_FIRSTLAST); 
            else
                return vec_type(
                    ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
                    issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
                return vec_type(
                    ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
                    ct() TMV_FIRSTLAST); 
            else 
                return vec_type(
                    ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
                    issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type diag() const
        { return vec_type(
                ptr(),size(),stepi()+stepj(),ct() TMV_FIRSTLAST); }

        inline vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            if (i>=0)
                if (uplo()==Upper)
                    return vec_type(
                        ptr()+i*stepj(),size()-i,
                        diagstep(),ct() TMV_FIRSTLAST); 
                else
                    return vec_type(
                        ptr()+i*stepi(),size()-i,
                        diagstep(),issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST);
            else
                if (uplo()==Upper)
                    return vec_type(
                        ptr()-i*stepj(),size()+i,
                        diagstep(),issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST);
                else
                    return vec_type(
                        ptr()-i*stepi(),size()+i,
                        diagstep(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            if (i>=0)
                if (uplo()==Upper)
                    return vec_type(
                        ptr()+i*stepj()+j1*diagstep(),j2-j1,
                        diagstep(),ct() TMV_FIRSTLAST); 
                else
                    return vec_type(
                        ptr()+i*stepi()+j1*diagstep(),j2-j1,
                        diagstep(),issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST);
            else
                if (uplo()==Upper)
                    return vec_type(
                        ptr()-i*stepj()+j1*diagstep(),j2-j1,
                        diagstep(),issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST);
                else
                    return vec_type(
                        ptr()-i*stepi()+j1*diagstep(),j2-j1,
                        diagstep(),ct() TMV_FIRSTLAST); 
        }

        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { upperBand().setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        {
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            upperBand().setAllTo(x); return *this; 
        }

        inline const type& addToAll(const T& x) const
        {
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            upperBand().addToAll(x); return *this; 
        }

        inline const type& clip(RT thresh) const
        { upperBand().clip(thresh); return *this; }

        inline const type& transposeSelf() const
        { if (!this->issym()) upperBand().conjugateSelf(); return *this; }

        inline const type& conjugateSelf() const
        { if (isComplex(T())) upperBand().conjugateSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        {
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            setZero(); diag().setAllTo(x); return *this; 
        }

        //
        // SubMatrix
        //

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            if ( (uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1) )
                return rec_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),stor(),ct() TMV_FIRSTLAST);
            else
                return rec_type(
                    ptr()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if ( (uplo()==Upper && i2-j1<=istep) || 
                 (uplo()==Lower && j2-i1<=jstep) )
                return rec_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,ct() TMV_FIRSTLAST);
            else
                return rec_type(
                    ptr()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),this->issym()?ct():TMV_ConjOf(T,ct()) 
                    TMV_FIRSTLAST);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,1,1));
            if ((i1+newnlo-j1<=0 && uplo()==Upper) || 
                (j1+newnhi-i1<=0 && uplo()==Lower))
                return band_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
                    stor(),ct() TMV_FIRSTLAST);
            else
                return band_type(
                    ptr()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
                    TMV_TransOf(stor()), this->issym()?ct():TMV_ConjOf(T,ct()) 
                    TMV_FIRSTLAST);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const int newnhi = TMV_MIN(nlo()+i1-j1,j2-j1-1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(base::hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            if ((i1+newnlo*istep<=j1 && uplo()==Upper) || 
                (j1+newnhi*jstep<=i1 && uplo()==Lower)) {
                const int newstepi = stepi()*istep;
                const int newstepj = stepj()*jstep;
                return band_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,newstor,ct() 
                    TMV_FIRSTLAST);
            } else {
                const int newstepi = stepj()*istep;
                const int newstepj = stepi()*jstep;
                return band_type(
                    ptr()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,
                    TMV_TransOf(newstor),this->issym()?ct():TMV_ConjOf(T,ct()) 
                    TMV_FIRSTLAST);
            }
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,n));
            if ((i-j<=0 && uplo()==Upper) || (j-i<=0 && uplo()==Lower))
                return vec_type(
                    ptr()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),ct() TMV_FIRSTLAST);
            else
                return vec_type(
                    ptr()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline sym_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(base::hasSubSymMatrix(i1,i2,1));
            return sym_type(
                ptr()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),sym(),uplo(),stor(),ct() TMV_FIRSTLAST);
        }

        inline sym_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(base::hasSubSymMatrix(i1,i2,istep));
            return sym_type(
                ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),sym(),uplo(),
                istep==1 ? stor() : NoMajor,ct() TMV_FIRSTLAST);
        }

        inline view_type subSymBandMatrix(int i1, int i2, int newnlo) const
        {
            TMVAssert(base::hasSubSymBandMatrix(i1,i2,newnlo,1));
            return view_type(
                ptr()+i1*diagstep(),i2-i1,newnlo,
                stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline view_type subSymBandMatrix(int i1, int i2) const
        { return subSymBandMatrix(i1,i2,TMV_MIN(this->nlo(),i2-i1-1)); }

        inline view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep) const
        {
            TMVAssert(base::hasSubSymBandMatrix(i1,i2,newnlo,istep));
            return view_type(
                ptr()+i1*diagstep(),(i2-i1)/istep,newnlo,
                istep*stepi(),istep*stepj(),istep*diagstep(),
                sym(),uplo(),istep==1 ? stor() : NoMajor,ct() TMV_FIRSTLAST);
        }

        inline band_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
            TMVAssert(k2 <= 1 || k1 >= 0);

            if (k1 >= 0) {
                const int newsize = size()-k1;
                const int newnhi = k2-k1-1;
                if (uplo() == Upper)
                    return band_type(
                        ptr()+k1*stepj(),
                        newsize, newsize, 0, newnhi,
                        stepi(), stepj(), diagstep(),
                        stor(), ct() TMV_FIRSTLAST );
                else
                    return band_type(
                        ptr()+k1*stepi(),
                        newsize, newsize, 0, newnhi, 
                        stepj(), stepi(), diagstep(),
                        TMV_TransOf(stor()), issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST );
            } else {
                const int newsize = size()+k2-1;
                const int newnlo = k2-k1-1;
                if (uplo() == Lower)
                    return band_type(
                        ptr()-k2*stepi(),
                        newsize, newsize, newnlo, 0,
                        stepi(), stepj(), diagstep(),
                        stor(), ct() TMV_FIRSTLAST );
                else
                    return band_type(
                        ptr()-k2*stepj(),
                        newsize, newsize, newnlo, 0,
                        stepj(), stepi(), diagstep(),
                        TMV_TransOf(stor()), issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST );
            }
        }

        inline view_type symDiagRange(int newnlo) const
        {
            TMVAssert(newnlo>=0 && newnlo <= nlo());
            return view_type(
                ptr(),size(),newnlo,
                stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline band_type upperBand() const
        {
            if (uplo() == Upper)
                return band_type(
                    ptr(),size(),size(),0,nlo(),
                    stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
            else
                return band_type(
                    ptr(),size(),size(),0,nlo(),
                    stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline band_type upperBandOff() const
        {
            TMVAssert(size()>0);
            TMVAssert(nlo()>0);
            if (uplo() == Upper)
                return band_type(
                    ptr()+stepj(),size()-1,size()-1,
                    0,nlo()-1,stepi(),stepj(),diagstep(),stor(),ct() 
                    TMV_FIRSTLAST);
            else
                return band_type(
                    ptr()+stepi(),size()-1,size()-1,
                    0,nlo()-1,stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline band_type lowerBand() const
        {
            if (uplo() == Lower)
                return band_type(
                    ptr(),size(),size(),nlo(),0,
                    stepi(),stepj(),diagstep(),stor(),ct() TMV_FIRSTLAST);
            else
                return band_type(
                    ptr(),size(),size(),nlo(),0,
                    stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline band_type lowerBandOff() const
        {
            TMVAssert(size()>0);
            TMVAssert(nlo()>0);
            if (uplo() == Lower)
                return band_type(
                    ptr()+stepi(),size()-1,size()-1,
                    nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct() 
                    TMV_FIRSTLAST);
            else
                return band_type(
                    ptr()+stepj(),size()-1,size()-1,
                    nlo()-1,0,stepj(),stepi(),diagstep(),TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline realpart_type realPart() const
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),size(),nlo(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                sym(),uplo(), isReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(this->issym());
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1,
                size(),nlo(),2*stepi(),2*stepj(),2*diagstep(),
                sym(),uplo(),NoMajor, NonConj
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
                ptr(),size(),nlo(),
                stepj(),stepi(),diagstep(),
                sym(),TMV_UTransOf(uplo()),TMV_TransOf(stor()),ct() 
                TMV_FIRSTLAST);
        }

        inline view_type conjugate() const
        {
            return view_type(
                ptr(),size(),nlo(),
                stepi(),stepj(),diagstep(),
                sym(),uplo(),stor(),TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline view_type adjoint() const
        {
            return view_type(
                ptr(),size(),nlo(),
                stepj(),stepi(),diagstep(),
                sym(),TMV_UTransOf(uplo()),TMV_TransOf(stor()),
                TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader) const;

        inline int size() const { return itss; }
        inline int nlo() const { return itslo; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline int diagstep() const { return itssd; }
        inline SymType sym() const { return itssym; }
        inline UpLoType uplo() const { return itsuplo; }
        inline StorageType stor() const { return itsstor; }
        inline ConjType ct() const { return itsct; }
        using base::issym;
        using base::iscm;
        using base::isrm;

        reference ref(int i, int j) const;

    protected :

        T*const itsm;
        const int itss;
        const int itslo;
        const int itssi;
        const int itssj;
        const int itssd;

        const SymType itssym;
        const UpLoType itsuplo;
        const StorageType itsstor;
        const ConjType itsct;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected:
#endif

        using base::okij;

    }; // SymBandMatrixView

    template <class T> 
    class SymBandMatrixView<T,FortranStyle> : 
        public SymBandMatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymBandMatrixView<T,FortranStyle> type;
        typedef GenSymBandMatrix<T> base;
        typedef ConstSymBandMatrixView<T,FortranStyle> const_type;
        typedef SymBandMatrixView<T,FortranStyle> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef SymBandMatrixView<T,CStyle> c_type;
        typedef VectorView<T,FortranStyle> vec_type;
        typedef MatrixView<T,FortranStyle> rec_type;
        typedef SymMatrixView<T,FortranStyle> sym_type;
        typedef BandMatrixView<T,FortranStyle> band_type;
        typedef SymBandMatrixView<RT,FortranStyle> realpart_type;

        //
        // Constructors
        //

        inline SymBandMatrixView(const type& rhs) : c_type(rhs) {}

        inline SymBandMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline SymBandMatrixView(
            T* _m, int _s, int _lo, int _si, int _sj, int _sd,
            SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
            TMV_PARAMFIRSTLAST(T) 
        ) :
            c_type(_m,_s,_lo,_si,_sj,_sd,_sym,_uplo,_stor,_ct 
                   TMV_FIRSTLAST1(_first,_last) ) {}

        virtual inline ~SymBandMatrixView() {} 

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

        inline const type& operator=(const GenSymBandMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenSymBandMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenSymBandMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(
            const AssignableToSymBandMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(
            const AssignableToSymBandMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenBandMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenBandMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        //
        // Access
        //

        inline TMV_RefType(T) operator()(int i,int j) const 
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j>0 && j<=this->size());
            TMVAssert(this->okij(i-1,j-1));
            return c_type::ref(i-1,j-1); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size());
            TMVAssert(this->okij(i-1,j1-1));
            TMVAssert(this->okij(i-1,j2-1));
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=this->size());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
            TMVAssert(this->okij(i1-1,j-1));
            TMVAssert(this->okij(i2-1,j-1));
            return c_type::col(j-1,i1-1,i2);
        }

        inline vec_type diag() const
        { return c_type::diag(); }

        inline vec_type diag(int i) const
        { return c_type::diag(i); }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-this->nlo() && i<=this->nlo());
            TMVAssert(i>=-this->size() && i<=this->size()); 
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size()-std::abs(i));
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

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& transposeSelf() const
        { c_type::transposeSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { c_type::setToIdentity(x); return *this; }

        //
        // SubMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return const_type(*this).hasSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline bool hasSubBandMatrix(
            int i1, int i2, int j1, int j2, int lo, int hi, 
            int istep, int jstep) const
        {
            return const_type(*this).hasSubBandMatrix(
                i1,i2,j1,j2,lo,hi,istep,jstep); 
        }

        inline bool hasSubVector (
            int i, int j, int istep, int jstep, int n) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,n); }

        inline bool hasSubSymMatrix(int i1, int i2, int istep) const
        { return const_type(*this).hasSubSymMatrix(i1,i2,istep); }

        inline bool hasSubSymBandMatrix(int i1, int i2, int lo, int istep) const
        { return const_type(*this).hasSubSymMatrix(i1,i2,lo,istep); }

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

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int lo, int hi) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,lo,hi,1,1));
            return c_type::subBandMatrix(i1-1,i2,j1-1,j2,lo,hi);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(this->nlo()+j1-i1,i2-i1);
            const int newnhi = TMV_MIN(this->nlo()+i1-j1,j2-j1);
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int lo, int hi, 
            int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,lo,hi,istep,jstep));
            return c_type::subBandMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,lo,hi,istep,jstep);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return c_type::subVector(i-1,j-1,istep,jstep,n);
        }

        inline sym_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return c_type::subSymMatrix(i1-1,i2);
        }

        inline sym_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return c_type::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline view_type subSymBandMatrix(int i1, int i2, int lo) const
        {
            TMVAssert(hasSubSymBandMatrix(i1,i2,lo,1));
            return c_type::subSymBandMatrix(i1-1,i2,lo);
        }

        inline view_type subSymBandMatrix(int i1, int i2) const
        { return subSymBandMatrix(i1,i2,TMV_MIN(this->nlo(),i2-i1)); }

        inline sym_type subSymMatrix(int i1, int i2, int lo, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,lo,istep));
            return c_type::subSymMatrix(i1-1,i2-1+istep,lo,istep);
        }

        inline band_type diagRange(int k1, int k2) const
        { return c_type::diagRange(k1,k2); }

        inline view_type symDiagRange(int newnlo) const
        { return c_type::symDiagRange(newnlo); }

        inline band_type upperBand() const
        { return c_type::upperBand(); }

        inline band_type upperBandOff() const
        { return c_type::upperBandOff(); }

        inline band_type lowerBand() const
        { return c_type::lowerBand(); }

        inline band_type lowerBandOff() const
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

    }; // FortranStyle SymBandMatrixView

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    class SymBandMatrix : public GenSymBandMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymBandMatrix<T,U,S,I> type;
        typedef type copy_type;
        typedef GenSymBandMatrix<T> base;
        typedef ConstSymBandMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstSymMatrixView<T,I> const_sym_type;
        typedef ConstBandMatrixView<T,I> const_band_type;
        typedef ConstSymBandMatrixView<RT,I> const_realpart_type;
        typedef SymBandMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef SymMatrixView<T,I> sym_type;
        typedef BandMatrixView<T,I> band_type;
        typedef SymBandMatrixView<RT,I> realpart_type;
        typedef T& reference;

        //
        // Constructors
        //

#define NEW_SIZE(s,lo) \
        linsize(BandStorageLength(S,s,s,lo,0)), \
        itsm1(linsize), itss(s), itslo(lo), \
        itssi(S==DiagMajor ? -s+1 : S==RowMajor ? lo : 1), \
        itssj(S==DiagMajor ? s : S==RowMajor ? 1 : lo), \
        itssd(S==DiagMajor ? 1 : lo+1), \
        itsm((S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize) 

#define NEW_SIZE2(ls,s,lo) \
        linsize(ls), \
        itsm1(linsize), itss(s), itslo(lo), \
        itssi(S==DiagMajor ? -s+1 : S==RowMajor ? lo : 1), \
        itssj(S==DiagMajor ? s : S==RowMajor ? 1 : lo), \
        itssd(S==DiagMajor ? 1 : lo+1), \
        itsm((S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize) 

        inline SymBandMatrix() : NEW_SIZE(0,0) 
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
        }

        inline SymBandMatrix(int s, int lo) : NEW_SIZE(s,lo) 
        {
            TMVAssert(s >= 0);
            TMVAssert(lo < s);
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline SymBandMatrix(int s, int lo, const T& x) : NEW_SIZE(s,lo) 
        {
            TMVAssert(s >= 0);
            TMVAssert(lo < s);
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            setAllTo(x);
        }

        inline SymBandMatrix(const type& rhs) : 
            NEW_SIZE2(rhs.linsize,rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
        }

        template <UpLoType U2, StorageType S2, IndexStyle I2> 
        inline SymBandMatrix(const SymBandMatrix<T,U2,S2,I2>& rhs) : 
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if ( (U==U2 && S==S2) ||
                 (S!=DiagMajor && (S2==RowMajor || S2==ColMajor) && 
                  S!=S2 && U!=U2) )
                std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
            else if (U==Upper) upperBand() = rhs.upperBand();
            else lowerBand() = rhs.lowerBand();
        }

        inline SymBandMatrix(const GenSymBandMatrix<RT>& rhs) :
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (rhs.issym()) {
                rhs.assignTosB(view());
            } else {
                if (U==Upper) 
                    upperBand() = rhs.upperBand();
                else 
                    lowerBand() = rhs.lowerBand();
            }
        }

        inline SymBandMatrix(const GenSymBandMatrix<CT>& rhs) :
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (rhs.issym()) {
                rhs.assignTosB(view());
            } else {
                if (U==Upper) 
                    upperBand() = rhs.upperBand();
                else 
                    lowerBand() = rhs.lowerBand();
            }
        }

        template <class T2> 
        inline SymBandMatrix(const GenSymBandMatrix<T2>& rhs) :
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) upperBand() = rhs.upperBand();
            else lowerBand() = rhs.lowerBand();
        }

        template <class T2> 
        inline SymBandMatrix(const GenSymBandMatrix<T2>& rhs, int newnlo) : 
            NEW_SIZE(rhs.size(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(newnlo <= rhs.nlo());
            if (U==Upper) upperBand() = rhs.upperBand().diagRange(0,newnlo+1);
            else lowerBand() = rhs.lowerBand().diagRange(-newnlo,1);
        }

        template <class T2> 
        inline SymBandMatrix(const GenBandMatrix<T2>& rhs) :
            NEW_SIZE(rhs.rowsize(),U==Upper?rhs.nhi():rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(rhs.rowsize() == rhs.colsize());
            if (U==Upper) upperBand() = rhs.upperBand();
            else lowerBand() = rhs.lowerBand();
        }

        template <class T2> 
        inline SymBandMatrix(const GenMatrix<T2>& rhs, int newnlo) :
            NEW_SIZE(rhs.rowsize(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) upperBand() = BandMatrixViewOf(rhs,0,nlo());
            else lowerBand() = BandMatrixViewOf(rhs,nlo(),0);
        }

        template <class T2> 
        inline SymBandMatrix(const GenBandMatrix<T2>& rhs, int newnlo) :
            NEW_SIZE(rhs.rowsize(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) upperBand() = BandMatrixViewOf(rhs,0,nlo());
            else lowerBand() = BandMatrixViewOf(rhs,nlo(),0);
        }

        template <class T2> 
        inline SymBandMatrix(const GenSymMatrix<T2>& rhs, int newnlo) :
            NEW_SIZE(rhs.rowsize(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) 
                upperBand() = BandMatrixViewOf(rhs.upperTri(),nlo());
            else 
                lowerBand() = BandMatrixViewOf(rhs.lowerTri(),nlo());
        }

        template <class T2> 
        inline SymBandMatrix(const GenDiagMatrix<T2>& m2) :
            NEW_SIZE(m2.size(),0)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            setZero();
            DiagMatrixViewOf(diag()) = m2;
        }

        inline SymBandMatrix(const AssignableToSymBandMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(m2.issym());
            m2.assignTosB(view());
        }

        inline SymBandMatrix(const AssignableToSymBandMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.nlo())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(m2.issym());
            m2.assignTosB(view());
        }

        inline SymBandMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),0)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            m2.assignToD(DiagMatrixViewOf(diag()));
            if (nlo()>0) upperBandOff().setZero();
        }

        inline SymBandMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),0)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(isComplex(T()));
            m2.assignToD(DiagMatrixViewOf(diag()));
            if (nlo()>0) upperBandOff().setZero();
        }

#undef NEW_SIZE
#undef NEW_SIZE2

        virtual inline ~SymBandMatrix() 
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
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            if (&m2 != this) {
                if (nlo() > m2.nlo()) 
                    lowerBand().diagRange(-nlo(),-m2.nlo()).setZero();
                if (S==DiagMajor) {
                    if (nlo() > m2.nlo()) {
                        std::copy(
                            m2.start_mem(),m2.start_mem()+m2.mem_used(),
                            itsm+m2.nlo()*stepi());
                    } else {
                        std::copy(
                            m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
                    }
                } else if (nlo()==m2.nlo()) {
                    std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
                } else {
                    lowerBand().diagRange(-m2.nlo(),1) = m2.lowerBand();
                }
            }
            return *this;
        }

        template <IndexStyle I2> 
        inline type& operator=(const SymBandMatrix<T,U,S,I2>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            if (nlo() > m2.nlo()) lowerBand().diagRange(-nlo(),-m2.nlo()).setZero();
            if (S==DiagMajor) {
                if (nlo() > m2.nlo()) {
                    std::copy(
                        m2.start_mem(),m2.start_mem()+m2.mem_used(),
                        itsm+m2.nlo()*stepi());
                } else {
                    std::copy(
                        m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
                }
            } else if (nlo()==m2.nlo()) {
                std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
            } else {
                lowerBand().diagRange(-m2.nlo(),1) = m2.lowerBand();
            }
            return *this;
        }

        inline type& operator=(const GenSymBandMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(m2.issym());
            m2.assignTosB(view());
            return *this; 
        }

        inline type& operator=(const GenSymBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(m2.issym());
            m2.assignTosB(view());
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenSymBandMatrix<T2>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T2()) || m2.issym());
            lowerBand() = m2.lowerBand();
            return *this; 
        }

        inline type& operator=(const T& x) 
        { return setToIdentity(x); }

        inline type& operator=(const AssignableToSymBandMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(m2.issym());
            m2.assignTosB(view());
            return *this;
        }

        inline type& operator=(const AssignableToSymBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.issym());
            m2.assignTosB(view());
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2) 
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(issym());
            view() = m2;
            return *this;
        }

        //
        // Access
        //

        inline T operator()(int i, int j) const
        {
            if (I==CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                return cref(i,j);
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                return cref(i-1,j-1);
            }
        }

        inline T& operator()(int i, int j) 
        {
            if (I==CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                TMVAssert(okij(i,j));
                return ref(i,j); 
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                TMVAssert(okij(i-1,j-1));
                return ref(i-1,j-1); 
            }
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        {
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return const_vec_type(
                    itsm+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj); 
            else
                return const_vec_type(
                    itsm+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),NonConj); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I==FortranStyle) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return const_vec_type(
                    itsm+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj); 
            else 
                return const_vec_type(
                    itsm+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),NonConj); 
        }

        inline const_vec_type diag() const
        { return const_vec_type(
                itsm,size(),diagstep(),NonConj); }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            i = std::abs(i);
            return const_vec_type(
                itsm+i*(U==Upper?stepj():stepi()),
                size()-i,diagstep(),NonConj); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            i = std::abs(i);
            return const_vec_type(
                itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
                j2-j1,diagstep(),NonConj);
        }

        inline vec_type row(int i, int j1, int j2)
        {
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return vec_type(
                    itsm+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj TMV_FIRSTLAST); 
            else
                return vec_type(
                    itsm+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I==FortranStyle) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return vec_type(
                    itsm+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj TMV_FIRSTLAST); 
            else 
                return vec_type(
                    itsm+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag()
        { return vec_type(
                itsm,size(),diagstep(),NonConj TMV_FIRSTLAST); }

        inline vec_type diag(int i) 
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            i = std::abs(i);
            return vec_type(
                itsm+i*(U==Upper?stepj():stepi()),
                size()-i,diagstep(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i, int j1, int j2) 
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            i = std::abs(i);
            return vec_type(
                itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
                j2-j1,diagstep(),NonConj TMV_FIRSTLAST);
        }


        //
        // Modifying Functions
        //

        inline type& setZero() 
        {
            std::fill_n(itsm1.get(),linsize,T(0));
            return *this;
        }

        inline type& setAllTo(const T& x) 
        { upperBand().setAllTo(x); return *this; }

        inline type& addToAll(const T& x) 
        { upperBand().addToAll(x); return *this; }

        inline type& clip(RT thresh) 
        { upperBand().clip(thresh); return *this; }

        inline type& conjugateSelf() 
        { if (isComplex(T())) upperBand().conjugateSelf(); return *this; }

        inline type& transposeSelf() 
        { return *this; }

        inline type& setToIdentity(const T& x=T(1)) 
        { setZero(); diag().setAllTo(x); return *this; }

        //
        // SubMatrix
        //

        inline const_rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ( (U==Upper && i2-j1<=1) || (U==Lower && j2-i1<=1) )
                return const_rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
            else
                return const_rec_type(
                    itsm+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),TMV_TransOf(S),NonConj);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            if ( (U==Upper && i2-j1<=istep) || (U==Lower && j2-i1<=jstep) )
                return const_rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj);
            else
                return const_rec_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),NonConj);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
                return const_band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
                    S,NonConj);
            else
                return const_band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),NonConj);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==CStyle?1:0));
            const int newnhi = TMV_MIN(nlo()+i1-j1,j2-j1-(I==CStyle?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(view().hasSubBandMatrix(
                    i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            if ((i1+newnlo*istep<=j1 && U==Upper) || 
                (j1+newnhi*jstep<=i1 && U==Lower)) {
                const int newstepi = stepi()*istep;
                const int newstepj = stepj()*jstep;
                return const_band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,newstor,NonConj);
            } else {
                const int newstepi = stepj()*istep;
                const int newstepj = stepi()*jstep;
                return const_band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,
                    TMV_TransOf(newstor),NonConj);
            }
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==FortranStyle) { --i; --j; }
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return const_vec_type(
                    itsm+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj);
            else
                return const_vec_type(
                    itsm+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),NonConj);
        }

        inline const_sym_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return const_sym_type(
                itsm+i1*diagstep(),i2-i1,
                stepi(),stepj(),Sym,U,S,NonConj);
        }

        inline const_sym_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return const_sym_type(
                itsm+i1*diagstep(),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Sym,U,
                istep==1 ? S : NoMajor, NonConj);
        }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo) const
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,1));
            if (I==FortranStyle) { --i1; }
            return const_view_type(
                itsm+i1*diagstep(),i2-i1,
                newnlo,stepi(),stepj(),diagstep(),Sym,U,S,NonConj);
        }

        inline const_view_type subSymBandMatrix(int i1, int i2) const
        {
            return subSymBandMatrix(
                i1,i2,TMV_MIN(nlo(),i2-i1-(I==CStyle?1:0))); 
        }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep) const
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return const_view_type(
                itsm+i1*diagstep(), (i2-i1)/istep, newnlo,
                istep*stepi(),istep*stepj(),istep*diagstep(),
                Sym,U,istep==1 ? S : NoMajor,NonConj);
        }

        inline const_band_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
            TMVAssert(k2 <= 1 || k1 >= 0);

            if (k1 >= 0) {
                const int newsize = size()-k1;
                const int newnhi = k2-k1-1;
                if (uplo() == Upper)
                    return const_band_type(
                        itsm+k1*stepj(), newsize, newsize, 0, newnhi,
                        stepi(), stepj(), diagstep(), stor(), ct());
                else
                    return const_band_type(
                        itsm+k1*stepi(), newsize, newsize, 0, newnhi,
                        stepj(), stepi(), diagstep(),
                        TMV_TransOf(stor()), issym()?ct():TMV_ConjOf(T,ct()));
            } else {
                const int newsize = size()+k2-1;
                const int newnlo = k2-k1-1;
                if (uplo() == Lower)
                    return const_band_type(
                        itsm-k2*stepi(), newsize, newsize, newnlo, 0,
                        stepi(), stepj(), diagstep(), stor(), ct());
                else
                    return const_band_type(
                        itsm-k2*stepj(),
                        newsize, newsize, newnlo, 0, stepj(), stepi(),
                        diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()));
            }
        }

        inline const_view_type symDiagRange(int newnlo) const
        {
            TMVAssert(newnlo>=0 && newnlo <= nlo());
            return const_view_type(
                itsm,size(),newnlo,
                stepi(),stepj(),diagstep(),Sym,U,S,NonConj);
        }

        inline const_band_type upperBand() const
        {
            if (uplo() == Upper)
                return const_band_type(
                    itsm,size(),size(),0,nlo(),
                    stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm,size(),size(),0,nlo(),
                    stepj(),stepi(),diagstep(),TMV_TransOf(S),NonConj);
        }

        inline const_band_type upperBandOff() const
        {
            if (uplo() == Upper)
                return const_band_type(
                    itsm+stepj(),size()-1,size()-1,
                    0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm+stepi(),size()-1,size()-1,0,nlo()-1,
                    stepj(),stepi(),diagstep(),TMV_TransOf(S),NonConj);
        }

        inline const_band_type lowerBand() const
        {
            if (uplo() == Lower)
                return const_band_type(
                    itsm,size(),size(),nlo(),0,
                    stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm,size(),size(),nlo(),0,stepj(),stepi(),
                    diagstep(),TMV_TransOf(S),NonConj);
        }

        inline const_band_type lowerBandOff() const
        {
            if (uplo() == Lower)
                return const_band_type(
                    itsm+stepi(),size()-1,size()-1,
                    nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm+stepj(),size()-1,size()-1,nlo()-1,0,stepj(),stepi(),
                    diagstep(),TMV_TransOf(S),NonConj);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm),size(),nlo(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                Sym,U,isReal(T())?S:NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm)+1,size(),nlo(),
                2*stepi(),2*stepj(),2*diagstep(),Sym,U,NoMajor,NonConj);
        } 

        inline const_view_type view() const
        {
            return const_view_type(
                itsm,size(),nlo(),stepi(),stepj(),diagstep(),Sym,U,S,NonConj);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm,size(),nlo(),stepj(),stepi(),diagstep(),
                Sym,TMV_UTransOf(U), TMV_TransOf(S),NonConj);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                itsm,size(),nlo(), stepi(),stepj(),diagstep(),
                Sym,U,S,TMV_ConjOf(T,NonConj));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ( (U==Upper && i2-j1<=1) || (U==Lower && j2-i1<=1) )
                return rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm+i1*stepj()+j1*stepi(),i2-i1,j2-j1,
                    stepj(),stepi(),TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ( (U==Upper && i2-j1<=istep) || (U==Lower && j2-i1<=jstep) )
                return rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),NonConj TMV_FIRSTLAST);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
                return band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
                    S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline band_type subBandMatrix(int i1, int i2, int j1, int j2)
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==CStyle?1:0));
            const int newnhi = TMV_MIN(nlo()+i1-j1,j2-j1-(I==CStyle?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
                                              istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const StorageType newstor=
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            if ((i1+newnlo*istep<=j1 && U==Upper) || 
                (j1+newnhi*jstep<=i1 && U==Lower)) {
                const int newstepi = stepi()*istep;
                const int newstepj = stepj()*jstep;
                return band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,newstor,NonConj 
                    TMV_FIRSTLAST);
            } else {
                const int newstepi = stepj()*istep;
                const int newstepj = stepi()*jstep;
                return band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,
                    TMV_TransOf(newstor),NonConj TMV_FIRSTLAST);
            }
        }

        inline vec_type subVector(int i, int j, int istep, int jstep, int n)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==FortranStyle) { --i; --j; }
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return vec_type(
                    itsm+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),NonConj TMV_FIRSTLAST);
        }

        inline sym_type subSymMatrix(int i1, int i2)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return sym_type(
                itsm+i1*diagstep(),i2-i1,
                stepi(),stepj(),Sym,U,S,NonConj TMV_FIRSTLAST);
        }

        inline sym_type subSymMatrix(int i1, int i2, int istep) 
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return sym_type(
                itsm+i1*diagstep(),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),Sym,U,
                istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymBandMatrix(int i1, int i2, int newnlo)
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,1));
            if (I==FortranStyle) { --i1; }
            return view_type(
                itsm+i1*diagstep(),i2-i1,newnlo,
                stepi(),stepj(),diagstep(),Sym,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymBandMatrix(int i1, int i2)
        {
            return subSymBandMatrix(
                i1,i2,TMV_MIN(nlo(),i2-i1-(I==CStyle?1:0))); 
        }

        inline view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep)
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return view_type(
                itsm+i1*diagstep(),(i2-i1)/istep,
                newnlo,istep*stepi(),istep*stepj(),istep*diagstep(),
                Sym,U,istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline band_type diagRange(int k1, int k2) 
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
            TMVAssert(k2 <= 1 || k1 >= 0);

            if (k1 >= 0) {
                const int newsize = size()-k1;
                const int newnhi = k2-k1-1;
                if (uplo() == Upper)
                    return band_type(
                        itsm+k1*stepj(),newsize, newsize, 0, newnhi,
                        stepi(), stepj(), diagstep(),
                        stor(), ct() TMV_FIRSTLAST );
                else
                    return band_type(
                        itsm+k1*stepi(),newsize, newsize, 0, newnhi,
                        stepj(), stepi(), diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST );
            } else {
                const int newsize = size()+k2-1;
                const int newnlo = k2-k1-1;
                if (uplo() == Lower)
                    return band_type(
                        itsm-k2*stepi(), newsize, newsize, newnlo, 0,
                        stepi(), stepj(), diagstep(),
                        stor(), ct() TMV_FIRSTLAST );
                else
                    return band_type(
                        itsm-k2*stepj(), newsize, newsize, newnlo, 0,
                        stepj(), stepi(), diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST );
            }
        }

        inline view_type symDiagRange(int newnlo)
        {
            TMVAssert(newnlo>=0 && newnlo <= nlo());
            return view_type(
                itsm,size(),newnlo,stepi(),stepj(),diagstep(),
                Sym,U,S,NonConj TMV_FIRSTLAST);
        }

        inline band_type upperBand()
        {
            if (uplo() == Upper)
                return band_type(
                    itsm,size(),size(),0,nlo(), stepi(),stepj(),diagstep(),
                    S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm,size(),size(),0,nlo(), stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline band_type upperBandOff()
        {
            if (uplo() == Upper)
                return band_type(
                    itsm+stepj(),size()-1,size()-1,
                    0,nlo()-1,stepi(),stepj(),diagstep(),
                    S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm+stepi(),size()-1,size()-1,
                    0,nlo()-1,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline band_type lowerBand()
        {
            if (uplo() == Lower)
                return band_type(
                    itsm,size(),size(),nlo(),0,
                    stepi(),stepj(),diagstep(),S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm,size(),size(),nlo(),0,
                    stepj(),stepi(),diagstep(),TMV_TransOf(S),
                    NonConj TMV_FIRSTLAST);
        }

        inline band_type lowerBandOff()
        {
            if (uplo() == Lower)
                return band_type(
                    itsm+stepi(),size()-1,size()-1,
                    nlo()-1,0,stepi(),stepj(),diagstep(),
                    S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm+stepj(),size()-1,size()-1,
                    nlo()-1,0,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm),size(),nlo(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                Sym,U,isReal(T())?S:NoMajor,NonConj
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
                reinterpret_cast<RT*>(itsm)+1,size(),nlo(),
                2*stepi(),2*stepj(),2*diagstep(),Sym,U,NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        } 

        inline view_type view() 
        {
            return view_type(
                itsm,size(),nlo(), stepi(),stepj(),diagstep(),
                Sym,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type transpose() 
        {
            return view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline view_type conjugate() 
        {
            return view_type(
                itsm,size(),nlo(), stepi(),stepj(),diagstep(),
                Sym,U,S,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline view_type adjoint() 
        {
            return view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj) 
                TMV_FIRSTLAST);
        }

        //
        // I/O
        //
        
        void read(const TMV_Reader& reader);

        inline int size() const { return itss; }
        inline int nlo() const { return itslo; }
        inline int mem_used() const { return linsize; }
        inline const T* start_mem() const { return itsm1.get(); }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline int diagstep() const { return itssd; }
        inline SymType sym() const { return Sym; }
        inline UpLoType uplo() const { return U; }
        inline StorageType stor() const { return S; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isdm() const { return S==DiagMajor; }
        inline bool isconj() const { return false; }
        inline bool isherm() const { return isReal(T()); }
        inline bool issym() const { return true; }

        inline T& ref(int i, int j)
        {
            if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                return itsm[i*itssi + j*itssj];
            else 
                return itsm[j*itssi + i*itssj];
        }

        inline T cref(int i, int j) const 
        {
            if (okij(i,j)) {
                if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                    return itsm[i*itssi + j*itssj];
                else 
                    return itsm[j*itssi + i*itssj];
            } else {
                return T(0);
            }
        }

        inline void resize(int s, int lo)
        {
            TMVAssert(s >= 0);
            TMVAssert(lo < s);
            linsize = BandStorageLength(S,s,s,lo,0);
            itsm1.resize(linsize);
            itss = s;
            itslo = lo;
            itssi = S==DiagMajor ? -s+1 : S==RowMajor ? lo : 1;
            itssj = S==DiagMajor ? s : S==RowMajor ? 1 : lo;
            itssd = S==DiagMajor ? 1 : lo+1;
            itsm = (S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi :
                itsm1.get();
#ifdef TMVFLDEBUG
            _first = itsm1.get();
            _last = _first+linsize;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

    protected :

        int linsize;
        AlignedArray<T> itsm1;
        int itss;
        int itslo;
        int itssi;
        int itssj;
        int itssd;
        T* itsm;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
    protected :
#endif

        inline bool okij(int i, int j) const
        { return (j+nlo() >= i && i+nlo() >= j); }

        template <IndexStyle I2>
        friend void Swap(
            SymBandMatrix<T,U,S,I>& m1, SymBandMatrix<T,U,S,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.nlo() == m2.nlo());
            m1.itsm1.swapWith(m2.itsm1);
            TMV_SWAP(m1.itsm,m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // SymBandMatrix

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    class HermBandMatrix : 
        public GenSymBandMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef HermBandMatrix<T,U,S,I> type;
        typedef type copy_type;
        typedef GenSymBandMatrix<T> base;
        typedef ConstSymBandMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstSymMatrixView<T,I> const_sym_type;
        typedef ConstBandMatrixView<T,I> const_band_type;
        typedef ConstSymBandMatrixView<RT,I> const_realpart_type;
        typedef SymBandMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef SymMatrixView<T,I> sym_type;
        typedef BandMatrixView<T,I> band_type;
        typedef SymBandMatrixView<RT,I> realpart_type;
        typedef TMV_RefType(T) reference;

        //
        // Constructors
        //

#define NEW_SIZE(s,lo) \
        linsize(BandStorageLength(S,s,s,lo,0)), \
        itsm1(linsize), itss(s), itslo(lo), \
        itssi(S==DiagMajor ? -s+1 : S==RowMajor ? lo : 1), \
        itssj(S==DiagMajor ? s : S==RowMajor ? 1 : lo), \
        itssd(S==DiagMajor ? 1 : lo+1), \
        itsm((S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize) 

#define NEW_SIZE2(ls,s,lo) \
        linsize(ls), \
        itsm1(linsize), itss(s), itslo(lo), \
        itssi(S==DiagMajor ? -s+1 : S==RowMajor ? lo : 1), \
        itssj(S==DiagMajor ? s : S==RowMajor ? 1 : lo), \
        itssd(S==DiagMajor ? 1 : lo+1), \
        itsm((S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi : itsm1.get()) \
        TMV_DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize) 

        inline HermBandMatrix() : NEW_SIZE(0,0) 
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
        }

        inline HermBandMatrix(int s, int lo) : NEW_SIZE(s,lo) 
        {
            TMVAssert(s >= 0);
            TMVAssert(lo < s);
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#else
            if (isComplex(T())) diag().imagPart().setZero();
#endif
        }

        inline HermBandMatrix(int s, int lo, const RT& x) : 
            NEW_SIZE(s,lo) 
        {
            TMVAssert(s >= 0);
            TMVAssert(lo < s);
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            setAllTo(x);
        }

        inline HermBandMatrix(const type& rhs) : 
            NEW_SIZE2(rhs.linsize,rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
        }

        template <UpLoType U2, StorageType S2, IndexStyle I2> 
        inline HermBandMatrix(const HermBandMatrix<T,U2,S2,I2>& rhs) : 
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==U2 && S==S2)
                std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
            else if (S!=DiagMajor && (S2==RowMajor || S2==ColMajor) && 
                     S!=S2 && U!=U2) {
                std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
                conjugateSelf();
            }
            else if (U==Upper) upperBand() = rhs.upperBand();
            else lowerBand() = rhs.lowerBand();
        }

        inline HermBandMatrix(const GenSymBandMatrix<RT>& rhs) :
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (rhs.isherm()) {
                rhs.assignTosB(view());
            } else {
                if (U==Upper) upperBand() = rhs.upperBand();
                else lowerBand() = rhs.lowerBand();
            }
        }

        inline HermBandMatrix(const GenSymBandMatrix<CT>& rhs) :
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (rhs.isherm()) {
                rhs.assignTosB(view());
            } else {
                if (U==Upper) upperBand() = rhs.upperBand();
                else lowerBand() = rhs.lowerBand();
                if (isComplex(T())) diag().imagPart().setZero();
            }
        }

        template <class T2> 
        inline HermBandMatrix(const GenSymBandMatrix<T2>& rhs) :
            NEW_SIZE(rhs.size(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) 
                upperBand() = rhs.upperBand();
            else 
                lowerBand() = rhs.lowerBand();
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermBandMatrix(const GenSymBandMatrix<T2>& rhs, int newnlo) : 
            NEW_SIZE(rhs.size(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(newnlo <= rhs.nlo());
            if (U==Upper) upperBand() = rhs.upperBand().diagRange(0,newnlo+1);
            else lowerBand() = rhs.lowerBand().diagRange(-newnlo,1);
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermBandMatrix(const GenBandMatrix<T2>& rhs) :
            NEW_SIZE(rhs.rowsize(),rhs.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) upperBand() = rhs.upperBand();
            else lowerBand() = rhs.lowerBand();
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermBandMatrix(const GenMatrix<T2>& rhs, int newnlo) :
            NEW_SIZE(rhs.rowsize(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) upperBand() = BandMatrixViewOf(rhs,0,nlo());
            else lowerBand() = BandMatrixViewOf(rhs,nlo(),0);
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermBandMatrix(const GenBandMatrix<T2>& rhs, int newnlo) :
            NEW_SIZE(rhs.rowsize(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) upperBand() = BandMatrixViewOf(rhs,0,nlo());
            else lowerBand() = BandMatrixViewOf(rhs,nlo(),0);
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermBandMatrix(const GenSymMatrix<T2>& rhs, int newnlo) :
            NEW_SIZE(rhs.rowsize(),newnlo)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            if (U==Upper) 
                upperBand() = BandMatrixViewOf(rhs.upperTri(),nlo());
            else 
                lowerBand() = BandMatrixViewOf(rhs.lowerTri(),nlo());
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermBandMatrix(const GenDiagMatrix<T2>& m2) :
            NEW_SIZE(m2.size(),0)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            DiagMatrixViewOf(diag()) = m2;
            if (nlo()>0) upperBandOff().setZero();
            if (isComplex(T())) diag().imagPart().setZero();
        }

        inline HermBandMatrix(const AssignableToSymBandMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),m2.nlo())
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(m2.isherm());
            m2.assignTosB(view());
        }

        inline HermBandMatrix(const AssignableToSymBandMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),m2.nlo())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            TMVAssert(m2.isherm());
            m2.assignTosB(view());
        }

        inline HermBandMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size(),0)
        {
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            m2.assignToD(DiagMatrixViewOf(diag()));
            if (nlo()>0) upperBandOff().setZero();
        }

        inline HermBandMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size(),0)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
            m2.assignToD(DiagMatrixViewOf(diag()));
            if (isComplex(T())) diag().imagPart().setZero();
            if (nlo()>0) upperBandOff().setZero();
        }

#undef NEW_SIZE

        virtual inline ~HermBandMatrix() 
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
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            if (&m2 != this) {
                if (nlo() > m2.nlo()) 
                    lowerBand().diagRange(-nlo(),-m2.nlo()).setZero();
                if (S==DiagMajor) {
                    if (nlo() > m2.nlo()) {
                        std::copy(
                            m2.start_mem(),m2.start_mem()+m2.mem_used(),
                            itsm+m2.nlo()*stepi());
                    } else {
                        std::copy(
                            m2.start_mem(),m2.start_mem()+m2.mem_used(),
                            itsm1.get());
                    }
                } else if (nlo()==m2.nlo()) {
                    std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
                } else {
                    lowerBand().diagRange(-m2.nlo(),1) = m2.lowerBand();
                }
            }
            return *this;
        }

        template <IndexStyle I2> 
        inline type& operator=(const HermBandMatrix<T,U,S,I2>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            if (nlo() > m2.nlo()) 
                lowerBand().diagRange(-nlo(),-m2.nlo()).setZero();
            if (S==DiagMajor)
                if (nlo() > m2.nlo())
                    std::copy(
                        m2.start_mem(),m2.start_mem()+m2.mem_used(),
                        itsm+m2.nlo()*stepi());
                else
                    std::copy(
                        m2.start_mem(),m2.start_mem()+m2.mem_used(),
                        itsm1.get());
            else if (nlo()==m2.nlo())
                std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
            else
                lowerBand().diagRange(-m2.nlo(),1) = m2.lowerBand();
            return *this;
        }

        inline type& operator=(const GenSymBandMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(m2.isherm());
            m2.assignTosB(view());
            return *this; 
        }

        inline type& operator=(const GenSymBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(m2.isherm());
            m2.assignTosB(view());
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenSymBandMatrix<T2>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(m2.isherm());
            lowerBand() = m2.lowerBand();
            return *this; 
        }

        inline type& operator=(const T& x) 
        { return setToIdentity(x); }

        inline type& operator=(const AssignableToSymBandMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(m2.isherm());
            m2.assignTosB(view());
            return *this;
        }

        inline type& operator=(const AssignableToSymBandMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.isherm());
            m2.assignTosB(view());
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2) 
        {
            TMVAssert(size() == m2.size());
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2) 
        { 
            TMVAssert(size() == m2.size());
            view() = m2;
            TMVAssert(this->isHermOK());
            return *this;
        }

        //
        // Access
        //

        inline T operator()(int i, int j) const
        {
            if (I==CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                return cref(i,j);
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                return cref(i-1,j-1);
            }
        }

        inline reference operator()(int i, int j) 
        {
            if (I==CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                TMVAssert(okij(i,j));
                return ref(i,j); 
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                TMVAssert(okij(i-1,j-1));
                return ref(i-1,j-1); 
            }
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        {
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return const_vec_type(
                    itsm+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj); 
            else
                return const_vec_type(
                    itsm+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),TMV_ConjOf(T,NonConj)); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I==FortranStyle) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return const_vec_type(
                    itsm+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj); 
            else 
                return const_vec_type(
                    itsm+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),TMV_ConjOf(T,NonConj)); 
        }

        inline const_vec_type diag() const
        { return const_vec_type(itsm,size(),diagstep(),NonConj); }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj); 
            i = std::abs(i);
            return const_vec_type(
                itsm+i*(U==Upper?stepj():stepi()),
                size()-i,diagstep(),newct); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj); 
            i = std::abs(i);
            return const_vec_type(
                itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
                j2-j1,diagstep(),newct);
        }

        inline vec_type row(int i, int j1, int j2)
        {
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            TMVAssert(j1==j2 || okij(i,j2-1));
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return vec_type(
                    itsm+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj TMV_FIRSTLAST); 
            else
                return vec_type(
                    itsm+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I==FortranStyle) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            TMVAssert(i1==i2 || okij(i2-1,j));
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return vec_type(
                    itsm+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj TMV_FIRSTLAST); 
            else 
                return vec_type(
                    itsm+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST); 
        }

        inline vec_type diag()
        { return vec_type(itsm,size(),diagstep(),NonConj TMV_FIRSTLAST); }

        inline vec_type diag(int i) 
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj); 
            i = std::abs(i);
            return vec_type(
                itsm+i*(U==Upper?stepj():stepi()),
                size()-i,diagstep(),newct TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i, int j1, int j2) 
        {
            TMVAssert(i>=-nlo() && i<=nlo());
            TMVAssert(i>=-size() && i<=size()); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj); 
            i = std::abs(i);
            return vec_type(
                itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
                j2-j1,diagstep(),newct TMV_FIRSTLAST);
        }


        //
        // Modifying Functions
        //

        inline type& setZero() 
        {
            std::fill_n(itsm1.get(),linsize,T(0));
            return *this; 
        }

        inline type& setAllTo(const T& x) 
        {
            TMVAssert(TMV_IMAG(x) == RT(0));
            upperBand().setAllTo(x); 
            return *this; 
        }

        inline type& addToAll(const T& x) 
        {
            TMVAssert(TMV_IMAG(x) == RT(0));
            upperBand().addToAll(x); 
            return *this; 
        }

        inline type& clip(RT thresh) 
        { upperBand().clip(thresh); return *this; }

        inline type& conjugateSelf() 
        { 
            if (isComplex(T()) && nlo()>0) upperBandOff().conjugateSelf();
            return *this; 
        }

        inline type& transposeSelf() 
        { return conjugateSelf(); }

        inline type& setToIdentity(const T& x=T(1)) 
        {
            TMVAssert(TMV_IMAG(x) == RT(0));
            setZero(); 
            diag().setAllTo(x); 
            return *this; 
        }

        //
        // SubMatrix
        //

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ( (U==Upper && i2-j1<=1) || (U==Lower && j2-i1<=1) )
                return const_rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
            else
                return const_rec_type(
                    itsm+i1*stepj()+j1*stepi(),i2-i1,j2-j1,
                    stepj(),stepi(),TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            if ( (U==Upper && i2-j1<=istep) || (U==Lower && j2-i1<=jstep) )
                return const_rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj);
            else
                return const_rec_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),TMV_ConjOf(T,NonConj));
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
                return const_band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
                    S,NonConj);
            else
                return const_band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2) const
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==CStyle?1:0));
            const int newnhi = TMV_MIN(nlo()+i1-j1,j2-j1-(I==CStyle?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline const_band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep) const
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
                                              istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const StorageType newstor=
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            if ((i1+newnlo*istep<=j1 && U==Upper) || 
                (j1+newnhi*jstep<=i1 && U==Lower)) {
                const int newstepi = stepi()*istep;
                const int newstepj = stepj()*jstep;
                return const_band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,newstor,NonConj);
            } else {
                const int newstepi = stepj()*istep;
                const int newstepj = stepi()*jstep;
                return const_band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,
                    TMV_TransOf(newstor),TMV_ConjOf(T,NonConj));
            }
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==FortranStyle) { --i; --j; }
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return const_vec_type(
                    itsm+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj);
            else
                return const_vec_type(
                    itsm+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),TMV_ConjOf(T,NonConj));
        }

        inline const_sym_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return const_sym_type(
                itsm+i1*diagstep(),i2-i1,
                stepi(),stepj(),Herm,U,S,NonConj);
        }

        inline const_sym_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return const_sym_type(
                itsm+i1*diagstep(),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,U,
                istep==1 ? S : NoMajor, NonConj);
        }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo) const
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,1));
            if (I==FortranStyle) { --i1; }
            return const_view_type(
                itsm+i1*diagstep(),i2-i1,
                newnlo,stepi(),stepj(),diagstep(),Herm,U,S,NonConj);
        }

        inline const_view_type subSymBandMatrix(int i1, int i2) const
        {
            return subSymBandMatrix(
                i1,i2,TMV_MIN(nlo(),i2-i1-(I==CStyle?1:0))); 
        }

        inline const_view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep) const
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return const_view_type(
                itsm+i1*diagstep(), (i2-i1)/istep, newnlo,
                istep*stepi(),istep*stepj(),istep*diagstep(),
                Herm,U,istep==1 ? S : NoMajor,NonConj);
        }

        inline const_band_type diagRange(int k1, int k2) const
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
            TMVAssert(k2 <= 1 || k1 >= 0);

            if (k1 >= 0) {
                const int newsize = size()-k1;
                const int newnhi = k2-k1-1;
                if (uplo() == Upper)
                    return const_band_type(
                        itsm+k1*stepj(), newsize, newsize, 0, newnhi,
                        stepi(), stepj(), diagstep(), stor(), ct());
                else
                    return const_band_type(
                        itsm+k1*stepi(), newsize, newsize, 0, newnhi,
                        stepj(), stepi(), diagstep(),
                        TMV_TransOf(stor()), issym()?ct():TMV_ConjOf(T,ct()));
            } else {
                const int newsize = size()+k2-1;
                const int newnlo = k2-k1-1;
                if (uplo() == Lower)
                    return const_band_type(
                        itsm-k2*stepi(), newsize, newsize, newnlo, 0,
                        stepi(), stepj(), diagstep(), stor(), ct());
                else
                    return const_band_type(
                        itsm-k2*stepj(), newsize, newsize, newnlo, 0,
                        stepj(), stepi(), diagstep(),
                        TMV_TransOf(stor()), issym()?ct():TMV_ConjOf(T,ct()));
            }
        }

        inline const_view_type symDiagRange(int newnlo) const
        {
            TMVAssert(newnlo>=0 && newnlo <= nlo());
            return const_view_type(
                itsm,size(),newnlo,
                stepi(),stepj(),diagstep(),Herm,U,S,NonConj);
        }

        inline const_band_type upperBand() const
        {
            if (uplo() == Upper)
                return const_band_type(
                    itsm,size(),size(),0,nlo(),
                    stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm,size(),size(),0,nlo(), stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_band_type upperBandOff() const
        {
            if (uplo() == Upper)
                return const_band_type(
                    itsm+stepj(),size()-1,size()-1,
                    0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm+stepi(),size()-1,size()-1,
                    0,nlo()-1,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_band_type lowerBand() const
        {
            if (uplo() == Lower)
                return const_band_type(
                    itsm,size(),size(),nlo(),0,
                    stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm,size(),size(),nlo(),0, stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_band_type lowerBandOff() const
        {
            if (uplo() == Lower)
                return const_band_type(
                    itsm+stepi(),size()-1,size()-1,
                    nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj);
            else
                return const_band_type(
                    itsm+stepj(),size()-1,size()-1,
                    nlo()-1,0,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm),size(),nlo(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                Herm,U,isReal(T())?S:NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm)+1,size(),nlo(),
                2*stepi(),2*stepj(),2*diagstep(),Herm,U,NoMajor,NonConj);
        } 

        inline const_view_type view() const
        {
            return const_view_type(
                itsm,size(),nlo(),
                stepi(),stepj(),diagstep(),Herm,U,S,NonConj);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),NonConj);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                itsm,size(),nlo(), stepi(),stepj(),diagstep(),
                Herm,U,S,TMV_ConjOf(T,NonConj));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ( (U==Upper && i2-j1<=1) || (U==Lower && j2-i1<=1) )
                return rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm+i1*stepj()+j1*stepi(), i2-i1,j2-j1,stepj(),stepi(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ( (U==Upper && i2-j1<=istep) || (U==Lower && j2-i1<=jstep) )
                return rec_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
                return band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
                    S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline band_type subBandMatrix(int i1, int i2, int j1, int j2)
        {
            const int newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-(I==CStyle?1:0));
            const int newnhi = TMV_MIN(nlo()+i1-j1,j2-j1-(I==CStyle?1:0));
            return subBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        inline band_type subBandMatrix(
            int i1, int i2, int j1, int j2, int newnlo, int newnhi,
            int istep, int jstep)
        {
            TMVAssert(view().hasSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
                                              istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            const StorageType newstor=
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            if ((i1+newnlo*istep<=j1 && U==Upper) || 
                (j1+newnhi*jstep<=i1 && U==Lower)) {
                const int newstepi = stepi()*istep;
                const int newstepj = stepj()*jstep;
                return band_type(
                    itsm+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,newstor,NonConj 
                    TMV_FIRSTLAST);
            } else {
                const int newstepi = stepj()*istep;
                const int newstepj = stepi()*jstep;
                return band_type(
                    itsm+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
                    newstepi,newstepj,newstepi+newstepj,
                    TMV_TransOf(newstor),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
            }
        }

        inline vec_type subVector(int i, int j, int istep, int jstep, int n)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==FortranStyle) { --i; --j; }
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return vec_type(
                    itsm+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),TMV_ConjOf(T,NonConj) 
                    TMV_FIRSTLAST);
        }

        inline sym_type subSymMatrix(int i1, int i2)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return sym_type(
                itsm+i1*diagstep(),i2-i1,
                stepi(),stepj(),Herm,U,S,NonConj TMV_FIRSTLAST);
        }

        inline sym_type subSymMatrix(int i1, int i2, int istep) 
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return sym_type(
                itsm+i1*diagstep(),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),Herm,U,
                istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymBandMatrix(int i1, int i2, int newnlo)
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,1));
            if (I==FortranStyle) { --i1; }
            return view_type(
                itsm+i1*diagstep(),i2-i1,newnlo,
                stepi(),stepj(),diagstep(),Herm,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymBandMatrix(int i1, int i2)
        {
            return subSymBandMatrix(
                i1,i2,TMV_MIN(nlo(),i2-i1-(I==CStyle?1:0))); 
        }

        inline view_type subSymBandMatrix(
            int i1, int i2, int newnlo, int istep)
        {
            TMVAssert(view().hasSubSymBandMatrix(i1,i2,newnlo,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return view_type(
                itsm+i1*diagstep(),(i2-i1)/istep,
                newnlo,istep*stepi(),istep*stepj(),istep*diagstep(),
                Herm,U,istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline band_type diagRange(int k1, int k2)
        {
            TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
            TMVAssert(k2 <= 1 || k1 >= 0);

            if (k1 >= 0) {
                const int newsize = size()-k1;
                const int newnhi = k2-k1-1;
                if (uplo() == Upper)
                    return band_type(
                        itsm+k1*stepj(), newsize, newsize, 0, newnhi,
                        stepi(), stepj(), diagstep(),
                        stor(), ct() TMV_FIRSTLAST );
                else
                    return band_type(
                        itsm+k1*stepi(), newsize, newsize, 0, newnhi,
                        stepj(), stepi(), diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST );
            } else {
                const int newsize = size()+k2-1;
                const int newnlo = k2-k1-1;
                if (uplo() == Lower)
                    return band_type(
                        itsm-k2*stepi(), newsize, newsize, newnlo, 0,
                        stepi(), stepj(), diagstep(),
                        stor(), ct() TMV_FIRSTLAST );
                else
                    return band_type(
                        itsm-k2*stepj(), newsize, newsize, newnlo, 0,
                        stepj(), stepi(), diagstep(), TMV_TransOf(stor()),
                        issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST );
            }
        }

        inline view_type symDiagRange(int newnlo)
        {
            TMVAssert(newnlo>=0 && newnlo <= nlo());
            return view_type(
                itsm,size(),newnlo,
                stepi(),stepj(),diagstep(),Herm,U,S,NonConj TMV_FIRSTLAST);
        }

        inline band_type upperBand()
        {
            if (uplo() == Upper)
                return band_type(
                    itsm,size(),size(),0,nlo(),
                    stepi(),stepj(),diagstep(),S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm,size(),size(),0,nlo(), stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline band_type upperBandOff()
        {
            if (uplo() == Upper)
                return band_type(
                    itsm+stepj(),size()-1,size()-1,
                    0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj 
                    TMV_FIRSTLAST);
            else
                return band_type(
                    itsm+stepi(),size()-1,size()-1,
                    0,nlo()-1,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline band_type lowerBand()
        {
            if (uplo() == Lower)
                return band_type(
                    itsm,size(),size(),nlo(),0,
                    stepi(),stepj(),diagstep(),S,NonConj TMV_FIRSTLAST);
            else
                return band_type(
                    itsm,size(),size(),nlo(),0, stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline band_type lowerBandOff()
        {
            if (uplo() == Lower)
                return band_type(
                    itsm+stepi(),size()-1,size()-1,
                    nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj 
                    TMV_FIRSTLAST);
            else
                return band_type(
                    itsm+stepj(),size()-1,size()-1,
                    nlo()-1,0,stepj(),stepi(),diagstep(),
                    TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm),size(),nlo(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? diagstep() : 2*diagstep(),
                Herm,U,isReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
            // The imaginary part of a Hermitian matrix is anti-symmetric
            // so this is illegal.
        {
            TMVAssert(TMV_FALSE);
            return realpart_type(
                0,0,0,0,0,0,Herm,U,S,NonConj TMV_FIRSTLAST1(0,0) );
        } 

        inline view_type view() 
        {
            return view_type(
                itsm,size(),nlo(), stepi(),stepj(),diagstep(),
                Herm,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type transpose() 
        {
            return view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline view_type conjugate() 
        {
            return view_type(
                itsm,size(),nlo(), stepi(),stepj(),diagstep(),
                Herm,U,S,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline view_type adjoint() 
        {
            return view_type(
                itsm,size(),nlo(), stepj(),stepi(),diagstep(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj) 
                TMV_FIRSTLAST);
        }

        //
        // I/O
        //
        
        void read(const TMV_Reader& reader);

        inline int size() const { return itss; }
        inline int nlo() const { return itslo; }
        inline int mem_used() const { return linsize; }
        inline const T* start_mem() const { return itsm1.get(); }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline int diagstep() const { return itssd; }
        inline SymType sym() const { return Herm; }
        inline UpLoType uplo() const { return U; }
        inline StorageType stor() const { return S; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isdm() const { return S==DiagMajor; }
        inline bool isconj() const { return false; }
        inline bool isherm() const { return true; }
        inline bool issym() const { return isReal(T()); }

        inline reference ref(int i, int j)
        {
            if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                return TMV_REF(itsm + i*stepi() + j*stepj(),NonConj);
            else 
                return TMV_REF(itsm + j*stepi() + i*stepj(),Conj);
        }

        inline T cref(int i, int j) const 
        {
            if (okij(i,j)) {
                if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                    return itsm[i*itssi + j*itssj];
                else 
                    return TMV_CONJ(itsm[j*itssi + i*itssj]);
            } else {
                return T(0);
            }
        }

        inline void resize(int s, int lo)
        {
            TMVAssert(s >= 0);
            TMVAssert(lo < s);
            linsize = BandStorageLength(S,s,s,lo,0);
            itsm1.resize(linsize);
            itss = s;
            itslo = lo;
            itssi = S==DiagMajor ? -s+1 : S==RowMajor ? lo : 1;
            itssj = S==DiagMajor ? s : S==RowMajor ? 1 : lo;
            itssd = S==DiagMajor ? 1 : lo+1;
            itsm = (S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi :
                itsm1.get();
#ifdef TMVFLDEBUG
            _first = itsm1.get();
            _last = _first+linsize;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#else
            if (isComplex(T())) diag().imagPart().setZero();
#endif
        }

    protected :

        int linsize;
        AlignedArray<T> itsm1;
        int itss;
        int itslo;
        int itssi;
        int itssj;
        int itssd;
        T* itsm;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
    protected :
#endif

        inline bool okij(int i, int j) const
        { return (j+nlo() >= i && i+nlo() >= j); }

        template <IndexStyle I2>
        friend void Swap(
            HermBandMatrix<T,U,S,I>& m1, HermBandMatrix<T,U,S,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.nlo() == m2.nlo());
            m1.itsm1.swapWith(m2.itsm1);
            TMV_SWAP(m1.itsm,m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // HermBandMatrix

    //---------------------------------------------------------------------------

    //
    // Special Creators: 
    //   SymTriDiagMatrix(v1,v2)
    //       Returns a SymBandMatrix with v1 as the diagonal, and v2 as
    //       the offdiagonal.
    //   HermTriDiagMatrix(v1,v2,uplo)
    //       Returns a HermBandMatrix with v1 as the diagonal, and v2 as
    //       the upper or lower offdiagonal, according to the value of uplo.
    //   SymBandMatrixViewOf(m,uplo,nlo)
    //   HermBandMatrixViewOf(m,uplo,nlo)
    //   SymBandMatrixViewOf(mptr,size,nlo,uplo,stor)
    //   HermBandMatrixViewOf(mptr,size,nlo,uplo,stor)
    //   SymBandMatrixViewOf(mptr,size,nlo,uplo,stepi,stepj)
    //   HermBandMatrixViewOf(mptr,size,nlo,uplo,stepi,stepj)
    //

    template <class T> 
    SymBandMatrix<T,Upper,DiagMajor> SymTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2);
    template <class T> 
    HermBandMatrix<T,Upper,DiagMajor> HermTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2, 
        UpLoType uplo=tmv::Upper);
    template <class T> 
    HermBandMatrix<std::complex<T>,Upper,DiagMajor> HermTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<std::complex<T> >& v2,
        UpLoType uplo);

    // From Matrix
    template <class T> 
    inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
        const GenMatrix<T>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymBandMatrixView<T>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const ConstMatrixView<T,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const Matrix<T,S,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        const MatrixView<T,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        Matrix<T,S,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const GenMatrix<T>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return ConstSymBandMatrixView<T>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const ConstMatrixView<T,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
        const Matrix<T,S,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T> HermBandMatrixViewOf(
        const MatrixView<T,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
        Matrix<T,S,I>& m, UpLoType uplo, int nlo)
    {
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    // From BandMatrix
    template <class T> 
    inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
        const GenBandMatrix<T>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymBandMatrixView<T>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const ConstBandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        const BandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const GenBandMatrix<T>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return ConstSymBandMatrixView<T>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const ConstBandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
        const BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T> HermBandMatrixViewOf(
        const BandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
        BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
    {
        if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
        TMVAssert(m.colsize()==m.rowsize());
#ifdef TMV_EXTRA_DEBUG
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
#endif
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
            TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    // From SymMatrix
    template <class T> 
    inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
        const GenSymMatrix<T>& m, int nlo)
    {
        return ConstSymBandMatrixView<T>(
            m.cptr(),m.size(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const ConstSymMatrixView<T,I>& m, int nlo)
    {
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const SymMatrix<T,U,S,I>& m, int nlo)
    {
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
        const HermMatrix<T,U,S,I>& m, int nlo)
    {
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        const SymMatrixView<T,I>& m, int nlo)
    {
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        SymMatrix<T,U,S,I>& m, int nlo)
    {
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
        HermMatrix<T,U,S,I>& m, int nlo)
    {
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const GenSymMatrix<T>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return ConstSymBandMatrixView<T>(
            m.cptr(),m.size(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
        const ConstSymMatrixView<T,I>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
        const SymMatrix<T,U,S,I>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
        const HermMatrix<T,U,S,I>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return ConstSymBandMatrixView<T,I>(
            m.cptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
        const SymMatrixView<T,I>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
        SymMatrix<T,U,S,I>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
        HermMatrix<T,U,S,I>& m, int nlo)
    {
        TMVAssert(m.isherm());
        return SymBandMatrixView<T,I>(
            m.ptr(),m.colsize(),nlo,
            m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
            m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last) ); 
    }

    // From ptr
    template <class T> 
    inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
        const T* m, int size, int nlo, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const int stepi = 
            stor == RowMajor ? nlo :
            stor == ColMajor ? 1 :
            -size + 1;
        const int stepj = 
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo :
            size;
        const T* m0 = (stor == DiagMajor && uplo == Lower) ? m-nlo*stepi : m;
        return ConstSymBandMatrixView<T>(
            m0,size,nlo,stepi,stepj,stepi+stepj,Sym,uplo,stor,NonConj);
    }

    template <class T> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const T* m, int size, int nlo, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const int stepi = 
            stor == RowMajor ? nlo :
            stor == ColMajor ? 1 :
            -size + 1;
        const int stepj = 
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo :
            size;
        const T* m0 = (stor == DiagMajor && uplo == Lower) ? m-nlo*stepi : m;
        return ConstSymBandMatrixView<T>(
            m0,size,nlo,stepi,stepj,stepi+stepj,Herm,uplo,stor,NonConj);
    }

    template <class T> 
    inline SymBandMatrixView<T> SymBandMatrixViewOf(
        T* m, int size, int nlo, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const int stepi = 
            stor == RowMajor ? nlo :
            stor == ColMajor ? 1 :
            -size + 1;
        const int stepj = 
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo :
            size;
        T* m0 = (stor == DiagMajor && uplo == Lower) ? m-nlo*stepi : m;
        return SymBandMatrixView<T>(
            m0,size,nlo,stepi,stepj,stepi+stepj,Sym,uplo,stor,NonConj
            TMV_FIRSTLAST1(m,m+BandStorageLength(stor,size,size,nlo,0)));
    }

    template <class T> 
    inline SymBandMatrixView<T> HermBandMatrixViewOf(
        T* m, int size, int nlo, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const int stepi = 
            stor == RowMajor ? nlo :
            stor == ColMajor ? 1 :
            -size + 1;
        const int stepj = 
            stor == RowMajor ? 1 :
            stor == ColMajor ? nlo :
            size;
        T* m0 = (stor == DiagMajor && uplo == Lower) ? m-nlo*stepi : m;
        return SymBandMatrixView<T>(
            m0,size,nlo,stepi,stepj,stepi+stepj,Herm,uplo,stor,NonConj
            TMV_FIRSTLAST1(m,m+BandStorageLength(stor,size,size,nlo,0)));
    }

    template <class T> 
    inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
        const T* m, int size, int nlo, UpLoType uplo, int stepi, int stepj)
    {
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const StorageType stor = 
            stepi == 1 ? ColMajor :
            stepj == 1 ? RowMajor :
            stepi+stepj == 1 ? DiagMajor :
            NoMajor;
        return ConstSymBandMatrixView<T>(
            m,size,nlo,stepi,stepj,stepi+stepj,Sym,uplo,stor,NonConj);
    }

    template <class T> 
    inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
        const T* m, int size, int nlo, UpLoType uplo, int stepi, int stepj)
    {
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const StorageType stor = 
            stepi == 1 ? ColMajor :
            stepj == 1 ? RowMajor :
            stepi+stepj == 1 ? DiagMajor :
            NoMajor;
        return SymBandMatrixView<T>(
            m,size,nlo,stepi,stepj,stepi+stepj,Herm,uplo,stor,NonConj);
    }

    template <class T> 
    inline SymBandMatrixView<T> SymBandMatrixViewOf(
        T* m, int size, int nlo, UpLoType uplo, int stepi, int stepj)
    {
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const StorageType stor = 
            stepi == 1 ? ColMajor :
            stepj == 1 ? RowMajor :
            stepi+stepj == 1 ? DiagMajor :
            NoMajor;
        return SymBandMatrixView<T>(
            m,size,nlo,stepi,stepj,stepi+stepj,Sym,uplo,stor,NonConj
            TMV_FIRSTLAST1(
                m + (stepi < 0 && uplo==Lower ? (stepi*nlo) : 
                     stepj < 0 && uplo==Upper ? (stepj*nhi) : 0),
                m + ((stepi < 0 && uplo==Lower ? (stepi*nlo) : 
                      stepj < 0 && uplo==Upper ? (stepj*nhi) : 0) +
                     BandStorageLength(stor,size,size,nlo,0))));
    }

    template <class T> 
    inline SymBandMatrixView<T> HermBandMatrixViewOf(
        T* m, int size, int nlo, UpLoType uplo, int stepi, int stepj)
    {
        TMVAssert(size>=0);
        TMVAssert(nlo<size);
        const StorageType stor = 
            stepi == 1 ? ColMajor :
            stepj == 1 ? RowMajor :
            stepi+stepj == 1 ? DiagMajor :
            NoMajor;
        return SymBandMatrixView<T>(
            m,size,nlo,stepi,stepj,stepi+stepj,Herm,uplo,stor,NonConj
            TMV_FIRSTLAST1(
                m + (stepi < 0 && uplo==Lower ? (stepi*nlo) : 
                     stepj < 0 && uplo==Upper ? (stepj*nhi) : 0),
                m + ((stepi < 0 && uplo==Lower ? (stepi*nlo) : 
                      stepj < 0 && uplo==Upper ? (stepj*nhi) : 0) +
                     BandStorageLength(stor,size,size,nlo,0))));
    }


    //
    // Swap Matrices
    //

    template <class T> 
    inline void Swap(
        const SymBandMatrixView<T>& m1, const SymBandMatrixView<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.nlo() == m2.nlo());
        TMVAssert(m1.issym() == m2.issym());
        TMVAssert(m1.isherm() == m2.isherm());
        Swap(m1.upperBand(),m2.upperBand()); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(
        const SymBandMatrixView<T>& m1, SymBandMatrix<T,U,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(
        SymBandMatrix<T,U,S,I>& m1, const SymBandMatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
    inline void Swap(
        SymBandMatrix<T,U1,S1,I1>& m1, SymBandMatrix<T,U2,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(
        const SymBandMatrixView<T>& m1, HermBandMatrix<T,U,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(
        HermBandMatrix<T,U,S,I>& m1, const SymBandMatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
    inline void Swap(
        HermBandMatrix<T,U1,S1,I1>& m1, HermBandMatrix<T,U2,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }



    //
    // Views:
    //

    template <class T> 
    inline ConstSymBandMatrixView<T> Transpose(const GenSymBandMatrix<T>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> Transpose(
        const ConstSymBandMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> Transpose(
        const SymBandMatrix<T,U,S,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> Transpose(const SymBandMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> Transpose(SymBandMatrix<T,U,S,I>& m)
    { return m.transpose(); }

    template <class T> 
    inline ConstSymBandMatrixView<T> Conjugate(const GenSymBandMatrix<T>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> Conjugate(
        const ConstSymBandMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> Conjugate(
        const SymBandMatrix<T,U,S,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> Conjugate(const SymBandMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> Conjugate(SymBandMatrix<T,U,S,I>& m)
    { return m.conjugate(); }

    template <class T> 
    inline ConstSymBandMatrixView<T> Adjoint(const GenSymBandMatrix<T>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> Adjoint(
        const ConstSymBandMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymBandMatrixView<T,I> Adjoint(const SymBandMatrix<T,U,S,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline SymBandMatrixView<T,I> Adjoint(const SymBandMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymBandMatrixView<T,I> Adjoint(SymBandMatrix<T,U,S,I>& m)
    { return m.adjoint(); }

    template <class T> 
    inline QuotXsB<T,T> Inverse(const GenSymBandMatrix<T>& m)
    { return m.inverse(); }


    //
    // SymBandMatrix ==, != SymBandMatrix
    //

    template <class T1, class T2> 
    inline bool operator==(
        const GenSymBandMatrix<T1>& m1, const GenSymBandMatrix<T2>& m2)
    { return m1.upperBand() == m2.upperBand(); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenSymBandMatrix<T1>& m1, const GenSymBandMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator==(
        const GenSymBandMatrix<T1>& m1, const GenMatrix<T2>& m2)
    {
        return 
            m1.upperBand() == BandMatrixViewOf(m2.upperTri()) &&
            m1.lowerBand() == BandMatrixViewOf(m2.lowerTri());
    }

    template <class T1, class T2> 
    inline bool operator==(
        const GenMatrix<T1>& m1, const GenSymBandMatrix<T2>& m2)
    { return m2 == m1; }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenSymBandMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenSymBandMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T>
    std::istream& operator>>(std::istream& is, const SymBandMatrixView<T>& m)
    { return is >> IOStyle() >> m; }

    template <class T, UpLoType U, StorageType S, IndexStyle I>
    std::istream& operator>>(std::istream& is, SymBandMatrix<T,U,S,I>& m)
    { return is >> IOStyle() >> m; }

    template <class T, UpLoType U, StorageType S, IndexStyle I>
    std::istream& operator>>(std::istream& is, HermBandMatrix<T,U,S,I>& m)
    { return is >> IOStyle() >> m; }

    template <class T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, const SymBandMatrixView<T>& m)
    { m.read(reader); return reader.getis(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::istream& operator>>(
        const TMV_Reader& reader, SymBandMatrix<T,U,S,I>& m)
    { m.read(reader); return reader.getis(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::istream& operator>>(
        const TMV_Reader& reader, HermBandMatrix<T,U,S,I>& m)
    { m.read(reader); return reader.getis(); }

} // namespace tmv

#endif
