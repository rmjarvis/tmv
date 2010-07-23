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
// This file defines the TMV SymMatrix and HermMatrix classes.
//
// Constructors:
//
//    SymMatrix is used for symmetric matrices, and HermMatrix is used
//    for Hermitian matrices.  For real matrices, these are the same thing:
//    A = A.transpose().
//    But for complex, they are different:
//    A_sym = A_sym.transpose()
//    A_herm = A_herm.adjoint()
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
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row(int i, int j1, int j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row
//        and be entirely in the upper or lower triangle.
//
//    Vector& col(int j, int i1, int i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column
//        and be entirely in the upper or lower triangle.
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
//        For HermMatrix, x must be real.
//    clip(RT thresh)
//    conjugateSelf()
//    transposeSelf()
//    setToIdentity(x = 1)
//    swapRowsCols(i1,i2)
//        Equivalent to swapping rows i1,i2 then swapping cols i1,i2.
//    permuteRowsCols(const int* p)
//    reversePermuteRowsCols(const int* p)
//        Perform a series of row/col swaps.
//    void Swap(SymMatrix& m1, SymMatrix& m2)
//        The SymMatrices must be the same size and Hermitianity.
//
// Views of a SymMatrix:
//
//    subMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the upper 
//        or lower triangle.
//
//    subVector(int i, int j, int istep, int jstep, int size)
//        Returns a VectorView which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.  The subvector must be completely with the upper or
//        lower triangle.
//
//    subSymMatrix(int i1, int i2, int istep)
//        Returns the SymMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with a SymMatrix of size 10, the x's below
//        are the original data, the O's are the subSymMatrix returned
//        with the command subSymMatrix(3,11,2), and the #'s are the 
//        subSymMatrix returned with subSymMatrix(0,3)
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
//    m.divideUsing(dt)
//    where dt is LU, CH, or SV
//     
//    m.lud(), m.chd(), m.svd(), and m.symsvd() return the
//        corresponding Divider classes.  
//
//    For SymMatrixes, LU actually does an LDLT decomposition.
//        ie. L(Lower Triangle) * Diagonal * L.transpose()
//        For non-positive definite matrices, D may have 2x2 blocks along
//        the diagonal, which can then be further decomposed via normal LU
//        decomposition.
//        
//    The option unique to hermitian matrixes is CH - Cholskey decomposition.
//        This is only appropriate if you know that your HermMatrix is 
//        positive definite (ie. all eigenvalues are positive).
//        (This is guaranteed, for example, if all the square 
//        submatrices have positive determinant.)
//        In this case, the SymMatrix can be decomposed into L*L.adjoint().
//
//    Finally, a difference for the SV Decomposition for Hermitian matrices
//        is that the decomposition can be done with V = Ut.  In this case,
//        the "singular values" are really the eigenvalues of A.
//        The only caveat is that they may be negative, whereas the usual
//        definition of singular values is that they are all positive.
//        Requiring positivity would destroy V = Ut, so it seemed more
//        useful to leave them as the actual eigenvalues.  Just keep that
//        in mind if you use the singular values for anything that expects
//        them to be positive.
//
//    If m is complex, symmetric (i.e. not hermitian), then the SVD
//        does not result in an eigenvalue decomposition, and we actually
//        just do a regular SVD.  In this case, the access is through
//        symsvd(), rather than the usual svd() method.


#ifndef TMV_SymMatrix_H
#define TMV_SymMatrix_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Array.h"
#include <vector>

namespace tmv {

    template <class T, class T1, class T2> 
    class OProdVV;
    template <class T, class T1, class T2> 
    class ProdMM;
    template <class T, class T1, class T2> 
    class ProdUL;
    template <class T, class T1, class T2> 
    class ProdLU;

    template <class T>
    struct SymCopyHelper // real
    { typedef SymMatrix<T> type; };
    template <class T>
    struct SymCopyHelper<std::complex<T> > // complex
    // Have to copy to matrix, since don't know whether herm or sym.
    { typedef Matrix<std::complex<T> > type; };

    template <class T> 
    class GenSymMatrix : 
        virtual public AssignableToSymMatrix<T>,
        public BaseMatrix<T>,
        private DivHelper<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenSymMatrix<T> type;
        typedef typename SymCopyHelper<T>::type copy_type;
        typedef ConstSymMatrixView<T> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef ConstMatrixView<T> const_rec_type;
        typedef ConstUpperTriMatrixView<T> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T> const_lowertri_type;
        typedef ConstSymMatrixView<RT> const_realpart_type;
        typedef SymMatrixView<T> nonconst_type;

        //
        // Constructors
        //

        inline GenSymMatrix() {}
        inline GenSymMatrix(const type&) {}
        virtual inline ~GenSymMatrix() {}

        //
        // Access Functions
        //

        using AssignableToSymMatrix<T>::size;
        inline size_t colsize() const { return size(); }
        inline size_t rowsize() const { return size(); }
        using AssignableToSymMatrix<T>::sym;

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>=0 && i<int(size()));
            TMVAssert(j>=0 && j<int(size()));
            return cref(i,j);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
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
            TMVAssert(j>=0 && j<int(size()));
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
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
        { return const_vec_type(
                cptr(),size(),stepi()+stepj(),ct()); }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            if (i>=0)
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()+i*stepj(),size()-i,
                        stepi()+stepj(),ct()); 
                else
                    return const_vec_type(
                        cptr()+i*stepi(),size()-i,
                        stepi()+stepj(),issym()?ct():TMV_ConjOf(T,ct()));
            else
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()-i*stepj(),size()+i,
                        stepi()+stepj(),issym()?ct():TMV_ConjOf(T,ct()));
                else
                    return const_vec_type(
                        cptr()-i*stepi(),size()+i,
                        stepi()+stepj(),ct()); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
            const int ds = stepi()+stepj();
            if (i>=0)
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()+i*stepj()+j1*ds,j2-j1,ds,ct()); 
                else
                    return const_vec_type(
                        cptr()+i*stepi()+j1*ds,j2-j1,ds,
                        issym()?ct():TMV_ConjOf(T,ct()));
            else
                if (uplo()==Upper)
                    return const_vec_type(
                        cptr()-i*stepj()+j1*ds,j2-j1,ds,
                        issym()?ct():TMV_ConjOf(T,ct()));
                else
                    return const_vec_type(
                        cptr()-i*stepi()+j1*ds,j2-j1,ds,ct()); 
        }

        template <class T2> 
        inline bool isSameAs(const BaseMatrix<T2>& ) const
        { return false; }

        inline bool isSameAs(const type& m2) const
        { 
            if (this == &m2) return true;
            else if (cptr()==m2.cptr() && size()==m2.size() && 
                     (isReal(T()) || sym() == m2.sym())) {
                if (uplo() == m2.uplo())
                    return (stepi() == m2.stepi() && stepj() == m2.stepj()
                            && ct() == m2.ct());
                else
                    return (stepi() == m2.stepj() && stepj() == m2.stepi()
                            && issym() == (ct()==m2.ct()));
            } else return false;
        }

        template <class T2>
        TMV_DEPRECATED(bool SameAs(const BaseMatrix<T2>& m2) const);
        TMV_DEPRECATED(bool SameAs(const type& m2) const)
        { return isSameAs(m2); }

        inline void assignToM(const MatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            assignToS(SymMatrixViewOf(m2,Upper));
            if (size() > 0)
                m2.lowerTri().offDiag() = 
                    m2.upperTri().offDiag().transpose();
        }

        inline void assignToM(const MatrixView<CT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            if (issym()) {
                assignToS(SymMatrixViewOf(m2,Upper));
                if (size() > 0)
                    m2.lowerTri().offDiag() = 
                        m2.upperTri().offDiag().transpose();
            } else {
                m2.diag().imagPart().setZero();
                assignToS(HermMatrixViewOf(m2,Upper));
                if (size() > 0)
                    m2.lowerTri().offDiag() = 
                        m2.upperTri().offDiag().adjoint();
            }
        }

        inline void assignToS(const SymMatrixView<RT>& m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == size());
            if (!isSameAs(m2)) m2.upperTri() = upperTri(); 
        }

        inline void assignToS(const SymMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(isReal(T()) || m2.sym() == sym());
            if (!isSameAs(m2)) m2.upperTri() = upperTri(); 
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
            if ((i2-1<=j1 && uplo()==Upper) || (j2-1<=i1 && uplo()==Lower))
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
            if ((i2-istep<=j1 && uplo()==Upper) || 
                (j2-jstep<=i1 && uplo()==Lower))
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

        bool hasSubVector(
            int i, int j, int istep, int jstep, int n) const;

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
                    cptr()+i*stepj()+j*stepi(),n,istep*stepj()+jstep*stepi(),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        bool hasSubSymMatrix(int i1, int i2, int istep) const;

        inline const_view_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return const_view_type(
                cptr()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),sym(),uplo(),stor(),ct());
        }

        inline const_view_type subSymMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return const_view_type(
                cptr()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),sym(),uplo(),
                istep==1 ? stor() : NoMajor,ct());
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            if (uplo() == Upper)
                return const_uppertri_type(
                    cptr(),size(),
                    stepi(),stepj(),dt,stor(),ct());
            else
                return const_uppertri_type(
                    cptr(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            if (uplo() == Lower)
                return const_lowertri_type(
                    cptr(),size(),
                    stepi(),stepj(),dt,stor(),ct());
            else
                return const_lowertri_type(
                    cptr(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(stor()),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                Sym, uplo(), isReal(T()) ? stor() : NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(issym());
            // The imaginary part of a Hermitian matrix is anti-symmetric
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                size(),2*stepi(),2*stepj(), Sym,uplo(),NoMajor,NonConj);
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
        TMV_DEPRECATED(const_view_type SubSymMatrix(int i1, int i2) const)
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubSymMatrix(
                int i1, int i2, int istep) const)
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_uppertri_type UpperTri(
                DiagType dt = NonUnitDiag) const)
        { return upperTri(); }
        TMV_DEPRECATED(const_lowertri_type LowerTri(
                DiagType dt = NonUnitDiag) const)
        { return lowerTri(); }
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
                cptr(),size(),stepi(),stepj(),
                sym(),uplo(),stor(),ct());
        }

        inline const_view_type transpose() const
        { 
            return const_view_type(
                cptr(),size(),stepj(),stepi(),
                sym(),TMV_UTransOf(uplo()),TMV_TransOf(stor()),ct());
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                cptr(),size(),stepi(),stepj(),
                sym(),uplo(),stor(),TMV_ConjOf(T,ct()));
        }

        inline const_view_type adjoint() const
        { 
            return const_view_type(
                cptr(),size(),stepj(),stepi(),
                sym(),TMV_UTransOf(uplo()),TMV_TransOf(stor()),
                TMV_ConjOf(T,ct()));
        }

        inline nonconst_type nonConst() const
        {
            return SymMatrixView<T>(
                const_cast<T*>(cptr()),size(),
                stepi(),stepj(),sym(),uplo(),stor(),ct()
                TMV_FIRSTLAST1(
                    cptr(),row(colsize()-1,0,colsize()).end().getP()));
        }

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


        //
        // Functions of Matrix
        //

        inline T det() const 
        { return DivHelper<T>::det(); }

        inline RT logDet(T* sign=0) const
        { return DivHelper<T>::logDet(sign); }

        inline T trace() const
        { return diag().sumElements(); }

        T sumElements() const;

        RT sumAbsElements() const;

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
            TMV_Warning("Calling SymMatrix::norm2 without previously "
                        "calling divideUsing(SV)");
            return doNorm2();
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
            TMV_Warning("Calling SymMatrix::condition without previously "
                        "calling divideUsing(SV)");
            return doCondition();
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

        template <class T1, UpLoType U, StorageType S, IndexStyle I> 
        inline void makeInverse(SymMatrix<T1,U,S,I>& sinv) const
        {
            TMVAssert(issym());
            makeInverse(sinv.view()); 
        }

        template <class T1, UpLoType U, StorageType S, IndexStyle I> 
        inline void makeInverse(HermMatrix<T1,U,S,I>& sinv) const
        {
            TMVAssert(isherm());
            makeInverse(sinv.view()); 
        }

        QuotXS<T,T> QInverse() const;
        inline QuotXS<T,T> inverse() const
        { return QInverse(); }

        auto_ptr<BaseMatrix<T> > newCopy() const;
        auto_ptr<BaseMatrix<T> > newView() const;
        auto_ptr<BaseMatrix<T> > newTranspose() const;
        auto_ptr<BaseMatrix<T> > newConjugate() const;
        auto_ptr<BaseMatrix<T> > newAdjoint() const;
        auto_ptr<BaseMatrix<T> > newInverse() const;

        typedef QuotXS<T,T> MyQuotXS;
        TMV_DEPRECATED(MyQuotXS Inverse() const)
        { return inverse(); }
        template <class T1>
        TMV_DEPRECATED(void Inverse(const MatrixView<T1>& minv) const);
        template <class T1, StorageType S, IndexStyle I>
        TMV_DEPRECATED(void Inverse(Matrix<T1,S,I>& minv) const);
        template <class T1>
        TMV_DEPRECATED(void Inverse(const SymMatrixView<T1>& minv) const);
        template <class T1, UpLoType U, StorageType S, IndexStyle I> 
        TMV_DEPRECATED(void Inverse(SymMatrix<T1,U,S,I>& minv) const);
        template <class T1, UpLoType U, StorageType S, IndexStyle I> 
        TMV_DEPRECATED(void Inverse(HermMatrix<T1,U,S,I>& minv) const);
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
            TMVAssert(dt == CH || dt == LU || dt == SV);
            TMVAssert(isherm() || dt != CH);
            DivHelper<T>::divideUsing(dt); 
        }

        inline const SymLDLDiv<T>& lud() const
        {
            divideUsing(LU);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(dynamic_cast<const SymLDLDiv<T>*>(this->getDiv()));
            return *dynamic_cast<const SymLDLDiv<T>*>(this->getDiv());
        }

        inline const HermCHDiv<T>& chd() const
        {
            TMVAssert(isherm());
            divideUsing(CH);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(dynamic_cast<const HermCHDiv<T>*>(this->getDiv()));
            return *dynamic_cast<const HermCHDiv<T>*>(this->getDiv());
        }

        inline const HermSVDiv<T>& svd() const
        {
            TMVAssert(isherm());
            divideUsing(SV);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(dynamic_cast<const HermSVDiv<T>*>(this->getDiv()));
            return *dynamic_cast<const HermSVDiv<T>*>(this->getDiv());
        }

        inline const SymSVDiv<T>& symsvd() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(issym());
            divideUsing(SV);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(dynamic_cast<const SymSVDiv<T>*>(this->getDiv()));
            return *dynamic_cast<const SymSVDiv<T>*>(this->getDiv());
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
        TMV_DEPRECATED(const SymLDLDiv<T>& LUD() const)
        { return lud(); }
        TMV_DEPRECATED(const HermCHDiv<T>& CHD() const)
        { return chd(); }
        TMV_DEPRECATED(const HermSVDiv<T>& SVD() const)
        { return svd(); }
        TMV_DEPRECATED(const SymSVDiv<T>& SymSVD() const)
        { return symsvd(); }

        //
        // I/O
        //

        void write(std::ostream& os) const;
        void writeCompact(std::ostream& os) const;
        void write(std::ostream& os, RT thresh) const;
        void writeCompact(std::ostream& os, RT thresh) const;

        TMV_DEPRECATED(void WriteCompact(std::ostream& fout) const)
        { writeCompact(fout); }
        TMV_DEPRECATED(void WriteCompact(std::ostream& fout, RT thresh) const)
        { writeCompact(fout,thresh); }

        virtual const T* cptr() const = 0;
        virtual int stepi() const = 0;
        virtual int stepj() const = 0;
        virtual UpLoType uplo() const = 0;
        virtual StorageType stor() const = 0;
        virtual ConjType ct() const = 0;
        virtual inline bool isrm() const { return stor() == RowMajor; }
        virtual inline bool iscm() const { return stor() == ColMajor; }
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

        void newDivider() const;
        inline const BaseMatrix<T>& getMatrix() const { return *this; }

    private :

        type& operator=(const type&);

    }; // GenSymMatrix

    template <class T> template <class T2>
    inline bool GenSymMatrix<T>::SameAs(const BaseMatrix<T2>& m2) const
    { return isSameAs(m2); }

    template <class T> template <class T1>
    inline void GenSymMatrix<T>::Inverse(const MatrixView<T1>& minv) const
    { makeInverse(minv); }

    template <class T> template <class T1, StorageType S, IndexStyle I>
    inline void GenSymMatrix<T>::Inverse(Matrix<T1,S,I>& minv) const
    { makeInverse(minv); }

    template <class T> template <StorageType S, IndexStyle I>
    inline void GenSymMatrix<T>::InverseATA(Matrix<T,S,I>& ata) const
    { makeInverseATA(ata); }

    template <class T> template <class T1>
    inline void GenSymMatrix<T>::Inverse(const SymMatrixView<T1>& minv) const
    { makeInverse(minv); }

    template <class T> 
    template <class T1, UpLoType U, StorageType S, IndexStyle I> 
    inline void GenSymMatrix<T>::Inverse(SymMatrix<T1,U,S,I>& minv) const
    { makeInverse(minv); }

    template <class T> 
    template <class T1, UpLoType U, StorageType S, IndexStyle I> 
    inline void GenSymMatrix<T>::Inverse(HermMatrix<T1,U,S,I>& minv) const
    { makeInverse(minv); }


    template <class T, IndexStyle I> 
    class ConstSymMatrixView : public GenSymMatrix<T>
    {
    public :

        typedef ConstSymMatrixView<T,I> type;
        typedef GenSymMatrix<T> base;

        inline ConstSymMatrixView(const type& rhs) :
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itssym(rhs.itssym), itsuplo(rhs.itsuplo), itsstor(rhs.itsstor),
            itsct(rhs.itsct) 
        {
        }

        inline ConstSymMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itss(rhs.size()), 
            itssi(rhs.stepi()), itssj(rhs.stepj()),
            itssym(rhs.sym()), itsuplo(rhs.uplo()), itsstor(rhs.stor()),
            itsct(rhs.ct()) 
        {
        }

        inline ConstSymMatrixView(
            const T* _m, size_t _s, int _si, int _sj,
            SymType _sym, UpLoType _uplo, StorageType _stor,
            ConjType _ct) : 
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
        { 
            TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
                      _si==1 : true);
        }

        virtual inline ~ConstSymMatrixView()
        {
#ifdef TMVDEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        inline size_t size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline SymType sym() const { return itssym; }
        inline UpLoType uplo() const { return itsuplo; }
        inline StorageType stor() const { return itsstor; }
        inline ConjType ct() const { return itsct; }

    protected :

        const T*const itsm;
        const size_t itss;
        const int itssi;
        const int itssj;

        const SymType itssym;
        const UpLoType itsuplo;
        const StorageType itsstor;
        const ConjType itsct;

    private :

        type& operator=(const type&);

    }; // ConstSymMatrixView

    template <class T> 
    class ConstSymMatrixView<T,FortranStyle> : 
        public ConstSymMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef ConstSymMatrixView<T,FortranStyle> type;
        typedef ConstSymMatrixView<T,CStyle> c_type;
        typedef GenSymMatrix<T> base;
        typedef ConstSymMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_lowertri_type;
        typedef ConstSymMatrixView<RT,FortranStyle> const_realpart_type;
        typedef SymMatrixView<T,FortranStyle> nonconst_type;

        inline ConstSymMatrixView(const type& rhs) : c_type(rhs) {}

        inline ConstSymMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline ConstSymMatrixView(const base& rhs) : c_type(rhs) {}

        inline ConstSymMatrixView(
            const T* _m, size_t _s, int _si, int _sj,
            SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct) : 
            c_type(_m,_s,_si,_sj,_sym,_uplo,_stor,_ct) {}

        virtual inline ~ConstSymMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>0 && i<=int(this->size()));
            TMVAssert(j>0 && j<=int(this->size()));
            return base::cref(i-1,j-1);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=int(this->size()));
            TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size()));
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=int(this->size()));
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
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
        // SubMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        bool hasSubVector(int i, int j, int istep, int jstep, int n) const;

        bool hasSubSymMatrix(int i1, int i2, int istep) const;

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

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return base::subVector(i-1,j-1,istep,jstep,n);
        }

        inline const_view_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return base::subSymMatrix(i1-1,i2);
        }

        inline const_view_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return base::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        { return base::upperTri(dt); }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        { return base::lowerTri(dt); }

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
        TMV_DEPRECATED(const_view_type SubSymMatrix(int i1, int i2) const)
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubSymMatrix(
                int i1, int i2, int istep) const)
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_uppertri_type UpperTri(
                DiagType dt = NonUnitDiag) const)
        { return upperTri(); }
        TMV_DEPRECATED(const_lowertri_type LowerTri(
                DiagType dt = NonUnitDiag) const)
        { return lowerTri(); }
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


    private :

        type& operator=(const type&);

    }; // FortranStyle ConstSymMatrixView

    template <class T, IndexStyle I> 
    class SymMatrixView : public GenSymMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymMatrixView<T,I> type;
        typedef GenSymMatrix<T> base;
        typedef SymMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef SymMatrixView<RT,I> realpart_type;
        typedef TMV_RefType(T) reference;

        //
        // Constructors
        //

        inline SymMatrixView(const type& rhs) : 
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itssym(rhs.sym()), itsuplo(rhs.uplo()), itsstor(rhs.stor()),
            itsct(rhs.ct()) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        inline SymMatrixView(
            T* _m, size_t _s, int _si, int _sj, 
            SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
            TMV_PARAMFIRSTLAST(T) 
        ) :
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(_stor==RowMajor ? _sj==1 : _stor==ColMajor ?
                      _si==1 : true); 
        }

        virtual inline ~SymMatrixView() 
        {
#ifdef TMVDEBUG
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
            if (!isSameAs(m2)) upperTri() = m2.upperTri(); 
            return *this; 
        }

        inline const type& operator=(const type& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T()) || m2.sym() == sym());
            if (!isSameAs(m2)) upperTri() = m2.upperTri(); 
            return *this; 
        }

        inline const type& operator=(const GenSymMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToS(*this);
            return *this; 
        }

        inline const type& operator=(const GenSymMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.sym() == sym());
            m2.assignToS(*this);
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenSymMatrix<T2>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(isReal(T2()) || m2.sym() == sym());
            if (!isSameAs(m2)) upperTri() = m2.upperTri(); 
            return *this; 
        }

        inline const type& operator=(const T& x) const 
        { 
            TMVAssert(issym() || TMV_IMAG(x) == RT(0));
            return setToIdentity(x); 
        }

        inline const type& operator=(const AssignableToSymMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToS(view());
            return *this;
        }

        inline const type& operator=(const AssignableToSymMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(sym() == m2.sym());
            m2.assignToS(view());
            return *this;
        }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToD(DiagMatrixViewOf(diag()));
            upperTri().offDiag().setZero();
            return *this;
        }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(issym());
            m2.assignToD(DiagMatrixViewOf(diag()));
            upperTri().offDiag().setZero();
            TMVAssert(this->isHermOK());
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const OProdVV<RT,T2,T2>& opvv) const
        {
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(
                    issym() ? opvv.getV2().view() : opvv.getV2().conjugate()));
            TMVAssert(issym() || TMV_IMAG(opvv.getX()) == RT(0));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const OProdVV<CT,T2,T2>& opvv) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(
                    issym() ? opvv.getV2().view() : opvv.getV2().conjugate()));
            TMVAssert(issym() || TMV_IMAG(opvv.getX()) == RT(0));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const ProdMM<RT,T2,T2>& pmm) const
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(issym() ? pmm.getM2().transpose() :
                                           pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const ProdMM<CT,T2,T2>& pmm) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const ProdUL<RT,T2,T2>& pmm) const
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const ProdUL<CT,T2,T2>& pmm) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const ProdLU<RT,T2,T2>& pmm) const
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <class T2> 
        inline const type& operator=(const ProdLU<CT,T2,T2>& pmm) const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }


        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        {
            TMVAssert(i>=0 && i<int(size()));
            TMVAssert(j>=0 && j<int(size()));
            return ref(i,j); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
            if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
                return vec_type(
                    ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
                    ct() TMV_FIRSTLAST); 
            else
                return vec_type(
                    ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<int(size()));
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
            if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
                return vec_type(
                    ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
                    ct() TMV_FIRSTLAST); 
            else 
                return vec_type(
                    ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type diag() const
        { return vec_type(ptr(),size(),stepi()+stepj(),ct() TMV_FIRSTLAST); }

        inline vec_type diag(int i) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            if (i>=0)
                if (uplo()==Upper)
                    return vec_type(
                        ptr()+i*stepj(),size()-i,
                        stepi()+stepj(),ct() TMV_FIRSTLAST); 
                else
                    return vec_type(
                        ptr()+i*stepi(),size()-i,
                        stepi()+stepj(),this->issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST);
            else
                if (uplo()==Upper)
                    return vec_type(
                        ptr()-i*stepj(),size()+i,
                        stepi()+stepj(),this->issym()?ct():TMV_ConjOf(T,ct()) 
                        TMV_FIRSTLAST);
                else
                    return vec_type(
                        ptr()-i*stepi(),size()+i,
                        stepi()+stepj(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
            const int ds = stepi()+stepj();
            if (i>=0)
                if (uplo()==Upper)
                    return vec_type(
                        ptr()+i*stepj()+j1*ds,j2-j1,ds,ct() 
                        TMV_FIRSTLAST); 
                else
                    return vec_type(
                        ptr()+i*stepi()+j1*ds,j2-j1,ds, 
                        this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
            else
                if (uplo()==Upper)
                    return vec_type(
                        ptr()-i*stepj()+j1*ds,j2-j1,ds, 
                        this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
                else
                    return vec_type(
                        ptr()-i*stepi()+j1*ds,j2-j1,ds,ct() 
                        TMV_FIRSTLAST); 
        }

        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { upperTri().setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        { 
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            upperTri().setAllTo(x); return *this; 
        }

        inline const type& clip(RT thresh) const
        { upperTri().clip(thresh); return *this; }

        inline const type& transposeSelf() const
        { if (!this->issym()) upperTri().conjugateSelf(); return *this; }

        inline const type& conjugateSelf() const
        { if (isComplex(T())) upperTri().conjugateSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { 
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            setZero(); diag().setAllTo(x); return *this; 
        }

        const type& swapRowsCols(int i1, int i2) const;

        const type& permuteRowsCols(const int* p, int i1, int i2) const;

        const type& reversePermuteRowsCols(const int* p, int i1, int i2) const;

        inline const type& permuteRowsCols(const int* p) const
        { return permuteRowsCols(p,0,size()); }

        inline const type& reversePermuteRowsCols(const int* p) const
        { return reversePermuteRowsCols(p,0,size()); }

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
        TMV_DEPRECATED(const type& SwapRowsCols(int i1, int i2) const)
        { return swapRowsCols(i1,i2); }
        TMV_DEPRECATED(const type& PermuteRowsCols(
                const int* p, int i1, int i2) const)
        { return permuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(const type& PermuteRowsCols(const int* p) const)
        { return permuteRowsCols(p); }
        TMV_DEPRECATED(const type& ReversePermuteRowsCols(
                const int* p, int i1, int i2) const)
        { return reversePermuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(const type& ReversePermuteRowsCols(const int* p) const)
        { return reversePermuteRowsCols(p); }


        //
        // subMatrix
        //

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            if ((i2-j1<=1 && uplo()==Upper) || (j2-i1<=1 && uplo()==Lower))
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
                this->iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                this->isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if ((i2-istep<=j1 && uplo()==Upper) || 
                (j2-jstep<=i1 && uplo()==Lower))
                return rec_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,ct() TMV_FIRSTLAST);
            else
                return rec_type(
                    ptr()+i1*stepj()+j1*stepi(),(i2-i1)/istep,(j2-j1)/jstep,
                    istep*stepj(),jstep*stepi(), TMV_TransOf(newstor),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type subVector(int i, int j,
                                  int istep, int jstep, int n) const
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,n));
            if ((i-j<=0 && uplo()==Upper) || (j-i<=0 && uplo()==Lower))
                return vec_type(
                    ptr()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),ct() TMV_FIRSTLAST);
            else
                return vec_type(
                    ptr()+i*stepj()+j*stepi(),n,istep*stepj()+jstep*stepi(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(base::hasSubSymMatrix(i1,i2,1));
            return view_type(
                ptr()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),sym(),uplo(),stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(base::hasSubSymMatrix(i1,i2,istep));
            return view_type(
                ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),sym(),uplo(),
                istep==1 ? stor() : NoMajor,ct() TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            if (uplo() == Upper)
                return uppertri_type(
                    ptr(),size(),
                    stepi(),stepj(),dt,stor(),ct() TMV_FIRSTLAST);
            else
                return uppertri_type(
                    ptr(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            if (uplo() == Lower)
                return lowertri_type(
                    ptr(),size(),
                    stepi(),stepj(),dt,stor(),ct() TMV_FIRSTLAST);
            else
                return lowertri_type(
                    ptr(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(stor()),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline realpart_type realPart() const
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
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
                size(),2*stepi(),2*stepj(),sym(),uplo(),NoMajor, NonConj
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
                ptr(),size(),stepj(),stepi(),
                sym(),TMV_UTransOf(uplo()),TMV_TransOf(stor()),ct() 
                TMV_FIRSTLAST);
        }

        inline view_type conjugate() const
        {
            return view_type(
                ptr(),size(),stepi(),stepj(),
                sym(),uplo(),stor(),TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline view_type adjoint() const
        {
            return view_type(
                ptr(),size(),stepj(),stepi(),
                sym(),TMV_UTransOf(uplo()),TMV_TransOf(stor()),
                TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
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
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2) const)
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubSymMatrix(
                int i1, int i2, int istep) const)
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(uppertri_type UpperTri(
                DiagType dt = NonUnitDiag) const)
        { return upperTri(); }
        TMV_DEPRECATED(lowertri_type LowerTri(
                DiagType dt = NonUnitDiag) const)
        { return lowerTri(); }
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


        //
        // I/O
        //

        void read(std::istream& is) const;

        inline size_t size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline SymType sym() const { return itssym; }
        inline UpLoType uplo() const { return itsuplo; }
        inline StorageType stor() const { return itsstor; }
        inline ConjType ct() const { return itsct; }
        using base::issym;
        using base::isherm;

        reference ref(int i, int j) const;

    protected :

        T*const itsm;
        const size_t itss;
        const int itssi;
        const int itssj;

        const SymType itssym;
        const UpLoType itsuplo;
        const StorageType itsstor;
        const ConjType itsct;

#ifdef TMVFLDEBUG
    public:
        const T*const _first;
        const T*const _last;
    protected:
#endif

    }; // SymMatrixView

    template <class T> 
    class SymMatrixView<T,FortranStyle> : public SymMatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymMatrixView<T,FortranStyle> type;
        typedef SymMatrixView<T,CStyle> c_type;
        typedef GenSymMatrix<T> base;
        typedef SymMatrixView<T,FortranStyle> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,FortranStyle> vec_type;
        typedef MatrixView<T,FortranStyle> rec_type;
        typedef UpperTriMatrixView<T,FortranStyle> uppertri_type;
        typedef LowerTriMatrixView<T,FortranStyle> lowertri_type;
        typedef SymMatrixView<RT,FortranStyle> realpart_type;
        typedef ConstSymMatrixView<T,FortranStyle> const_type;

        //
        // Constructors
        //

        inline SymMatrixView(const type& rhs) : c_type(rhs) {}

        inline SymMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline SymMatrixView(
            T* _m, size_t _s, int _si, int _sj, 
            SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_s,_si,_sj,_sym,_uplo,_stor,_ct
                   TMV_FIRSTLAST1(_first,_last) ) {}

        virtual inline ~SymMatrixView() {} 

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

        inline const type& operator=(const GenSymMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenSymMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenSymMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(const AssignableToSymMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const AssignableToSymMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2>
        inline const type& operator=(const OProdVV<T,T2,T2>& opvv) const
        { c_type::operator=(opvv); return *this; }

        template <class T2>
        inline const type& operator=(const ProdMM<T,T2,T2>& pmm) const
        { c_type::operator=(pmm); return *this; }

        template <class T2>
        inline const type& operator=(const ProdLU<T,T2,T2>& pmm) const
        { c_type::operator=(pmm); return *this; }

        template <class T2>
        inline const type& operator=(const ProdUL<T,T2,T2>& pmm) const
        { c_type::operator=(pmm); return *this; }


        //
        // Access
        //

        inline TMV_RefType(T) operator()(int i,int j) const 
        {
            TMVAssert(i>0 && i<=int(this->size()));
            TMVAssert(j>0 && j<=int(this->size()));
            return c_type::ref(i-1,j-1); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=int(this->size()));
            TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size()));
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=int(this->size()));
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
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

        inline const type& clip(RT thresh) const
        { c_type::clip(thresh); return *this; }

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& transposeSelf() const
        { c_type::transposeSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { c_type::setToIdentity(x); return *this; }

        inline const type& swapRowsCols(int i1, int i2) const
        {
            TMVAssert(i1>0 && i1<=int(this->size()));
            TMVAssert(i2>0 && i2<=int(this->size()));
            c_type::swapRowsCols(i1-1,i2-1); 
            return *this; 
        }

        inline const type& permuteRowsCols(
            const int* p, int i1, int i2) const
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
            c_type::permuteRowsCols(p,i1-1,i2); 
            return *this; 
        }

        inline const type& reversePermuteRowsCols(
            const int* p, int i1, int i2) const
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
            c_type::reversePermuteRowsCols(p,i1-1,i2); 
            return *this; 
        }

        inline const type& permuteRowsCols(const int* p) const
        { c_type::permuteRowsCols(p); return *this; }

        inline const type& reversePermuteRowsCols(const int* p) const
        { c_type::reversePermuteRowsCols(p); return *this; }

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
        TMV_DEPRECATED(const type& SwapRowsCols(int i1, int i2) const)
        { return swapRowsCols(i1,i2); }
        TMV_DEPRECATED(const type& PermuteRowsCols(
                const int* p, int i1, int i2) const)
        { return permuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(const type& PermuteRowsCols(const int* p) const)
        { return permuteRowsCols(p); }
        TMV_DEPRECATED(const type& ReversePermuteRowsCols(
                const int* p, int i1, int i2) const)
        { return reversePermuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(const type& ReversePermuteRowsCols(const int* p) const)
        { return reversePermuteRowsCols(p); }

        //
        // subMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return const_type(*this).hasSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int n) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,n); }


        inline bool hasSubSymMatrix(int i1, int i2, int istep) const
        { return const_type(*this).hasSubSymMatrix(i1,i2,istep); }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
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
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return c_type::subVector(i-1,j-1,istep,jstep,n);
        }

        inline view_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return c_type::subSymMatrix(i1-1,i2);
        }

        inline view_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return c_type::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        { return c_type::upperTri(dt); }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        { return c_type::lowerTri(dt); }

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

        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2) const)
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep) const)
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s) const)
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2) const)
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2, int istep) const)
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(uppertri_type UpperTri(DiagType dt = NonUnitDiag) const)
        { return upperTri(); }
        TMV_DEPRECATED(lowertri_type LowerTri(DiagType dt = NonUnitDiag) const)
        { return lowerTri(); }
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

    }; // FortranStyle SymMatrixView

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    class SymMatrix : public GenSymMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymMatrix<T,U,S,I> type;
        typedef type copy_type;
        typedef GenSymMatrix<T> base;
        typedef ConstSymMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstUpperTriMatrixView<T,I> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,I> const_lowertri_type;
        typedef ConstSymMatrixView<RT,I> const_realpart_type;
        typedef SymMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef SymMatrixView<RT,I> realpart_type;
        typedef T& reference;

        //
        // Constructors
        //

#define NEW_SIZE(s) itslen((s)*(s)), itsm(itslen), itss(s)  \
        TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen) 

        explicit inline SymMatrix(size_t _size) : NEW_SIZE(_size) 
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        inline SymMatrix(size_t _size, const T& x) : NEW_SIZE(_size)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            setAllTo(x);
        }

        inline SymMatrix(size_t _size, const T* vv) : NEW_SIZE(_size)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(vv,vv+itslen,itsm.get());
        }

        inline SymMatrix(size_t _size, const std::vector<T>& vv) :
            NEW_SIZE(_size)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == itslen);
            std::copy(vv.begin(),vv.end(),itsm.get());
        }

        inline SymMatrix(const type& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        template <IndexStyle I2> 
        inline SymMatrix(const SymMatrix<T,U,S,I2>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        inline SymMatrix(const GenSymMatrix<RT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            if (rhs.issym()) 
                rhs.assignToS(view());
            else {
                if (U == Upper)
                    upperTri() = rhs.upperTri();
                else
                    lowerTri() = rhs.lowerTri();
            }
        }

        inline SymMatrix(const GenSymMatrix<CT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor);
            if (rhs.issym()) 
                rhs.assignToS(view());
            else {
                if (U == Upper)
                    upperTri() = rhs.upperTri();
                else
                    lowerTri() = rhs.lowerTri();
            }
        }

        template <class T2> 
        inline SymMatrix(const GenSymMatrix<T2>& rhs) : NEW_SIZE(rhs.size())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            if (U == Upper)
                upperTri() = rhs.upperTri();
            else
                lowerTri() = rhs.lowerTri();
        }

        template <IndexStyle I2> 
        inline explicit SymMatrix(const Matrix<T,S,I2>& rhs) :
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(rhs.isSquare());
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        template <class T2> 
        inline SymMatrix(const GenMatrix<T2>& rhs) :
            NEW_SIZE(rhs.rowsize())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(rhs.isSquare());
            if (U == Upper)
                upperTri() = rhs.upperTri();
            else
                lowerTri() = rhs.lowerTri();
        }

        inline SymMatrix(const AssignableToSymMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.issym());
            m2.assignToS(view());
        }

        inline SymMatrix(const AssignableToSymMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.issym());
            m2.assignToS(view());
        }

        inline explicit SymMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(size() == m2.size());
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        inline explicit SymMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        template <class T2> 
        inline SymMatrix(const OProdVV<T,T2,T2>& opvv) :
            NEW_SIZE(opvv.colsize())
        {
            TMVAssert(opvv.colsize() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().view()));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
        }

        template <class T2> 
        inline SymMatrix(const ProdMM<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <class T2> 
        inline SymMatrix(const ProdUL<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <class T2> 
        inline SymMatrix(const ProdLU<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

#undef NEW_SIZE

        virtual inline ~SymMatrix() 
        {
#ifdef TMVDEBUG
            setAllTo(T(999));
#endif
        }

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { 
            TMVAssert(size() == m2.size());
            if (&m2 != this) 
                std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
            return *this;
        }

        template <IndexStyle I2>
        inline type& operator=(const SymMatrix<T,U,S,I2>& m2)
        { 
            TMVAssert(size() == m2.size());
            if (&m2 != this) 
                std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
            return *this;
        }

        inline type& operator=(const GenSymMatrix<RT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(m2.issym());
            m2.assignToS(view());
            return *this; 
        }

        inline type& operator=(const GenSymMatrix<CT>& m2)
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.issym());
            m2.assignToS(view());
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenSymMatrix<T2>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T2()) || m2.issym());
            upperTri() = m2.upperTri(); 
            return *this; 
        }

        inline type& operator=(const T& x) 
        { return setToIdentity(x); }

        inline type& operator=(const AssignableToSymMatrix<RT>& m2) 
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(m2.issym());
            m2.assignToS(view());
            return *this;
        }

        inline type& operator=(const AssignableToSymMatrix<CT>& m2) 
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.issym());
            m2.assignToS(view());
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
            view() = m2;
            return *this;
        }

        template <class T2> 
        inline type& operator=(const OProdVV<T,T2,T2>& opvv)
        {
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().view()));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const ProdMM<T,T2,T2>& pmm) 
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const ProdUL<T,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const ProdLU<T,T2,T2>& pmm) 
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }


        //
        // Access
        //

        inline T operator()(int i, int j) const
        { 
            if (I==CStyle) {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j>=0 && j<int(size()));
                return cref(i,j); 
            } else {
                TMVAssert(i>0 && i<=int(size()));
                TMVAssert(j>0 && j<=int(size()));
                return cref(i-1,j-1); 
            }
        }

        inline T& operator()(int i, int j) 
        { 
            if (I==CStyle) {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j>=0 && j<int(size()));
                return ref(i,j); 
            } else {
                TMVAssert(i>0 && i<=int(size()));
                TMVAssert(j>0 && j<=int(size()));
                return ref(i-1,j-1); 
            }
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=int(size()));
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size()));
                --j1;
            } else {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
            }
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return const_vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj); 
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),NonConj); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=int(size()));
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size()));
                --i1;
            } else {
                TMVAssert(j>=0 && j<int(size()));
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
            }
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return const_vec_type(
                    itsm.get()+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj); 
            else 
                return const_vec_type(
                    itsm.get()+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),NonConj); 
        }

        inline const_vec_type diag() const
        { return const_vec_type(
                itsm.get(),size(),stepi()+stepj(),NonConj); }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            i = std::abs(i);
            return const_vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi()),
                size()-i,stepi()+stepj(),NonConj); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
            }
            i = std::abs(i);
            const int ds = stepi()+stepj();
            return const_vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,NonConj);
        }

        inline vec_type row(int i, int j1, int j2)
        { 
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=int(size()));
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size()));
                --j1;
            } else {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
            }
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj TMV_FIRSTLAST); 
            else
                return vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=int(size()));
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size()));
                --i1;
            } else {
                TMVAssert(j>=0 && j<int(size()));
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
            }
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return vec_type(
                    itsm.get()+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj TMV_FIRSTLAST); 
            else 
                return vec_type(
                    itsm.get()+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag()
        {
            return vec_type(
                itsm.get(),size(),stepi()+stepj(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i) 
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            i = std::abs(i);
            TMVAssert(i<=int(size())); 
            return vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi()),
                size()-i,stepi()+stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i, int j1, int j2) 
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
            }
            i = std::abs(i);
            const int ds = stepi()+stepj();
            return vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,NonConj TMV_FIRSTLAST);
        }


        //
        // Modifying Functions
        //

        inline type& setZero() 
        { std::fill_n(itsm.get(),itslen,T(0)); return *this; }

        inline type& setAllTo(const T& x) 
        { upperTri().setAllTo(x); return *this; }

        inline type& clip(RT thresh) 
        { upperTri().clip(thresh); return *this; }

        inline type& conjugateSelf() 
        { if (isComplex(T())) upperTri().conjugateSelf(); return *this; }

        inline type& transposeSelf() 
        { return *this; }

        inline type& setToIdentity(const T& x=T(1)) 
        { setZero(); diag().setAllTo(x); return *this; }

        inline type& swapRowsCols(int i1, int i2)
        { view().swapRowsCols(i1,i2); return *this; }

        inline type& permuteRowsCols(const int* p, int i1, int i2)
        { view().permuteRowsCols(p,i1,i2); return *this; }

        inline type& reversePermuteRowsCols(const int* p, int i1, int i2)
        { view().reversePermuteRowsCols(p,i1,i2); return *this; }

        inline type& permuteRowsCols(const int* p)
        { view().permuteRowsCols(p); return *this; }

        inline type& reversePermuteRowsCols(const int* p)
        { view().reversePermuteRowsCols(p); return *this; }

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
        TMV_DEPRECATED(type& SwapRowsCols(int i1, int i2))
        { return swapRowsCols(i1,i2); }
        TMV_DEPRECATED(type& PermuteRowsCols(
                const int* p, int i1, int i2))
        { return permuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(type& PermuteRowsCols(const int* p))
        { return permuteRowsCols(p); }
        TMV_DEPRECATED(type& ReversePermuteRowsCols(
                const int* p, int i1, int i2))
        { return reversePermuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(type& ReversePermuteRowsCols(const int* p))
        { return reversePermuteRowsCols(p); }


        //
        // subMatrix
        //

        inline const_rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
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
            if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),NonConj);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==FortranStyle) { --i; --j; }
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return const_vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj);
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),NonConj);
        }

        inline const_view_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),Sym,U,S,NonConj);
        }

        inline const_view_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Sym,U,
                istep==1 ? S : NoMajor, NonConj);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            return U==Upper ? 
                const_uppertri_type(
                    itsm.get(),size(),
                    stepi(),stepj(),dt,S,NonConj) :
                const_uppertri_type(
                    itsm.get(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(S),NonConj);
        }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            return U==Lower ? 
                const_lowertri_type(
                    itsm.get(),size(),
                    stepi(),stepj(),dt,S,NonConj) :
                const_lowertri_type(
                    itsm.get(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(S),NonConj);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                Sym,U,isReal(T())?S:NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get())+1,size(),
                2*stepi(),2*stepj(),Sym,U,NoMajor,NonConj);
        } 

        inline const_view_type view() const
        { 
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),
                Sym,U,S,NonConj);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),NonConj);
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),
                Sym,U,S,TMV_ConjOf(T,NonConj));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),TMV_TransOf(S),NonConj 
                    TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),NonConj TMV_FIRSTLAST);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int n)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==FortranStyle) { --i; --j; }
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(int i1, int i2)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return view_type(
                itsm.get()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),Sym,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(int i1, int i2, int istep) 
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return view_type(
                itsm.get()+i1*(stepi()+stepj()),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),Sym,U,
                istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag)
        {
            return U==Upper ? 
                uppertri_type(
                    itsm.get(),size(),
                    stepi(),stepj(),dt,S,NonConj TMV_FIRSTLAST) :
                uppertri_type(
                    itsm.get(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag)
        {
            return U==Lower ? 
                lowertri_type(
                    itsm.get(),size(),
                    stepi(),stepj(),dt,S,NonConj TMV_FIRSTLAST) :
                lowertri_type(
                    itsm.get(),size(),
                    stepj(),stepi(),dt,TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
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
                reinterpret_cast<RT*>(itsm.get())+1,size(),
                2*stepi(),2*stepj(),Sym,U,NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        } 

        inline view_type view() 
        { 
            return view_type(
                itsm.get(),size(),stepi(),stepj(),
                Sym,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type transpose() 
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline view_type conjugate() 
        { 
            return view_type(
                itsm.get(),size(),stepi(),stepj(),
                Sym,U,S,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline view_type adjoint() 
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj) 
                TMV_FIRSTLAST);
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
        TMV_DEPRECATED(const_view_type SubSymMatrix(int i1, int i2) const)
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubSymMatrix(
                int i1, int i2, int istep) const)
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_uppertri_type UpperTri(
                DiagType dt = NonUnitDiag) const)
        { return upperTri(); }
        TMV_DEPRECATED(const_lowertri_type LowerTri(
                DiagType dt = NonUnitDiag) const)
        { return lowerTri(); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }

        TMV_DEPRECATED(rec_type SubMatrix(int i1, int i2, int j1, int j2))
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep))
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s))
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2))
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2, int istep))
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(uppertri_type UpperTri(DiagType dt = NonUnitDiag))
        { return upperTri(); }
        TMV_DEPRECATED(lowertri_type LowerTri(DiagType dt = NonUnitDiag))
        { return lowerTri(); }
        TMV_DEPRECATED(realpart_type Real())
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag())
        { return imagPart(); }
        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Transpose())
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint())
        { return adjoint(); }


        inline size_t size() const { return itss; }
        inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        inline int stepi() const { return S==RowMajor ? itss : 1; }
        inline int stepj() const { return S==RowMajor ? 1 : itss; }
        inline SymType sym() const { return Sym; }
        inline UpLoType uplo() const { return U; }
        inline StorageType stor() const { return S; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isconj() const { return false; }
        inline bool isherm() const { return isReal(T()); }
        inline bool issym() const { return true; }
        inline bool isupper() const { return U == Upper; }

        inline T& ref(int i, int j)
        {
            if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                if (S == RowMajor) return itsm.get()[i*itss + j];
                else return itsm.get()[j*itss + i];
            else 
                if (S == RowMajor) return itsm.get()[j*itss + i];
                else return itsm.get()[i*itss + j];
        }

        inline T cref(int i, int j) const 
        {
            if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                if (S == RowMajor) return itsm.get()[i*itss + j];
                else return itsm.get()[j*itss + i];
            else 
                if (S == RowMajor) return itsm.get()[j*itss + i];
                else return itsm.get()[i*itss + j];
        }

    protected :

        const size_t itslen;
        AlignedArray<T> itsm;
        const size_t itss;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

        template <IndexStyle I2>
        friend void Swap(SymMatrix<T,U,S,I>& m1, SymMatrix<T,U,S,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            m1.itsm.swapWith(m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // SymMatrix

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    class HermMatrix : public GenSymMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef HermMatrix<T,U,S,I> type;
        typedef type copy_type;
        typedef GenSymMatrix<T> base;
        typedef ConstSymMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstUpperTriMatrixView<T,I> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,I> const_lowertri_type;
        typedef ConstSymMatrixView<RT,I> const_realpart_type;
        typedef SymMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef SymMatrixView<RT,I> realpart_type;
        typedef TMV_RefType(T) reference;

        //
        // Constructors
        //

#define NEW_SIZE(s) itslen((s)*(s)), itsm(itslen), itss(s)  \
        TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

        explicit inline HermMatrix(size_t _size) : NEW_SIZE(_size) 
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        HermMatrix(size_t _size, const RT& x) : NEW_SIZE(_size)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            setAllTo(x);
        }

        HermMatrix(size_t _size, const T* vv) : NEW_SIZE(_size)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(vv,vv+itslen,itsm.get());
            TMVAssert(this->isHermOK());
        }

        HermMatrix(size_t _size, const std::vector<T>& vv) : NEW_SIZE(_size)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == itslen);
            std::copy(vv.begin(),vv.end(),itsm.get());
            TMVAssert(this->isHermOK());
        }

        HermMatrix(const type& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        template <IndexStyle I2> 
        HermMatrix(const HermMatrix<T,U,S,I2>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        HermMatrix(const GenSymMatrix<RT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            if (rhs.isherm()) rhs.assignToS(view());
            else if (U == Upper) upperTri() = rhs.upperTri();
            else lowerTri() = rhs.lowerTri();
        }

        HermMatrix(const GenSymMatrix<CT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor);
            if (rhs.isherm())
                rhs.assignToS(view());
            else {
                if (U == Upper) upperTri() = rhs.upperTri();
                else lowerTri() = rhs.lowerTri();
                diag().imagPart().setZero();
            }
        }

        template <class T2> 
        inline HermMatrix(const GenSymMatrix<T2>& rhs) : NEW_SIZE(rhs.size())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            if (U == Upper) upperTri() = rhs.upperTri();
            else lowerTri() = rhs.lowerTri();
            if (isComplex(T()) && isComplex(T2()) && rhs.issym()) 
                diag().imagPart().setZero();
        }

        template <IndexStyle I2> 
        inline explicit HermMatrix(const Matrix<T,S,I2>& rhs) : 
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(rhs.isSquare());
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermMatrix(const GenMatrix<T2>& rhs) : NEW_SIZE(rhs.rowsize())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(rhs.isSquare());
            if (U == Upper)
                upperTri() = rhs.upperTri();
            else
                lowerTri() = rhs.lowerTri();
            if (isComplex(T()) && isComplex(T2())) diag().imagPart().setZero();
        }

        inline HermMatrix(const AssignableToSymMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.isherm());
            m2.assignToS(view());
        }

        inline HermMatrix(const AssignableToSymMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(m2.isherm());
            m2.assignToS(view());
        }

        inline explicit HermMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            m2.assignToD(DiagMatrixViewOf(diag()));
            upperTri().offDiag().setZero();
        }

        inline explicit HermMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        { 
            TMVAssert(S==RowMajor || S==ColMajor);
            m2.assignToD(DiagMatrixViewOf(diag().realPart()));
            upperTri().offDiag().setZero();
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <class T2> 
        inline HermMatrix(const OProdVV<T,T2,T2>& opvv) :
            NEW_SIZE(opvv.colsize())
        {
            TMVAssert(opvv.colsize() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().conjugate()));
            TMVAssert(TMV_IMAG(opvv.getX()) == RT(0));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
        }

        template <class T2> 
        inline HermMatrix(const ProdMM<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            TMVAssert(TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <class T2> 
        inline HermMatrix(const ProdUL<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            TMVAssert(TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <class T2> 
        inline HermMatrix(const ProdLU<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            TMVAssert(TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }


#undef NEW_SIZE

        virtual inline ~HermMatrix() 
        {
#ifdef TMVDEBUG
            setAllTo(T(999));
#endif
        }

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { 
            TMVAssert(m2.size() == size());
            if (&m2 != this) 
                std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
            return *this;
        }

        template <IndexStyle I2>
        inline type& operator=(const HermMatrix<T,U,S,I2>& m2)
        { 
            TMVAssert(m2.size() == size());
            if (&m2 != this) 
                std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
            return *this;
        }

        inline type& operator=(const GenSymMatrix<RT>& m2)
        {
            TMVAssert(m2.size() == size());
            TMVAssert(m2.isherm());
            m2.assignToS(view());
            return *this; 
        }

        inline type& operator=(const GenSymMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(m2.size() == size());
            TMVAssert(m2.isherm());
            m2.assignToS(view());
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenSymMatrix<T2>& m2)
        {
            TMVAssert(m2.size() == size());
            TMVAssert(m2.isherm());
            upperTri() = m2.upperTri(); 
            return *this; 
        }

        inline type& operator=(const T& x) 
        {
            TMVAssert(TMV_IMAG(x) == RT(0));
            return setToIdentity(x); 
        }

        inline type& operator=(
            const AssignableToSymMatrix<RT>& m2) 
        { 
            TMVAssert(size() == m2.colsize());
            TMVAssert(m2.isherm());
            m2.assignToS(view());
            return *this;
        }

        inline type& operator=(
            const AssignableToSymMatrix<CT>& m2) 
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.colsize());
            TMVAssert(m2.isherm());
            m2.assignToS(view());
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2) 
        { 
            TMVAssert(size() == m2.colsize());
            view() = m2;
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2) 
        { 
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.colsize());
            view() = m2;
            TMVAssert(this->isHermOK());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const OProdVV<T,T2,T2>& opvv)
        {
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().conjugate()));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const ProdMM<T,T2,T2>& pmm) 
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const ProdUL<T,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const ProdLU<T,T2,T2>& pmm) 
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }


        //
        // Access
        //

        inline T operator()(int i, int j) const
        {
            if (I==CStyle) {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j>=0 && j<int(size()));
                return cref(i,j); 
            } else {
                TMVAssert(i>0 && i<=int(size()));
                TMVAssert(j>0 && j<=int(size()));
                return cref(i-1,j-1); 
            }
        }

        inline reference operator()(int i, int j) 
        { 
            if (I==CStyle) {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j>=0 && j<int(size()));
                return ref(i,j); 
            } else {
                TMVAssert(i>0 && i<=int(size()));
                TMVAssert(j>0 && j<=int(size()));
                return ref(i-1,j-1); 
            }
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=int(size()));
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size()));
                --j1;
            } else {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
            }
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return const_vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj); 
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),TMV_ConjOf(T,NonConj)); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=int(size()));
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size()));
                --i1;
            } else {
                TMVAssert(j>=0 && j<int(size()));
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
            }
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return const_vec_type(
                    itsm.get()+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj); 
            else 
                return const_vec_type(
                    itsm.get()+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),TMV_ConjOf(T,NonConj)); 
        }

        inline const_vec_type diag() const
        { return const_vec_type(
                itsm.get(),size(),stepi()+stepj(),NonConj); }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            TMVAssert(i<=int(size())); 
            return const_vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi()),
                size()-i,stepi()+stepj(),newct);
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
            }
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            const int ds = stepi()+stepj();
            return const_vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,newct);
        }

        inline vec_type row(int i, int j1, int j2)
        { 
            if (I==FortranStyle) {
                TMVAssert(i>0 && i<=int(size()));
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size()));
                --j1;
            } else {
                TMVAssert(i>=0 && i<int(size()));
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
            }
            if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
                return vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj TMV_FIRSTLAST); 
            else
                return vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=int(size()));
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size()));
                --i1;
            } else {
                TMVAssert(j>=0 && j<int(size()));
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
            }
            if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
                return vec_type(
                    itsm.get()+i1*stepi()+j*stepj(),
                    i2-i1,stepi(),NonConj TMV_FIRSTLAST); 
            else 
                return vec_type(
                    itsm.get()+i1*stepj()+j*stepi(),
                    i2-i1,stepj(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST); 
        }

        inline vec_type diag()
        {
            return vec_type(
                itsm.get(),size(),stepi()+stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type diag(int i) 
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            TMVAssert(i<=int(size())); 
            return vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi()),size()-i,
                stepi()+stepj(),newct TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i, int j1, int j2) 
        {
            TMVAssert(i>=-int(size()) && i<=int(size())); 
            if (I==FortranStyle) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
            }
            ConjType newct = 
                ((i>0) == (U==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            const int ds = stepi()+stepj();
            return vec_type(
                itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,newct TMV_FIRSTLAST);
        }

        //
        // Modifying Functions
        //

        inline type& setZero() 
        {
            std::fill_n(itsm.get(),itslen,T(0));
            return *this;
        }

        inline type& setAllTo(const T& x) 
        { 
            TMVAssert(TMV_IMAG(x) == RT(0));
            upperTri().setAllTo(x); 
            return *this; 
        }

        inline type& clip(RT thresh) 
        { upperTri().clip(thresh); return *this; }

        inline type& conjugateSelf() 
        { if (isComplex(T())) upperTri().conjugateSelf(); return *this; }

        inline type& transposeSelf() 
        { conjugateSelf(); }

        inline type& setToIdentity(const T& x=T(1)) 
        { 
            TMVAssert(TMV_IMAG(x) == RT(0));
            setZero(); 
            diag().setAllTo(x); 
            return *this; 
        }

        inline type& swapRowsCols(int i1, int i2)
        { view().swapRowsCols(i1,i2); return *this; }

        inline type& permuteRowsCols(const int* p, int i1, int i2)
        { view().permuteRowsCols(p,i1,i2); return *this; }

        inline type& reversePermuteRowsCols(const int* p, int i1, int i2)
        { view().reversePermuteRowsCols(p,i1,i2); return *this; }

        inline type& permuteRowsCols(const int* p)
        { view().permuteRowsCols(p); return *this; }

        inline type& reversePermuteRowsCols(const int* p)
        { view().reversePermuteRowsCols(p); return *this; }

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
        TMV_DEPRECATED(type& SwapRowsCols(int i1, int i2))
        { return swapRowsCols(i1,i2); }
        TMV_DEPRECATED(type& PermuteRowsCols(
                const int* p, int i1, int i2))
        { return permuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(type& PermuteRowsCols(const int* p))
        { return permuteRowsCols(p); }
        TMV_DEPRECATED(type& ReversePermuteRowsCols(
                const int* p, int i1, int i2))
        { return reversePermuteRowsCols(p,i1,i2); }
        TMV_DEPRECATED(type& ReversePermuteRowsCols(const int* p))
        { return reversePermuteRowsCols(p); }


        //
        // subMatrix
        //

        inline const_rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),i2-i1,j2-j1,
                    stepj(),stepi(),TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),TMV_ConjOf(T,NonConj));
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int n) const
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,n));
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return const_vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj);
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),TMV_ConjOf(T,NonConj));
        }

        inline const_view_type subSymMatrix(int i1, int i2) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),Herm,U,S,NonConj);
        }

        inline const_view_type subSymMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,U,
                istep==1 ? S : NoMajor, NonConj);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            return U==Upper ? 
                const_uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),
                    dt,S,NonConj) :
                const_uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            return U==Lower ? 
                const_lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),
                    dt,S,NonConj) :
                const_lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                Herm,U,isReal(T())?S:NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        { 
            // The imaginary part of a Hermitian matrix is anti-symmetric
            // so this is illegal.
            TMVAssert(TMV_FALSE);
            return const_realpart_type(0,0,0,0,Herm,U,S,NonConj);
        }

        inline const_view_type view() const
        { 
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),
                Herm,U,S,NonConj);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),NonConj);
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),
                Herm,U,S,TMV_ConjOf(T,NonConj));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),S,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),TMV_TransOf(S),
                    TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            const StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    newstor,NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_TransOf(newstor),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int n)
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,n));
            if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
                return vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),TMV_ConjOf(T,NonConj) 
                    TMV_FIRSTLAST);
        }

        inline SymMatrixView<T,I> subSymMatrix(int i1, int i2)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return SymMatrixView<T,I>(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),Herm,U,S,NonConj TMV_FIRSTLAST);
        }

        inline SymMatrixView<T,I> subSymMatrix(int i1, int i2, int istep) 
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return SymMatrixView<T,I>(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,U,
                istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag)
        {
            return U==Upper ? 
                uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),
                    dt,S,NonConj TMV_FIRSTLAST) :
                uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag)
        {
            return U==Lower ? 
                lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),
                    dt,S,NonConj TMV_FIRSTLAST) :
                lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(
                    itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                Herm,U,isReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        { 
            // The imaginary part of a Hermitian matrix is anti-symmetric
            // so this is illegal.
            TMVAssert(TMV_FALSE);
            return realpart_type(0,0,0,0,Herm,U,S,NonConj TMV_FIRSTLAST1(0,0) );
        }

        inline view_type view() 
        { 
            return view_type(
                itsm.get(),size(),stepi(),stepj(),
                Herm,U,S,NonConj TMV_FIRSTLAST);
        }

        inline view_type transpose() 
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline view_type conjugate() 
        { 
            return view_type(
                itsm.get(),size(),stepi(),stepj(),
                Herm,U,S,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline view_type adjoint() 
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),TMV_TransOf(S),TMV_ConjOf(T,NonConj) 
                TMV_FIRSTLAST);
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
        TMV_DEPRECATED(const_view_type SubSymMatrix(int i1, int i2) const)
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubSymMatrix(
                int i1, int i2, int istep) const)
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_uppertri_type UpperTri(
                DiagType dt = NonUnitDiag) const)
        { return upperTri(); }
        TMV_DEPRECATED(const_lowertri_type LowerTri(
                DiagType dt = NonUnitDiag) const)
        { return lowerTri(); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }

        TMV_DEPRECATED(rec_type SubMatrix(int i1, int i2, int j1, int j2))
        { return subMatrix(i1,i2,j1,j2); }
        TMV_DEPRECATED(rec_type SubMatrix(
                int i1, int i2, int j1, int j2, int istep, int jstep))
        { return subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_DEPRECATED(vec_type SubVector(
                int i, int j, int istep, int jstep, int s))
        { return subVector(i,j,istep,jstep,s); }
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2))
        { return subSymMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubSymMatrix(int i1, int i2, int istep))
        { return subSymMatrix(i1,i2,istep); }
        TMV_DEPRECATED(uppertri_type UpperTri(DiagType dt = NonUnitDiag))
        { return upperTri(); }
        TMV_DEPRECATED(lowertri_type LowerTri(DiagType dt = NonUnitDiag))
        { return lowerTri(); }
        TMV_DEPRECATED(realpart_type Real())
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag())
        { return imagPart(); }
        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Transpose())
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint())
        { return adjoint(); }


        inline size_t size() const { return itss; }
        inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        inline int stepi() const { return S==RowMajor ? itss : 1; }
        inline int stepj() const { return S==RowMajor ? 1 : itss; }
        inline SymType sym() const { return Herm; }
        inline UpLoType uplo() const { return U; }
        inline StorageType stor() const { return S; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isconj() const { return false; }
        inline bool isherm() const { return true; }
        inline bool issym() const { return isReal(T()); }
        inline bool isupper() const { return U == Upper; }

        inline reference ref(int i, int j)
        {
            if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                if (S == RowMajor) 
                    return TMV_REF(itsm.get() + i*itss + j, NonConj);
                else 
                    return TMV_REF(itsm.get() + j*itss + i, NonConj);
            else 
                if (S == RowMajor) 
                    return TMV_REF(itsm.get() + j*itss + i, Conj);
                else 
                    return TMV_REF(itsm.get() + i*itss + j, Conj);
        }

        inline T cref(int i, int j) const 
        {
            if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
                if (S == RowMajor) 
                    return itsm.get()[i*itss + j];
                else 
                    return itsm.get()[j*itss + i];
            else 
                if (S == RowMajor) 
                    return TMV_CONJ(itsm.get()[j*itss + i]);
                else 
                    return TMV_CONJ(itsm.get()[i*itss + j]);
        }

    protected :

        const size_t itslen;
        AlignedArray<T> itsm;
        const size_t itss;

#ifdef TMVFLDEBUG
    public:
        const T*const _first;
        const T*const _last;
#endif

        template <IndexStyle I2>
        friend void Swap(HermMatrix<T,U,S,I>& m1, HermMatrix<T,U,S,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            m1.itsm.swapWith(m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // HermMatrix

    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   SymMatrixViewOf(m,uplo)
    //   HermMatrixViewOf(m,uplo)
    //   SymMatrixViewOf(mptr,size,uplo,stor)
    //   HermMatrixViewOf(mptr,size,uplo,stor)
    //

    template <class T> 
    inline ConstSymMatrixView<T> SymMatrixViewOf(
        const GenMatrix<T>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymMatrixView<T>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
            Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T,I> SymMatrixViewOf(
        const ConstMatrixView<T,I>& m, UpLoType uplo)
    { 
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymMatrixView<T,I>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
            Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> SymMatrixViewOf(
        const Matrix<T,S,I>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymMatrixView<T,I>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
            Sym,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymMatrixView<T,I> SymMatrixViewOf(
        const MatrixView<T,I>& m, UpLoType uplo)
    { 
        TMVAssert(m.colsize()==m.rowsize());
        return SymMatrixView<T,I>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
            Sym,uplo,m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last)); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> SymMatrixViewOf(
        Matrix<T,S,I>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return SymMatrixView<T,I>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
            Sym,uplo,m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last)); 
    }

    template <class T> 
    inline ConstSymMatrixView<T> HermMatrixViewOf(
        const GenMatrix<T>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return ConstSymMatrixView<T>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
            Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T> HermMatrixViewOf(
        const ConstMatrixView<T,I>& m, UpLoType uplo)
    { 
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return ConstSymMatrixView<T,I>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
            Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> HermMatrixViewOf(
        const Matrix<T,S,I>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return ConstSymMatrixView<T,I>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
            Herm,uplo,m.stor(),m.ct()); 
    }

    template <class T, IndexStyle I> 
    inline SymMatrixView<T> HermMatrixViewOf(
        const MatrixView<T,I>& m, UpLoType uplo)
    { 
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return SymMatrixView<T,I>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
            Herm,uplo,m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last)); 
    }

    template <class T, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> HermMatrixViewOf(
        Matrix<T,S,I>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) || 
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return SymMatrixView<T,I>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
            Herm,uplo,m.stor(),m.ct() TMV_FIRSTLAST1(m._first,m._last)); 
    }

    template <class T> 
    ConstSymMatrixView<T> SymMatrixViewOf(
        const T* vv, size_t size, UpLoType uplo, StorageType stor);

    template <class T> 
    ConstSymMatrixView<T> HermMatrixViewOf(
        const T* vv, size_t size, UpLoType uplo, StorageType stor);

    template <class T> 
    SymMatrixView<T> SymMatrixViewOf(
        T* vv, size_t size, UpLoType uplo, StorageType stor);

    template <class T> 
    SymMatrixView<T> HermMatrixViewOf(
        T* vv, size_t size, UpLoType uplo, StorageType stor);

    //
    // Swap Matrices
    //

    template <class T> 
    inline void Swap(const SymMatrixView<T>& m1, const SymMatrixView<T>& m2)
    { 
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.issym() == m2.issym());
        TMVAssert(m1.isherm() == m2.isherm());
        Swap(m1.upperTri(),m2.upperTri()); 
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(const SymMatrixView<T>& m1, SymMatrix<T,U,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(SymMatrix<T,U,S,I>& m1, const SymMatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
    inline void Swap(SymMatrix<T,U1,S1,I1>& m1, SymMatrix<T,U2,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(const SymMatrixView<T>& m1, HermMatrix<T,U,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline void Swap(HermMatrix<T,U,S,I>& m1, const SymMatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
    inline void Swap(HermMatrix<T,U1,S1,I1>& m1, HermMatrix<T,U2,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }



    //
    // Views:
    //

    template <class T> 
    inline ConstSymMatrixView<T> Transpose(const GenSymMatrix<T>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Transpose(const ConstSymMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Transpose(const SymMatrix<T,U,S,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline SymMatrixView<T,I> Transpose(const SymMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> Transpose(SymMatrix<T,U,S,I>& m)
    { return m.transpose(); }

    template <class T> 
    inline ConstSymMatrixView<T> Conjugate(const GenSymMatrix<T>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Conjugate(const ConstSymMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Conjugate(const SymMatrix<T,U,S,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline SymMatrixView<T,I> Conjugate(const SymMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> Conjugate(SymMatrix<T,U,S,I>& m)
    { return m.conjugate(); }

    template <class T> 
    inline ConstSymMatrixView<T> Adjoint(const GenSymMatrix<T>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Adjoint(const ConstSymMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline ConstSymMatrixView<T,I> Adjoint(const SymMatrix<T,U,S,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline SymMatrixView<T,I> Adjoint(const SymMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline SymMatrixView<T,I> Adjoint(SymMatrix<T,U,S,I>& m)
    { return m.adjoint(); }

    template <class T> 
    inline QuotXS<T,T> Inverse(const GenSymMatrix<T>& m)
    { return m.inverse(); }

    //
    // SymMatrix ==, != SymMatrix
    //

    template <class T1, class T2> 
    inline bool operator==(
        const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
    { return m1.upperTri() == m2.upperTri(); }
    template <class T1, class T2> 
    inline bool operator!=(
        const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& is, auto_ptr<SymMatrix<T,U,S,I> >& m);

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& is, auto_ptr<HermMatrix<T,U,S,I> >& m);

    template <class T> 
    std::istream& operator>>(std::istream& is, const SymMatrixView<T>& m);

    template <class T, UpLoType U, StorageType S, IndexStyle I>
    inline std::istream& operator>>(std::istream& is, SymMatrix<T,U,S,I>& m)
    { return is>>m.view(); }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(std::istream& is, HermMatrix<T,U,S,I>& m)
    { return is>>m.view(); }

} // namespace tmv

#endif
