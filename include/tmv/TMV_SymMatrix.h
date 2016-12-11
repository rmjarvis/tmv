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
// This file defines the TMV SymMatrix and HermMatrix classes.
//
// SymMatrix is used for symmetric matrices, and
// HermMatrix is used for Hermitian matrices.
// For real matrices, these are the same thing:
//     A = A.transpose().
// But for complex, they are different:
//     A_sym = A_sym.transpose()
//     A_herm = A_herm.adjoint()
//
// For these notes, I will always write SymMatrix, but (except where
// otherwise indicated) everything applies the same for Sym and Herm.
// Also, the Views keep track of sym/herm difference with a parameter,
// so it is always a GenSymMatrix, ConstSymMatrixView, or
// SymMatrixView - never Herm in any of these.
//
// Caveat: Complex Hermitian matrices are such that A = At, which
// implies that their diagonal elements are real.  Many routines
// involving HermMatrixes assume the reality of the diagonal.
// However, it is possible to assign a non-real value to a diagonal
// element.  If the user does this, the results are undefined.
//
// As usual, the first template parameter is the type of the data,
// and the optional second template parameter specifies the known
// attributes.  The valid attributs for a SymMatrix are:
// - ColMajor or RoMajor
// - CStyle or FortranStyle
// - Lower or Upper
// The defaults are (ColMajor,CStyle,Lower) if you do not specify
// otherwise.
//
// The storage and index style follow the same meaning as for regular
// Matrices.
// Lower or Upper refers to which triangular half stores the actual data.
//
// Constructors:
//
//    SymMatrix<T,A>(int n)
//        Makes a Symmetric Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SymMatrix<T,A>(int n, T x)
//        Makes a Symmetric Matrix with values all initialized to x
//        For Hermitian matrixces, x must be real.
//
//    SymMatrix<T,A>(const Matrix<T>& m)
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
//    int colsize() const
//    int rowsize() const
//    int size() const
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
//    addToAll(T x)
//        For HermMatrix, x must be real in both setAllTo and addToAll.
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
//          size
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(size-1,0) ... m(size-1,size-1) )
//
//    is >> m
//    is >> CompactIO() >> m
//        Reads m from istream is in either format
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

    template <typename T, typename T1, typename T2>
    class OProdVV;
    template <typename T, typename T1, typename T2>
    class ProdMM;
    template <typename T, typename T1, typename T2>
    class ProdUL;
    template <typename T, typename T1, typename T2>
    class ProdLU;

    template <typename T>
    struct SymCopyHelper // real
    { typedef SymMatrix<T> type; };
    template <typename T>
    struct SymCopyHelper<std::complex<T> > // complex
    // Have to copy to matrix, since don't know whether herm or sym.
    { typedef Matrix<std::complex<T> > type; };

    template <typename T>
    class GenSymMatrix :
        virtual public AssignableToSymMatrix<T>,
        public BaseMatrix<T>,
        public DivHelper<T>
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
        inline ptrdiff_t colsize() const { return size(); }
        inline ptrdiff_t rowsize() const { return size(); }
        using AssignableToSymMatrix<T>::sym;

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            return cref(i,j);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
                return const_vec_type(
                    cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
                    ct());
            else
                return const_vec_type(
                    cptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
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

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-size() && i<=size());
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

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=-size() && i<=size());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            const ptrdiff_t ds = stepi()+stepj();
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

        template <typename T2>
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

        inline void assignToM(MatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            assignToS(SymMatrixViewOf(m2,Upper));
            if (size() > 0)
                m2.lowerTri().offDiag() =
                    m2.upperTri().offDiag().transpose();
        }

        inline void assignToM(MatrixView<CT> m2) const
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

        inline void assignToS(SymMatrixView<RT> m2) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(m2.size() == size());
            if (!isSameAs(m2)) m2.upperTri() = upperTri();
        }

        inline void assignToS(SymMatrixView<CT> m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(isReal(T()) || m2.sym() == sym());
            if (!isSameAs(m2)) m2.upperTri() = upperTri();
        }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            if ( (uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1) )
                return const_rec_type(
                    cptr()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),ct());
            else
                return const_rec_type(
                    cptr()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if ( (uplo()==Upper && i2-j1<=istep) ||
                 (uplo()==Lower && j2-i1<=istep) )
                return const_rec_type(
                    cptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    ct());
            else
                return const_rec_type(
                    cptr()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        bool hasSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const;

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
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

        bool hasSubSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const;

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return const_view_type(
                cptr()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),sym(),uplo(),ct());
        }

        inline const_view_type subSymMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return const_view_type(
                cptr()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),sym(),uplo(),ct());
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            if (uplo() == Upper)
                return const_uppertri_type(
                    cptr(),size(),stepi(),stepj(),dt,ct());
            else
                return const_uppertri_type(
                    cptr(),size(),stepj(),stepi(),dt,
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_uppertri_type unitUpperTri() const
        {
            if (uplo() == Upper)
                return const_uppertri_type(
                    cptr(),size(),stepi(),stepj(),UnitDiag,ct());
            else
                return const_uppertri_type(
                    cptr(),size(),stepj(),stepi(),UnitDiag,
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            if (uplo() == Lower)
                return const_lowertri_type(
                    cptr(),size(),stepi(),stepj(),dt,ct());
            else
                return const_lowertri_type(
                    cptr(),size(),stepj(),stepi(),dt,
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_lowertri_type unitLowerTri() const
        {
            if (uplo() == Lower)
                return const_lowertri_type(
                    cptr(),size(),
                    stepi(),stepj(),UnitDiag,ct());
            else
                return const_lowertri_type(
                    cptr(),size(),stepj(),stepi(),UnitDiag,
                    issym()?ct():TMV_ConjOf(T,ct()));
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), Sym, uplo(), NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(issym());
            // The imaginary part of a Hermitian matrix is anti-symmetric
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                size(),2*stepi(),2*stepj(), Sym,uplo(),NonConj);
        }

        //
        // Views
        //

        inline const_view_type view() const
        {
            return const_view_type(
                cptr(),size(),stepi(),stepj(),sym(),uplo(),ct());
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                cptr(),size(),stepj(),stepi(),sym(),TMV_UTransOf(uplo()),ct());
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                cptr(),size(),stepi(),stepj(),sym(),uplo(),TMV_ConjOf(T,ct()));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                cptr(),size(),stepj(),stepi(),sym(),TMV_UTransOf(uplo()),
                TMV_ConjOf(T,ct()));
        }

        inline nonconst_type nonConst() const
        {
            return SymMatrixView<T>(
                const_cast<T*>(cptr()),size(),stepi(),stepj(),sym(),uplo(),ct()
                TMV_FIRSTLAST1(
                    cptr(),row(colsize()-1,0,colsize()).end().getP()));
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
        { return norm1(); }

        RT maxAbsElement() const;
        RT maxAbs2Element() const;

        RT doNorm2() const;
        inline RT norm2() const
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

        template <typename T1>
        void doMakeInverse(SymMatrixView<T1> sinv) const;

        template <typename T1>
        inline void makeInverse(SymMatrixView<T1> minv) const
        {
            TMVAssert(minv.size() == size());
            TMVAssert(isherm() == minv.isherm());
            TMVAssert(issym() == minv.issym());
            doMakeInverse(minv);
        }

        inline void makeInverse(MatrixView<T> minv) const
        { DivHelper<T>::makeInverse(minv); }

        template <typename T1>
        inline void makeInverse(MatrixView<T1> minv) const
        { DivHelper<T>::makeInverse(minv); }

        template <typename T1, int A>
        inline void makeInverse(Matrix<T1,A>& minv) const
        { DivHelper<T>::makeInverse(minv); }

        template <typename T1, int A>
        inline void makeInverse(SymMatrix<T1,A>& sinv) const
        {
            TMVAssert(issym());
            makeInverse(sinv.view());
        }

        template <typename T1, int A>
        inline void makeInverse(HermMatrix<T1,A>& sinv) const
        {
            TMVAssert(isherm());
            makeInverse(sinv.view());
        }

        QuotXS<T,T> QInverse() const;
        inline QuotXS<T,T> inverse() const
        { return QInverse(); }

        //
        // Division Control
        //

        void setDiv() const;

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
            TMVAssert(divIsLUDiv());
            return static_cast<const SymLDLDiv<T>&>(*this->getDiv());
        }

        inline const HermCHDiv<T>& chd() const
        {
            TMVAssert(isherm());
            divideUsing(CH);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsCHDiv());
            return static_cast<const HermCHDiv<T>&>(*this->getDiv());
        }

        inline const HermSVDiv<T>& svd() const
        {
            TMVAssert(isherm());
            divideUsing(SV);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsHermSVDiv());
            return static_cast<const HermSVDiv<T>&>(*this->getDiv());
        }

        inline const SymSVDiv<T>& symsvd() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(issym());
            divideUsing(SV);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsSymSVDiv());
            return static_cast<const SymSVDiv<T>&>(*this->getDiv());
        }


        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        virtual const T* cptr() const = 0;
        virtual ptrdiff_t stepi() const = 0;
        virtual ptrdiff_t stepj() const = 0;
        virtual UpLoType uplo() const = 0;
        virtual ConjType ct() const = 0;
        virtual inline bool isrm() const { return stepj() == 1; }
        virtual inline bool iscm() const { return stepi() == 1; }
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

        virtual T cref(ptrdiff_t i, ptrdiff_t j) const;

    protected :

        inline const BaseMatrix<T>& getMatrix() const { return *this; }

    private :

        type& operator=(const type&);

        bool divIsLUDiv() const;
        bool divIsCHDiv() const;
        bool divIsHermSVDiv() const;
        bool divIsSymSVDiv() const;

    }; // GenSymMatrix

    template <typename T, int A>
    class ConstSymMatrixView : public GenSymMatrix<T>
    {
    public :

        typedef ConstSymMatrixView<T,A> type;
        typedef GenSymMatrix<T> base;

        inline ConstSymMatrixView(const type& rhs) :
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itssym(rhs.itssym), itsuplo(rhs.itsuplo),
            itsct(rhs.itsct)
        { TMVAssert(Attrib<A>::viewok);  }

        inline ConstSymMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itss(rhs.size()),
            itssi(rhs.stepi()), itssj(rhs.stepj()),
            itssym(rhs.sym()), itsuplo(rhs.uplo()),
            itsct(rhs.ct())
        { TMVAssert(Attrib<A>::viewok);  }

        inline ConstSymMatrixView(
            const T* _m, ptrdiff_t _s, ptrdiff_t _si, ptrdiff_t _sj,
            SymType _sym, UpLoType _uplo, ConjType _ct) :
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itssym(_sym), itsuplo(_uplo), itsct(_ct)
        { TMVAssert(Attrib<A>::viewok);  }

        virtual inline ~ConstSymMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        inline ptrdiff_t size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline ptrdiff_t stepi() const { return itssi; }
        inline ptrdiff_t stepj() const { return itssj; }
        inline SymType sym() const { return itssym; }
        inline UpLoType uplo() const { return itsuplo; }
        inline ConjType ct() const { return itsct; }

    protected :

        const T*const itsm;
        const ptrdiff_t itss;
        const ptrdiff_t itssi;
        const ptrdiff_t itssj;

        const SymType itssym;
        const UpLoType itsuplo;
        const ConjType itsct;

    private :

        type& operator=(const type&);

    }; // ConstSymMatrixView

    template <typename T>
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
            const T* _m, ptrdiff_t _s, ptrdiff_t _si, ptrdiff_t _sj,
            SymType _sym, UpLoType _uplo, ConjType _ct) :
            c_type(_m,_s,_si,_sj,_sym,_uplo,_ct) {}

        virtual inline ~ConstSymMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j>0 && j<=this->size());
            return base::cref(i-1,j-1);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size());
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>0 && j<=this->size());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
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
        // SubMatrix
        //

        bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const;

        bool hasSubVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const;

        bool hasSubSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const;

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

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return base::subVector(i-1,j-1,istep,jstep,n);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return base::subSymMatrix(i1-1,i2);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return base::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        { return base::upperTri(dt); }

        inline const_uppertri_type unitUpperTri() const
        { return base::unitUpperTri(); }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        { return base::lowerTri(dt); }

        inline const_lowertri_type unitLowerTri() const
        { return base::unitLowerTri(); }

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

    }; // FortranStyle ConstSymMatrixView

    template <typename T, int A>
    class SymMatrixView : public GenSymMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymMatrixView<T,A> type;
        typedef GenSymMatrix<T> base;
        typedef SymMatrixView<T,A> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef VectorView<T,A> vec_type;
        typedef MatrixView<T,A> rec_type;
        typedef UpperTriMatrixView<T,A> uppertri_type;
        typedef LowerTriMatrixView<T,A> lowertri_type;
        typedef SymMatrixView<RT,A> realpart_type;
        typedef ConstSymMatrixView<T,A> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,A> const_vec_type;
        typedef ConstMatrixView<T,A> const_rec_type;
        typedef ConstUpperTriMatrixView<T,A> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,A> const_lowertri_type;
        typedef ConstSymMatrixView<RT,A> const_realpart_type;
        typedef typename RefHelper<T>::reference reference;

        //
        // Constructors
        //

        inline SymMatrixView(const type& rhs) :
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itssym(rhs.sym()), itsuplo(rhs.uplo()),
            itsct(rhs.ct()) TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        { TMVAssert(Attrib<A>::viewok);  }

        inline SymMatrixView(
            T* _m, ptrdiff_t _s, ptrdiff_t _si, ptrdiff_t _sj,
            SymType _sym, UpLoType _uplo, ConjType _ct
            TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itssym(_sym), itsuplo(_uplo), itsct(_ct)
            TMV_DEFFIRSTLAST(_first,_last)
        { TMVAssert(Attrib<A>::viewok);  }

        virtual inline ~SymMatrixView()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T()) || m2.sym() == sym());
            if (!this->isSameAs(m2)) upperTri() = m2.upperTri();
            return *this;
        }

        inline type& operator=(const GenSymMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            m2.assignToS(*this);
            return *this;
        }

        inline type& operator=(const GenSymMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(m2.sym() == sym());
            m2.assignToS(*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const GenSymMatrix<T2>& m2)
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(isReal(T2()) || m2.sym() == sym());
            if (!this->isSameAs(m2)) upperTri() = m2.upperTri();
            return *this;
        }

        inline type& operator=(const T& x)
        {
            TMVAssert(issym() || TMV_IMAG(x) == RT(0));
            return setToIdentity(x);
        }

        inline type& operator=(const AssignableToSymMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            m2.assignToS(view());
            return *this;
        }

        inline type& operator=(const AssignableToSymMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(sym() == m2.sym());
            m2.assignToS(view());
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        {
            TMVAssert(size() == m2.size());
            m2.assignToD(DiagMatrixViewOf(diag()));
            upperTri().offDiag().setZero();
            return *this;
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2)
        {
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            TMVAssert(issym());
            m2.assignToD(DiagMatrixViewOf(diag()));
            upperTri().offDiag().setZero();
            TMVAssert(this->isHermOK());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const OProdVV<RT,T2,T2>& opvv)
        {
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(
                    issym() ? opvv.getV2().view() : opvv.getV2().conjugate()));
            TMVAssert(issym() || TMV_IMAG(opvv.getX()) == RT(0));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const OProdVV<CT,T2,T2>& opvv)
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

        template <typename T2>
        inline type& operator=(const ProdMM<RT,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(issym() ? pmm.getM2().transpose() :
                                           pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdMM<CT,T2,T2>& pmm)
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

        template <typename T2>
        inline type& operator=(const ProdUL<RT,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdUL<CT,T2,T2>& pmm)
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

        template <typename T2>
        inline type& operator=(const ProdLU<RT,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(
                    issym() ? pmm.getM2().transpose() : pmm.getM2().adjoint()));
            TMVAssert(issym() || TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdLU<CT,T2,T2>& pmm)
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

        inline reference operator()(ptrdiff_t i,ptrdiff_t j)
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            return ref(i,j);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
                return vec_type(
                    ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
                    ct() TMV_FIRSTLAST);
            else
                return vec_type(
                    ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
                return vec_type(
                    ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
                    ct() TMV_FIRSTLAST);
            else
                return vec_type(
                    ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type diag()
        { return vec_type(ptr(),size(),stepi()+stepj(),ct() TMV_FIRSTLAST); }

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-size() && i<=size());
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

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=-size() && i<=size());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            const ptrdiff_t ds = stepi()+stepj();
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


        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
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

        inline type& setZero()
        { upperTri().setZero(); return *this; }

        inline type& setAllTo(const T& x)
        {
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            upperTri().setAllTo(x); return *this;
        }

        inline type& addToAll(const T& x)
        {
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            upperTri().addToAll(x); return *this;
        }

        inline type& clip(RT thresh)
        { upperTri().clip(thresh); return *this; }

        inline type& transposeSelf()
        { if (!this->issym()) upperTri().conjugateSelf(); return *this; }

        inline type& conjugateSelf()
        { if (isComplex(T())) upperTri().conjugateSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1))
        {
            TMVAssert(TMV_IMAG(x)==RT(0) || this->issym());
            setZero(); diag().setAllTo(x); return *this;
        }

        type& swapRowsCols(ptrdiff_t i1, ptrdiff_t i2);

        type& permuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2);

        type& reversePermuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2);

        inline type& permuteRowsCols(const ptrdiff_t* p)
        { return permuteRowsCols(p,0,size()); }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p)
        { return reversePermuteRowsCols(p,0,size()); }

        //
        // subMatrix
        //

        inline rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            if ( (uplo()==Upper && i2-j1<=1) ||
                 (uplo()==Lower && j2-i1<=1) )
                return rec_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),ct() TMV_FIRSTLAST);
            else
                return rec_type(
                    ptr()+i1*stepj()+j1*stepi(),i2-i1,j2-j1,stepj(),stepi(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if ( (uplo()==Upper && i2-j1<=istep) ||
                 (uplo()==Lower && j2-i1<=istep) )
                return rec_type(
                    ptr()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    ct() TMV_FIRSTLAST);
            else
                return rec_type(
                    ptr()+i1*stepj()+j1*stepi(),(i2-i1)/istep,(j2-j1)/jstep,
                    istep*stepj(),jstep*stepi(),
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline vec_type subVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n)
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

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(base::hasSubSymMatrix(i1,i2,1));
            return view_type(
                ptr()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),sym(),uplo(),ct() TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(base::hasSubSymMatrix(i1,i2,istep));
            return view_type(
                ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),sym(),uplo(), ct() TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag)
        {
            if (uplo() == Upper)
                return uppertri_type(
                    ptr(),size(),stepi(),stepj(),dt,ct() TMV_FIRSTLAST);
            else
                return uppertri_type(
                    ptr(),size(),stepj(),stepi(),dt,
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline uppertri_type unitUpperTri()
        {
            if (uplo() == Upper)
                return uppertri_type(
                    ptr(),size(),stepi(),stepj(),UnitDiag,ct() TMV_FIRSTLAST);
            else
                return uppertri_type(
                    ptr(),size(),stepj(),stepi(),UnitDiag,
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag)
        {
            if (uplo() == Lower)
                return lowertri_type(
                    ptr(),size(),stepi(),stepj(),dt,ct() TMV_FIRSTLAST);
            else
                return lowertri_type(
                    ptr(),size(),stepj(),stepi(),dt,
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri()
        {
            if (uplo() == Lower)
                return lowertri_type(
                    ptr(),size(),stepi(),stepj(),UnitDiag,ct() TMV_FIRSTLAST);
            else
                return lowertri_type(
                    ptr(),size(),stepj(),stepi(),UnitDiag,
                    this->issym()?ct():TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), sym(),uplo(),NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        {
            TMVAssert(isComplex(T()));
            TMVAssert(this->issym());
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1,
                size(),2*stepi(),2*stepj(),sym(),uplo(),NonConj
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
                ptr(),size(),stepj(),stepi(),sym(),TMV_UTransOf(uplo()),ct()
                TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return view_type(
                ptr(),size(),stepi(),stepj(),sym(),uplo(),TMV_ConjOf(T,ct())
                TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return view_type(
                ptr(),size(),stepj(),stepi(),sym(),TMV_UTransOf(uplo()),
                TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline const_rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::subMatrix(i1,i2,j1,j2); }
        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::subMatrix(i1,i2,j1,j2,istep,jstep); }
        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
        { return base::subVector(i,j,istep,jstep,n); }
        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::subSymMatrix(i1,i2); }
        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return base::subSymMatrix(i1,i2,istep); }
        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        { return base::upperTri(dt); }
        inline const_uppertri_type unitUpperTri() const
        { return base::unitUpperTri(); }
        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        { return base::lowerTri(dt); }
        inline const_lowertri_type unitLowerTri() const
        { return base::unitLowerTri(); }
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

        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        inline ptrdiff_t size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }
        inline ptrdiff_t stepi() const { return itssi; }
        inline ptrdiff_t stepj() const { return itssj; }
        inline SymType sym() const { return itssym; }
        inline UpLoType uplo() const { return itsuplo; }
        inline ConjType ct() const { return itsct; }
        using base::issym;
        using base::isherm;

        reference ref(ptrdiff_t i, ptrdiff_t j);

    protected :

        T*const itsm;
        const ptrdiff_t itss;
        const ptrdiff_t itssi;
        const ptrdiff_t itssj;

        const SymType itssym;
        const UpLoType itsuplo;
        const ConjType itsct;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected:
#endif

    }; // SymMatrixView

    template <typename T>
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
        typedef ConstSymMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_lowertri_type;
        typedef ConstSymMatrixView<RT,FortranStyle> const_realpart_type;
        typedef ConstSymMatrixView<T,FortranStyle> const_type;
        typedef typename RefHelper<T>::reference reference;

        //
        // Constructors
        //

        inline SymMatrixView(const type& rhs) : c_type(rhs) {}

        inline SymMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline SymMatrixView(
            T* _m, ptrdiff_t _s, ptrdiff_t _si, ptrdiff_t _sj,
            SymType _sym, UpLoType _uplo, ConjType _ct
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_s,_si,_sj,_sym,_uplo,_ct
                   TMV_FIRSTLAST1(_first,_last) ) {}

        virtual inline ~SymMatrixView() {}

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const c_type& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenSymMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenSymMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        template <typename T2>
        inline type& operator=(const GenSymMatrix<T2>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const T& x)
        { c_type::operator=(x); return *this; }

        inline type& operator=(const AssignableToSymMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const AssignableToSymMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenDiagMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        template <typename T2>
        inline type& operator=(const OProdVV<T,T2,T2>& opvv)
        { c_type::operator=(opvv); return *this; }

        template <typename T2>
        inline type& operator=(const ProdMM<T,T2,T2>& pmm)
        { c_type::operator=(pmm); return *this; }

        template <typename T2>
        inline type& operator=(const ProdLU<T,T2,T2>& pmm)
        { c_type::operator=(pmm); return *this; }

        template <typename T2>
        inline type& operator=(const ProdUL<T,T2,T2>& pmm)
        { c_type::operator=(pmm); return *this; }


        //
        // Access
        //

        inline reference operator()(ptrdiff_t i,ptrdiff_t j)
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j>0 && j<=this->size());
            return c_type::ref(i-1,j-1);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size());
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(j>0 && j<=this->size());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
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
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j>0 && j<=this->size());
            return c_type::cref(i-1,j-1);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>0 && i<=this->size());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->size());
            return c_type::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>0 && j<=this->size());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
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

        inline type& conjugateSelf()
        { c_type::conjugateSelf(); return *this; }

        inline type& transposeSelf()
        { c_type::transposeSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1))
        { c_type::setToIdentity(x); return *this; }

        inline type& swapRowsCols(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1<=this->size());
            TMVAssert(i2>0 && i2<=this->size());
            c_type::swapRowsCols(i1-1,i2-1);
            return *this;
        }

        inline type& permuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
            c_type::permuteRowsCols(p,i1-1,i2);
            return *this;
        }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->size());
            c_type::reversePermuteRowsCols(p,i1-1,i2);
            return *this;
        }

        inline type& permuteRowsCols(const ptrdiff_t* p)
        { c_type::permuteRowsCols(p); return *this; }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p)
        { c_type::reversePermuteRowsCols(p); return *this; }

        //
        // subMatrix
        //

        inline bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return const_type(*this).hasSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline bool hasSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,n); }


        inline bool hasSubSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return const_type(*this).hasSubSymMatrix(i1,i2,istep); }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
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

        inline vec_type subVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n)
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return c_type::subVector(i-1,j-1,istep,jstep,n);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return c_type::subSymMatrix(i1-1,i2);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return c_type::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag)
        { return c_type::upperTri(dt); }

        inline uppertri_type unitUpperTri()
        { return c_type::unitUpperTri(); }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag)
        { return c_type::lowerTri(dt); }

        inline lowertri_type unitLowerTri()
        { return c_type::unitLowerTri(); }

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

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::subMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,n));
            return c_type::subVector(i-1,j-1,istep,jstep,n);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,1));
            return c_type::subSymMatrix(i1-1,i2);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(hasSubSymMatrix(i1,i2,istep));
            return c_type::subSymMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        { return c_type::upperTri(dt); }

        inline const_uppertri_type unitUpperTri() const
        { return c_type::unitUpperTri(); }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        { return c_type::lowerTri(dt); }

        inline const_lowertri_type unitLowerTri() const
        { return c_type::unitLowerTri(); }

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

    }; // FortranStyle SymMatrixView

    template <typename T, int A>
    class SymMatrix : public GenSymMatrix<T>
    {
    public:

        enum { S = A & AllStorageType };
        enum { U = A & Upper };
        enum { I = A & FortranStyle };

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SymMatrix<T,A> type;
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

        explicit inline SymMatrix(ptrdiff_t _size=0) : NEW_SIZE(_size)
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(_size >= 0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline SymMatrix(ptrdiff_t _size, const T& x) : NEW_SIZE(_size)
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(_size >= 0);
            setAllTo(x);
        }

        inline SymMatrix(const type& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        inline SymMatrix(const GenSymMatrix<RT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            if (rhs.issym())
                rhs.assignToS(view());
            else {
                if (uplo() == Upper)
                    upperTri() = rhs.upperTri();
                else
                    lowerTri() = rhs.lowerTri();
            }
        }

        inline SymMatrix(const GenSymMatrix<CT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(isComplex(T()));
            if (rhs.issym())
                rhs.assignToS(view());
            else {
                if (uplo() == Upper)
                    upperTri() = rhs.upperTri();
                else
                    lowerTri() = rhs.lowerTri();
            }
        }

        template <typename T2>
        inline SymMatrix(const GenSymMatrix<T2>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            if (uplo() == Upper)
                upperTri() = rhs.upperTri();
            else
                lowerTri() = rhs.lowerTri();
        }

        template <typename T2>
        inline SymMatrix(const GenMatrix<T2>& rhs) :
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(rhs.isSquare());
            if (uplo() == Upper)
                upperTri() = rhs.upperTri();
            else
                lowerTri() = rhs.lowerTri();
        }

        inline SymMatrix(const AssignableToSymMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(m2.issym());
            m2.assignToS(view());
        }

        inline SymMatrix(const AssignableToSymMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(isComplex(T()));
            TMVAssert(m2.issym());
            m2.assignToS(view());
        }

        inline explicit SymMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(size() == m2.size());
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        inline explicit SymMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(isComplex(T()));
            TMVAssert(size() == m2.size());
            setZero();
            m2.assignToD(DiagMatrixViewOf(diag()));
        }

        template <typename T2>
        inline SymMatrix(const OProdVV<T,T2,T2>& opvv) :
            NEW_SIZE(opvv.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(opvv.colsize() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().view()));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
        }

        template <typename T2>
        inline SymMatrix(const ProdMM<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <typename T2>
        inline SymMatrix(const ProdUL<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <typename T2>
        inline SymMatrix(const ProdLU<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

#undef NEW_SIZE

        virtual inline ~SymMatrix()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
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

        template <typename T2>
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

        template <typename T2>
        inline type& operator=(const OProdVV<T,T2,T2>& opvv)
        {
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().view()));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdMM<T,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdUL<T,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().transpose()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <typename T2>
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

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            if (I==int(CStyle)) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                return cref(i,j);
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                return cref(i-1,j-1);
            }
        }

        inline T& operator()(ptrdiff_t i, ptrdiff_t j)
        {
            if (I==int(CStyle)) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                return ref(i,j);
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            if ((uplo()==Upper && i-j1<=0) || (uplo()==Lower && j2-i<=1))
                return const_vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj);
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),NonConj);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            if ((uplo()==Upper && i2-j<=1) || (uplo()==Lower && j-i1<=0))
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

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-size() && i<=size());
            i = std::abs(i);
            return const_vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi()),
                size()-i,stepi()+stepj(),NonConj);
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=-size() && i<=size());
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            i = std::abs(i);
            const ptrdiff_t ds = stepi()+stepj();
            return const_vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,NonConj);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            if ((uplo()==Upper && i-j1<=0) || (uplo()==Lower && j2-i<=1))
                return vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            if ((uplo()==Upper && i2-j<=1) || (uplo()==Lower && j-i1<=0))
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

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-size() && i<=size());
            i = std::abs(i);
            TMVAssert(i<=size());
            return vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi()),
                size()-i,stepi()+stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=-size() && i<=size());
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            i = std::abs(i);
            const ptrdiff_t ds = stepi()+stepj();
            return vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,NonConj TMV_FIRSTLAST);
        }


        //
        // Modifying Functions
        //

        inline type& setZero()
        { std::fill_n(itsm.get(),itslen,T(0)); return *this; }

        inline type& setAllTo(const T& x)
        { upperTri().setAllTo(x); return *this; }

        inline type& addToAll(const T& x)
        { upperTri().addToAll(x); return *this; }

        inline type& clip(RT thresh)
        { upperTri().clip(thresh); return *this; }

        inline type& conjugateSelf()
        { if (isComplex(T())) upperTri().conjugateSelf(); return *this; }

        inline type& transposeSelf()
        { return *this; }

        inline type& setToIdentity(const T& x=T(1))
        { setZero(); diag().setAllTo(x); return *this; }

        inline type& swapRowsCols(ptrdiff_t i1, ptrdiff_t i2)
        { view().swapRowsCols(i1,i2); return *this; }

        inline type& permuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        { view().permuteRowsCols(p,i1,i2); return *this; }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        { view().reversePermuteRowsCols(p,i1,i2); return *this; }

        inline type& permuteRowsCols(const ptrdiff_t* p)
        { view().permuteRowsCols(p); return *this; }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p)
        { view().reversePermuteRowsCols(p); return *this; }

        //
        // subMatrix
        //

        inline const_rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==int(FortranStyle)) { --i1; --j1; }
            if ( (uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1) )
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),NonConj);
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ( (uplo()==Upper && i2-j1<=istep) ||
                 (uplo()==Lower && j2-i1<=jstep) )
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    NonConj);
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==int(FortranStyle)) { --i; --j; }
            if ((uplo()==Upper && i-j<=0) || (uplo()==Lower && j-i<=0))
                return const_vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj);
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),NonConj);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==int(FortranStyle)) { --i1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),Sym,uplo(),NonConj);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==int(FortranStyle)) { --i1; i2+=istep-1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Sym,uplo(),
                NonConj);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            return uplo()==Upper ?
                const_uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj) :
                const_uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),dt,NonConj);
        }

        inline const_uppertri_type unitUpperTri() const
        {
            return uplo()==Upper ?
                const_uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj) :
                const_uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),UnitDiag,NonConj);
        }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            return uplo()==Lower ?
                const_lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj) :
                const_lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),dt,NonConj);
        }

        inline const_lowertri_type unitLowerTri() const
        {
            return uplo()==Lower ?
                const_lowertri_type(
                    itsm.get(),size(),
                    stepi(),stepj(),UnitDiag,NonConj) :
                const_lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),UnitDiag,NonConj);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), Sym,uplo(),NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get())+1,size(),
                2*stepi(),2*stepj(),Sym,uplo(),NonConj);
        }

        inline const_view_type view() const
        {
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),Sym,uplo(),NonConj);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),Sym,TMV_UTransOf(U),NonConj);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),Sym,uplo(),
                TMV_ConjOf(T,NonConj));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),TMV_ConjOf(T,NonConj));
        }

        inline rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==int(FortranStyle)) { --i1; --j1; }
            if ( (uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1) )
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),NonConj TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ( (uplo()==Upper && i2-j1<=istep) ||
                 (uplo()==Lower && j2-i1<=jstep) )
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    NonConj TMV_FIRSTLAST);
        }

        inline vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,n));
            if (I==int(FortranStyle)) { --i; --j; }
            if ((uplo()==Upper && i-j<=0) || (uplo()==Lower && j-i<=0))
                return vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==int(FortranStyle)) { --i1; }
            return view_type(
                itsm.get()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),Sym,uplo(),NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==int(FortranStyle)) { --i1; i2+=istep-1; }
            return view_type(
                itsm.get()+i1*(stepi()+stepj()),(i2-i1)/istep,
                istep*stepi(),istep*stepj(),Sym,uplo(), NonConj TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag)
        {
            return uplo()==Upper ?
                uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj
                    TMV_FIRSTLAST) :
                uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),dt,NonConj
                    TMV_FIRSTLAST);
        }

        inline uppertri_type unitUpperTri()
        {
            return uplo()==Upper ?
                uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj
                    TMV_FIRSTLAST) :
                uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),UnitDiag,NonConj
                    TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag)
        {
            return uplo()==Lower ?
                lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj
                    TMV_FIRSTLAST) :
                lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),dt,NonConj
                    TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri()
        {
            return uplo()==Lower ?
                lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj
                    TMV_FIRSTLAST) :
                lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),UnitDiag,NonConj
                    TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), Sym,uplo(),NonConj
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
                2*stepi(),2*stepj(),Sym,uplo(),NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline view_type view()
        {
            return view_type(
                itsm.get(),size(),stepi(),stepj(),Sym,uplo(),NonConj
                TMV_FIRSTLAST);
        }

        inline view_type transpose()
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),NonConj TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return view_type(
                itsm.get(),size(),stepi(),stepj(),
                Sym,uplo(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Sym,TMV_UTransOf(U),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        inline ptrdiff_t size() const { return itss; }
        inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        inline ptrdiff_t stepi() const { return S==int(RowMajor) ? itss : 1; }
        inline ptrdiff_t stepj() const { return S==int(RowMajor) ? 1 : itss; }
        inline SymType sym() const { return Sym; }
        inline UpLoType uplo() const { return static_cast<UpLoType>(U); }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==int(RowMajor); }
        inline bool iscm() const { return S==int(ColMajor); }
        inline bool isconj() const { return false; }
        inline bool isherm() const { return isReal(T()); }
        inline bool issym() const { return true; }
        inline bool isupper() const { return U==int(Upper); }

        inline T& ref(ptrdiff_t i, ptrdiff_t j)
        {
            if ((uplo()==Upper && i <= j) || (uplo()==Lower && i>=j))
                return itsm.get()[i*stepi() + j*stepj()];
            else
                return itsm.get()[j*stepi() + i*stepj()];
        }

        inline T cref(ptrdiff_t i, ptrdiff_t j) const
        {
            if ((uplo()==Upper && i <= j) || (uplo()==Lower && i>=j))
                return itsm.get()[i*stepi() + j*stepj()];
            else
                return itsm.get()[j*stepi() + i*stepj()];
        }

        inline void resize(ptrdiff_t s)
        {
            TMVAssert(s >= 0);
            itslen = s*s;
            itsm.resize(itslen);
            itss = s;
            DivHelper<T>::resetDivType();
#ifdef TMVFLDEBUG
            _first = itsm.get();
            _last = _first+itslen;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

    protected :

        ptrdiff_t itslen;
        AlignedArray<T> itsm;
        ptrdiff_t itss;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

        friend void Swap(SymMatrix<T,A>& m1, SymMatrix<T,A>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            m1.itsm.swapWith(m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // SymMatrix

    template <typename T, int A>
    class HermMatrix : public GenSymMatrix<T>
    {
    public:

        enum { S = A & AllStorageType };
        enum { U = A & Upper };
        enum { I = A & FortranStyle };

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef HermMatrix<T,A> type;
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
        typedef typename RefHelper<T>::reference reference;

        //
        // Constructors
        //

#define NEW_SIZE(s) itslen((s)*(s)), itsm(itslen), itss(s)  \
        TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

        explicit inline HermMatrix(ptrdiff_t _size=0) : NEW_SIZE(_size)
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(_size >= 0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#else
            if (isComplex(T())) diag().imagPart().setZero();
#endif
        }

        HermMatrix(ptrdiff_t _size, const RT& x) : NEW_SIZE(_size)
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(_size >= 0);
            setAllTo(x);
        }

        HermMatrix(const type& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        HermMatrix(const GenSymMatrix<RT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            if (rhs.isherm()) rhs.assignToS(view());
            else if (uplo() == Upper) upperTri() = rhs.upperTri();
            else lowerTri() = rhs.lowerTri();
        }

        HermMatrix(const GenSymMatrix<CT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(isComplex(T()));
            if (rhs.isherm())
                rhs.assignToS(view());
            else {
                if (uplo() == Upper) upperTri() = rhs.upperTri();
                else lowerTri() = rhs.lowerTri();
                diag().imagPart().setZero();
            }
        }

        template <typename T2>
        inline HermMatrix(const GenSymMatrix<T2>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            if (uplo() == Upper) upperTri() = rhs.upperTri();
            else lowerTri() = rhs.lowerTri();
            if (isComplex(T()) && isComplex(T2()) && rhs.issym())
                diag().imagPart().setZero();
        }

        template <typename T2>
        inline HermMatrix(const GenMatrix<T2>& rhs) : NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(rhs.isSquare());
            if (uplo() == Upper)
                upperTri() = rhs.upperTri();
            else
                lowerTri() = rhs.lowerTri();
            if (isComplex(T()) && isComplex(T2())) diag().imagPart().setZero();
        }

        inline HermMatrix(const AssignableToSymMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(m2.isherm());
            m2.assignToS(view());
        }

        inline HermMatrix(const AssignableToSymMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(isComplex(T()));
            TMVAssert(m2.isherm());
            m2.assignToS(view());
        }

        inline explicit HermMatrix(const GenDiagMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            m2.assignToD(DiagMatrixViewOf(diag()));
            upperTri().offDiag().setZero();
        }

        inline explicit HermMatrix(const GenDiagMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            m2.assignToD(DiagMatrixViewOf(diag().realPart()));
            upperTri().offDiag().setZero();
            if (isComplex(T())) diag().imagPart().setZero();
        }

        template <typename T2>
        inline HermMatrix(const OProdVV<T,T2,T2>& opvv) :
            NEW_SIZE(opvv.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(opvv.colsize() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().conjugate()));
            TMVAssert(TMV_IMAG(opvv.getX()) == RT(0));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
        }

        template <typename T2>
        inline HermMatrix(const ProdMM<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            TMVAssert(TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <typename T2>
        inline HermMatrix(const ProdUL<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            TMVAssert(TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

        template <typename T2>
        inline HermMatrix(const ProdLU<T,T2,T2>& pmm) :
            NEW_SIZE(pmm.colsize())
        {
            TMVAssert(Attrib<A>::symmatrixok);
            TMVAssert(pmm.colsize() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            TMVAssert(TMV_IMAG(pmm.getX()) == RT(0));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
        }

#undef NEW_SIZE

        virtual inline ~HermMatrix()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
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

        template <typename T2>
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

        template <typename T2>
        inline type& operator=(const OProdVV<T,T2,T2>& opvv)
        {
            TMVAssert(size() == opvv.colsize());
            TMVAssert(size() == opvv.rowsize());
            TMVAssert(opvv.getV1().isSameAs(opvv.getV2().conjugate()));
            Rank1Update<false>(opvv.getX(),opvv.getV1(),view());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdMM<T,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const ProdUL<T,T2,T2>& pmm)
        {
            TMVAssert(size() == pmm.colsize());
            TMVAssert(size() == pmm.rowsize());
            TMVAssert(pmm.getM1().isSameAs(pmm.getM2().adjoint()));
            RankKUpdate<false>(pmm.getX(),pmm.getM1(),view());
            return *this;
        }

        template <typename T2>
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

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            if (I==int(CStyle)) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                return cref(i,j);
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                return cref(i-1,j-1);
            }
        }

        inline reference operator()(ptrdiff_t i, ptrdiff_t j)
        {
            if (I==int(CStyle)) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                return ref(i,j);
            } else {
                TMVAssert(i>0 && i<=size());
                TMVAssert(j>0 && j<=size());
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            if ((uplo()==Upper && i-j1<=0) || (uplo()==Lower && j2-i<=1))
                return const_vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj);
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),TMV_ConjOf(T,NonConj));
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            if ((uplo()==Upper && i2-j<=1) || (uplo()==Lower && j-i1<=0))
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

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-size() && i<=size());
            ConjType newct =
                ((i>0) == (uplo()==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            TMVAssert(i<=size());
            return const_vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi()),
                size()-i,stepi()+stepj(),newct);
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=-size() && i<=size());
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            ConjType newct =
                ((i>0) == (uplo()==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            const ptrdiff_t ds = stepi()+stepj();
            return const_vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi())+j1*ds,
                j2-j1,ds,newct);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(i>0 && i<=size());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size());
                --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size());
            }
            if ((uplo()==Upper && i-j1<=0) || (uplo()==Lower && j2-i<=1))
                return vec_type(
                    itsm.get()+i*stepi()+j1*stepj(),
                    j2-j1,stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm.get()+i*stepj()+j1*stepi(),
                    j2-j1,stepi(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j>0 && j<=size());
                --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=size());
                --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=size());
            }
            if ((uplo()==Upper && i2-j<=1) || (uplo()==Lower && j-i1<=0))
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

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-size() && i<=size());
            ConjType newct =
                ((i>0) == (uplo()==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            TMVAssert(i<=size());
            return vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi()),size()-i,
                stepi()+stepj(),newct TMV_FIRSTLAST);
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=-size() && i<=size());
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=size()-std::abs(i));
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=size()-std::abs(i));
            }
            ConjType newct =
                ((i>0) == (uplo()==Upper)) ? NonConj : TMV_ConjOf(T,NonConj);
            i = std::abs(i);
            const ptrdiff_t ds = stepi()+stepj();
            return vec_type(
                itsm.get()+i*(uplo()==Upper?stepj():stepi())+j1*ds,
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

        inline type& addToAll(const T& x)
        {
            TMVAssert(TMV_IMAG(x) == RT(0));
            upperTri().addToAll(x);
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

        inline type& swapRowsCols(ptrdiff_t i1, ptrdiff_t i2)
        { view().swapRowsCols(i1,i2); return *this; }

        inline type& permuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        { view().permuteRowsCols(p,i1,i2); return *this; }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        { view().reversePermuteRowsCols(p,i1,i2); return *this; }

        inline type& permuteRowsCols(const ptrdiff_t* p)
        { view().permuteRowsCols(p); return *this; }

        inline type& reversePermuteRowsCols(const ptrdiff_t* p)
        { view().reversePermuteRowsCols(p); return *this; }

        //
        // subMatrix
        //

        inline const_rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==int(FortranStyle)) { --i1; --j1; }
            if ((uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1))
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),i2-i1,j2-j1,
                    stepj(),stepi(),TMV_ConjOf(T,NonConj));
        }

        inline const_rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ((uplo()==Upper && i2-j1<=istep) ||
                (uplo()==Lower && j2-i1<=jstep))
                return const_rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    NonConj);
            else
                return const_rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_ConjOf(T,NonConj));
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,n));
            if ((uplo()==Upper && i-j<=0) || (uplo()==Lower && j-i<=0))
                return const_vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj);
            else
                return const_vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),TMV_ConjOf(T,NonConj));
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==int(FortranStyle)) { --i1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),Herm,uplo(),NonConj);
        }

        inline const_view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==int(FortranStyle)) { --i1; i2+=istep-1; }
            return const_view_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,uplo(), NonConj);
        }

        inline const_uppertri_type upperTri(DiagType dt = NonUnitDiag) const
        {
            return uplo()==Upper ?
                const_uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj) :
                const_uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_ConjOf(T,NonConj));
        }

        inline const_uppertri_type unitUpperTri() const
        {
            return uplo()==Upper ?
                const_uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj) :
                const_uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    UnitDiag,TMV_ConjOf(T,NonConj));
        }

        inline const_lowertri_type lowerTri(DiagType dt = NonUnitDiag) const
        {
            return uplo()==Lower ?
                const_lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj) :
                const_lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_ConjOf(T,NonConj));
        }

        inline const_lowertri_type unitLowerTri() const
        {
            return uplo()==Lower ?
                const_lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj) :
                const_lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    UnitDiag,TMV_ConjOf(T,NonConj));
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), Herm,uplo(),NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            // The imaginary part of a Hermitian matrix is anti-symmetric
            // so this is illegal.
            TMVAssert(TMV_FALSE);
            return const_realpart_type(0,0,0,0,Herm,uplo(),NonConj);
        }

        inline const_view_type view() const
        {
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),Herm,uplo(),NonConj);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),NonConj);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                itsm.get(),size(),stepi(),stepj(),
                Herm,uplo(),TMV_ConjOf(T,NonConj));
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),TMV_ConjOf(T,NonConj));
        }

        inline rec_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==int(FortranStyle)) { --i1; --j1; }
            if ((uplo()==Upper && i2-j1<=1) || (uplo()==Lower && j2-i1<=1))
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    i2-i1,j2-j1,stepi(),stepj(),NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    i2-i1,j2-j1,stepj(),stepi(),
                    TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==int(FortranStyle)) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            if ((uplo()==Upper && i2-j1<=istep) ||
                (uplo()==Lower && j2-i1<=jstep))
                return rec_type(
                    itsm.get()+i1*stepi()+j1*stepj(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                    NonConj TMV_FIRSTLAST);
            else
                return rec_type(
                    itsm.get()+i1*stepj()+j1*stepi(),
                    (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
                    TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n)
        {
            TMVAssert(base::hasSubVector(i,j,istep,jstep,n));
            if ((uplo()==Upper && i-j<=0) || (uplo()==Lower && j-i<=0))
                return vec_type(
                    itsm.get()+i*stepi()+j*stepj(),n,
                    istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
            else
                return vec_type(
                    itsm.get()+i*stepj()+j*stepi(),n,
                    istep*stepj()+jstep*stepi(),TMV_ConjOf(T,NonConj)
                    TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,1));
            if (I==int(FortranStyle)) { --i1; }
            return view_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),Herm,uplo(),NonConj TMV_FIRSTLAST);
        }

        inline view_type subSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(view().hasSubSymMatrix(i1,i2,istep));
            if (I==int(FortranStyle)) { --i1; i2+=istep-1; }
            return view_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,uplo(), NonConj
                TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt = NonUnitDiag)
        {
            return uplo()==Upper ?
                uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj
                    TMV_FIRSTLAST) :
                uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),dt,TMV_ConjOf(T,NonConj)
                    TMV_FIRSTLAST);
        }

        inline uppertri_type unitUpperTri()
        {
            return uplo()==Upper ?
                uppertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj
                    TMV_FIRSTLAST) :
                uppertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    UnitDiag,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt = NonUnitDiag)
        {
            return uplo()==Lower ?
                lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),dt,NonConj
                    TMV_FIRSTLAST) :
                lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    dt,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri()
        {
            return uplo()==Lower ?
                lowertri_type(
                    itsm.get(),size(),stepi(),stepj(),UnitDiag,NonConj
                    TMV_FIRSTLAST) :
                lowertri_type(
                    itsm.get(),size(),stepj(),stepi(),
                    UnitDiag,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(
                    itsm.get()),size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), Herm,uplo(),NonConj
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
            return realpart_type(0,0,0,0,Herm,uplo(),NonConj
                                 TMV_FIRSTLAST1(0,0) );
        }

        inline view_type view()
        {
            return view_type(
                itsm.get(),size(),stepi(),stepj(),Herm,uplo(),NonConj
                TMV_FIRSTLAST);
        }

        inline view_type transpose()
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),Herm,TMV_UTransOf(U),NonConj
                TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return view_type(
                itsm.get(),size(),stepi(),stepj(),
                Herm,uplo(),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return view_type(
                itsm.get(),size(),stepj(),stepi(),
                Herm,TMV_UTransOf(U),TMV_ConjOf(T,NonConj)
                TMV_FIRSTLAST);
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        inline ptrdiff_t size() const { return itss; }
        inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        inline ptrdiff_t stepi() const { return S==int(RowMajor) ? itss : 1; }
        inline ptrdiff_t stepj() const { return S==int(RowMajor) ? 1 : itss; }
        inline SymType sym() const { return Herm; }
        inline UpLoType uplo() const { return static_cast<UpLoType>(U); }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==int(RowMajor); }
        inline bool iscm() const { return S==int(ColMajor); }
        inline bool isconj() const { return false; }
        inline bool isherm() const { return true; }
        inline bool issym() const { return isReal(T()); }
        inline bool isupper() const { return U==int(Upper); }

        inline reference ref(ptrdiff_t i, ptrdiff_t j)
        {
            if ((uplo()==Upper && i <= j) || (uplo()==Lower && i>=j))
                return RefHelper<T>::makeRef(
                    itsm.get() + i*stepi() + j*stepj(),NonConj);
            else
                return RefHelper<T>::makeRef(
                    itsm.get() + j*stepi() + i*stepj(),Conj);
        }

        inline T cref(ptrdiff_t i, ptrdiff_t j) const
        {
            if ((uplo()==Upper && i <= j) || (uplo()==Lower && i>=j))
                return itsm.get()[i*stepi() + j*stepj()];
            else
                return TMV_CONJ(itsm.get()[j*stepi() + i*stepj()]);
        }

        inline void resize(ptrdiff_t s)
        {
            TMVAssert(s >= 0);
            itslen = s*s;
            itsm.resize(itslen);
            itss = s;
            DivHelper<T>::resetDivType();
#ifdef TMVFLDEBUG
            _first = itsm.get();
            _last = _first+itslen;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#else
            if (isComplex(T())) diag().imagPart().setZero();
#endif
        }

    protected :

        ptrdiff_t itslen;
        AlignedArray<T> itsm;
        ptrdiff_t itss;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

        friend void Swap(HermMatrix<T,A>& m1, HermMatrix<T,A>& m2)
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
    //   SymMatrixViewOf(mptr,size,uplo,stepi,stepj)
    //   HermMatrixViewOf(mptr,size,uplo,stepi,stepj)
    //

    template <typename T>
    inline ConstSymMatrixView<T> SymMatrixViewOf(
        const GenMatrix<T>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymMatrixView<T>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),Sym,uplo,m.ct());
    }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A> SymMatrixViewOf(
        const ConstMatrixView<T,A>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymMatrixView<T,A>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),Sym,uplo,m.ct());
    }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> SymMatrixViewOf(
        const Matrix<T,A>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return ConstSymMatrixView<T,A&FortranStyle>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),Sym,uplo,m.ct());
    }

    template <typename T, int A>
    inline SymMatrixView<T,A> SymMatrixViewOf(
        MatrixView<T,A> m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return SymMatrixView<T,A>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),Sym,uplo,m.ct()
            TMV_FIRSTLAST1(m._first,m._last));
    }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> SymMatrixViewOf(
        Matrix<T,A>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        return SymMatrixView<T,A&FortranStyle>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),Sym,uplo,m.ct()
            TMV_FIRSTLAST1(m._first,m._last));
    }

    template <typename T>
    inline ConstSymMatrixView<T> HermMatrixViewOf(
        const GenMatrix<T>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) ||
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return ConstSymMatrixView<T>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),Herm,uplo,m.ct());
    }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A> HermMatrixViewOf(
        const ConstMatrixView<T,A>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) ||
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return ConstSymMatrixView<T,A>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),Herm,uplo,m.ct());
    }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> HermMatrixViewOf(
        const Matrix<T,A>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) ||
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return ConstSymMatrixView<T,A&FortranStyle>(
            m.cptr(),m.rowsize(),m.stepi(),m.stepj(),Herm,uplo,m.ct());
    }

    template <typename T, int A>
    inline SymMatrixView<T,A> HermMatrixViewOf(
        MatrixView<T,A> m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) ||
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return SymMatrixView<T,A>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),Herm,uplo,m.ct()
            TMV_FIRSTLAST1(m._first,m._last));
    }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> HermMatrixViewOf(
        Matrix<T,A>& m, UpLoType uplo)
    {
        TMVAssert(m.colsize()==m.rowsize());
        TMVAssert(isReal(T()) ||
                  m.diag().imagPart().normInf() == TMV_RealType(T)(0));
        return SymMatrixView<T,A&FortranStyle>(
            m.ptr(),m.rowsize(),m.stepi(),m.stepj(),Herm,uplo,m.ct()
            TMV_FIRSTLAST1(m._first,m._last));
    }

    template <typename T>
    inline ConstSymMatrixView<T> SymMatrixViewOf(
        const T* m, ptrdiff_t size, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(size>=0);
        const ptrdiff_t stepi = stor == RowMajor ? size : 1;
        const ptrdiff_t stepj = stor == RowMajor ? 1 : size;
        return ConstSymMatrixView<T>(
            m,size,stepi,stepj,Sym,uplo,NonConj);
    }

    template <typename T>
    inline ConstSymMatrixView<T> HermMatrixViewOf(
        const T* m, ptrdiff_t size, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(size>=0);
        const ptrdiff_t stepi = stor == RowMajor ? size : 1;
        const ptrdiff_t stepj = stor == RowMajor ? 1 : size;
        return ConstSymMatrixView<T>(m,size,stepi,stepj,Herm,uplo,NonConj);
    }

    template <typename T>
    inline SymMatrixView<T> SymMatrixViewOf(
        T* m, ptrdiff_t size, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(size>=0);
        const ptrdiff_t stepi = stor == RowMajor ? size : 1;
        const ptrdiff_t stepj = stor == RowMajor ? 1 : size;
        return SymMatrixView<T>(
            m,size,stepi,stepj,Sym,uplo,NonConj TMV_FIRSTLAST1(m,m+size*size));
    }

    template <typename T>
    inline SymMatrixView<T> HermMatrixViewOf(
        T* m, ptrdiff_t size, UpLoType uplo, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(size>=0);
        const ptrdiff_t stepi = stor == RowMajor ? size : 1;
        const ptrdiff_t stepj = stor == RowMajor ? 1 : size;
        return SymMatrixView<T>(
            m,size,stepi,stepj,Herm,uplo,NonConj TMV_FIRSTLAST1(m,m+size*size));
    }

    template <typename T>
    inline ConstSymMatrixView<T> SymMatrixViewOf(
        const T* m, ptrdiff_t size, UpLoType uplo, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(size>=0);
        return ConstSymMatrixView<T>(
            m,size,stepi,stepj,Sym,uplo,NonConj);
    }

    template <typename T>
    inline ConstSymMatrixView<T> HermMatrixViewOf(
        const T* m, ptrdiff_t size, UpLoType uplo, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(size>=0);
        return ConstSymMatrixView<T>(
            m,size,stepi,stepj,Herm,uplo,NonConj);
    }

    template <typename T>
    inline SymMatrixView<T> SymMatrixViewOf(
        T* m, ptrdiff_t size, UpLoType uplo, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(size>=0);
        return SymMatrixView<T>(
            m,size,stepi,stepj,Sym,uplo,NonConj TMV_FIRSTLAST1(m,m+size*size));
    }

    template <typename T>
    inline SymMatrixView<T> HermMatrixViewOf(
        T* m, ptrdiff_t size, UpLoType uplo, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(size>=0);
        return SymMatrixView<T>(
            m,size,stepi,stepj,Herm,uplo,NonConj TMV_FIRSTLAST1(m,m+size*size));
    }

    //
    // Swap Matrices
    //

    template <typename T>
    inline void Swap(SymMatrixView<T> m1, SymMatrixView<T> m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.issym() == m2.issym());
        TMVAssert(m1.isherm() == m2.isherm());
        Swap(m1.upperTri(),m2.upperTri());
    }

    template <typename T, int A>
    inline void Swap(SymMatrixView<T> m1, SymMatrix<T,A>& m2)
    { Swap(m1,m2.view()); }

    template <typename T, int A>
    inline void Swap(SymMatrix<T,A>& m1, SymMatrixView<T> m2)
    { Swap(m1.view(),m2); }

    template <typename T, int A1, int A2>
    inline void Swap(SymMatrix<T,A1>& m1, SymMatrix<T,A2>& m2)
    { Swap(m1.view(),m2.view()); }

    template <typename T, int A>
    inline void Swap(SymMatrixView<T> m1, HermMatrix<T,A>& m2)
    { Swap(m1,m2.view()); }

    template <typename T, int A>
    inline void Swap(HermMatrix<T,A>& m1, SymMatrixView<T> m2)
    { Swap(m1.view(),m2); }

    template <typename T, int A1, int A2>
    inline void Swap(HermMatrix<T,A1>& m1, HermMatrix<T,A2>& m2)
    { Swap(m1.view(),m2.view()); }



    //
    // Views:
    //

    template <typename T>
    inline ConstSymMatrixView<T> Transpose(const GenSymMatrix<T>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A> Transpose(const ConstSymMatrixView<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> Transpose(
        const SymMatrix<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> Transpose(
        const HermMatrix<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline SymMatrixView<T,A> Transpose(SymMatrixView<T,A> m)
    { return m.transpose(); }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> Transpose(SymMatrix<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> Transpose(HermMatrix<T,A>& m)
    { return m.transpose(); }

    template <typename T>
    inline ConstSymMatrixView<T> Conjugate(const GenSymMatrix<T>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A> Conjugate(const ConstSymMatrixView<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> Conjugate(
        const SymMatrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> Conjugate(
        const HermMatrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline SymMatrixView<T,A> Conjugate(SymMatrixView<T,A> m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> Conjugate(SymMatrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> Conjugate(HermMatrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T>
    inline ConstSymMatrixView<T> Adjoint(const GenSymMatrix<T>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A> Adjoint(const ConstSymMatrixView<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> Adjoint(
        const SymMatrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstSymMatrixView<T,A&FortranStyle> Adjoint(
        const HermMatrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline SymMatrixView<T,A> Adjoint(SymMatrixView<T,A> m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> Adjoint(SymMatrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline SymMatrixView<T,A&FortranStyle> Adjoint(HermMatrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T>
    inline QuotXS<T,T> Inverse(const GenSymMatrix<T>& m)
    { return m.inverse(); }

    //
    // SymMatrix ==, != SymMatrix
    //

    template <typename T1, typename T2>
    inline bool operator==(
        const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
    { return m1.upperTri() == m2.upperTri(); }

    template <typename T1, typename T2>
    inline bool operator!=(
        const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <typename T1, typename T2>
    inline bool operator==(
        const GenSymMatrix<T1>& m1, const GenMatrix<T2>& m2)
    {
        return
            m1.upperTri() == m2.upperTri() &&
            m1.lowerTri() == m2.lowerTri();
    }

    template <typename T1, typename T2>
    inline bool operator==(
        const GenMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
    { return m2 == m1; }

    template <typename T1, typename T2>
    inline bool operator!=(
        const GenSymMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <typename T1, typename T2>
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <typename T>
    inline std::istream& operator>>(std::istream& is, SymMatrixView<T> m)
    { return is >> IOStyle() >> m; }

    template <typename T, int A>
    inline std::istream& operator>>(std::istream& is, SymMatrix<T,A>& m)
    { return is >> IOStyle() >> m; }

    template <typename T, int A>
    inline std::istream& operator>>(std::istream& is, HermMatrix<T,A>& m)
    { return is >> IOStyle() >> m; }

    template <typename T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, SymMatrixView<T> m)
    { m.read(reader); return reader.getis(); }

    template <typename T, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, SymMatrix<T,A>& m)
    { m.read(reader); return reader.getis(); }

    template <typename T, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, HermMatrix<T,A>& m)
    { m.read(reader); return reader.getis(); }

} // namespace tmv

#endif
