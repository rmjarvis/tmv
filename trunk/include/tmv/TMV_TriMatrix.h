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
// This file defines the TMV TriMatrix class.
//
// Constructors:
//
//    There are two TriMatrix classes: UpperTriMatrix<T> and LowerTriMatrix<T>
//    For these notes, I will just write TriMatrix, but for all uses,
//    you need to write "Upper" or "Lower" before the "Tri".
//
//    In addition to the type template parameter (T), TriMatrixes have three
//    additional template parameters:
//        DiagType dt = UnitDiag || NonUnitDiag 
//        StorageType stor = RowMajor || ColMajor
//        IndexType I = CStyle || FortranStyle
//
//        They each have default values, so you can omit all three,
//        stor and I, or just stor.  
//        The default values are: {NonUnitDiag, RowMajor, CStyle}
//
//        If dt is UnitDiag, then the diagonal elements are not
//        actually stored or referenced.  The are all taken to be = 1.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    TriMatrix<T,dt,stor,I>(size_t n)
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    TriMatrix<T,dt,stor,I>(size_t n, T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    TriMatrix<T,dt,stor,I>(size_t n, T* vv)
//    TriMatrix<T,dt,stor,I>(size_t n, const std::vector<T>& vv)
//        Makes a Triangular Matrix with column size = row size = n
//        which copies the values from vv.
//
//    TriMatrix<T,dt,stor,I>(const Matrix<T>& m)
//    TriMatrix<T,dt,stor,I>(const TriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//
// Special Creators:
//
//    ConstUpperTriMatrixView UpperTriMatrixViewOf(
//            const T* m, size_t size, StorageType stor)
//    ConstUpperTriMatrixView UnitUpperTriMatrixViewOf(
//            const T* m, size_t size, StorageType stor)
//    UpperTriMatrixView UnitUpperTriMatrixViewOf(
//            T* m, size_t size, StorageType stor)
//    UpperTriMatrixView UnitUnitUpperTriMatrixViewOf(
//            T* m, size_t size, StorageType stor)
//        Returns a TriMatrixView of the elements in m, using the 
//        actual elements m for the storage.  The Unit versions return
//        views with dt = UnitDiag, the non-Unit versions return views
//        with dt = NonUnitDiag.
//        There are also corresponding LowerTriMatrix versions of these.
//        
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the TriMatrix
//
//    value_type operator()(int i, int j) const
//    value_type cref(int i, int j) const
//        Return the (i,j) element of the TriMatrix
//        The first one respects the index-style of the underlying matrix.
//        The second, cref, always uses CStyle indexing and does not 
//        do any checking of the valididty of i,j.
//
//    reference operator()(int i, int j)
//    reference ref(int i, int j)
//        Return a reference to the (i,j) element of the matrix
//        The first one respects the index-style of the underlying matrix.
//        The second, ref, always uses CStyle indexing and does not 
//        do any checking of the valididty of i,j.
//   
//    row_sub_type row(int i, int j1, int j2)
//    const_row_sub_type row(int i, int j1, int j2) const
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row.
//
//    col_sub_type col(int j, int i1, int i2)
//    const_col_sub_type col(int j, int i1, int i2) const
//        Return a portion of the jth column
//        This range must be a valid range for the requested column.
//
//    diag_type diag()
//    const_diag_type diag() const
//        Return the main diagonal
//        The TriMatrix must be NonUnitDiag.
//
//    diag_sub_type diag(int i)
//    diag_sub_type diag(int i, int j1, int j2)
//    const_diag_sub_type diag(int i) const
//    const_diag_sub_type diag(int i, int j1, int j2) const
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal SubVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//        i>0 will give an error for a LowerTriMatrix
//        i<0 will give an error for an UpperTriMatrix
//        i=0 will give an error for a UnitDiag TriMatrix
//
// Functions of Matrices:
//
//    Most of these are the same as for a regular matrix, so I 
//    only give the full description for the new functionality.
//
//    value_type m.det() const    or Det(m)
//    real_type m.logDet(value_type* sign=0) const   or LogDet(m,sign)
//    value_type m.trace() const    or Trace(m)
//    real_type m.norm() const    or Norm(m)
//    real_type m.normF() const    or NormF(m)
//    real_type m.normSq() const    or NormSq()
//    real_type m.normSq(real_type scale) const
//    real_type m.norm1() const    or Norm1(m)
//    real_type m.norm2() const    or Norm2(m)
//    real_type m.normInf() const    or NormInf(m)
//    value_type m.sumElements() const    or SumElements(m) 
//    real_type m.sumAbsElements() const    or SumAbsElements(m) 
//    value_type m.maxElement() const    or MaxElement(m) 
//    value_type m.minElement() const    or MinElement(m) 
//    real_type m.maxAbsElement() const    or MaxAbsElement(m) 
//    real_type m.minAbsElement() const    or MinAbsElement(m) 
//
//    void m.makeInverse(minv) const
//        This function allows minv to be either a regular Matrix
//        or a TriMatrix (of the same Upper or Lower as m).
//    void m.makeInverseATA(invata) const
//    inverse_type m.inverse() const    or Inverse(m)
//
//
// Modifying Functions
//
//    type& setZero()
//    type& setAllTo(value_type x)
//    type& addToAll(value_type x)
//    type& clip(real_type thresh)
//    type& applyToAll(const F& f)
//    type& conjugateSelf()
//    type& setToIdentity(value_type x = 1)
//    Swap(TriMatrix& m1, TriMatrix& m2)
//
//    type& invertSelf()
//        Change the TriMatrix into its own inverse.  This can be done
//        efficiently in place without requiring extra storage, so 
//        this function provides that functionality.
//
//
// Views of a TriMatrix:
//
//    (As usual, all of these have a const_version as well.)
//
//    submatrix_type subMatrix(int i1, int i2, int j1, int j2,
//            int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the TriMatrix.
//
//    subvector_type subVector(int i, int j, int istep, int jstep, int size)
//        Returns a SubVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//
//    subtrimatrix_type subTriMatrix(int i1, int i2, int istep)
//        Returns the TriMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with an UpperTriMatrix of size 10, the x's below
//        are the original data, the O's are the SubTriMatrix returned
//        with the command subTriMatrix(3,11,2), and the #'s are the 
//        SubTriMatrix returned with subTriMatrix(0,3)
//
//        ###xxxxxxx
//         ##xxxxxxx
//          #xxxxxxx
//           OxOxOxO
//            xxxxxx
//             OxOxO
//              xxxx
//               OxO
//                xx
//                 O
//
//    offdiag_type offDiag()
//        Returns the (NonUnitDiag) TriMatrix of all the off-diagonal
//        elements of a TriMatrix.
//
//    unitdiag_type viewAsUnitDiag()
//        Re-view a NonUnitDiag (or UknownDiag) TriMatrix with a view that
//        takes the diagonal elements to all be equal to 1.
//
//    nonunitdiag_type viewAsNonUnitDiag()
//        Re-view an UnknownUnitDiag TriMatrix as NonUnitDiag.
//
//    unknowndiag_type viewAsUnknownDiag()
//        Re-view a UnitDiag or NonUnitDiag TriMatrix as UnknownDiag.
//
//    realpart_type realPart()
//    imagpart_type imagPart()
//        For a complex TriMatrix, returns the real or imaginary part
//        as a real TriMatrix.
//
//    view_type view()
//    conjugate_type conjugate()
//    transpose_type transpose()
//    adjoint_type adjoint()
//        Note that the Transpose or Adjoint of an UpperTriMatrix returns 
//        a view which is a LowerTriMatrix, and vice versa.
//
//    nonconj_type nonConj()
//        Returns a mutable view of the (const) TriMatrix
//
//    nonconst_type nonConst()
//        Returns a mutable view of a const Matrix.
//
//    const_view_type constView()
//        Returns a const view of a mutable Matrix.
//
//
// I/O: 
//
//    os << m 
//        Writes m to ostream os in the usual Matrix format
//
//    m.writeCompact(os)
//        Writes m to ostream os in the following compact format:
//        For an UpperTriMatrix:
//          size 
//          ( m(0,0) m(0,1) ... m(0,size) )
//          ( m(1,1) .. m(1,size) )
//          ...
//          ( m(size,size) )
//
//        For a LowerTriMatrix:
//          size 
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(size,0) ... m(size,size) )
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the TriMatrix to be read, you can
//        use this form where mptr is an auto_ptr to an undefined TriMatrix.
//        (Note: if the DiagType for the TriMatrix is UnitDiag, then
//        all of the diagonals read in must be = 1.)
//
//
// Division Control Functions:
//
//    Most of the point of using TriMatrixes is that they are easy
//    to divide using either forward substitution or back substitution.
//    Therefore, the only division available for TriMatrixes is 
//    this variety.  To do something else (like SVD), you need to 
//    copy it to a regular matrix.
//


#ifndef TMV_TriMatrix_H
#define TMV_TriMatrix_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Matrix.h"
#include <vector>

namespace tmv {

    // A helper struct to provide a valid reference when the DiagType
    // is UnitDiag (or UnknownDiag)
    // It is very similar to ConjRef, but also has a boolean isunit
    // which indicates whether the reference is on the diagonal of a 
    // UnitDiag TriMatrix.  If so, it can be an rhs with the value 1, 
    // but not a lhs.
    template <class T, bool C>
    class TriRef 
    {
    public:
        typedef typename AuxRef<T,C>::reference reference;
        explicit inline TriRef(bool _u, T& _v) : isunit(_u), itsref(_v) {}
        inline TriRef(const TriRef<T,C>& rhs) : 
            isunit(rhs.isunit), itsref(rhs.itsref) {}
        inline ~TriRef() {}

        inline operator T() const { return val(); }
        inline reference getRef() { return ref(); }
        inline T operator-() const { return -val(); }

        template <class T2>
        inline TriRef<T,C>& operator=(const TriRef<T2,C>& rhs)
        { assign(rhs.val()); return *this; }
        template <class T2>
        inline TriRef<T,C>& operator=(T2 rhs)
        { assign(rhs); return *this; }

        template <class T2>
        inline TriRef<T,C>& operator+=(const TriRef<T2,C>& x2)
        { assign(val() + x2.val()); return *this; }
        template <class T2>
        inline TriRef<T,C>& operator+=(T2 x2)
        { assign(val() + x2); return *this; }
        template <class T2>
        inline typename Traits2<T,T2>::type operator+(const TriRef<T2,C>& x2)
        { return val() + x2.val(); }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator+(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()+x2; }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator+(
            T2 x1, const TriRef<T,C>& x2)
        { return x1+x2.val(); }
        template <class T2>
        inline friend T2& operator+=(T2& x1, const TriRef<T,C>& x2)
        { return x1 += x2.val(); }

        template <class T2>
        inline TriRef<T,C>& operator-=(const TriRef<T2,C>& x2) 
        { assign(val() - x2.val()); return *this; }
        template <class T2>
        inline TriRef<T,C>& operator-=(T2 x2)
        { assign(val() - x2); return *this; }
        template <class T2>
        inline typename Traits2<T,T2>::type operator-(const TriRef<T2,C>& x2)
        { return val()-x2.val(); }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator-(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()-x2; }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator-(
            T2 x1, const TriRef<T,C>& x2)
        { return x1-x2.val(); }
        template <class T2>
        inline friend T2& operator-=(T2& x1, const TriRef<T,C>& x2)
        { return x1 -= x2.val(); }

        template <class T2>
        inline TriRef<T,C>& operator*=(const TriRef<T2,C>& x2) 
        { assign(x2.val() * val()); return *this; }
        template <class T2>
        inline TriRef<T,C>& operator*=(T2 x2) 
        { assign(x2 * val()); return *this; }
        template <class T2>
        inline typename Traits2<T,T2>::type operator*(const TriRef<T2,C> x2)
        { return val()*x2.val(); }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator*(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()*x2; }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator*(
            T2 x1, const TriRef<T,C>& x2)
        { return x1*x2.val(); }
        template <class T2>
        inline friend T2& operator*=(T2& x1, const TriRef<T,C>& x2)
        {
            if (x2.isunit) return x1;
            else return x1 *= x2.val(); 
        }

        template <class T2>
        inline TriRef<T,C>& operator/=(const TriRef<T2,C>& x2) 
        { assign(val() / x2.val()); return *this; }
        template <class T2>
        inline TriRef<T,C>& operator/=(T2 x2) 
        { assign(val() / x2); return *this; }
        template <class T2>
        inline typename Traits2<T,T2>::type operator/(const TriRef<T2,C>& x2)
        { return val()/x2.val(); }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator/(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()/x2; }
        template <class T2>
        inline friend typename Traits2<T,T2>::type operator/(
            T2 x1, const TriRef<T,C>& x2)
        { return x1/x2.val(); }
        template <class T2>
        inline friend T2& operator/=(T2& x1, const TriRef<T,C>& x2)
        {
            if (x2.isunit) return x1;
            else return x1 /= x2.val(); 
        }

        template <class T2>
        inline bool operator==(const TriRef<T2,C>& x2) const
        { return val() == x2.val(); }
        template <class T2>
        inline bool operator==(T2 x2) const 
        { return val() == x2; }
        template <class T2>
        inline friend bool operator==(T2 x1, const TriRef<T,C>& x2)
        { return x1 == val(); }
        template <class T2>
        inline bool operator!=(const TriRef<T2,C>& x2) const
        { return !(operator==(x2)); }
        template <class T2>
        inline bool operator!=(T2 x2) const 
        { return !(operator==(x2)); }
        template <class T2>
        inline friend bool operator!=(T2 x1, const TriRef<T,C>& x2)
        { return !(x2==x1); }

        inline friend std::ostream& operator<<(std::ostream& os, TriRef<T,C> x)
        { os << x.val(); return os; }
        inline friend std::istream& operator>>(std::istream& is, TriRef<T,C> x)
        { is >> x.ref(); return is; }

    private:

        template <class T2>
        void check(T2 x) const
        {
            TMVAssert(
                (!isunit || x == T2(1)) && 
                "Trying to assign to the diagonal of a UnitDiag TriMatrix.");
        }
        void check() const
        {
            TMVAssert(
                !isunit && 
                "Trying to assign to the diagonal of a UnitDiag TriMatrix.");
        }
        reference ref() { check(); return itsref; }
        template <class T2>
        void assign(T2 x) { check(x); if (!isunit) itsref = x; }
        T val() const { return isunit ? T(1) : T(itsref); }

        const bool isunit;
        reference itsref;
    };

    // Overload some functions to work with TriRef
    template <class T, bool C> 
    inline T TMV_CONJ(const TriRef<T,C>& x) 
    { return TMV_CONJ(T(x)); }
    template <class T, bool C> 
    inline typename Traits<T>::real_type TMV_NORM(const TriRef<T,C>& x) 
    { return TMV_NORM(T(x)); }
    template <class T, bool C> 
    inline typename Traits<T>::real_type TMV_ABS(const TriRef<T,C>& x) 
    { return TMV_ABS(T(x)); }
    template <class T, bool C> 
    inline T TMV_SQR(const TriRef<T,C>& x) 
    { return TMV_SQR(T(x)); }
    template <class T, bool C> 
    inline T TMV_SQRT(const TriRef<T,C>& x) 
    { return TMV_SQRT(T(x)); }
    template <class T, bool C> 
    inline typename Traits<T>::real_type TMV_REAL(const TriRef<T,C>& x) 
    { return TMV_REAL(T(x)); }
    template <class T, bool C> 
    inline typename Traits<T>::real_type TMV_IMAG(const TriRef<T,C>& x) 
    { return TMV_IMAG(T(x)); }

    template <class T, bool C1, bool C2> 
    inline void TMV_SWAP(TriRef<T,C1> x1, TriRef<T,C2> x2)
    { return TMV_SWAP(x1.GetRef(),x2.GetRef()); }
    template <class T, bool C2> 
    inline void TMV_SWAP(T& x1, TriRef<T,C2> x2)
    { return TMV_SWAP(x1,x2.GetRef()); }
    template <class T, bool C1> 
    inline void TMV_SWAP(TriRef<T,C1> x1, T& x2)
    { return TMV_SWAP(x1.GetRef(),x2); }


    template <class T, DiagType D, StorageType S, IndexStyle I> 
    struct Traits<UpperTriMatrix<T,D,S,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef UpperTriMatrix<T,D,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==ColMajor ? 1 : UNKNOWN) };
        enum { _stepj = (S==RowMajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = false };
        enum { _unit = (D == UnitDiag) };
        enum { _unknowndiag = false };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _hasdivider = false };

        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { notC = iscomplex };

        typedef ConstVectorView<T,_stepi,false,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,false,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_view_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,CStyle> 
            const_cview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,FortranStyle> 
            const_fview_type;
        typedef ConstUpperTriMatrixView<T,D> const_xview_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag> const_xdview_type;
        typedef ConstUpperTriMatrixView<T,D,1,_stepj,false,I> 
            const_cmview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,1,false,I> 
            const_rmview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstLowerTriMatrixView<T,D,_stepj,_stepi,false,I> 
            const_transpose_type;
        typedef ConstLowerTriMatrixView<T,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_offdiag_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            const_unitdiag_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_nonunitdiag_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            const_unknowndiag_type;
        typedef ConstUpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,I> nonconst_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        typedef TriRef<T,false> reference;

        typedef VectorView<T,_stepi,false,I> col_sub_type;
        typedef VectorView<T,_stepj,false,I> row_sub_type;
        typedef VectorView<T,_diagstep,false,I> diag_type;
        typedef VectorView<T,_diagstep,false,I> diag_sub_type;

        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            subtrimatrix_type;
        typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> 
            subtrimatrix_step_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_type;

        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,I> view_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,CStyle> cview_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,FortranStyle> 
            fview_type;
        typedef UpperTriMatrixView<T,D> xview_type;
        typedef UpperTriMatrixView<T,UnknownDiag> xdview_type;
        typedef UpperTriMatrixView<T,D,1,_stepj,false,I> cmview_type;
        typedef UpperTriMatrixView<T,D,_stepi,1,false,I> rmview_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,notC,I> conjugate_type;
        typedef LowerTriMatrixView<T,D,_stepj,_stepi,false,I> transpose_type;
        typedef LowerTriMatrixView<T,D,_stepj,_stepi,notC,I> adjoint_type;

        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            offdiag_type;
        typedef UpperTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            unitdiag_type;
        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            nonunitdiag_type;
        typedef UpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            unknowndiag_type;
        typedef UpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,I> nonconj_type;
    };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

    // A Helper class to make a TriMatrix from a BaseMatrix
    // If the BaseMatrix is assignable to the Tri then we do the 
    // assign.  But if not, then we allow the construction - we
    // just make the TriMatrix from the correct portion of the matrix.
    template <bool assignable_to_tri>
    struct TriCopy // assignable
    {
        template <class M1, class M2>
        static void copy(
            const BaseMatrix<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
        { m1.newAssignTo(m2); }
    };
    template <>
    struct TriCopy<false>
    {
        template <class M1, class M2>
        static void copy(
            const BaseMatrix<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
        {
            const bool up = M2::_upper;
            if (m2.isunit()) {
                typename M2::offdiag_type m2o = m2.offDiag();
                Maybe<up>::uppertri(m1.calc()).offDiag().newAssignTo(m2o);
            } else {
                typename M2::nonunitdiag_type m2nu = m2.viewAsNonUnitDiag();
                Maybe<up>::uppertri(m1.calc()).newAssignTo(m2nu);
            }
        }
    };

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    class UpperTriMatrix : 
        public BaseMatrix_Tri_Mutable<UpperTriMatrix<T,D,S,I> >
    {
    public:
        typedef UpperTriMatrix<T,D,S,I> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename Traits<type>::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _unit = Traits<type>::_unit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _conj = Traits<type>::_conj };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };

        //
        // Constructors
        //

        explicit inline UpperTriMatrix(size_t n=0) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
        }

        inline UpperTriMatrix(size_t n, T x) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            Maybe<_unit>::offdiag(*this).setAllTo(x);
        }

        inline UpperTriMatrix(size_t n, const T* vv) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            VectorView<T,1> lv(ptr(),n*n);
            ConstVectorView<T,1>(vv,n*n).newAssignTo(lv);
        }

        inline UpperTriMatrix(size_t n, const std::vector<T>& vv) :
            itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TMVAssert(vv.size() == n*n);
            VectorView<T,1> lv(ptr(),n*n);
            ConstVectorView<T,1>(&vv[0],n*n).newAssignTo(lv);
        }

        inline UpperTriMatrix(const type& m2) :
            itss(m2.size()), itsm(itss*itss)
        { 
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            m2.newAssignTo(*this);
        }

        template <class M2> 
        inline UpperTriMatrix(const BaseMatrix<M2>& m2) :
            itss(m2.rowsize()), itsm(itss*itss)
        { 
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert((
                (M2::_calc && ShapeTraits<M2::_shape>::upper) || assignable));
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());

            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        inline UpperTriMatrix(const BaseMatrix_Tri<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(M2::_upper);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            Maybe<_unit && !M2::_unit>::unitview(m2).newAssignTo(*this);
        }

        template <class M2>
        inline UpperTriMatrix(const BaseMatrix_Diag<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            typename type::diag_type d = this->diag();
            this->setZero();
            m2.calc().diag().newAssignTo(d);
        }

        inline ~UpperTriMatrix()
        {
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(999));
#endif
        }


        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { 
            if (&m2 != this) base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        inline type& operator=(T x)
        {
            base_mut::operator=(x);
            return *this;
        }


        // 
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]);
        }

        inline reference ref(int i, int j)
        {
            return reference(
                isunit() && i==j,
                itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()] ); 
        }

        inline void swapWith(type& m2)
        {
            TMVAssert(m2.size() == size());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }

        inline void resize(const size_t s)
        {
            itss = s;
            itsm.resize(s*s);
        }

        inline size_t size() const { return itss; }
        inline int stepi() const { return S==RowMajor ? itss : 1; }
        inline int stepj() const { return S==RowMajor ? 1 : itss; }
        inline DiagType dt() const { return D; }
        inline bool isunit() const { return D == UnitDiag; }
        inline bool isconj() const { return false; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline StorageType stor() const { return S; }

    protected :

        size_t itss;
        AlignedArray<T> itsm;

    }; // UpperTriMatrix

    template <class T, DiagType D, StorageType S>
    class UpperTriMatrixF : public UpperTriMatrix<T,D,S,FortranStyle>
    {
    public:

        typedef UpperTriMatrixF<T,D,S> type;
        typedef UpperTriMatrix<T,D,S,FortranStyle> mtype;

        explicit inline UpperTriMatrixF(size_t s) : mtype(s) {}
        inline UpperTriMatrixF(size_t s, T x) : mtype(s,x) {}
        inline UpperTriMatrixF(size_t s, const T* vv) : mtype(s,vv) {}
        inline UpperTriMatrixF(size_t s, const std::vector<T>& vv) :
            mtype(s,vv) {}
        template <class M2> 
        inline UpperTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        inline UpperTriMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        inline UpperTriMatrixF(const BaseMatrix_Tri<M2>& m2) : mtype(m2) {}
        template <class M2>
        inline UpperTriMatrixF(const BaseMatrix_Rec<M2>& m2) : mtype(m2) {}
        inline ~UpperTriMatrixF() {}

        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(T x)
        { mtype::operator=(x); return *this; }

    }; // UpperTriMatrixF

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstUpperTriMatrixView<T,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstUpperTriMatrixView<T,D,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        enum { DD = (D == UnknownDiag ? NonUnitDiag : D) };
        typedef UpperTriMatrix<T,DiagType(DD),Sj==1?RowMajor:ColMajor,I> 
            copy_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _unit = (D == UnitDiag) };
        enum { _unknowndiag = (D == UnknownDiag) };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _hasdivider = false };

        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };

        typedef ConstVectorView<T,_stepi,C,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,I> 
            const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,I> const_view_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,FortranStyle> 
            const_fview_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstUpperTriMatrixView<T,D,1,_stepj,C,I> const_cmview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,1,C,I> const_rmview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstLowerTriMatrixView<T,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstLowerTriMatrixView<T,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_offdiag_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            const_unitdiag_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_nonunitdiag_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknowndiag_type;
        typedef ConstUpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,C,I> nonconst_type;
    };

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class ConstUpperTriMatrixView :
        public BaseMatrix_Tri<ConstUpperTriMatrixView<T,D,Si,Sj,C,I> >
    {
    public:
        typedef ConstUpperTriMatrixView<T,D,Si,Sj,C,I> type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _unit = Traits<type>::_unit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _conj = Traits<type>::_conj };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };

        //
        // Constructors
        //
        inline ConstUpperTriMatrixView(
            const T* m, size_t s, bool u, int si, int sj) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline ConstUpperTriMatrixView(const T* m, size_t s, bool u, int si) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstUpperTriMatrixView(const T* m, size_t s, bool u) :
            itsm(m), itss(s), itsu(u), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstUpperTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixView(
            const ConstUpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixView(
            const UpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixView(
            const ConstSmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixView(
            const SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~ConstUpperTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }

    private :
        inline void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }

        inline T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                DoConj<C>(itsm[i*stepi() + j*stepj()]));
        }

        inline size_t colsize() const { return itss; }
        inline size_t rowsize() const { return itss; }
        inline size_t size() const { return itss; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline bool isconj() const { return C; }
        inline bool isunit() const { return itsu; }
        inline bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        inline bool iscm() const
        { 
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

#ifdef TMV_DEBUG
        const T* itsm;
#else
        const T*const itsm;
#endif
        const size_t itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstUpperTriMatrixView

    template <class T, DiagType D, int Si, int Sj, bool C>
    class ConstUpperTriMatrixViewF :
        public ConstUpperTriMatrixView<T,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef ConstUpperTriMatrixViewF<T,D,Si,Sj,C> type;
        typedef ConstUpperTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

        inline ConstUpperTriMatrixViewF(const T* m, size_t s, int si, int sj) :
            mtype(m,s,si,sj) {}
        inline ConstUpperTriMatrixViewF(const T* m, size_t s, int si) :
            mtype(m,s,si) {}
        inline ConstUpperTriMatrixViewF(const T* m, size_t s) :
            mtype(m,s) {}
        inline ConstUpperTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixViewF(
            const ConstUpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixViewF(
            const UpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixViewF(
            const ConstSmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstUpperTriMatrixViewF(
            const SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        inline ~ConstUpperTriMatrixViewF() {}

    private :
        inline void operator=(const type& m2);

    }; // ConstUpperTriMatrixViewF

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<UpperTriMatrixView<T,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef UpperTriMatrixView<T,D,Si,Sj,C,I> type;
        typedef const ConstUpperTriMatrixView<T,D,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;
        enum { DD = (D == UnknownDiag ? NonUnitDiag : D) };
        typedef UpperTriMatrix<T,DiagType(DD),Sj==1?RowMajor:ColMajor,I> 
            copy_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _unit = (D == UnitDiag) };
        enum { _unknowndiag = (D == UnknownDiag) };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _hasdivider = false };

        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };

        typedef ConstVectorView<T,_stepi,C,I> const_col_type;
        typedef ConstVectorView<T,_stepi,C,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,I> 
            const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,I> const_view_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,C,FortranStyle> 
            const_fview_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstUpperTriMatrixView<T,D,1,_stepj,C,I> const_cmview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,1,C,I> const_rmview_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstLowerTriMatrixView<T,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstLowerTriMatrixView<T,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,1,C,I> 
            const_offdiag_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,1,C,I> 
            const_unitdiag_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,1,C,I> 
            const_nonunitdiag_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,1,C,I> 
            const_unknowndiag_type;
        typedef ConstUpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,C,I> nonconst_type;

        typedef TriRef<T,C> reference;

        typedef VectorView<T,_stepi,C,I> col_sub_type;
        typedef VectorView<T,_stepj,C,I> row_sub_type;
        typedef VectorView<T,_diagstep,C,I> diag_type;
        typedef VectorView<T,_diagstep,C,I> diag_sub_type;

        typedef UpperTriMatrixView<T,D,_stepi,_stepj,C,I> subtrimatrix_type;
        typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> 
            subtrimatrix_step_type;
        typedef MatrixView<T,_stepi,_stepj,C,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,C,I> subvector_type;

        typedef UpperTriMatrixView<T,D,_stepi,_stepj,C,I> view_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,C,CStyle> cview_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,C,FortranStyle> 
            fview_type;
        typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
        typedef UpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            xdview_type;
        typedef UpperTriMatrixView<T,D,1,_stepj,C,I> cmview_type;
        typedef UpperTriMatrixView<T,D,_stepi,1,C,I> rmview_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,notC,I> conjugate_type;
        typedef LowerTriMatrixView<T,D,_stepj,_stepi,C,I> transpose_type;
        typedef LowerTriMatrixView<T,D,_stepj,_stepi,notC,I> adjoint_type;

        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            offdiag_type;
        typedef UpperTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            unitdiag_type;
        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            nonunitdiag_type;
        typedef UpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            unknowndiag_type;
        typedef UpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef UpperTriMatrixView<T,D,_stepi,_stepj,false,I> nonconj_type;
    };

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class UpperTriMatrixView :
        public BaseMatrix_Tri_Mutable<UpperTriMatrixView<T,D,Si,Sj,C,I> >
    {
    public:
        typedef UpperTriMatrixView<T,D,Si,Sj,C,I> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename Traits<type>::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _unit = Traits<type>::_unit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _conj = Traits<type>::_conj };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };

        //
        // Constructors
        //

        inline UpperTriMatrixView(T* m, size_t s, bool u, int si, int sj) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline UpperTriMatrixView(T* m, size_t s, bool u, int si) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline UpperTriMatrixView(T* m, size_t s, bool u) :
            itsm(m), itss(s), itsu(u), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        inline UpperTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), itsu(m2.isunit()), 
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline UpperTriMatrixView(UpperTriMatrixView<T,D2,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline UpperTriMatrixView(
            SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2
        ) :
            itsm(m2.ptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~UpperTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }


        //
        // Op = 
        //

        inline type& operator=(const type& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        inline type& operator=(const T x)
        {
            base_mut::operator=(x);
            return *this;
        }


        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                DoConj<C>(itsm[i*stepi() + j*stepj()]));
        }

        inline reference ref(int i, int j)
        { return reference(isunit() && i==j,itsm[i*stepi()+j*stepj()]); }

        inline size_t colsize() const { return itss; }
        inline size_t rowsize() const { return itss; }
        inline size_t size() const { return itss; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline bool isconj() const { return C; }
        inline bool isunit() const { return itsu; }
        inline bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        inline bool iscm() const
        { 
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

#ifdef TMV_DEBUG
        T* itsm;
#else
        T*const itsm;
#endif
        const size_t itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // UpperTriMatrixView

    template <class T, DiagType D, int Si, int Sj, bool C>
    class UpperTriMatrixViewF :
        public UpperTriMatrixView<T,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef UpperTriMatrixViewF<T,D,Si,Sj,C> type;
        typedef UpperTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

        inline UpperTriMatrixViewF(T* m, size_t s, bool u, int si, int sj) :
            mtype(m,s,si,sj) {}
        inline UpperTriMatrixViewF(T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline UpperTriMatrixViewF(T* m, size_t s, bool u) : mtype(m,s) {}
        inline UpperTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline UpperTriMatrixViewF(UpperTriMatrixView<T,D2,Si2,Sj2,C,I2> m2) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline UpperTriMatrixViewF(
            SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2) : mtype(m2) {}
        inline ~UpperTriMatrixViewF() {}

        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(const T x)
        { mtype::operator=(x); return *this; }

    }; // UpperTriMatrixViewF


    template <class T, DiagType D, StorageType S, IndexStyle I> 
    struct Traits<LowerTriMatrix<T,D,S,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef LowerTriMatrix<T,D,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==ColMajor ? 1 : UNKNOWN) };
        enum { _stepj = (S==RowMajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = false };
        enum { _unit = (D == UnitDiag) };
        enum { _unknowndiag = false };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _hasdivider = false };

        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { notC = iscomplex };

        typedef ConstVectorView<T,_stepi,false,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,false,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_view_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,CStyle> 
            const_cview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,FortranStyle> 
            const_fview_type;
        typedef ConstLowerTriMatrixView<T,D> const_xview_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag> const_xdview_type;
        typedef ConstLowerTriMatrixView<T,D,1,_stepj,false,I> 
            const_cmview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,1,false,I> 
            const_rmview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstUpperTriMatrixView<T,D,_stepj,_stepi,false,I> 
            const_transpose_type;
        typedef ConstUpperTriMatrixView<T,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_offdiag_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            const_unitdiag_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_nonunitdiag_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            const_unknowndiag_type;
        typedef ConstLowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,I> nonconst_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        typedef TriRef<T,false> reference;

        typedef VectorView<T,_stepi,false,I> col_sub_type;
        typedef VectorView<T,_stepj,false,I> row_sub_type;
        typedef VectorView<T,_diagstep,false,I> diag_type;
        typedef VectorView<T,_diagstep,false,I> diag_sub_type;

        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            subtrimatrix_type;
        typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> 
            subtrimatrix_step_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_type;

        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,I> view_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,CStyle> cview_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,FortranStyle> 
            fview_type;
        typedef LowerTriMatrixView<T,D> xview_type;
        typedef LowerTriMatrixView<T,UnknownDiag> xdview_type;
        typedef LowerTriMatrixView<T,D,1,_stepj,false,I> cmview_type;
        typedef LowerTriMatrixView<T,D,_stepi,1,false,I> rmview_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,notC,I> conjugate_type;
        typedef UpperTriMatrixView<T,D,_stepj,_stepi,false,I> transpose_type;
        typedef UpperTriMatrixView<T,D,_stepj,_stepi,notC,I> adjoint_type;

        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            offdiag_type;
        typedef LowerTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            unitdiag_type;
        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            nonunitdiag_type;
        typedef LowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            unknowndiag_type;
        typedef LowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,I> nonconj_type;
    };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    class LowerTriMatrix : 
        public BaseMatrix_Tri_Mutable<LowerTriMatrix<T,D,S,I> >
    {
    public:
        typedef LowerTriMatrix<T,D,S,I> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename Traits<type>::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _unit = Traits<type>::_unit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _conj = Traits<type>::_conj };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };

        //
        // Constructors
        //

        explicit inline LowerTriMatrix(size_t n=0) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
        }

        inline LowerTriMatrix(size_t n, T x) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            Maybe<_unit>::offdiag(*this).setAllTo(x);
        }

        inline LowerTriMatrix(size_t n, const T* vv) : itss(n), itsm(n*n)
        {
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            VectorView<T,1> lv(ptr(),n*n);
            ConstVectorView<T,1>(vv,n*n).newAssignTo(lv);
        }

        inline LowerTriMatrix(size_t n, const std::vector<T>& vv) :
            itss(n), itsm(n*n)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TMVAssert(vv.size() == n*n);
            VectorView<T,1> lv(ptr(),n*n);
            ConstVectorView<T,1>(&vv[0],n*n).newAssignTo(lv);
        }

        inline LowerTriMatrix(const type& m2) :
            itss(m2.size()), itsm(itss*itss)
        { 
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            m2.newAssignTo(*this);
        }

        template <class M2> 
        inline LowerTriMatrix(const BaseMatrix<M2>& m2) :
            itss(m2.colsize()), itsm(itss*itss)
        { 
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert((
                (M2::_calc && ShapeTraits<M2::_shape>::lower) || assignable));
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        inline LowerTriMatrix(const BaseMatrix_Tri<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(M2::_lower);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            Maybe<_unit && !M2::_unit>::unitview(m2).newAssignTo(*this);
        }

        template <class M2>
        inline LowerTriMatrix(const BaseMatrix_Diag<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            typename type::diag_type d = this->diag();
            this->setZero();
            m2.calc().diag().newAssignTo(d);
        }

        inline ~LowerTriMatrix()
        {
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(999));
#endif
        }


        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { 
            if (&m2 != this) base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        inline type& operator=(T x)
        {
            base_mut::operator=(x);
            return *this;
        }


        // 
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]);
        }

        inline reference ref(int i, int j)
        {
            return reference(
                isunit() && i==j,
                itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()] ); 
        }

        inline void swapWith(type& m2)
        {
            TMVAssert(m2.size() == size());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }

        inline void resize(const size_t s)
        {
            itss = s;
            itsm.resize(s*s);
        }

        inline size_t size() const { return itss; }
        inline int stepi() const { return S==RowMajor ? itss : 1; }
        inline int stepj() const { return S==RowMajor ? 1 : itss; }
        inline DiagType dt() const { return D; }
        inline bool isunit() const { return D == UnitDiag; }
        inline bool isconj() const { return false; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline StorageType stor() const { return S; }

    protected :

        size_t itss;
        AlignedArray<T> itsm;

    }; // LowerTriMatrix

    template <class T, DiagType D, StorageType S>
    class LowerTriMatrixF : public LowerTriMatrix<T,D,S,FortranStyle>
    {
    public:

        typedef LowerTriMatrixF<T,D,S> type;
        typedef LowerTriMatrix<T,D,S,FortranStyle> mtype;

        explicit inline LowerTriMatrixF(size_t s) : mtype(s) {}
        inline LowerTriMatrixF(size_t s, T x) : mtype(s,x) {}
        inline LowerTriMatrixF(size_t s, const T* vv) : mtype(s,vv) {}
        inline LowerTriMatrixF(size_t s, const std::vector<T>& vv) :
            mtype(s,vv) {}
        template <class M2> 
        inline LowerTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        inline LowerTriMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        inline LowerTriMatrixF(const BaseMatrix_Tri<M2>& m2) : mtype(m2) {}
        template <class M2>
        inline LowerTriMatrixF(const BaseMatrix_Rec<M2>& m2) : mtype(m2) {}
        inline ~LowerTriMatrixF() {}

        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(T x)
        { mtype::operator=(x); return *this; }

    }; // LowerTriMatrixF

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstLowerTriMatrixView<T,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstLowerTriMatrixView<T,D,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        enum { DD = (D == UnknownDiag ? NonUnitDiag : D) };
        typedef LowerTriMatrix<T,DiagType(DD),Sj==1?RowMajor:ColMajor,I> 
            copy_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _unit = (D == UnitDiag) };
        enum { _unknowndiag = (D == UnknownDiag) };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _hasdivider = false };

        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };

        typedef ConstVectorView<T,_stepi,C,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,I> 
            const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,I> const_view_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,FortranStyle> 
            const_fview_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstLowerTriMatrixView<T,D,1,_stepj,C,I> const_cmview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,1,C,I> const_rmview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstUpperTriMatrixView<T,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstUpperTriMatrixView<T,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_offdiag_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            const_unitdiag_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_nonunitdiag_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknowndiag_type;
        typedef ConstLowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,C,I> nonconst_type;
    };

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class ConstLowerTriMatrixView :
        public BaseMatrix_Tri<ConstLowerTriMatrixView<T,D,Si,Sj,C,I> >
    {
    public:
        typedef ConstLowerTriMatrixView<T,D,Si,Sj,C,I> type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _unit = Traits<type>::_unit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _conj = Traits<type>::_conj };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };

        //
        // Constructors
        //
        inline ConstLowerTriMatrixView(
            const T* m, size_t s, bool u, int si, int sj) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline ConstLowerTriMatrixView(const T* m, size_t s, bool u, int si) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstLowerTriMatrixView(const T* m, size_t s, bool u) :
            itsm(m), itss(s), itsu(u), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstLowerTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixView(
            const ConstLowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixView(
            const LowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixView(
            const ConstSmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixView(
            const SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~ConstLowerTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }

    private :
        inline void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }

        inline T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                DoConj<C>(itsm[i*stepi() + j*stepj()]));
        }

        inline size_t colsize() const { return itss; }
        inline size_t rowsize() const { return itss; }
        inline size_t size() const { return itss; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline bool isconj() const { return C; }
        inline bool isunit() const { return itsu; }
        inline bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        inline bool iscm() const
        { 
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

#ifdef TMV_DEBUG
        const T* itsm;
#else
        const T*const itsm;
#endif
        const size_t itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstLowerTriMatrixView

    template <class T, DiagType D, int Si, int Sj, bool C>
    class ConstLowerTriMatrixViewF :
        public ConstLowerTriMatrixView<T,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef ConstLowerTriMatrixViewF<T,D,Si,Sj,C> type;
        typedef ConstLowerTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

        inline ConstLowerTriMatrixViewF(
            const T* m, size_t s, bool u, int si, int sj) :
            mtype(m,s,si,sj) {}
        inline ConstLowerTriMatrixViewF(const T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline ConstLowerTriMatrixViewF(const T* m, size_t s, bool u) :
            mtype(m,s) {}
        inline ConstLowerTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixViewF(
            const ConstLowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixViewF(
            const LowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixViewF(
            const ConstSmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstLowerTriMatrixViewF(
            const SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        inline ~ConstLowerTriMatrixViewF() {}

    private :
        inline void operator=(const type& m2);

    }; // ConstLowerTriMatrixViewF

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<LowerTriMatrixView<T,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef LowerTriMatrixView<T,D,Si,Sj,C,I> type;
        typedef const ConstLowerTriMatrixView<T,D,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;
        enum { DD = (D == UnknownDiag ? NonUnitDiag : D) };
        typedef LowerTriMatrix<T,DiagType(DD),Sj==1?RowMajor:ColMajor,I> 
            copy_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _unit = (D == UnitDiag) };
        enum { _unknowndiag = (D == UnknownDiag) };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _hasdivider = false };

        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };

        typedef ConstVectorView<T,_stepi,C,I> const_col_type;
        typedef ConstVectorView<T,_stepi,C,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,I> 
            const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,I> const_view_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,C,FortranStyle> 
            const_fview_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstLowerTriMatrixView<T,D,1,_stepj,C,I> const_cmview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,1,C,I> const_rmview_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstUpperTriMatrixView<T,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstUpperTriMatrixView<T,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,1,C,I> 
            const_offdiag_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,1,C,I> 
            const_unitdiag_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,1,C,I> 
            const_nonunitdiag_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,1,C,I> 
            const_unknowndiag_type;
        typedef ConstLowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,C,I> nonconst_type;

        typedef TriRef<T,C> reference;

        typedef VectorView<T,_stepi,C,I> col_sub_type;
        typedef VectorView<T,_stepj,C,I> row_sub_type;
        typedef VectorView<T,_diagstep,C,I> diag_type;
        typedef VectorView<T,_diagstep,C,I> diag_sub_type;

        typedef LowerTriMatrixView<T,D,_stepi,_stepj,C,I> subtrimatrix_type;
        typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> 
            subtrimatrix_step_type;
        typedef MatrixView<T,_stepi,_stepj,C,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,C,I> subvector_type;

        typedef LowerTriMatrixView<T,D,_stepi,_stepj,C,I> view_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,C,CStyle> cview_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,C,FortranStyle> 
            fview_type;
        typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
        typedef LowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            xdview_type;
        typedef LowerTriMatrixView<T,D,1,_stepj,C,I> cmview_type;
        typedef LowerTriMatrixView<T,D,_stepi,1,C,I> rmview_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,notC,I> conjugate_type;
        typedef UpperTriMatrixView<T,D,_stepj,_stepi,C,I> transpose_type;
        typedef UpperTriMatrixView<T,D,_stepj,_stepi,notC,I> adjoint_type;

        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            offdiag_type;
        typedef LowerTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            unitdiag_type;
        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            nonunitdiag_type;
        typedef LowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            unknowndiag_type;
        typedef LowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef LowerTriMatrixView<T,D,_stepi,_stepj,false,I> nonconj_type;
    };

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class LowerTriMatrixView :
        public BaseMatrix_Tri_Mutable<LowerTriMatrixView<T,D,Si,Sj,C,I> >
    {
    public:
        typedef LowerTriMatrixView<T,D,Si,Sj,C,I> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _unit = Traits<type>::_unit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _conj = Traits<type>::_conj };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };

        //
        // Constructors
        //

        inline LowerTriMatrixView(T* m, size_t s, bool u, int si, int sj) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline LowerTriMatrixView(T* m, size_t s, bool u, int si) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline LowerTriMatrixView(T* m, size_t s, bool u) :
            itsm(m), itss(s), itsu(u), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        inline LowerTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), itsu(m2.isunit()), 
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline LowerTriMatrixView(LowerTriMatrixView<T,D2,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline LowerTriMatrixView(
            SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2
        ) :
            itsm(m2.ptr()), itss(m2.size()), itsu(m2.isunit()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~LowerTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }


        //
        // Op = 
        //

        inline type& operator=(const type& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        inline type& operator=(const T x)
        {
            base_mut::operator=(x);
            return *this;
        }


        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                DoConj<C>(itsm[i*stepi() + j*stepj()]));
        }

        inline reference ref(int i, int j)
        { return reference(isunit() && i==j,itsm[i*stepi()+j*stepj()]); }

        inline size_t colsize() const { return itss; }
        inline size_t rowsize() const { return itss; }
        inline size_t size() const { return itss; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline bool isconj() const { return C; }
        inline bool isunit() const { return itsu; }
        inline bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        inline bool iscm() const
        { 
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

#ifdef TMV_DEBUG
        T* itsm;
#else
        T*const itsm;
#endif
        const size_t itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // LowerTriMatrixView

    template <class T, DiagType D, int Si, int Sj, bool C>
    class LowerTriMatrixViewF :
        public LowerTriMatrixView<T,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef LowerTriMatrixViewF<T,D,Si,Sj,C> type;
        typedef LowerTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

        inline LowerTriMatrixViewF(T* m, size_t s, bool u, int si, int sj) :
            mtype(m,s,si,sj) {}
        inline LowerTriMatrixViewF(T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline LowerTriMatrixViewF(T* m, size_t s, bool u) : mtype(m,s) {}
        inline LowerTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline LowerTriMatrixViewF(LowerTriMatrixView<T,D2,Si2,Sj2,C,I2> m2) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline LowerTriMatrixViewF(
            SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2) : mtype(m2) {}
        inline ~LowerTriMatrixViewF() {}

        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Tri<M2>& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix_Diag<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(const T x)
        { mtype::operator=(x); return *this; }

    }; // LowerTriMatrixViewF



    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   UpperTriMatrixViewOf(T* m, n, S)
    //   UpperTriMatrixViewOf(T* m, n, si, sj)
    //   UnitUpperTriMatrixViewOf(T* m, n, S)
    //   UnitUpperTriMatrixViewOf(T* m, n, si, sj)
    //   LowerTriMatrixViewOf(T* m, n, S)
    //   LowerTriMatrixViewOf(T* m, n, si, sj)
    //   UnitLowerTriMatrixViewOf(T* m, n, S)
    //   UnitLowerTriMatrixViewOf(T* m, n, si, sj)
    //

    template <class T> 
    inline UpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return UpperTriMatrixView<T,NonUnitDiag>(m,size,false,size,1);
        else
            return UpperTriMatrixView<T,NonUnitDiag>(m,size,false,1,size);
    }

    template <class T> 
    inline ConstUpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,false,size,1);
        else
            return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,false,1,size);
    }

    template <class T> 
    inline UpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return UpperTriMatrixView<T,NonUnitDiag>(m,size,false,stepi,stepj); }

    template <class T> 
    inline ConstUpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,false,stepi,stepj); }

    template <class T> 
    inline UpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return UpperTriMatrixView<T,UnitDiag>(m,size,true,size,1);
        else
            return UpperTriMatrixView<T,UnitDiag>(m,size,true,1,size);
    }

    template <class T> 
    inline ConstUpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T,UnitDiag>(m,size,true,size,1);
        else
            return ConstUpperTriMatrixView<T,UnitDiag>(m,size,true,1,size);
    }

    template <class T> 
    inline UpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return UpperTriMatrixView<T,UnitDiag>(m,size,true,stepi,stepj); }

    template <class T> 
    inline ConstUpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstUpperTriMatrixView<T,UnitDiag>(m,size,true,stepi,stepj); }

    template <class T> 
    inline LowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return LowerTriMatrixView<T,NonUnitDiag>(m,size,false,size,1);
        else
            return LowerTriMatrixView<T,NonUnitDiag>(m,size,false,1,size);
    }

    template <class T> 
    inline ConstLowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,false,size,1);
        else
            return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,false,1,size);
    }

    template <class T> 
    inline LowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return LowerTriMatrixView<T,NonUnitDiag>(m,size,false,stepi,stepj); }

    template <class T> 
    inline ConstLowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,false,stepi,stepj); }

    template <class T> 
    inline LowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return LowerTriMatrixView<T,UnitDiag>(m,size,true,size,1);
        else
            return LowerTriMatrixView<T,UnitDiag>(m,size,true,1,size);
    }

    template <class T> 
    inline ConstLowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T,UnitDiag>(m,size,true,size,1);
        else
            return ConstLowerTriMatrixView<T,UnitDiag>(m,size,true,1,size);
    }

    template <class T> 
    inline LowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return LowerTriMatrixView<T,UnitDiag>(m,size,true,stepi,stepj); }

    template <class T> 
    inline ConstLowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstLowerTriMatrixView<T,UnitDiag>(m,size,true,stepi,stepj); }



    //
    // Swap
    //

    template <class T, DiagType D, StorageType S, IndexStyle I>
    inline void Swap(UpperTriMatrix<T,D,S,I>& m1, UpperTriMatrix<T,D,S,I>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, DiagType D, 
              int Si, int Sj, bool C, IndexStyle I>
    inline void Swap(
        BaseMatrix_Tri<M>& m1, UpperTriMatrixView<T,D,Si,Sj,C,I> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, DiagType D, 
              int Si, int Sj, bool C, IndexStyle I>
    inline void Swap(
        UpperTriMatrixView<T,D,Si,Sj,C,I> m1, BaseMatrix_Tri<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, 
              int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        UpperTriMatrixView<T,D,Si1,Sj1,C1,I1> m1,
        UpperTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }

    template <class T, DiagType D, StorageType S, IndexStyle I>
    inline void Swap(LowerTriMatrix<T,D,S,I>& m1, LowerTriMatrix<T,D,S,I>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, DiagType D, 
              int Si, int Sj, bool C, IndexStyle I>
    inline void Swap(
        BaseMatrix_Tri<M>& m1, LowerTriMatrixView<T,D,Si,Sj,C,I> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, DiagType D, 
              int Si, int Sj, bool C, IndexStyle I>
    inline void Swap(
        LowerTriMatrixView<T,D,Si,Sj,C,I> m1, BaseMatrix_Tri<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, 
              int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        LowerTriMatrixView<T,D,Si1,Sj1,C1,I1> m1, 
        LowerTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }


    //
    // TMV_Text 
    //

    template <class T, DiagType D, StorageType S, IndexStyle I>
    inline std::string TMV_Text(const UpperTriMatrix<T,D,S,I>& m)
    {
        std::ostringstream s;
        s << "UpperTriMatrix<"<<TMV_Text(T())<<","<<TMV_Text(D);
        if (D == UnknownDiag) s << "("<<(m.isunit()?"Unit":"NonUnit")<<")";
        s << ","<<TMV_Text(S)<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
        const ConstUpperTriMatrixView<T,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstUpperTriMatrixView<"<<TMV_Text(T())<<","<<TMV_Text(D);
        if (D == UnknownDiag) s << "("<<(m.isunit()?"Unit":"NonUnit")<<")";
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(const UpperTriMatrixView<T,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "UpperTriMatrixView<"<<TMV_Text(T())<<","<<TMV_Text(D);
        if (D == UnknownDiag) s << "("<<(m.isunit()?"Unit":"NonUnit")<<")";
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, DiagType D, StorageType S, IndexStyle I>
    inline std::string TMV_Text(const LowerTriMatrix<T,D,S,I>& m)
    {
        std::ostringstream s;
        s << "LowerTriMatrix<"<<TMV_Text(T())<<","<<TMV_Text(D);
        if (D == UnknownDiag) s << "("<<(m.isunit()?"Unit":"NonUnit")<<")";
        s << ","<<TMV_Text(S)<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
        const ConstLowerTriMatrixView<T,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstLowerTriMatrixView<"<<TMV_Text(T())<<","<<TMV_Text(D);
        if (D == UnknownDiag) s << "("<<(m.isunit()?"Unit":"NonUnit")<<")";
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(const LowerTriMatrixView<T,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "LowerTriMatrixView<"<<TMV_Text(T())<<","<<TMV_Text(D);
        if (D == UnknownDiag) s << "("<<(m.isunit()?"Unit":"NonUnit")<<")";
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }


} // namespace tmv

#endif
