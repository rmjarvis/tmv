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
//    In addition to the type template parameter (T), TriMatrixes have two
//    additional template parameters:
//        DiagType dt = UnitDiag || NonUnitDiag 
//        StorageType stor = RowMajor || ColMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {NonUnitDiag, RowMajor}
//
//        If dt is UnitDiag, then the diagonal elements are not
//        actually stored or referenced.  The are all taken to be = 1.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    TriMatrix<T,dt,stor>(int n)
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    TriMatrix<T,dt,stor>(int n, T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    TriMatrix<T,dt,stor>(const Matrix<T>& m)
//    TriMatrix<T,dt,stor>(const TriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//    ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
//        const T* m, size, stor)
//    UpperTriMatrixView<T> UpperTriMatrixViewOf(T* m, size, stor)
//    ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
//        const T* m, size, stor)
//    LowerTriMatrixView<T> LowerTriMatrixViewOf(T* m, size, stor)
//        Make a TriMatrix view of the specific location in memory
//        specified by m.  The size and storage are parameters.
//
//    ConstUpperTriMatrixView<T> UnitUpperTriMatrixViewOf(
//        const T* m, size, stor)
//    UpperTriMatrixView<T> UnitUpperTriMatrixViewOf(T* m, size, stor)
//    ConstLowerTriMatrixView<T> UnitLowerTriMatrixViewOf(
//        const T* m, size, stor)
//    LowerTriMatrixView<T> UnitLowerTriMatrixViewOf(T* m, size, stor)
//        Same as above, but with unit-diagonal.
//
// Access Functions
//
//    int colsize() const
//    int rowsize() const
//    int size() const
//        Return the dimensions of the TriMatrix
//
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the TriMatrix
//
//    Vector& row(int i, int j1, int j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row.
//
//    Vector& col(int j, int i1, int i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column.
//
//    Vector& diag()
//        Return the main diagonal
//        The TriMatrix must be NonUnitDiag.
//
//    Vector& diag(int i, int j1, int j2)
//    Vector& diag(int i)
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal subVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//        i>0 will give an error for a LowerTriMatrix
//        i<0 will give an error for an UpperTriMatrix
//        i=0 will give an error for a UnitDiag TriMatrix
//
// Modifying Functions:
//
//    setZero()
//    setAllTo(T x)
//    addToAll(T x)
//    conjugateSelf()
//    setToIdentity(x = 1)
//    void Swap(TriMatrix& m1, TriMatrix& m2)
//        The TriMatrices must be the same size and shape (Upper or Lower).
//
// Views of a TriMatrix:
//
//    subMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the TriMatrix.
//
//    subVector(int i, int j, int istep, int jstep, int size)
//        Returns a subVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//
//    subTriMatrix(int i1, int i2, int istep)
//        Returns the TriMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with an UpperTriMatrix of size 10, the x's below
//        are the original data, the O's are the subTriMatrix returned
//        with the command subTriMatrix(3,11,2), and the #'s are the 
//        subTriMatrix returned with subTriMatrix(0,3)
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
//    offDiag()
//        Returns the (NonUnitDiag) TriMatrix of all the off-diagonal
//        elements of a TriMatrix.
//
//    realPart(), imagPart()
//        For a complex TriMatrix, returns the real or imaginary part
//        as a real TriMatrix.
//
//    view()
//    transpose()
//    adjoint()
//        Note that the transpose or adjoint of an UpperTriMatrix returns 
//        a view which is a LowerTriMatrix, and vice versa.
//    conjugate()
//
//    viewAsUnitDiag()
//        Returns a UnitDiag view of a TriMatrix.
//        i.e. a TriMatrix of the same elements, but treating the diagonl
//        as all 1's.
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
//    m.makeInverse(minv) // Takes either a TriMatrix or Matrix argument
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
//    is >> CompactIO() >> m
//        reads m from istream in either format
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

#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Array.h"
#include <vector>

namespace tmv {

    template <class T, bool C>
    struct TriRefHelper // C = false
    {
        typedef T& reference;
        TriRefHelper(T& _r, bool ) : itsref(_r) {}
        TriRefHelper(const TriRefHelper& rhs) : itsref(rhs.itsref) {}
        reference itsref;
    };

    template <class T>
    struct TriRefHelper<T,true> // T is real
    {
        typedef T& reference;
        TriRefHelper(T& _r, bool ) : itsref(_r) {}
        TriRefHelper(const TriRefHelper& rhs) : itsref(rhs.itsref) {}
        reference itsref;
    };

    template <class T>
    struct TriRefHelper<std::complex<T>,true> // complex and maybe conj
    {
        typedef VarConjRef<std::complex<T> > reference;
        TriRefHelper(std::complex<T>& _r, bool _c) : 
            itsref(_r,_c?Conj:NonConj) {}
        TriRefHelper(const TriRefHelper& rhs) : itsref(rhs.itsref) {}
        reference itsref;
    };

    // A helper struct to provide a valid reference when the DiagType
    // is UnitDiag (or UnknownDiag)
    // It is very similar to VarConjRef, but also has a boolean isunit
    // which indicates whether the reference is on the diagonal of a 
    // UnitDiag TriMatrix.  If so, it can be an rhs with the value 1, 
    // but not a lhs.
    // If C is true, then the underlying reference is a VarCojRef.
    // Otherwise it is a simple T&.
    template <class T, bool C>
    class TriRef
    {
    public:
        typedef typename TriRefHelper<T,C>::reference reference;
        explicit inline TriRef(bool _u, T& _v, bool _c) :
            isunit(_u), helper(_v,_c) {}
        inline TriRef(const TriRef<T,C>& rhs) :
            isunit(rhs.isunit), helper(rhs.helper) {}
        inline ~TriRef() {}

        inline operator T() const { return val(); }
        inline reference getRef() { return ref(); }
        inline T operator-() const { return -val(); }

        inline TriRef<T,C>& operator=(const TriRef<T,C>& rhs)
        { assign(rhs.val()); return *this; }
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
        { return x1 == x2.val(); }
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
        void check(T2 TMV_DEBUGPARAM(x)) const
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
        reference ref() 
        {
            check();
            return helper.itsref; 
        }
        template <class T2>
        void assign(T2 x) 
        {
            check(x);
            if (!isunit) helper.itsref = x; 
        }
        T val() const 
        { return isunit ? T(1) : T(helper.itsref); }

        const bool isunit;
        TriRefHelper<T,C> helper;
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
    { TMV_SWAP(x1.getRef(),x2.getRef()); }
    template <class T, bool C2>
    inline void TMV_SWAP(T& x1, TriRef<T,C2> x2)
    { TMV_SWAP(x1,x2.getRef()); }
    template <class T, bool C1>
    inline void TMV_SWAP(TriRef<T,C1> x1, T& x2)
    { TMV_SWAP(x1.getRef(),x2); }

    // Another helper class to deal with the case of the regular 
    // UpperTriMatrix and LowerTriMatrix classes (i.e. not views)
    // Here, we can use a simple T& for the reference if D is NonUnitDiag
    // Also, here C is always false.
    template <class T, DiagType D>
    struct TriRefHelper2 // D = NonUnitDiag
    {
        typedef T& reference;
        static reference makeRef(bool, T& r) { return r; }
    };
    template <class T>
    struct TriRefHelper2<T,UnitDiag>
    {
        typedef TriRef<T,false> reference;
        static reference makeRef(bool ondiag, T& r) 
        { return reference(ondiag,r,false); }
    };

    template <class T> 
    class GenUpperTriMatrix : 
        virtual public AssignableToUpperTriMatrix<T>,
        public BaseMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenUpperTriMatrix<T> type;
        typedef UpperTriMatrix<T> copy_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef ConstMatrixView<T> const_rec_type;
        typedef ConstUpperTriMatrixView<T> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T> const_lowertri_type;
        typedef const_uppertri_type const_view_type;
        typedef const_lowertri_type const_transpose_type;
        typedef const_uppertri_type const_conjugate_type;
        typedef const_lowertri_type const_adjoint_type;
        typedef ConstUpperTriMatrixView<RT> const_realpart_type;
        typedef UpperTriMatrixView<T> nonconst_type;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;

        //
        // Constructors
        //

        inline GenUpperTriMatrix() {}
        inline GenUpperTriMatrix(const type&) {}
        virtual inline ~GenUpperTriMatrix() {}

        //
        // Access Functions
        //

        using AssignableToUpperTriMatrix<T>::size;
        using AssignableToUpperTriMatrix<T>::dt;
        inline int colsize() const { return size(); }
        inline int rowsize() const { return size(); }

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            return cref(i,j);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            TMVAssert(j1==j2 || okij(i,j1));
            return const_vec_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct()); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1<=i2 && i2<=size());
            TMVAssert(i1==i2 || okij(i2-1,j));
            return const_vec_type(
                cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct()); 
        }

        inline const_vec_type diag() const
        {
            TMVAssert(!isunit());
            return const_vec_type(cptr(),size(),stepi()+stepj(),ct()); 
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            return const_vec_type(
                cptr()+i*stepj(),size()-i,stepi()+stepj(),ct()); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            TMVAssert(j1>=0 && j1 <= j2 && j2 <= size()-i);
            const int ds = stepi()+stepj();
            return const_vec_type(cptr()+i*stepj()+j1*ds,j2-j1,ds,ct());
        }

        template <class T2> 
        inline bool isSameAs(const BaseMatrix<T2>& ) const
        { return false; }

        inline bool isSameAs(const GenUpperTriMatrix<T>& m2) const
        { 
            return 
                this == &m2 || 
                ( cptr()==m2.cptr() && size()==m2.size() && 
                  dt() == m2.dt() && ct() == m2.ct() &&
                  stepi()==m2.stepi() && stepj()==m2.stepj() 
                );
        }

        inline void assignToM(const MatrixView<RT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TMVAssert(isReal(T()));
            assignToU(m2.upperTri(dt()));
            if (isunit()) m2.diag().setAllTo(RT(1));
            if (size() > 0) m2.lowerTri().offDiag().setZero();
        }

        inline void assignToM(const MatrixView<CT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            assignToU(m2.upperTri(dt()));
            if (isunit()) m2.diag().setAllTo(T(1));
            if (size() > 0) m2.lowerTri().offDiag().setZero();
        }

        inline void assignToU(const UpperTriMatrixView<RT>& m2) const
        { 
            TMVAssert(m2.size() == size());
            TMVAssert(isunit() || !m2.isunit());
            TMVAssert(isReal(T()));
            if (!isSameAs(m2)) Copy(*this,m2); 
        }

        inline void assignToU(const UpperTriMatrixView<CT>& m2) const
        { 
            TMVAssert(m2.size() == size());
            TMVAssert(isunit() || !m2.isunit());
            if (!isSameAs(m2)) Copy(*this,m2); 
        }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const;

        bool hasSubTriMatrix(int i1, int i2, int istep) const;

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        {
            return const_rec_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(), stepj(), stor(), ct());
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return const_rec_type(
                cptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(), 
                newstor, ct());
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            return const_vec_type(
                cptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),ct());
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return cSubVector(i,j,istep,jstep,size);
        }

        inline const_uppertri_type cSubTriMatrix(int i1, int i2) const
        {
            return const_uppertri_type(
                cptr()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),dt(),stor(),ct());
        }

        inline const_uppertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return cSubTriMatrix(i1,i2);
        }

        inline const_uppertri_type cSubTriMatrix(
            int i1, int i2, int istep) const
        {
            return const_uppertri_type(
                cptr()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
                istep==1 ? stor() : NoMajor,ct());
        }

        inline const_uppertri_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return cSubTriMatrix(i1,i2,istep);
        }

        inline const_uppertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return const_uppertri_type(
                cptr()+noff*stepj(),size()-noff,
                stepi(),stepj(),NonUnitDiag,stor(),ct());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                dt(), isReal(T()) ? stor() : NoMajor, NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            // Since Imag of a UnitDiag TriMatrix has 0's on diagonal.
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1, size(),
                2*stepi(), 2*stepj(), dt(), NoMajor, NonConj);
        }

        inline const_uppertri_type view() const
        { 
            return const_uppertri_type(
                cptr(),size(),stepi(),stepj(),dt(),stor(),ct());
        }

        inline const_uppertri_type viewAsUnitDiag() const
        { 
            return const_uppertri_type(
                cptr(),size(),stepi(),stepj(),UnitDiag,stor(),ct());
        }

        inline const_lowertri_type transpose() const
        { 
            return const_lowertri_type(
                cptr(),size(),stepj(),stepi(),dt(),TMV_TransOf(stor()),ct());
        }

        inline const_uppertri_type conjugate() const
        { 
            return const_uppertri_type(
                cptr(),size(),stepi(),stepj(),dt(),stor(),TMV_ConjOf(T,ct()));
        }

        inline const_lowertri_type adjoint() const
        { 
            return const_lowertri_type(
                cptr(),size(),stepj(),stepi(),dt(),
                TMV_TransOf(stor()),TMV_ConjOf(T,ct()));
        }

        inline nonconst_type nonConst() const
        {
            const int n=size();
            return nonconst_type(
                const_cast<T*>(cptr()),n,
                stepi(),dt(),stor(),ct()
                TMV_FIRSTLAST1(cptr(),row(n-1,n-1,n).end().getP()));
        }


        //
        // Functions of Matrix
        //

        T det() const;

        RT logDet(T* sign=0) const;

        inline T trace() const
        { return isunit() ? T(RT(size())) : diag().sumElements(); }

        T sumElements() const;

        RT sumAbsElements() const;

        RT sumAbs2Elements() const;

        RT norm() const
        { return normF(); }

        RT normF() const;

        RT normSq(const RT scale = RT(1)) const;

        RT norm1() const;

        RT doNorm2() const;
        inline RT norm2() const
        { return doNorm2(); }

        RT doCondition() const;
        inline RT condition() const
        { return doCondition(); }

        RT normInf() const;

        RT maxAbsElement() const;
        RT maxAbs2Element() const;

        bool isSingular() const { return det() == T(0); }

        template <class T1> 
        void doMakeInverse(const UpperTriMatrixView<T1>& minv) const;
        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;
        void doMakeInverseATA(const MatrixView<T>& ata) const;

        inline void makeInverse(const MatrixView<T>& minv) const
        {
            TMVAssert(minv.colsize() == size());
            TMVAssert(minv.rowsize() == size());
            doMakeInverse(minv);
        }

        template <class T1> 
        inline void makeInverse(const MatrixView<T1>& minv) const
        {
            TMVAssert(minv.colsize() == size());
            TMVAssert(minv.rowsize() == size());
            doMakeInverse(minv);
        }

        template <class T1> 
        inline void makeInverse(const UpperTriMatrixView<T1>& minv) const
        {
            TMVAssert(minv.size() == size());
            doMakeInverse(minv);
        }

        QuotXU<T,T> QInverse() const;
        inline QuotXU<T,T> inverse() const
        { return QInverse(); }

        inline void makeInverseATA(const MatrixView<T>& ata) const
        {
            TMVAssert(ata.colsize() == size());
            TMVAssert(ata.rowsize() == size());
            doMakeInverseATA(ata);
        }

        template <class T1, DiagType D, StorageType S, IndexStyle I> 
        inline void makeInverse(UpperTriMatrix<T1,D,S,I>& minv) const
        { 
            TMVAssert(D==NonUnitDiag || isunit());
            makeInverse(minv.view()); 
        }

        template <class T1, StorageType S, IndexStyle I> 
        inline void makeInverse(Matrix<T1,S,I>& minv) const
        { makeInverse(minv.view()); }

        template <StorageType S, IndexStyle I> 
        inline void makeInverseATA(Matrix<T,S,I>& minv) const
        { makeInverseATA(minv.view()); }

        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        //
        // Arithmetic Helpers
        //

        template <class T1> 
        void doLDivEq(const VectorView<T1>& v) const;
        template <class T1, class T0> 
        void doLDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const;
        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;
        template <class T1, class T0> 
        void doLDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;
        template <class T1> 
        void doLDivEq(const UpperTriMatrixView<T1>& m) const;
        template <class T1, class T0> 
        void doLDiv(
            const GenUpperTriMatrix<T1>& m1,
            const UpperTriMatrixView<T0>& m0) const;

        template <class T1> 
        inline void LDivEq(const VectorView<T1>& v) const
        { 
            TMVAssert(v.size() == size());
            doLDivEq(v);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { 
            TMVAssert(v0.size() == size());
            TMVAssert(v1.size() == size());
            doLDiv(v1,v0);
        }
        template <class T1> 
        inline void RDivEq(const VectorView<T1>& v) const
        { transpose().LDivEq(v); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { transpose().LDiv(v1,v0); }

        template <class T1> 
        inline void LDivEq(const MatrixView<T1>& m) const
        { 
            TMVAssert(m.colsize() == size());
            doLDivEq(m);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { 
            TMVAssert(m0.colsize() == size());
            TMVAssert(m1.colsize() == size());
            TMVAssert(m0.rowsize() == m1.rowsize());
            doLDiv(m1,m0);
        }
        template <class T1> 
        inline void RDivEq(const MatrixView<T1>& m) const
        { transpose().LDivEq(m.transpose()); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { transpose().LDiv(m1.transpose(),m0.transpose()); }

        template <class T1> 
        inline void LDivEq(const UpperTriMatrixView<T1>& m) const
        { 
            TMVAssert(m.colsize() == size());
            doLDivEq(m);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenUpperTriMatrix<T1>& m1,
            const UpperTriMatrixView<T0>& m0) const
        { 
            TMVAssert(m0.size() == size());
            TMVAssert(m1.size() == size());
            doLDiv(m1,m0);
        }
        template <class T1> 
        inline void RDivEq(const UpperTriMatrixView<T1>& m) const
        { transpose().LDivEq(m.transpose()); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenUpperTriMatrix<T1>& m1,
            const UpperTriMatrixView<T0>& m0) const
        { transpose().LDiv(m1.transpose(),m0.transpose()); }

        // For easier compatibility with regular matrices:
        inline void divideInPlace() const {}
        inline void saveDiv() const {}
        inline void setDiv() const {}
        inline void unsetDiv() const {}
        inline void resetDiv() const {}
        inline bool divIsSet() const { return true; }
        inline void divideUsing(DivType TMV_DEBUGPARAM(dt)) const
        { TMVAssert(dt == LU); }

        virtual const T* cptr() const = 0;
        virtual int stepi() const = 0;
        virtual int stepj() const = 0;
        virtual StorageType stor() const = 0;
        virtual ConjType ct() const = 0;
        inline bool isrm() const { return stor() == RowMajor; }
        inline bool iscm() const { return stor() == ColMajor; }
        inline bool isunit() const { return dt() == UnitDiag; }
        inline bool isconj() const
        {
            TMVAssert(isComplex(T()) || ct()==NonConj);
            return isComplex(T()) && ct()==Conj;
        }

        virtual inline T cref(int i, int j) const
        {
            if (isunit() && i==j) return T(1);
            else if (i>j) return T(0);
            else {
                const T* mi = cptr() + i*stepi() + j*stepj();
                return (isconj() ? TMV_CONJ(*mi) : *mi);
            }
        }

        inline int rowstart(int i) const { return i; }
        inline int rowend(int ) const { return size(); }

        inline int colstart(int ) const { return 0; }
        inline int colend(int j) const { return j+1; }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        { return const_rowmajor_iterator(this,size(),size()); }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        { return const_colmajor_iterator(this,0,size()); }

    protected :

        inline bool okij(int i, int j) const
        {
            TMVAssert(i>=0 && i < size());
            TMVAssert(j>=0 && j < size());
            if (isunit()) return i<j; else return i<=j;
        }

    private :

        type& operator=(const type&);

    }; // GenUpperTriMatrix

    template <class T> 
    class GenLowerTriMatrix : 
        virtual public AssignableToLowerTriMatrix<T>,
        public BaseMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenLowerTriMatrix<T> type;
        typedef LowerTriMatrix<T> copy_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef ConstMatrixView<T> const_rec_type;
        typedef ConstUpperTriMatrixView<T> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T> const_lowertri_type;
        typedef const_lowertri_type const_view_type;
        typedef const_uppertri_type const_transpose_type;
        typedef const_lowertri_type const_conjugate_type;
        typedef const_uppertri_type const_adjoint_type;
        typedef ConstLowerTriMatrixView<RT> const_realpart_type;
        typedef LowerTriMatrixView<T> nonconst_type;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;

        //
        // Constructors
        //

        inline GenLowerTriMatrix() {}
        inline GenLowerTriMatrix(const type&) {}
        virtual inline ~GenLowerTriMatrix() {}

        //
        // Access Functions
        //

        using AssignableToLowerTriMatrix<T>::size;
        using AssignableToLowerTriMatrix<T>::dt;
        inline int colsize() const { return size(); }
        inline int rowsize() const { return size(); }

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            return cref(i,j);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            TMVAssert(j1==j2 || okij(i,j2-1));
            return const_vec_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct()); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1<=i2 && i2<=size());
            TMVAssert(i1==i2 || okij(i1,j));
            return const_vec_type(
                cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct()); 
        }

        inline const_vec_type diag() const
        {
            TMVAssert(!isunit());
            return const_vec_type(cptr(),size(),stepi()+stepj(),ct()); 
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            return const_vec_type(
                cptr()-i*stepi(),size()-i,stepi()+stepj(),ct());
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            TMVAssert(j1>=0 && j1 <= j2 && j2 <= size()+i);
            const int ds = stepi()+stepj();
            return const_vec_type(cptr()-i*stepi()+j1*ds,j2-j1,ds,ct());
        }

        template <class T2> 
        inline bool isSameAs(const BaseMatrix<T2>& ) const
        { return false; }

        inline bool isSameAs(const GenLowerTriMatrix<T>& m2) const
        { 
            if (this == &m2) return true;
            else return (cptr()==m2.cptr() && size()==m2.size() && 
                         dt() == m2.dt() && ct() == m2.ct() &&
                         stepi()==m2.stepi() && stepj()==m2.stepj());
        }

        inline void assignToM(const MatrixView<RT>& m2) const
        { transpose().assignToM(m2.transpose()); }
        inline void assignToM(const MatrixView<CT>& m2) const
        { transpose().assignToM(m2.transpose()); }
        inline void assignToL(const LowerTriMatrixView<RT>& m2) const
        { transpose().assignToU(m2.transpose()); }
        inline void assignToL(const LowerTriMatrixView<CT>& m2) const
        { transpose().assignToU(m2.transpose()); }

        //
        // subMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return transpose().hasSubMatrix(j1,j2,i1,i2,jstep,istep); }

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const
        { return transpose().hasSubVector(j,i,jstep,istep,size); }

        inline bool hasSubTriMatrix(int i1, int i2, int istep) const
        { return transpose().hasSubTriMatrix(i1,i2,istep); }

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        {
            return const_rec_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(), stepj(), stor(), ct());
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return const_rec_type(
                cptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(), 
                newstor,ct());
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            return const_vec_type(
                cptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),ct());
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return cSubVector(i,j,istep,jstep,size);
        }

        inline const_lowertri_type cSubTriMatrix(int i1, int i2) const
        {
            return const_lowertri_type(
                cptr()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),dt(),stor(),ct());
        }

        inline const_lowertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return cSubTriMatrix(i1,i2);
        }

        inline const_lowertri_type cSubTriMatrix(
            int i1, int i2, int istep) const
        {
            return const_lowertri_type(
                cptr()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
                istep==1 ? stor() : NoMajor,ct());
        }

        inline const_lowertri_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return cSubTriMatrix(i1,i2,istep);
        }

        inline const_lowertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return const_lowertri_type(
                cptr()+noff*stepi(),size()-noff,
                stepi(),stepj(),NonUnitDiag,stor(),ct());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                dt(), isReal(T()) ? stor() : NoMajor, NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1, size(),
                2*stepi(), 2*stepj(), dt(), NoMajor, NonConj);
        }

        inline const_lowertri_type view() const
        { 
            return const_lowertri_type(
                cptr(),size(),
                stepi(),stepj(),dt(),stor(),ct());
        }

        inline const_lowertri_type viewAsUnitDiag() const
        { 
            return const_lowertri_type(
                cptr(),size(),
                stepi(),stepj(),UnitDiag,stor(),ct());
        }

        inline const_uppertri_type transpose() const
        { 
            return const_uppertri_type(
                cptr(),size(),
                stepj(),stepi(),dt(),TMV_TransOf(stor()),ct());
        }

        inline const_lowertri_type conjugate() const
        { 
            return const_lowertri_type(
                cptr(),size(),
                stepi(),stepj(),dt(),stor(),TMV_ConjOf(T,ct()));
        }

        inline const_uppertri_type adjoint() const
        { 
            return const_uppertri_type(
                cptr(),size(),
                stepj(),stepi(),dt(),TMV_TransOf(stor()),TMV_ConjOf(T,ct()));
        }

        inline nonconst_type nonConst() const
        {
            const int n=size();
            return nonconst_type(
                const_cast<T*>(cptr()),n,
                stepi(),dt(),stor(),ct()
                TMV_FIRSTLAST1(cptr(),row(n-1,n-1,n).end().getP()));
        }

        //
        // Functions of Matrix
        //

        T det() const
        { return transpose().det(); }

        RT logDet(T* sign=0) const
        { return transpose().logDet(sign); }

        inline T trace() const
        { return isunit() ? T(RT(size())) : diag().sumElements(); }

        inline T sumElements() const
        { return transpose().sumElements(); }

        inline RT sumAbsElements() const
        { return transpose().sumAbsElements(); }

        inline RT sumAbs2Elements() const
        { return transpose().sumAbs2Elements(); }

        inline RT norm() const 
        { return transpose().normF(); }

        inline RT normF() const 
        { return transpose().normF(); }

        inline RT normSq(const RT scale = RT(1)) const
        { return transpose().normSq(scale); }

        inline RT norm1() const
        { return transpose().normInf(); }

        inline RT doNorm2() const
        { return transpose().doNorm2(); }

        inline RT doCondition() const
        { return transpose().doCondition(); }

        inline RT norm2() const
        { return transpose().norm2(); }

        inline RT condition() const
        { return transpose().condition(); }

        inline RT normInf() const
        { return transpose().norm1(); }

        inline RT maxAbsElement() const
        { return transpose().maxAbsElement(); }

        inline RT maxAbs2Element() const
        { return transpose().maxAbs2Element(); }

        bool isSingular() const { return det() == T(0); }

        QuotXL<T,T> QInverse() const;
        inline QuotXL<T,T> inverse() const
        { return QInverse(); }

        template <class T1> 
        inline void makeInverse(const LowerTriMatrixView<T1>& minv) const
        {
            TMVAssert(minv.size() == size());
            transpose().makeInverse(minv.transpose());
        }

        inline void makeInverse(const MatrixView<T>& minv) const
        {
            TMVAssert(minv.colsize() == size());
            TMVAssert(minv.rowsize() == size());
            transpose().makeInverse(minv.transpose());
        }

        template <class T1> 
        inline void makeInverse(const MatrixView<T1>& minv) const
        {
            TMVAssert(minv.colsize() == size());
            TMVAssert(minv.rowsize() == size());
            transpose().makeInverse(minv.transpose());
        }

        void doMakeInverseATA(const MatrixView<T>& minv) const;

        inline void makeInverseATA(const MatrixView<T>& ata) const
        {
            TMVAssert(ata.colsize() == size());
            TMVAssert(ata.rowsize() == size());
            doMakeInverseATA(ata);
        }

        template <class T1, DiagType D, StorageType S, IndexStyle I> 
        inline void makeInverse(LowerTriMatrix<T1,D,S,I>& minv) const
        { 
            TMVAssert(D==NonUnitDiag || isunit());
            makeInverse(minv.view()); 
        }

        template <class T1, StorageType S, IndexStyle I> 
        inline void makeInverse(Matrix<T1,S,I>& minv) const
        { makeInverse(minv.view()); }

        template <StorageType S, IndexStyle I> 
        inline void makeInverseATA(Matrix<T,S,I>& minv) const
        { makeInverseATA(minv.view()); }

        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        //
        // Arithmetic Helpers
        //

        template <class T1> 
        void doLDivEq(const VectorView<T1>& v) const;
        template <class T1, class T0> 
        void doLDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const;
        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;
        template <class T1, class T0> 
        void doLDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;
        template <class T1> 
        void doLDivEq(const LowerTriMatrixView<T1>& m) const;
        template <class T1, class T0> 
        void doLDiv(
            const GenLowerTriMatrix<T1>& m1,
            const LowerTriMatrixView<T0>& m0) const;

        template <class T1> 
        inline void LDivEq(const VectorView<T1>& v) const
        { 
            TMVAssert(v.size() == size());
            doLDivEq(v);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { 
            TMVAssert(v0.size() == size());
            TMVAssert(v1.size() == size());
            doLDiv(v1,v0);
        }
        template <class T1> 
        inline void RDivEq(const VectorView<T1>& v) const
        { transpose().LDivEq(v); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { transpose().LDiv(v1,v0); }

        template <class T1> 
        inline void LDivEq(const MatrixView<T1>& m) const
        { 
            TMVAssert(m.colsize() == size());
            doLDivEq(m);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { 
            TMVAssert(m0.colsize() == size());
            TMVAssert(m1.colsize() == size());
            TMVAssert(m0.rowsize() == m1.rowsize());
            doLDiv(m1,m0);
        }
        template <class T1> 
        inline void RDivEq(const MatrixView<T1>& m) const
        { transpose().LDivEq(m.transpose()); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { transpose().LDiv(m1.transpose(),m0.transpose()); }

        template <class T1> 
        inline void LDivEq(const LowerTriMatrixView<T1>& m) const
        { 
            TMVAssert(m.colsize() == size());
            doLDivEq(m);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenLowerTriMatrix<T1>& m1,
            const LowerTriMatrixView<T0>& m0) const
        { 
            TMVAssert(m0.size() == size());
            TMVAssert(m1.size() == size());
            doLDiv(m1,m0);
        }
        template <class T1> 
        inline void RDivEq(const LowerTriMatrixView<T1>& m) const
        { transpose().LDivEq(m.transpose()); }
        template <class T1, class T0> 
        inline void RDiv(
            const GenLowerTriMatrix<T1>& m1,
            const LowerTriMatrixView<T0>& m0) const
        { transpose().LDiv(m1.transpose(),m0.transpose()); }


        // For easier compatibility with regular matrices:
        inline void divideInPlace() const {}
        inline void saveDiv() const {}
        inline void setDiv() const {}
        inline void unsetDiv() const {}
        inline void resetDiv() const {}
        inline bool divIsSet() const { return true; }
        inline void divideUsing(DivType TMV_DEBUGPARAM(dt)) const
        { TMVAssert(dt == LU); }
        inline bool checkDecomp(std::ostream* fout=0) const { return true; }
        inline bool checkDecomp(
                const BaseMatrix<T>& m2, std::ostream* fout=0) const
        { return true; }

        virtual const T* cptr() const = 0;
        virtual int stepi() const = 0;
        virtual int stepj() const = 0;
        virtual StorageType stor() const = 0;
        virtual ConjType ct() const = 0;
        inline bool isrm() const { return stor() == RowMajor; }
        inline bool iscm() const { return stor() == ColMajor; }
        inline bool isunit() const { return dt() == UnitDiag; }
        inline bool isconj() const
        {
            TMVAssert(isComplex(T()) || ct()==NonConj);
            return isComplex(T()) && ct()==Conj;
        }

        virtual inline T cref(int i, int j) const
        {
            if (isunit() && i==j) return T(1);
            else if (i<j) return T(0);
            else {
                const T* mi = cptr() + i*stepi() + j*stepj();
                return (isconj() ? TMV_CONJ(*mi) : *mi);
            }
        }

        inline int rowstart(int ) const { return 0; }
        inline int rowend(int i) const { return i+1; }

        inline int colstart(int j) const { return j; }
        inline int colend(int ) const { return size(); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        { return const_rowmajor_iterator(this,size(),0); }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        { return const_colmajor_iterator(this,size(),size()); }


    protected :

        inline bool okij(int i, int j) const
        {
            TMVAssert(i>=0 && i < size());
            TMVAssert(j>=0 && j < size());
            if (isunit()) return i>j; else return i>=j;
        }

    private :

        type& operator=(const type&);

    }; // GenLowerTriMatrix

    template <class T, IndexStyle I> 
    class ConstUpperTriMatrixView : public GenUpperTriMatrix<T>
    {
    public :

        typedef GenUpperTriMatrix<T> base;
        typedef ConstUpperTriMatrixView<T,I> type;

        inline ConstUpperTriMatrixView(const type& rhs) :
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct) {}

        inline ConstUpperTriMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itss(rhs.size()), 
            itssi(rhs.stepi()), itssj(rhs.stepj()),
            itsdiag(rhs.dt()), itsstor(rhs.stor()), itsct(rhs.ct()) {}

        inline ConstUpperTriMatrixView(
            const T* _m, int _s, int _si, int _sj,
            DiagType _dt, StorageType _stor, ConjType _ct) : 
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itsdiag(_dt), itsstor(_stor), itsct(_ct)
        { 
            TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
                      _si==1 : true);
        }

        virtual inline ~ConstUpperTriMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        inline int size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline StorageType stor() const { return itsstor; }
        inline DiagType dt() const { return itsdiag; }
        inline ConjType ct() const { return itsct; }

    protected :

        const T*const itsm;
        const int itss;
        const int itssi;
        const int itssj;
        DiagType itsdiag;
        StorageType itsstor;
        ConjType itsct;

    private :

        type& operator=(const type&);

    }; // ConstUpperTriMatrixView

    template <class T, IndexStyle I> 
    class ConstLowerTriMatrixView : public GenLowerTriMatrix<T>
    {
    public :

        typedef GenLowerTriMatrix<T> base;
        typedef ConstLowerTriMatrixView<T,I> type;

        inline ConstLowerTriMatrixView(const type& rhs) :
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct) {}

        inline ConstLowerTriMatrixView(const GenLowerTriMatrix<T>& rhs) :
            itsm(rhs.cptr()), itss(rhs.size()), 
            itssi(rhs.stepi()), itssj(rhs.stepj()),
            itsdiag(rhs.dt()), itsstor(rhs.stor()), itsct(rhs.ct()) {}

        inline ConstLowerTriMatrixView(
            const T* _m, int _s, int _si, int _sj,
            DiagType _dt, StorageType _stor, ConjType _ct) : 
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itsdiag(_dt), itsstor(_stor), itsct(_ct)
        { 
            TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
                      _si==1 : true);
        }

        virtual inline ~ConstLowerTriMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        inline int size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline StorageType stor() const { return itsstor; }
        inline DiagType dt() const { return itsdiag; }
        inline ConjType ct() const { return itsct; }

    protected :

        const T*const itsm;
        const int itss;
        const int itssi;
        const int itssj;
        DiagType itsdiag;
        StorageType itsstor;
        ConjType itsct;

    private :

        type& operator=(const type&);

    }; // ConstLowerTriMatrixView

    template <class T> 
    class ConstUpperTriMatrixView<T,FortranStyle> : 
        public ConstUpperTriMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef GenUpperTriMatrix<T> base;
        typedef ConstUpperTriMatrixView<T,FortranStyle> type;
        typedef ConstUpperTriMatrixView<T,CStyle> c_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_lowertri_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_uppertri_type;
        typedef const_uppertri_type const_view_type;
        typedef const_lowertri_type const_transpose_type;
        typedef const_uppertri_type const_conjugate_type;
        typedef const_lowertri_type const_adjoint_type;
        typedef ConstUpperTriMatrixView<RT,FortranStyle> const_realpart_type;

        inline ConstUpperTriMatrixView(const type& rhs) : c_type(rhs) {}

        inline ConstUpperTriMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline ConstUpperTriMatrixView(const base& rhs) : c_type(rhs) {}

        inline ConstUpperTriMatrixView(
            const T* _m, int _s, int _si, int _sj,
            DiagType indt, StorageType instor, ConjType inct) : 
            c_type(_m,_s,_si,_sj,indt,instor,inct) 
        {}

        virtual inline ~ConstUpperTriMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>0 && i<=size());
            TMVAssert(j>0 && j<=size());
            return cref(i-1,j-1);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=size());
            TMVAssert(j1>0 && j1<=j2 && j2<=size());
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=size());
            TMVAssert(i1>0 && i1<=i2 && i2<=size());
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
        // subMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const;

        bool hasSubTriMatrix(int i1, int i2, int istep) const;

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return base::cSubVector(i-1,j-1,istep,jstep,size);
        }

        inline const_uppertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return base::cSubTriMatrix(i1-1,i2);
        }

        inline const_uppertri_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return base::cSubTriMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_uppertri_type offDiag(int noff=1) const
        { return base::offDiag(noff); }

        inline const_realpart_type realPart() const
        { return base::realPart(); }

        inline const_realpart_type imagPart() const
        { return base::imagPart(); }

        inline const_uppertri_type view() const
        { return base::view(); }

        inline const_uppertri_type viewAsUnitDiag() const
        { return base::viewAsUnitDiag(); }

        inline const_lowertri_type transpose() const
        { return base::transpose(); }

        inline const_uppertri_type conjugate() const
        { return base::conjugate(); }

        inline const_lowertri_type adjoint() const
        { return base::adjoint(); }

        using c_type::size;
        using base::isunit;

        using c_type::cref;

    protected :

        using base::okij;

    private :

        type& operator=(const type&);

    }; // FortranStyle ConstUpperTriMatrixView

    template <class T> 
    class ConstLowerTriMatrixView<T,FortranStyle> : 
        public ConstLowerTriMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef GenLowerTriMatrix<T> base;
        typedef ConstLowerTriMatrixView<T,FortranStyle> type;
        typedef ConstLowerTriMatrixView<T,CStyle> c_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<T,FortranStyle> const_rec_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_lowertri_type;
        typedef const_lowertri_type const_view_type;
        typedef const_uppertri_type const_transpose_type;
        typedef const_lowertri_type const_conjugate_type;
        typedef const_uppertri_type const_adjoint_type;
        typedef ConstLowerTriMatrixView<RT,FortranStyle> const_realpart_type;

        inline ConstLowerTriMatrixView(const type& rhs) : c_type(rhs) {}

        inline ConstLowerTriMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline ConstLowerTriMatrixView(const base& rhs) : c_type(rhs) {}

        inline ConstLowerTriMatrixView(
            const T* _m, int _s, int _si, int _sj,
            DiagType indt, StorageType instor, ConjType inct) : 
            c_type(_m,_s,_si,_sj,indt,instor,inct) 
        {}

        virtual inline ~ConstLowerTriMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(int i, int j) const
        {
            TMVAssert(i>0 && i<=size());
            TMVAssert(j>0 && j<=size());
            return cref(i-1,j-1);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=size());
            TMVAssert(j1>0 && j1<=j2 && j2<=size());
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=size());
            TMVAssert(i1>0 && i1<=i2 && i2<=size());
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
        // subMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return transpose().hasSubMatrix(j1,j2,i1,i2,jstep,istep); }

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int size) const
        { return transpose().hasSubVector(j,i,jstep,istep,size); }

        inline bool hasSubTriMatrix(int i1, int i2, int istep) const
        { return transpose().hasSubTriMatrix(i1,i2,istep); }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return base::cSubVector(i-1,j-1,istep,jstep,size);
        }

        inline const_lowertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return base::cSubTriMatrix(i1-1,i2);
        }

        inline const_lowertri_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return base::cSubTriMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_lowertri_type offDiag(int noff=1) const
        { return base::offDiag(noff); }

        inline const_realpart_type realPart() const
        { return base::realPart(); }

        inline const_realpart_type imagPart() const
        { return base::imagPart(); }

        inline const_lowertri_type view() const
        { return base::view(); }

        inline const_lowertri_type viewAsUnitDiag() const
        { return base::viewAsUnitDiag(); }

        inline const_uppertri_type transpose() const
        { return base::transpose(); }

        inline const_lowertri_type conjugate() const
        { return base::conjugate(); }

        inline const_uppertri_type adjoint() const
        { return base::adjoint(); }

        using c_type::size;
        using base::isunit;

        using c_type::cref;

    protected :

        using base::okij;

    private :

        type& operator=(const type&);

    }; // FortranStyle ConstLowerTriMatrixView

    template <class T, IndexStyle I> 
    class UpperTriMatrixView : public GenUpperTriMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenUpperTriMatrix<T> base;
        typedef UpperTriMatrixView<T,I> type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef uppertri_type view_type;
        typedef lowertri_type transpose_type;
        typedef uppertri_type conjugate_type;
        typedef lowertri_type adjoint_type;
        typedef MatrixView<T,I> rec_type;
        typedef VectorView<T,I> vec_type;
        typedef UpperTriMatrixView<RT,I> realpart_type;
        typedef TriRef<T,true> reference;
        typedef RMIt<const type> rowmajor_iterator;
        typedef CMIt<const type> colmajor_iterator;

        //
        // Constructors
        //

        inline UpperTriMatrixView(const type& rhs) : 
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        inline UpperTriMatrixView(
            T* _m, int _s, int _si, int _sj,
            DiagType _dt, StorageType _stor, ConjType _ct
            TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itsdiag(_dt), itsstor(_stor), itsct(_ct) 
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(_stor==RowMajor ? _sj == 1 :
                      _stor==ColMajor ? _si==1 : true);
        }

        virtual inline ~UpperTriMatrixView() 
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
            m2.assignToU(*this);
            return *this; 
        }

        inline const type& operator=(const type& m2)
        {
            TMVAssert(size() == m2.size());
            m2.assignToU(*this);
            return *this; 
        }

        inline const type& operator=(const GenUpperTriMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToU(*this);
            return *this; 
        }

        inline const type& operator=(const GenUpperTriMatrix<CT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            m2.assignToU(*this);
            return *this; 
        }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToD(DiagMatrixViewOf(diag()));
            offDiag().setZero();
            return *this; 
        }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            m2.assignToD(DiagMatrixViewOf(diag()));
            offDiag().setZero();
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenUpperTriMatrix<T2>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(!isunit() || m2.isunit());
            Copy(m2,*this); 
            return *this; 
        }

        inline const type& operator=(const T& x) const 
        { TMVAssert(!isunit() || x==T(1)); return setToIdentity(x); }

        inline const type& operator=(
            const AssignableToUpperTriMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
            m2.assignToU(view());
            return *this;
        }

        inline const type& operator=(
            const AssignableToUpperTriMatrix<CT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
            m2.assignToU(view());
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return MyListAssigner(rowmajor_begin(),size()*(size()+1)/2,x); }


        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            TMVAssert(i<=j);
            return ref(i,j); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            TMVAssert(j1==j2 || okij(i,j1));
            return vec_type(
                ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1<=i2 && i2<=size());
            TMVAssert(i1==i2 || okij(i2-1,j));
            return vec_type(
                ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type diag() const
        {
            TMVAssert(!isunit());
            return vec_type(
                ptr(),size(),stepi()+stepj(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i) const
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            return vec_type(
                ptr()+i*stepj(),size()-i,(stepi()+stepj()),ct() 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            TMVAssert(j1>=0 && j1 <= j2 && j2 <= size()-i);
            const int ds = stepi()+stepj();
            return vec_type(
                ptr()+i*stepj()+j1*ds,j2-j1,ds,ct() TMV_FIRSTLAST);
        }

        //
        // Modifying Functions
        //

        const type& setZero() const;

        const type& setAllTo(const T& x) const;

        const type& addToAll(const T& x) const;

        const type& clip(RT thresh) const;

        const type& conjugateSelf() const;

        const type& invertSelf() const;

        const type& setToIdentity(const T& x=T(1)) const;

        //
        // subMatrix
        //

        using base::hasSubMatrix;
        using base::hasSubVector;
        using base::hasSubTriMatrix;

        inline rec_type cSubMatrix(int i1, int i2, int j1, int j2) const
        {
            return rec_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(),stepj(),stor(),ct() TMV_FIRSTLAST );
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return rec_type(
                ptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,ct() TMV_FIRSTLAST );
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            return vec_type(
                ptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),ct() TMV_FIRSTLAST );
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return cSubVector(i,j,istep,jstep,size);
        }

        inline uppertri_type cSubTriMatrix(int i1, int i2) const
        {
            return uppertri_type(
                ptr()+i1*(stepi()+stepj()),i2-i1,
                stepi(),stepj(),dt(),stor(),ct() TMV_FIRSTLAST);
        }

        inline uppertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return cSubTriMatrix(i1,i2);
        }

        inline uppertri_type cSubTriMatrix(int i1, int i2, int istep) const
        {
            return uppertri_type(
                ptr()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
                istep==1 ? stor() : NoMajor,ct() TMV_FIRSTLAST);
        }

        inline uppertri_type subTriMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return cSubTriMatrix(i1,i2,istep);
        }

        inline uppertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return uppertri_type(
                ptr()+noff*stepj(),size()-noff,
                stepi(),stepj(),NonUnitDiag,stor(),ct() TMV_FIRSTLAST);
        }

        inline realpart_type realPart() const
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                dt(), isReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1, size(),
                2*stepi(), 2*stepj(), dt(), NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline uppertri_type view() const
        { return *this; }

        inline uppertri_type viewAsUnitDiag() const
        { 
            return uppertri_type(
                ptr(),size(),stepi(),stepj(),UnitDiag,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline lowertri_type transpose() const
        {
            return lowertri_type(
                ptr(),size(),stepj(),stepi(),dt(),TMV_TransOf(stor()),ct() 
                TMV_FIRSTLAST);
        }

        inline uppertri_type conjugate() const
        {
            return uppertri_type(
                ptr(),size(),stepi(),stepj(),dt(),stor(),TMV_ConjOf(T,ct()) 
                TMV_FIRSTLAST);
        }

        inline lowertri_type adjoint() const
        {
            return lowertri_type(
                ptr(),size(),stepj(),stepi(),dt(),
                TMV_TransOf(stor()),TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader) const;

        inline int size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        using base::isconj;
        using base::isrm;
        using base::iscm;
        using base::isunit;
        inline StorageType stor() const { return itsstor; }
        inline DiagType dt() const { return itsdiag; }
        inline ConjType ct() const { return itsct; }

        inline reference ref(int i, int j) const
        {
            T* mi = ptr() + i*stepi() + j*stepj();
            return reference(isunit() && i==j,*mi,ct());
        }

        inline rowmajor_iterator rowmajor_begin() const
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end() const
        { return rowmajor_iterator(this,size(),size()); }

        inline colmajor_iterator colmajor_begin() const
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end() const
        { return colmajor_iterator(this,0,size()); }


    protected :

        T*const itsm;
        const int itss;
        const int itssi;
        const int itssj;
        DiagType itsdiag;
        StorageType itsstor;
        ConjType itsct;

        using base::okij;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

    }; // UpperTriMatrixView

    template <class T, IndexStyle I> 
    class LowerTriMatrixView : public GenLowerTriMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenLowerTriMatrix<T> base;
        typedef LowerTriMatrixView<T,I> type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef lowertri_type view_type;
        typedef uppertri_type transpose_type;
        typedef lowertri_type conjugate_type;
        typedef uppertri_type adjoint_type;
        typedef MatrixView<T,I> rec_type;
        typedef VectorView<T,I> vec_type;
        typedef LowerTriMatrixView<RT,I> realpart_type;
        typedef TriRef<T,true> reference;
        typedef RMIt<const type> rowmajor_iterator;
        typedef CMIt<const type> colmajor_iterator;

        //
        // Constructors
        //

        inline LowerTriMatrixView(const type& rhs) : 
            itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
            itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        inline LowerTriMatrixView(
            T* _m, int _s, int _si, int _sj,
            DiagType _dt, StorageType _stor, ConjType _ct
            TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itss(_s), itssi(_si), itssj(_sj),
            itsdiag(_dt), itsstor(_stor), itsct(_ct) 
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
                      _si==1 : true);
        }

        virtual inline ~LowerTriMatrixView()
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
            m2.assignToL(*this);
            return *this; 
        }

        inline const type& operator=(const type& m2) 
        {
            TMVAssert(size() == m2.size());
            m2.assignToL(*this);
            return *this; 
        }

        inline const type& operator=(const GenLowerTriMatrix<RT>& m2) const
        {
            TMVAssert(size() == m2.size());
            m2.assignToL(*this);
            return *this; 
        }

        inline const type& operator=(const GenLowerTriMatrix<CT>& m2) const
        {
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            m2.assignToL(*this);
            return *this; 
        }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            transpose() = m2;
            return *this; 
        }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            transpose() = m2;
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenLowerTriMatrix<T2>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!isunit() || m2.isunit());
            Copy(m2.transpose(),transpose()); 
            return *this; 
        }

        inline const type& operator=(const T& x) const 
        { 
            TMVAssert(!isunit()); 
            return setToIdentity(x); 
        }

        inline const type& operator=(
            const AssignableToLowerTriMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
            m2.assignToL(view());
            return *this;
        }

        inline const type& operator=(
            const AssignableToLowerTriMatrix<CT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
            m2.assignToL(view());
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return MyListAssigner(rowmajor_begin(),size()*(size()+1)/2,x); }

        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        {
            TMVAssert(i>=0 && i<size());
            TMVAssert(j>=0 && j<size());
            TMVAssert(i>=j);
            return ref(i,j); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<size());
            TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            TMVAssert(j1==j2 || okij(i,j2-1));
            return vec_type(
                ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<size());
            TMVAssert(i1>=0 && i1<=i2 && i2<=size());
            TMVAssert(i1==i2 || okij(i1,j));
            return vec_type(
                ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type diag() const
        {
            TMVAssert(!isunit());
            return vec_type(
                ptr(),size(),stepi()+stepj(),ct() TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i) const
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            return vec_type(
                ptr()-i*stepi(),size()+i,stepi()+stepj(),ct() TMV_FIRSTLAST);
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            TMVAssert(j1>=0 && j1 <= j2 && j2 <= size()+i);
            const int ds = stepi()+stepj();
            return vec_type(
                ptr()-i*stepi()+j1*ds,j2-j1,ds,ct() TMV_FIRSTLAST);
        }

        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { transpose().setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        { transpose().setAllTo(x); return *this; }

        inline const type& addToAll(const T& x) const
        { transpose().addToAll(x); return *this; }

        inline const type& clip(RT thresh) const
        { transpose().clip(thresh); return *this; }

        inline const type& conjugateSelf() const
        { transpose().conjugateSelf(); return *this; }

        inline const type& invertSelf() const
        { transpose().invertSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { transpose().setToIdentity(x); return *this; }

        //
        // subMatrix
        //

        using base::hasSubMatrix;
        using base::hasSubVector;
        using base::hasSubTriMatrix;

        inline rec_type cSubMatrix(int i1, int i2, int j1, int j2) const
        {
            return rec_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(),stepj(),stor(),ct() TMV_FIRSTLAST );
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            const StorageType newstor =
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return rec_type(
                ptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,ct() TMV_FIRSTLAST );
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            return vec_type(
                ptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),ct() TMV_FIRSTLAST );
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(hasSubVector(i,j,istep,jstep,size));
            return cSubVector(i,j,istep,jstep,size);
        }

        inline lowertri_type cSubTriMatrix(int i1, int i2) const
        {
            return lowertri_type(
                ptr()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),dt(),stor(),ct() TMV_FIRSTLAST);
        }

        inline lowertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return cSubTriMatrix(i1,i2);
        }

        inline lowertri_type cSubTriMatrix(int i1, int i2, int istep) const
        {
            return lowertri_type(
                ptr()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
                istep==1 ? stor() : NoMajor,ct() TMV_FIRSTLAST);
        }

        inline lowertri_type subTriMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return cSubTriMatrix(i1,i2,istep);
        }

        inline lowertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return lowertri_type(
                ptr()+noff*stepi(),size()-noff,
                stepi(),stepj(),NonUnitDiag,stor(),ct() TMV_FIRSTLAST);
        }

        inline realpart_type realPart() const
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                dt(), isReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1, size(),
                2*stepi(), 2*stepj(), dt(), NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline lowertri_type view() const
        { return *this; }

        inline lowertri_type viewAsUnitDiag() const
        {
            return lowertri_type(
                ptr(),size(),
                stepi(),stepj(),UnitDiag,stor(),ct() TMV_FIRSTLAST);
        }

        inline uppertri_type transpose() const
        {
            return uppertri_type(
                ptr(),size(),
                stepj(),stepi(),dt(),TMV_TransOf(stor()),ct() TMV_FIRSTLAST);
        }

        inline lowertri_type conjugate() const
        {
            return lowertri_type(
                ptr(),size(),
                stepi(),stepj(),dt(),stor(),TMV_ConjOf(T,ct()) TMV_FIRSTLAST);
        }

        inline uppertri_type adjoint() const
        {
            return uppertri_type(
                ptr(),size(),
                stepj(),stepi(),dt(),TMV_TransOf(stor()),TMV_ConjOf(T,ct()) 
                TMV_FIRSTLAST);
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader) const;

        inline int size() const { return itss; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        using base::isconj;
        using base::isrm;
        using base::iscm;
        using base::isunit;
        inline StorageType stor() const { return itsstor; }
        inline DiagType dt() const { return itsdiag; }
        inline ConjType ct() const { return itsct; }

        inline reference ref(int i, int j) const
        {
            T* mi = ptr() + i*stepi() + j*stepj();
            return reference(isunit() && i==j,*mi,ct());
        }

        inline rowmajor_iterator rowmajor_begin() const
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end() const
        { return rowmajor_iterator(this,size(),0); }

        inline colmajor_iterator colmajor_begin() const
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end() const
        { return colmajor_iterator(this,size(),size()); }


    protected :

        T*const itsm;
        const int itss;
        const int itssi;
        const int itssj;
        DiagType itsdiag;
        StorageType itsstor;
        ConjType itsct;

        using base::okij;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

    }; // LowerTriMatrixView

    template <class T> 
    class UpperTriMatrixView<T,FortranStyle> : 
        public UpperTriMatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenUpperTriMatrix<T> base;
        typedef UpperTriMatrixView<T,FortranStyle> type;
        typedef UpperTriMatrixView<T,CStyle> c_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_type;
        typedef LowerTriMatrixView<T,FortranStyle> lowertri_type;
        typedef UpperTriMatrixView<T,FortranStyle> uppertri_type;
        typedef uppertri_type view_type;
        typedef lowertri_type transpose_type;
        typedef uppertri_type conjugate_type;
        typedef lowertri_type adjoint_type;
        typedef MatrixView<T,FortranStyle> rec_type;
        typedef VectorView<T,FortranStyle> vec_type;
        typedef UpperTriMatrixView<RT,FortranStyle> realpart_type;
        typedef TriRef<T,true> reference;

        //
        // Constructors
        //

        inline UpperTriMatrixView(const type& rhs) : c_type(rhs) {}

        inline UpperTriMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline UpperTriMatrixView(
            T* _m, int _s, int _si, int _sj,
            DiagType indt, StorageType instor, ConjType inct 
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_s,_si,_sj,indt,instor,inct 
                   TMV_FIRSTLAST1(_first,_last) ) 
        {}

        virtual inline ~UpperTriMatrixView() {} 

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

        inline const type& operator=(const GenUpperTriMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenUpperTriMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenUpperTriMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(
            const AssignableToUpperTriMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(
            const AssignableToUpperTriMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        typedef typename c_type::MyListAssigner MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return c_type::operator<<(x); }

        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        { 
            TMVAssert(i>0 && i <= size());
            TMVAssert(j>0 && j <= size());
            TMVAssert(i<=j);
            return ref(i-1,j-1); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=size());
            TMVAssert(j1>0 && j1<=j2 && j2<=size());
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=size());
            TMVAssert(i1>0 && i1<=i2 && i2<=size());
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

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& invertSelf() const
        { c_type::invertSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { c_type::setToIdentity(x); return *this; }

        //
        // SubMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return const_type(*this).hasSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,s); }

        inline bool hasSubTriMatrix(int i1, int i2, int istep) const
        { return const_type(*this).hasSubTriMatrix(i1,i2,istep); }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return c_type::cSubVector(i-1,j-1,istep,jstep,s);
        }

        inline uppertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return c_type::cSubTriMatrix(i1-1,i2);
        }

        inline uppertri_type subTriMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return c_type::cSubTriMatrix(i1-1,i2-1+istep,istep);
        }

        inline uppertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return c_type::offDiag(noff);
        }

        inline realpart_type realPart() const
        { return c_type::realPart(); }

        inline realpart_type imagPart() const
        { return c_type::imagPart(); }

        inline uppertri_type view() const
        { return *this; }

        inline uppertri_type viewAsUnitDiag() const
        { return c_type::viewAsUnitDiag(); }

        inline lowertri_type transpose() const
        { return c_type::transpose(); }

        inline uppertri_type conjugate() const
        { return c_type::conjugate(); }

        inline lowertri_type adjoint() const
        { return c_type::adjoint(); }

        using c_type::size;
        using c_type::ref;

    protected :
        using base::okij;

    }; // FortranStyle UpperTriMatrixView

    template <class T> 
    class LowerTriMatrixView<T,FortranStyle> : 
        public LowerTriMatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenLowerTriMatrix<T> base;
        typedef LowerTriMatrixView<T,FortranStyle> type;
        typedef LowerTriMatrixView<T,CStyle> c_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_type;
        typedef UpperTriMatrixView<T,FortranStyle> uppertri_type;
        typedef LowerTriMatrixView<T,FortranStyle> lowertri_type;
        typedef lowertri_type view_type;
        typedef uppertri_type transpose_type;
        typedef lowertri_type conjugate_type;
        typedef uppertri_type adjoint_type;
        typedef MatrixView<T,FortranStyle> rec_type;
        typedef VectorView<T,FortranStyle> vec_type;
        typedef LowerTriMatrixView<RT,FortranStyle> realpart_type;
        typedef TriRef<T,true> reference;

        //
        // Constructors
        //

        inline LowerTriMatrixView(const type& rhs) : c_type(rhs) {}

        inline LowerTriMatrixView(const c_type& rhs) : c_type(rhs) {}

        inline LowerTriMatrixView(
            T* _m, int _s, int _si, int _sj,
            DiagType indt, StorageType instor, ConjType inct 
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_s,_si,_sj,indt,instor,inct 
                   TMV_FIRSTLAST1(_first,_last) ) 
        {}

        virtual inline ~LowerTriMatrixView() {} 

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

        inline const type& operator=(const GenLowerTriMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenLowerTriMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenLowerTriMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(
            const AssignableToLowerTriMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(
            const AssignableToLowerTriMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        typedef typename c_type::MyListAssigner MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return c_type::operator<<(x); }

        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        { 
            TMVAssert(i>0 && i <= size());
            TMVAssert(j>0 && j <= size());
            TMVAssert(i>=j);
            return ref(i-1,j-1); 
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=size());
            TMVAssert(j1>0 && j1<=j2 && j2<=size());
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=size());
            TMVAssert(i1>0 && i1<=i2 && i2<=size());
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

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& invertSelf() const
        { c_type::invertSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { c_type::setToIdentity(x); return *this; }

        //
        // subMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return const_type(*this).hasSubMatrix(i1,i2,j1,j2,istep,jstep); }

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,s); }

        inline bool hasSubTriMatrix(int i1, int i2, int istep) const
        { return const_type(*this).hasSubTriMatrix(i1,i2,istep); }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return c_type::cSubVector(i-1,j-1,istep,jstep,s);
        }

        inline lowertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,1));
            return c_type::cSubTriMatrix(i1-1,i2);
        }

        inline lowertri_type subTriMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubTriMatrix(i1,i2,istep));
            return c_type::cSubTriMatrix(i1-1,i2-1+istep,istep);
        }

        inline lowertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return c_type::offDiag(noff);
        }

        inline realpart_type realPart() const
        { return c_type::realPart(); }

        inline realpart_type imagPart() const
        { return c_type::imagPart(); }

        inline lowertri_type view() const
        { return *this; }

        inline lowertri_type viewAsUnitDiag() const
        { return c_type::viewAsUnitDiag(); }

        inline uppertri_type transpose() const
        { return c_type::transpose(); }

        inline lowertri_type conjugate() const
        { return c_type::conjugate(); }

        inline uppertri_type adjoint() const
        { return c_type::adjoint(); }

        using c_type::size;
        using c_type::ref;

    protected :
        using base::okij;

    }; // FortranStyle LowerTriMatrixView


    template <class T, DiagType D, StorageType S, IndexStyle I> 
    class UpperTriMatrix : public GenUpperTriMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenUpperTriMatrix<T> base;
        typedef UpperTriMatrix<T,D,S,I> type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstUpperTriMatrixView<T,I> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,I> const_lowertri_type;
        typedef const_uppertri_type const_view_type;
        typedef const_lowertri_type const_transpose_type;
        typedef const_uppertri_type const_conjugate_type;
        typedef const_lowertri_type const_adjoint_type;
        typedef ConstUpperTriMatrixView<RT,I> const_realpart_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef uppertri_type view_type;
        typedef lowertri_type transpose_type;
        typedef uppertri_type conjugate_type;
        typedef lowertri_type adjoint_type;
        typedef UpperTriMatrixView<RT,I> realpart_type;
        typedef typename TriRefHelper2<T,D>::reference reference;
        typedef RMIt<type> rowmajor_iterator;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;

        //
        // Constructors
        //

#define NEW_SIZE(s) \
        itslen((s)*(s)), itsm(itslen), itss(s) \
        TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

        explicit inline UpperTriMatrix(int _size=0) : NEW_SIZE(_size) 
        {
            TMVAssert(_size >= 0);
            TMVAssert(S==RowMajor || S==ColMajor); 
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline UpperTriMatrix(int _size, const T& x) : NEW_SIZE(_size)
        {
            TMVAssert(_size >= 0);
            TMVAssert(S==RowMajor || S==ColMajor);
            setAllTo(x);
        }

        template <class T2> 
        inline UpperTriMatrix(const GenMatrix<T2>& rhs) :
            NEW_SIZE(rhs.rowsize())
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor); 
            Copy(rhs.upperTri(D),view()); 
        }

        template <class T2> 
        inline UpperTriMatrix(const GenUpperTriMatrix<T2>& rhs) :
            NEW_SIZE(rhs.size())
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor); 
            if (isunit() && !rhs.isunit()) {
                if (rhs.size() > 0)
                    Copy(rhs.offDiag(),offDiag());
            } else {
                Copy(rhs,view());
            }
        }

        inline UpperTriMatrix(const type& rhs) :
            itslen(rhs.itslen), itsm(itslen), itss(rhs.itss)
            TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        template <DiagType D2, IndexStyle I2> 
        inline UpperTriMatrix(const UpperTriMatrix<T,D2,S,I2>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
            if (D==NonUnitDiag && D2==UnitDiag) diag().setAllTo(T(1));
        }

        template <IndexStyle I2> 
        inline UpperTriMatrix(const Matrix<T,S,I2>& rhs) :
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            if (rhs.isSquare())
                std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
            else
                Copy(rhs.upperTri(D),view());
        }

        inline UpperTriMatrix(const GenMatrix<T>& rhs) :
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            Copy(rhs.upperTri(D),view());
        }

        inline UpperTriMatrix(const GenUpperTriMatrix<RT>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            if (isunit() && !rhs.isunit()) {
                if (rhs.size() > 0) offDiag() = rhs.offDiag();
            } else 
                rhs.assignToU(view());
        }

        inline UpperTriMatrix(const GenUpperTriMatrix<CT>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()));
            if (isunit() && !rhs.isunit()) {
                if (rhs.size() > 0) offDiag() = rhs.offDiag();
            } else 
                rhs.assignToU(view());
        }

        inline UpperTriMatrix(const AssignableToUpperTriMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToU(view());
        }

        inline UpperTriMatrix(const AssignableToUpperTriMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToU(view());
        }

#undef NEW_SIZE

        virtual inline ~UpperTriMatrix()
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

        template <IndexStyle I2> 
        inline type& operator=(const UpperTriMatrix<T,D,S,I2>& m2)
        { 
            TMVAssert(size() == m2.size());
            if (&m2 != this) 
                std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<RT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToU(view());
            return *this;
        }

        inline type& operator=(const GenUpperTriMatrix<CT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            TMVAssert(isComplex(T()));
            m2.assignToU(view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const GenUpperTriMatrix<T2>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T2()) || sComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            Copy(m2,view());
            return *this;
        }

        inline type& operator=(const T& x) 
        { 
            TMVAssert(!this->isunit() || x==T(1));
            return setToIdentity(x); 
        }

        inline type& operator=(const AssignableToUpperTriMatrix<RT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToU(view());
            return *this;
        }

        inline type& operator=(const AssignableToUpperTriMatrix<CT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            TMVAssert(isComplex(T()));
            m2.assignToU(view());
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(rowmajor_begin(),size()*(size()+1)/2,x); }

        //
        // Access
        //

        inline T operator()(int i, int j) const
        {
            if (I == CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
            } else {
                TMVAssert(i>0 && i<=size()); --i;
                TMVAssert(j>0 && j<=size()); --j;
            }
            return cref(i,j);
        }

        inline reference operator()(int i, int j) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                TMVAssert(i<=j);
            } else {
                TMVAssert(i>0 && i<= size()); --i;
                TMVAssert(j>0 && j<= size()); --j;
                TMVAssert(i<=j);
            }
            return ref(i,j);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            if (I==FortranStyle) { 
                TMVAssert(i>0 && i<=size()); --i;
                TMVAssert(j1>0 && j1<=j2 && j2<=size()); --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            return const_vec_type(
                itsm.get()+i*stepi()+j1*stepj(),j2-j1,stepj(),NonConj);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=size()); --j;
                TMVAssert(i1>0 && i1<=i2 && i2<=size()); --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1<=i2 && i2<=size()); 
            }
            TMVAssert(i1==i2 || okij(i2-1,j));
            return const_vec_type(
                itsm.get()+i1*stepi()+j*stepj(),i2-i1,stepi(),NonConj);
        }

        inline const_vec_type diag() const
        {
            TMVAssert(!isunit());
            return const_vec_type(itsm.get(),size(),stepi()+stepj(),NonConj); 
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            return const_vec_type(
                itsm.get()+i*stepj(),size()-i,stepi()+stepj(),NonConj);
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            if (I == FortranStyle) {
                TMVAssert(j1 > 0 && j1<=j2 && j2<=size()-i); --j1; 
            } else {
                TMVAssert(j1>=0 && j1<=j2 && j2<=size()-i); 
            }
            const int ds = stepi()+stepj();
            return const_vec_type(
                itsm.get()+i*stepj()+j1*ds,j2-j1,ds,NonConj);
        }

        inline vec_type row(int i, int j1, int j2)
        { 
            if (I==FortranStyle) { 
                TMVAssert(i>0 && i<=size()); --i;
                TMVAssert(j1>0 && j1<=j2 && j2<=size()); --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j1));
            return vec_type(
                itsm.get()+i*stepi()+j1*stepj(),
                j2-j1,stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=size()); --j;
                TMVAssert(i1>0 && i1<=i2 && i2<=size()); --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1<=i2 && i2<=size()); 
            }
            TMVAssert(i1==i2 || okij(i2-1,j));
            return vec_type(
                itsm.get()+i1*stepi()+j*stepj(),
                i2-i1,stepi(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type diag()
        {
            TMVAssert(!isunit());
            return vec_type(
                itsm.get(),size(),stepi()+stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i)
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            return vec_type(
                itsm.get()+i*stepj(),size()-i,stepi()+stepj(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i, int j1, int j2)
        {
            TMVAssert(isunit() ? i>0 : i>=0);
            TMVAssert(i<=size()); 
            if (I == FortranStyle) {
                TMVAssert(j1 > 0 && j1<=j2 && j2<=size()-i); --j1; 
            } else {
                TMVAssert(j1>=0 && j1<=j2 && j2<=size()-i); 
            }
            const int ds = stepi()+stepj();
            return vec_type(
                itsm.get()+i*stepj()+j1*ds,j2-j1,ds,NonConj TMV_FIRSTLAST);
        }

        //
        // Modifying Functions
        //

        inline type& setZero() 
        { fill_n(itsm.get(),itslen,T(0)); return *this; }

        inline type& setAllTo(const T& x) 
        { VectorViewOf(itsm.get(),itslen).setAllTo(x); return *this; }

        inline type& addToAll(const T& x) 
        { 
            TMVAssert(!isunit());
            VectorViewOf(itsm.get(),itslen).addToAll(x);
            return *this; 
        }

        inline type& clip(RT thresh)
        { VectorViewOf(itsm.get(),itslen).clip(thresh); return *this; }

        inline type& conjugateSelf() 
        { VectorViewOf(itsm.get(),itslen).conjugateSelf(); return *this; }

        inline type& invertSelf()
        { view().invertSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1)) 
        {
            TMVAssert(!isunit() || x==T(1));
            setZero(); if (!isunit()) diag().setAllTo(x);
            return *this;
        }

        //
        // subMatrix
        //

        inline const_rec_type cSubMatrix(int i1, int i2, int j1, int j2) const
        {
            return const_rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1,stepi(),stepj(),S,NonConj);
        }

        inline const_rec_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            return const_rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor, NonConj);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            return const_vec_type(
                itsm.get()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),NonConj);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I==FortranStyle) { --i; --j; }
            return cSubVector(i,j,istep,jstep,size);
        }

        inline const_uppertri_type cSubTriMatrix(int i1, int i2) const
        {
            return const_uppertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),D,S,NonConj);
        }

        inline const_uppertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return cSubTriMatrix(i1,i2);
        }

        inline const_uppertri_type cSubTriMatrix(
            int i1, int i2, int istep) const
        {
            return const_uppertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep, istep*stepi(),istep*stepj(), D,
                istep==1 ? S : NoMajor, NonConj);
        }

        inline const_uppertri_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return cSubTriMatrix(i1,i2,istep);
        }

        inline const_uppertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return const_uppertri_type(
                itsm.get()+noff*stepj(),
                size()-noff,stepi(),stepj(),NonUnitDiag,S,NonConj);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<RT*>(itsm.get()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                D, isReal(T()) ? S : NoMajor, NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return const_realpart_type(
                reinterpret_cast<RT*>(itsm.get())+1, size(),
                2*stepi(), 2*stepj(), NonUnitDiag, NoMajor, NonConj);
        }

        inline rec_type cSubMatrix(int i1, int i2, int j1, int j2) 
        {
            return rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(),stepj(),S,NonConj TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) 
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            StorageType newstor = S == RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            return rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,NonConj TMV_FIRSTLAST);
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(int i, int j, int istep, int jstep, int size)
        {
            return vec_type(
                itsm.get()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type subVector(int i, int j, int istep, int jstep, int size)
        {
            TMVAssert(size >= 0);
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I==FortranStyle) { --i; --j; }
            return cSubVector(i,j,istep,jstep,size);
        }

        inline uppertri_type cSubTriMatrix(int i1, int i2)
        {
            return uppertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),D,S,NonConj TMV_FIRSTLAST);
        }

        inline uppertri_type subTriMatrix(int i1, int i2)
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return cSubTriMatrix(i1,i2);
        }

        inline uppertri_type cSubTriMatrix(int i1, int i2, int istep) 
        {
            return uppertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),D,
                istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline uppertri_type subTriMatrix(int i1, int i2, int istep) 
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return cSubTriMatrix(i1,i2,istep);
        }

        inline uppertri_type offDiag(int noff=1)
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return uppertri_type(
                itsm.get()+noff*stepj(),size()-noff,
                stepi(),stepj(),NonUnitDiag,S,NonConj TMV_FIRSTLAST);
        }

        inline realpart_type realPart() 
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm.get()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                D, isReal(T()) ? S : NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return realpart_type(
                reinterpret_cast<RT*>(itsm.get())+1, size(),
                2*stepi(), 2*stepj(), NonUnitDiag, NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline const_uppertri_type view() const
        { 
            return const_uppertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,NonConj);
        }

        inline const_uppertri_type viewAsUnitDiag() const
        { 
            return const_uppertri_type(
                itsm.get(),size(),stepi(),stepj(),UnitDiag,S,NonConj);
        }

        inline const_lowertri_type transpose() const
        { 
            return const_lowertri_type(
                itsm.get(),size(),stepj(),stepi(),D,TMV_TransOf(S),NonConj);
        }

        inline const_uppertri_type conjugate() const
        { 
            return const_uppertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,TMV_ConjOf(T,NonConj));
        }

        inline const_lowertri_type adjoint() const
        { 
            return const_lowertri_type(
                itsm.get(),size(),stepj(),stepi(),D,
                TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline uppertri_type view() 
        { 
            return uppertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,NonConj TMV_FIRSTLAST);
        }

        inline uppertri_type viewAsUnitDiag() 
        { 
            return uppertri_type(
                itsm.get(),size(),stepi(),stepj(),UnitDiag,S,NonConj 
                TMV_FIRSTLAST);
        }

        inline lowertri_type transpose() 
        { 
            return lowertri_type(
                itsm.get(),size(),stepj(),stepi(),D,TMV_TransOf(S),NonConj 
                TMV_FIRSTLAST);
        }

        inline uppertri_type conjugate() 
        { 
            return uppertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,TMV_ConjOf(T,NonConj) 
                TMV_FIRSTLAST);
        }

        inline lowertri_type adjoint() 
        { 
            return lowertri_type(
                itsm.get(),size(),stepj(),stepi(),D,
                TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        //
        // I/O
        //
        
        void read(const TMV_Reader& reader);

        inline int size() const { return itss; }
        inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        inline int stepi() const { return S==RowMajor ? itss : 1; }
        inline int stepj() const { return S==RowMajor ? 1 : itss; }
        inline DiagType dt() const { return D; }
        inline StorageType stor() const { return S; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isunit() const { return D == UnitDiag; }
        inline bool isconj() const { return false; }

        inline reference ref(int i, int j)
        {
            return TriRefHelper2<T,D>::makeRef(
                i==j,
                itsm.get()[S==RowMajor ? i*itss + j : j*itss + i]); 
        }

        inline T cref(int i, int j) const 
        { 
            return (
                (isunit() && i==j) ? T(1) :
                (i>j) ? T(0) :
                itsm.get()[S==RowMajor ? i*itss + j : j*itss + i]); 
        }

        inline void resize(int s)
        {
            TMVAssert(s >= 0);
            itslen = s*s;
            itsm.resize(itslen);
            itss = s;
#ifdef TMVFLDEBUG
            _first = itsm.get();
            _last = _first+itslen;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline rowmajor_iterator rowmajor_begin() 
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end()
        { return rowmajor_iterator(this,size(),size()); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        { return const_rowmajor_iterator(this,size(),size()); }

        inline colmajor_iterator colmajor_begin() 
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end() 
        { return colmajor_iterator(this,0,size()); }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        { return const_colmajor_iterator(this,0,size()); }

    protected :

        int itslen;
        AlignedArray<T> itsm;
        int itss;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
    protected :
#endif

        inline bool okij(int i, int j) const
        {
            TMVAssert(i>=0 && i < size());
            TMVAssert(j>=0 && j < size());
            if (isunit()) return i<j; else return i<=j;
        }

        template <IndexStyle I2>
        friend void Swap(
            UpperTriMatrix<T,D,S,I>& m1, UpperTriMatrix<T,D,S,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            m1.itsm.swapWith(m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // UpperTriMatrix

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    class LowerTriMatrix : public GenLowerTriMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenLowerTriMatrix<T> base;
        typedef LowerTriMatrix<T,D,S,I> type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstMatrixView<T,I> const_rec_type;
        typedef ConstUpperTriMatrixView<T,I> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,I> const_lowertri_type;
        typedef const_lowertri_type const_view_type;
        typedef const_uppertri_type const_transpose_type;
        typedef const_lowertri_type const_conjugate_type;
        typedef const_uppertri_type const_adjoint_type;
        typedef ConstLowerTriMatrixView<RT,I> const_realpart_type;
        typedef VectorView<T,I> vec_type;
        typedef MatrixView<T,I> rec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef lowertri_type view_type;
        typedef uppertri_type transpose_type;
        typedef lowertri_type conjugate_type;
        typedef uppertri_type adjoint_type;
        typedef LowerTriMatrixView<RT,I> realpart_type;
        typedef typename TriRefHelper2<T,D>::reference reference;
        typedef RMIt<type> rowmajor_iterator;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;

        //
        // Constructors
        //

#define NEW_SIZE(s) \
        itslen((s)*(s)), itsm(itslen), itss(s) \
        TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

        explicit inline LowerTriMatrix(int _size=0) : NEW_SIZE(_size) 
        {
            TMVAssert(_size >= 0);
            TMVAssert(S==RowMajor || S==ColMajor); 
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline LowerTriMatrix(int _size, const T& x) : NEW_SIZE(_size)
        {
            TMVAssert(_size >= 0);
            TMVAssert(S==RowMajor || S==ColMajor);
            setAllTo(x);
        }

        template <class T2> 
        inline LowerTriMatrix(const GenMatrix<T2>& rhs) :
            NEW_SIZE(rhs.colsize())
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor); 
            Copy(rhs.lowerTri(D).transpose(),transpose()); 
        }

        template <class T2> 
        inline LowerTriMatrix(const GenLowerTriMatrix<T2>& rhs) :
            NEW_SIZE(rhs.size())
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(S==RowMajor || S==ColMajor); 
            if (isunit() && !rhs.isunit()) {
                if (rhs.size() > 0) 
                    Copy(rhs.offDiag().transpose(),offDiag().transpose());
            } else {
                Copy(rhs.transpose(),transpose());
            }
        }

        inline LowerTriMatrix(const type& rhs) :
            itslen(rhs.itslen), itsm(itslen), itss(rhs.itss)
            TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
        }

        template <DiagType D2, IndexStyle I2> 
        inline LowerTriMatrix(const LowerTriMatrix<T,D2,S,I2>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
            if (D==NonUnitDiag && D2==UnitDiag) diag().setAllTo(T(1));
        }

        template <IndexStyle I2> 
        inline LowerTriMatrix(const Matrix<T,S,I2>& rhs) :
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            if (rhs.isSquare())
                std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
            else
                Copy(rhs.lowerTri(D).transpose(),transpose());
        }

        inline LowerTriMatrix(const GenMatrix<T>& rhs) :
            NEW_SIZE(rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            Copy(rhs.lowerTri(D).transpose(),transpose());
        }

        inline LowerTriMatrix(const GenLowerTriMatrix<RT>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            if (isunit() && !rhs.isunit()) {
                if (rhs.size() > 0) offDiag() = rhs.offDiag();
            } else 
                rhs.assignToL(view());
        }

        inline LowerTriMatrix(const GenLowerTriMatrix<CT>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()));
            if (isunit() && !rhs.isunit()) {
                if (rhs.size() > 0) offDiag() = rhs.offDiag();
            } else 
                rhs.assignToL(view());
        }

        inline LowerTriMatrix(const AssignableToLowerTriMatrix<RT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToL(view());
        }

        inline LowerTriMatrix(const AssignableToLowerTriMatrix<CT>& m2) :
            NEW_SIZE(m2.size())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToL(view());
        }

#undef NEW_SIZE

        virtual inline ~LowerTriMatrix() 
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

        template <IndexStyle I2> 
        inline type& operator=(const LowerTriMatrix<T,D,S,I2>& m2)
        { 
            TMVAssert(size() == m2.size());
            if (&m2 != this) 
                std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<RT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToL(view());
            return *this;
        }

        inline type& operator=(const GenLowerTriMatrix<CT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToL(view());
            return *this;
        }

        template <class T2> 
        inline type& operator=(const GenLowerTriMatrix<T2>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            Copy(m2,view());
            return *this;
        }

        inline type& operator=(const T& x) 
        { 
            TMVAssert(!this->isunit() || x==T(1));
            return setToIdentity(x); 
        }

        inline type& operator=(const AssignableToLowerTriMatrix<RT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            m2.assignToL(view());
            return *this;
        }

        inline type& operator=(const AssignableToLowerTriMatrix<CT>& m2)
        { 
            TMVAssert(size() == m2.size());
            TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
            TMVAssert(isComplex(T()));
            m2.assignToL(view());
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(rowmajor_begin(),size()*(size()+1)/2,x); }

        //
        // Access
        //

        inline T operator()(int i, int j) const
        {
            if (I == CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
            } else {
                TMVAssert(i>0 && i<=size()); --i;
                TMVAssert(j>0 && j<=size()); --j;
            }
            return cref(i,j);
        }

        inline reference operator()(int i, int j) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j>=0 && j<size());
                TMVAssert(i>=j);
            } else {
                TMVAssert(i>0 && i<= size()); --i;
                TMVAssert(j>0 && j<= size()); --j;
                TMVAssert(i>=j);
            }
            return ref(i,j);
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            if (I==FortranStyle) { 
                TMVAssert(i>0 && i<=size()); --i;
                TMVAssert(j1>0 && j1<=j2 && j2<=size()); --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j2-1));
            return const_vec_type(
                itsm.get()+i*stepi()+j1*stepj(),j2-j1,stepj(),NonConj);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=size()); --j;
                TMVAssert(i1>0 && i1<=i2 && i2<=size()); --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1<=i2 && i2<=size());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            return const_vec_type(
                itsm.get()+i1*stepi()+j*stepj(),i2-i1,stepi(),NonConj);
        }

        inline const_vec_type diag() const
        {
            TMVAssert(!isunit());
            return const_vec_type(itsm.get(),size(),stepi()+stepj(),NonConj); 
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            return const_vec_type(
                itsm.get()-i*stepi(),size()+i,stepi()+stepj(),NonConj);
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            if (I == FortranStyle) { 
                TMVAssert(j1>0 && j1<=j2 && j2<=size()+i); --j1;
            } else {
                TMVAssert(j1>=0 && j1<=j2 && j2<=size()+i);
            }
            const int ds = stepi()+stepj();
            return const_vec_type(
                itsm.get()-i*stepi()+j1*ds,j2-j1,ds,NonConj);
        }

        inline vec_type row(int i, int j1, int j2)
        { 
            if (I==FortranStyle) { 
                TMVAssert(i>0 && i<=size()); --i;
                TMVAssert(j1>0 && j1<=j2 && j2<=size()); --j1;
            } else {
                TMVAssert(i>=0 && i<size());
                TMVAssert(j1>=0 && j1<=j2 && j2<=size());
            }
            TMVAssert(j1==j2 || okij(i,j2-1));
            return vec_type(
                itsm.get()+i*stepi()+j1*stepj(),
                j2-j1,stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I==FortranStyle) { 
                TMVAssert(j>0 && j<=size()); --j;
                TMVAssert(i1>0 && i1<=i2 && i2<=size()); --i1;
            } else {
                TMVAssert(j>=0 && j<size());
                TMVAssert(i1>=0 && i1<=i2 && i2<=size());
            }
            TMVAssert(i1==i2 || okij(i1,j));
            return vec_type(
                itsm.get()+i1*stepi()+j*stepj(),i2-i1,stepi(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type diag()
        {
            TMVAssert(!isunit());
            return vec_type(
                itsm.get(),size(),stepi()+stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag(int i)
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            return vec_type(
                itsm.get()-i*stepi(),size()+i,stepi()+stepj(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i, int j1, int j2)
        {
            TMVAssert(i>=-size()); 
            TMVAssert(isunit() ? i<0 : i<=0);
            if (I == FortranStyle) { 
                TMVAssert(j1>0 && j1<=j2 && j2<=size()+i); --j1;
            } else {
                TMVAssert(j1>=0 && j1<=j2 && j2<=size()+i);
            }
            const int ds = stepi()+stepj();
            return vec_type(
                itsm.get()-i*stepi()+j1*ds,j2-j1,ds,NonConj TMV_FIRSTLAST);
        }

        //
        // Modifying Functions
        //

        inline type& setZero() 
        { fill_n(itsm.get(),itslen,T(0)); return *this; }

        inline type& setAllTo(const T& x) 
        { VectorViewOf(itsm.get(),itslen).setAllTo(x); return *this; }

        inline type& addToAll(const T& x) 
        {
            TMVAssert(!isunit());
            VectorViewOf(itsm.get(),itslen).addToAll(x);
            return *this; 
        }

        inline type& clip(RT thresh)
        { VectorViewOf(itsm.get(),itslen).clip(thresh); return *this; }

        inline type& conjugateSelf() 
        { VectorViewOf(itsm.get(),itslen).conjugateSelf(); return *this; }

        inline type& invertSelf()
        { view().invertSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1)) 
        { 
            TMVAssert(!isunit() || x == T(1));
            setZero(); if (!isunit()) diag().setAllTo(x);
            return *this;
        }

        //
        // subMatrix
        //

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        {
            return const_rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1,stepi(),stepj(),S,NonConj);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            return const_rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,NonConj);
        }

        inline const_rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            return const_vec_type(
                itsm.get()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),NonConj);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I==FortranStyle) { --i; --j; }
            return cSubVector(i,j,istep,jstep,size);
        }

        inline const_lowertri_type cSubTriMatrix(int i1, int i2) const
        {
            return const_lowertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),D,S,NonConj);
        }

        inline const_lowertri_type subTriMatrix(int i1, int i2) const
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return cSubTriMatrix(i1,i2);
        }

        inline const_lowertri_type cSubTriMatrix(
            int i1, int i2, int istep) const
        {
            return const_lowertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),D,
                istep==1 ? S : NoMajor, NonConj);
        }

        inline const_lowertri_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return cSubTriMatrix(i1,i2,istep);
        }

        inline const_lowertri_type offDiag(int noff=1) const
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return const_lowertri_type(
                itsm.get()+noff*stepi(),size()-noff,
                stepi(),stepj(),NonUnitDiag,S,NonConj);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<RT*>(itsm.get()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                D, isReal(T()) ? S : NoMajor, NonConj
            );
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return const_realpart_type(
                reinterpret_cast<RT*>(itsm.get())+1, size(),
                2*stepi(), 2*stepj(), NonUnitDiag, NoMajor, NonConj
            );
        }

        inline rec_type cSubMatrix(int i1, int i2, int j1, int j2) 
        {
            return rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                i2-i1, j2-j1, stepi(),stepj(),S,NonConj TMV_FIRSTLAST );
        }

        inline rec_type subMatrix(int i1, int i2, int j1, int j2) 
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline rec_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            StorageType newstor = S==RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            return rec_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,NonConj TMV_FIRSTLAST );
        }

        inline rec_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(
            int i, int j, int istep, int jstep, int size)
        {
            TMVAssert(size >= 0);
            return vec_type(
                itsm.get()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),NonConj TMV_FIRSTLAST );
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int size)
        {
            TMVAssert(size >= 0);
            TMVAssert(view().hasSubVector(i,j,istep,jstep,size));
            if (I==FortranStyle) { --i; --j; }
            return cSubVector(i,j,istep,jstep,size);
        }

        inline lowertri_type cSubTriMatrix(int i1, int i2)
        {
            return lowertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                i2-i1,stepi(),stepj(),D,S,NonConj TMV_FIRSTLAST);
        }

        inline lowertri_type subTriMatrix(int i1, int i2)
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,1));
            if (I==FortranStyle) { --i1; }
            return cSubTriMatrix(i1,i2);
        }

        inline lowertri_type cSubTriMatrix(int i1, int i2, int istep) 
        {
            return lowertri_type(
                itsm.get()+i1*(stepi()+stepj()),
                (i2-i1)/istep,istep*stepi(),istep*stepj(),D,
                istep==1 ? S : NoMajor,NonConj TMV_FIRSTLAST);
        }

        inline lowertri_type subTriMatrix(int i1, int i2, int istep) 
        {
            TMVAssert(view().hasSubTriMatrix(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2+=istep-1; }
            return cSubTriMatrix(i1,i2,istep);
        }

        inline lowertri_type offDiag(int noff=1)
        {
            TMVAssert(noff >= 1);
            TMVAssert(noff <= size());
            return lowertri_type(
                itsm.get()+noff*stepi(),size()-noff,
                stepi(),stepj(),NonUnitDiag,S,NonConj TMV_FIRSTLAST);
        }

        inline realpart_type realPart() 
        {
            return realpart_type(
                reinterpret_cast<RT*>(itsm.get()), size(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                D, isReal(T()) ? S : NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        {
            TMVAssert(isComplex(T()));
            TMVAssert(!isunit());
            return realpart_type(
                reinterpret_cast<RT*>(itsm.get())+1,
                size(), 2*stepi(), 2*stepj(), NonUnitDiag, NoMajor, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline const_lowertri_type view() const
        { 
            return const_lowertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,NonConj);
        }

        inline const_lowertri_type viewAsUnitDiag() const
        { 
            return const_lowertri_type(
                itsm.get(),size(),stepi(),stepj(),UnitDiag,S,NonConj);
        }

        inline const_uppertri_type transpose() const
        { 
            return const_uppertri_type(
                itsm.get(),size(),stepj(),stepi(),D,TMV_TransOf(S),NonConj);
        }

        inline const_lowertri_type conjugate() const
        { 
            return const_lowertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,TMV_ConjOf(T,NonConj));
        }

        inline const_uppertri_type adjoint() const
        { 
            return const_uppertri_type(
                itsm.get(),size(),stepj(),stepi(),D,
                TMV_TransOf(S),TMV_ConjOf(T,NonConj));
        }

        inline lowertri_type view() 
        { 
            return lowertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,NonConj TMV_FIRSTLAST);
        }

        inline lowertri_type viewAsUnitDiag() 
        { 
            return lowertri_type(
                itsm.get(),size(),stepi(),stepj(),UnitDiag,S,NonConj 
                TMV_FIRSTLAST);
        }

        inline uppertri_type transpose() 
        { 
            return uppertri_type(
                itsm.get(),size(),stepj(),stepi(),D,
                TMV_TransOf(S),NonConj TMV_FIRSTLAST);
        }

        inline lowertri_type conjugate() 
        { 
            return lowertri_type(
                itsm.get(),size(),stepi(),stepj(),D,S,TMV_ConjOf(T,NonConj) 
                TMV_FIRSTLAST);
        }

        inline uppertri_type adjoint() 
        { 
            return uppertri_type(
                itsm.get(),size(),stepj(),stepi(),D,
                TMV_TransOf(S),TMV_ConjOf(T,NonConj) TMV_FIRSTLAST);
        }

        //
        // I/O
        //
        
        void read(const TMV_Reader& reader);

        inline int size() const { return itss; }
        inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        inline int stepi() const { return S==RowMajor ? itss : 1; }
        inline int stepj() const { return S==RowMajor ? 1 : itss; }
        inline DiagType dt() const { return D; }
        inline StorageType stor() const { return S; }
        inline ConjType ct() const { return NonConj; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isunit() const { return D == UnitDiag; }
        inline bool isconj() const { return false; }

        inline reference ref(int i, int j)
        { 
            return TriRefHelper2<T,D>::makeRef(
                i==j,
                itsm.get()[S==RowMajor ? i*itss + j : j*itss + i]); 
        }

        inline T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i<j) ? T(0) :
                itsm.get()[S==RowMajor ? i*itss + j : j*itss + i]); 
        }

        inline void resize(int s)
        {
            TMVAssert(s >= 0);
            itslen = s*s;
            itsm.resize(itslen);
            itss = s;
#ifdef TMVFLDEBUG
            _first = itsm.get();
            _last = _first+itslen;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline rowmajor_iterator rowmajor_begin() 
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end() 
        { return rowmajor_iterator(this,size(),0); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        { return const_rowmajor_iterator(this,size(),0); }

        inline colmajor_iterator colmajor_begin() 
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end() 
        { return colmajor_iterator(this,size(),size()); }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        { return const_colmajor_iterator(this,size(),size()); }

    protected :

        int itslen;
        AlignedArray<T> itsm;
        int itss;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected :
#endif

        inline bool okij(int i, int j) const
        {
            TMVAssert(i>=0 && i < size());
            TMVAssert(j>=0 && j < size());
            return isunit() ? i>j : i>=j;
        }

        template <IndexStyle I2>
        friend void Swap(
            LowerTriMatrix<T,D,S,I>& m1, LowerTriMatrix<T,D,S,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            m1.itsm.swapWith(m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // LowerTriMatrix

    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   UpperTriMatrixViewOf(T* m, n, S, D)
    //   UnitUpperTriMatrixViewOf(T* m, n, S)
    //   UpperTriMatrixViewOf(T* m, n, Si, Sj, D)
    //   UnitUpperTriMatrixViewOf(T* m, n, Si, Sj)
    //

    template <class T> 
    inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
        const T* m, int size, StorageType stor, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T>(
                m,size,size,1,dt,RowMajor,NonConj);
        else
            return ConstUpperTriMatrixView<T>(
                m,size,1,size,dt,ColMajor,NonConj);
    }

    template <class T> 
    inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
        T* m, int size, StorageType stor, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return UpperTriMatrixView<T>(
                m,size,size,1,dt,RowMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
        else
            return UpperTriMatrixView<T>(
                m,size,1,size,dt,ColMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
    }

    template <class T> 
    inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
        const T* m, int size, StorageType stor, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T>(
                m,size,size,1,dt,RowMajor,NonConj);
        else
            return ConstLowerTriMatrixView<T>(
                m,size,1,size,dt,ColMajor,NonConj);
    }

    template <class T> 
    inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
        T* m, int size, StorageType stor, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return LowerTriMatrixView<T>(
                m,size,size,1,dt,RowMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
        else
            return LowerTriMatrixView<T>(
                m,size,1,size,dt,ColMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
    }

    template <class T> 
    inline ConstUpperTriMatrixView<T> UnitUpperTriMatrixViewOf(
        const T* m, int size, StorageType stor)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T>(
                m,size,size,1,UnitDiag,RowMajor,NonConj);
        else
            return ConstUpperTriMatrixView<T>(
                m,size,1,size,UnitDiag,ColMajor,NonConj);
    }

    template <class T> 
    inline UpperTriMatrixView<T> UnitUpperTriMatrixViewOf(
        T* m, int size, StorageType stor)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return UpperTriMatrixView<T>(
                m,size,size,1,UnitDiag,RowMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
        else
            return UpperTriMatrixView<T>(
                m,size,1,size,UnitDiag,ColMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
    }

    template <class T> 
    inline ConstLowerTriMatrixView<T> UnitLowerTriMatrixViewOf(
        const T* m, int size, StorageType stor)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T>(
                m,size,size,1,UnitDiag,RowMajor,NonConj);
        else
            return ConstLowerTriMatrixView<T>(
                m,size,1,size,UnitDiag,ColMajor,NonConj);
    }

    template <class T> 
    inline LowerTriMatrixView<T> UnitLowerTriMatrixViewOf(
        T* m, int size, StorageType stor)
    {
        TMVAssert(size >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return LowerTriMatrixView<T>(
                m,size,size,1,UnitDiag,RowMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
        else
            return LowerTriMatrixView<T>(
                m,size,1,size,UnitDiag,ColMajor,NonConj 
                TMV_FIRSTLAST1(m,m+size*size));
    }


    template <class T> 
    inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
        const T* m, int size, int stepi, int stepj, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return ConstUpperTriMatrixView<T>(
            m,size,stepi,stepj,dt,stor,NonConj);
    }

    template <class T> 
    inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
        T* m, int size, int stepi, int stepj, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return UpperTriMatrixView<T>(
            m,size,stepi,stepj,dt,stor,NonConj
            TMV_FIRSTLAST1(m,m+size*size));
    }

    template <class T> 
    inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
        const T* m, int size, int stepi, int stepj, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return ConstLowerTriMatrixView<T>(
            m,size,stepi,stepj,dt,stor,NonConj);
    }

    template <class T> 
    inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
        T* m, int size, int stepi, int stepj, DiagType dt=NonUnitDiag)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return LowerTriMatrixView<T>(
            m,size,stepi,stepj,dt,stor,NonConj 
            TMV_FIRSTLAST1(m,m+size*size));
    }

    template <class T> 
    inline ConstUpperTriMatrixView<T> UnitUpperTriMatrixViewOf(
        const T* m, int size, int stepi, int stepj)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return ConstUpperTriMatrixView<T>(
            m,size,stepi,stepj,UnitDiag,stor,NonConj);
    }

    template <class T> 
    inline UpperTriMatrixView<T> UnitUpperTriMatrixViewOf(
        T* m, int size, int stepi, int stepj)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return UpperTriMatrixView<T>(
            m,size,size,1,UnitDiag,stor,NonConj 
            TMV_FIRSTLAST1(m,m+size*size));
    }

    template <class T> 
    inline ConstLowerTriMatrixView<T> UnitLowerTriMatrixViewOf(
        const T* m, int size, int stepi, int stepj)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return ConstLowerTriMatrixView<T>(
            m,size,size,1,UnitDiag,stor,NonConj);
    }

    template <class T> 
    inline LowerTriMatrixView<T> UnitLowerTriMatrixViewOf(
        T* m, int size, int stepi, int stepj)
    {
        TMVAssert(size >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor : stepj==1 ? RowMajor : NoMajor );
        return LowerTriMatrixView<T>(
            m,size,size,1,UnitDiag,stor,NonConj 
            TMV_FIRSTLAST1(m,m+size*size));
    }

    //
    // Copy
    //

    template <class T1, class T2> 
    inline void nonUnitDiagCopy(
        const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2)
    {
        TMVAssert(isReal(T1()) || isComplex(T2()));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.dt() == NonUnitDiag);
        TMVAssert(m2.dt() == NonUnitDiag);
        const int N = m1.size();

        if (!m1.isSameAs(m2) && m1.size() > 0) {
            if (m1.iscm() && m2.iscm()) 
                for(int j=0;j<N;++j) m2.col(j,0,j+1) = m1.col(j,0,j+1);
            else 
                for(int i=0;i<N;++i) m2.row(i,i,N) = m1.row(i,i,N);
        }
    }

    template <class T1, class T2> 
    inline void Copy(
        const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2)
    {
        TMVAssert(isReal(T1()) || isComplex(T2()));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());

        if (m1.isunit()) {
            if (m1.size() > 0)
                nonUnitDiagCopy(m1.offDiag(),m2.offDiag());
            if (!m2.isunit())
                m2.diag().setAllTo(T2(1));
        } else {
            nonUnitDiagCopy(m1,m2);
        }
    }


    //
    // Swap Matrices
    //

    template <class T> 
    void Swap(
        const UpperTriMatrixView<T>& m1, const UpperTriMatrixView<T>& m2);

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline void Swap(
        const UpperTriMatrixView<T>& m1, UpperTriMatrix<T,D,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline void Swap(
        UpperTriMatrix<T,D,S,I>& m1, const UpperTriMatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, DiagType D, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Swap(
        UpperTriMatrix<T,D,S1,I1>& m1, UpperTriMatrix<T,D,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }

    template <class T> 
    inline void Swap(
        const LowerTriMatrixView<T>& m1, const LowerTriMatrixView<T>& m2)
    { Swap(m1.transpose(),m2.transpose()); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline void Swap(
        const LowerTriMatrixView<T>& m1, LowerTriMatrix<T,D,S,I>& m2)
    { Swap(m1.transpose(),m2.transpose()); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline void Swap(
        LowerTriMatrix<T,D,S,I>& m1, const LowerTriMatrixView<T>& m2)
    { Swap(m1.transpose(),m2.transpose()); }

    template <class T, DiagType D, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Swap(
        LowerTriMatrix<T,D,S1,I1>& m1, LowerTriMatrix<T,D,S2,I2>& m2)
    { Swap(m1.transpose(),m2.transpose()); }


    //
    // Views:
    //

    template <class T> 
    inline ConstLowerTriMatrixView<T> Transpose(const GenUpperTriMatrix<T>& m)
    { return m.transpose(); }
    template <class T> 
    inline ConstUpperTriMatrixView<T> Transpose(const GenLowerTriMatrix<T>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Transpose(
        const ConstUpperTriMatrixView<T,I>& m)
    { return m.transpose(); }
    template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Transpose(
        const ConstLowerTriMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Transpose(
        const UpperTriMatrix<T,D,S,I>& m)
    { return m.transpose(); }
    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Transpose(
        const LowerTriMatrix<T,D,S,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Transpose(const UpperTriMatrixView<T,I>& m)
    { return m.transpose(); }
    template <class T, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Transpose(const LowerTriMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Transpose(UpperTriMatrix<T,D,S,I>& m)
    { return m.transpose(); }
    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Transpose(LowerTriMatrix<T,D,S,I>& m)
    { return m.transpose(); }

    template <class T> 
    inline ConstUpperTriMatrixView<T> Conjugate(const GenUpperTriMatrix<T>& m)
    { return m.conjugate(); }
    template <class T> 
    inline ConstLowerTriMatrixView<T> Conjugate(const GenLowerTriMatrix<T>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Conjugate(
        const ConstUpperTriMatrixView<T,I>& m)
    { return m.conjugate(); }
    template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Conjugate(
        const ConstLowerTriMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Conjugate(
        const UpperTriMatrix<T,D,S,I>& m)
    { return m.conjugate(); }
    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Conjugate(
        const LowerTriMatrix<T,D,S,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Conjugate(const UpperTriMatrixView<T,I>& m)
    { return m.conjugate(); }
    template <class T, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Conjugate(const LowerTriMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Conjugate(UpperTriMatrix<T,D,S,I>& m)
    { return m.conjugate(); }
    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Conjugate(LowerTriMatrix<T,D,S,I>& m)
    { return m.conjugate(); }

    template <class T> 
    inline ConstLowerTriMatrixView<T> Adjoint(const GenUpperTriMatrix<T>& m)
    { return m.adjoint(); }
    template <class T> 
    inline ConstUpperTriMatrixView<T> Adjoint(const GenLowerTriMatrix<T>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Adjoint(
        const ConstUpperTriMatrixView<T,I>& m)
    { return m.adjoint(); }
    template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Adjoint(
        const ConstLowerTriMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Adjoint(
        const UpperTriMatrix<T,D,S,I>& m)
    { return m.adjoint(); }
    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Adjoint(
        const LowerTriMatrix<T,D,S,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Adjoint(const UpperTriMatrixView<T,I>& m)
    { return m.adjoint(); }
    template <class T, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Adjoint(const LowerTriMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Adjoint(UpperTriMatrix<T,D,S,I>& m)
    { return m.adjoint(); }
    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Adjoint(LowerTriMatrix<T,D,S,I>& m)
    { return m.adjoint(); }

    template <class T> 
    inline QuotXU<T,T> Inverse(const GenUpperTriMatrix<T>& m)
    { return m.inverse(); }
    template <class T> 
    inline QuotXL<T,T> Inverse(const GenLowerTriMatrix<T>& m)
    { return m.inverse(); }


    //
    // TriMatrix ==, != TriMatrix
    //

    template <class T1, class T2> 
    bool operator==(
        const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2);

    template <class T1, class T2> 
    inline bool operator==(
        const GenLowerTriMatrix<T1>& m1, const GenLowerTriMatrix<T2>& m2)
    { return m1.transpose() == m2.transpose(); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenLowerTriMatrix<T1>& m1, const GenLowerTriMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator==(
        const GenUpperTriMatrix<T1>& m1, const GenMatrix<T2>& m2)
    {
        return 
            m1 == m2.upperTri() &&
            m2.lowerTri().offDiag().maxAbs2Element() == T2(0);
    }

    template <class T1, class T2> 
    inline bool operator==(
        const GenLowerTriMatrix<T1>& m1, const GenMatrix<T2>& m2)
    {
        return 
            m1 == m2.lowerTri() &&
            m2.upperTri().offDiag().maxAbs2Element() == T2(0);
    }

    template <class T1, class T2> 
    inline bool operator==(
        const GenMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
    { return m2 == m1; }

    template <class T1, class T2> 
    inline bool operator==(
        const GenMatrix<T1>& m1, const GenLowerTriMatrix<T2>& m2)
    { return m2 == m1; }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenUpperTriMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenLowerTriMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenLowerTriMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T>
    std::istream& operator>>(std::istream& is, const UpperTriMatrixView<T>& m)
    { return is >> IOStyle() >> m; }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    std::istream& operator>>(std::istream& is, UpperTriMatrix<T,D,S,I>& m)
    { return is >> IOStyle() >> m; }

    template <class T>
    std::istream& operator>>(std::istream& is, const LowerTriMatrixView<T>& m)
    { return is >> IOStyle() >> m; }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    std::istream& operator>>(std::istream& is, LowerTriMatrix<T,D,S,I>& m)
    { return is >> IOStyle() >> m; }

    template <class T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, const UpperTriMatrixView<T>& m)
    { m.read(reader); return reader.getis(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(
        const TMV_Reader& reader, UpperTriMatrix<T,D,S,I>& m)
    { m.read(reader); return reader.getis(); }

    template <class T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, const LowerTriMatrixView<T>& m)
    { m.read(reader); return reader.getis(); }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(
        const TMV_Reader& reader, LowerTriMatrix<T,D,S,I>& m)
    { m.read(reader); return reader.getis(); }


    template <class T, bool C>
    inline std::string TMV_Text(const TriRef<T,C>& ref)
    { 
        return std::string("TriRef<") + TMV_Text(T()) + "," +
            TMV_Text(C ? Conj : NonConj) + ">";
    }

} // namespace tmv

#endif
