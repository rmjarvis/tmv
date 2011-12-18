

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
//    In addition to the type template parameter (T), TriMatrixes have an
//    attributes parameter that accepts the following values:
//        NonUnitDiag || UnitDiag 
//        ColMajor || RowMajor
//        CStyle || FortranStyle
//
//    The defaults are NonUnitDiag, ColMajor, CStyle, and you only need
//    to specify the ones you want to change from the default.
//    So UpperTriMatrix<T,RowMajor> means NonUnitDiag, RowMajor, CStyle.
//    Multiple parameters are combined with the bitwise | operator.
//    e.g. UpperTriMatrix<T,UnitDiag | FortranStyle>.
//
//    If the matrix is UnitDiag, then the diagonal elements are not
//    actually stored or referenced.  The are all taken to be = 1.
//
//    The storage options follow the same meaning as for regular Matrices.
//
//    TriMatrix<T,A>(size_t n)
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    TriMatrix<T,A>(size_t n, T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    TriMatrix<T,A>(const Matrix<T>& m)
//    TriMatrix<T,A>(const TriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//
// Special Creators:
//
//    ConstUpperTriMatrixView UpperTriMatrixViewOf(
//            const T* m, size_t size, StorageType stor)
//    ConstUpperTriMatrixView UnitUpperTriMatrixViewOf(
//            const T* m, size_t size, StorageType stor)
//    UpperTriMatrixView UpperTriMatrixViewOf(
//            T* m, size_t size, StorageType stor)
//    UpperTriMatrixView UnitUpperTriMatrixViewOf(
//            T* m, size_t size, StorageType stor)
//        Returns a TriMatrixView of the elements in m, using the 
//        actual elements m for the storage.  The Unit versions return
//        views with dt = UnitDiag, the non-Unit versions return views
//        with dt = NonUnitDiag.
//        There are also corresponding LowerTriMatrix versions of these.
//        And there are versions with a final dt parameter that 
//        returns an UnknownDiag TriMatrixView with the value set by dt.
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
//    real_type m.maxAbsElement() const    or MaxAbsElement(m) 
//    real_type m.maxAbs2Element() const    or MaxAbs2Element(m) 
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
//    unknowndiag_type viewAsUnknownDiag(dt)
//        Re-view a UnitDiag or NonUnitDiag TriMatrix as UnknownDiag.
//        If dt is omitted, just use the current value.
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
//    nonconst_type nonConst()
//    noalias_type noAlias()
//    alias_type alias()
//    const_view_type constView()
//        Just like the regular Matrix versions
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
#include "TMV_VIt.h"
#include "TMV_MIt.h"
#include "TMV_Array.h"

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
        explicit TriRef(bool _u, T& _v) : isunit(_u), itsref(_v) {}
        TriRef(const TriRef<T,C>& rhs) : 
            isunit(rhs.isunit), itsref(rhs.itsref) {}
        ~TriRef() {}

        TMV_INLINE operator T() const { return val(); }
        TMV_INLINE reference getRef() { return ref(); }
        T operator-() const { return -val(); }

        template <class T2>
        TriRef<T,C>& operator=(const TriRef<T2,C>& rhs)
        { assign(rhs.val()); return *this; }
        template <class T2>
        TriRef<T,C>& operator=(T2 rhs)
        { assign(rhs); return *this; }

        template <class T2>
        TriRef<T,C>& operator+=(const TriRef<T2,C>& x2)
        { assign(val() + x2.val()); return *this; }
        template <class T2>
        TriRef<T,C>& operator+=(T2 x2)
        { assign(val() + x2); return *this; }
        template <class T2>
        typename Traits2<T,T2>::type operator+(const TriRef<T2,C>& x2)
        { return val() + x2.val(); }
        template <class T2>
        friend typename Traits2<T,T2>::type operator+(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()+x2; }
        template <class T2>
        friend typename Traits2<T,T2>::type operator+(
            T2 x1, const TriRef<T,C>& x2)
        { return x1+x2.val(); }
        template <class T2>
        friend T2& operator+=(T2& x1, const TriRef<T,C>& x2)
        { return x1 += x2.val(); }

        template <class T2>
        TriRef<T,C>& operator-=(const TriRef<T2,C>& x2) 
        { assign(val() - x2.val()); return *this; }
        template <class T2>
        TriRef<T,C>& operator-=(T2 x2)
        { assign(val() - x2); return *this; }
        template <class T2>
        typename Traits2<T,T2>::type operator-(const TriRef<T2,C>& x2)
        { return val()-x2.val(); }
        template <class T2>
        friend typename Traits2<T,T2>::type operator-(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()-x2; }
        template <class T2>
        friend typename Traits2<T,T2>::type operator-(
            T2 x1, const TriRef<T,C>& x2)
        { return x1-x2.val(); }
        template <class T2>
        friend T2& operator-=(T2& x1, const TriRef<T,C>& x2)
        { return x1 -= x2.val(); }

        template <class T2>
        TriRef<T,C>& operator*=(const TriRef<T2,C>& x2) 
        { assign(x2.val() * val()); return *this; }
        template <class T2>
        TriRef<T,C>& operator*=(T2 x2) 
        { assign(x2 * val()); return *this; }
        template <class T2>
        typename Traits2<T,T2>::type operator*(const TriRef<T2,C> x2)
        { return val()*x2.val(); }
        template <class T2>
        friend typename Traits2<T,T2>::type operator*(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()*x2; }
        template <class T2>
        friend typename Traits2<T,T2>::type operator*(
            T2 x1, const TriRef<T,C>& x2)
        { return x1*x2.val(); }
        template <class T2>
        friend T2& operator*=(T2& x1, const TriRef<T,C>& x2)
        {
            if (x2.isunit) return x1;
            else return x1 *= x2.val(); 
        }

        template <class T2>
        TriRef<T,C>& operator/=(const TriRef<T2,C>& x2) 
        { assign(val() / x2.val()); return *this; }
        template <class T2>
        TriRef<T,C>& operator/=(T2 x2) 
        { assign(val() / x2); return *this; }
        template <class T2>
        typename Traits2<T,T2>::type operator/(const TriRef<T2,C>& x2)
        { return val()/x2.val(); }
        template <class T2>
        friend typename Traits2<T,T2>::type operator/(
            const TriRef<T,C>& x1, T2 x2)
        { return x1.val()/x2; }
        template <class T2>
        friend typename Traits2<T,T2>::type operator/(
            T2 x1, const TriRef<T,C>& x2)
        { return x1/x2.val(); }
        template <class T2>
        friend T2& operator/=(T2& x1, const TriRef<T,C>& x2)
        {
            if (x2.isunit) return x1;
            else return x1 /= x2.val(); 
        }

        template <class T2>
        bool operator==(const TriRef<T2,C>& x2) const
        { return val() == x2.val(); }
        template <class T2>
        bool operator==(T2 x2) const 
        { return val() == x2; }
        template <class T2>
        friend bool operator==(T2 x1, const TriRef<T,C>& x2)
        { return x1 == x2.val(); }
        template <class T2>
        bool operator!=(const TriRef<T2,C>& x2) const
        { return !(operator==(x2)); }
        template <class T2>
        bool operator!=(T2 x2) const 
        { return !(operator==(x2)); }
        template <class T2>
        friend bool operator!=(T2 x1, const TriRef<T,C>& x2)
        { return !(x2==x1); }

        friend std::ostream& operator<<(std::ostream& os, TriRef<T,C> x)
        { os << x.val(); return os; }
        friend std::istream& operator>>(std::istream& is, TriRef<T,C> x)
        { is >> x.ref(); return is; }

    private:

        template <class T2>
        TMV_INLINE_ND void check(T2 TMV_DEBUG_PARAM(x)) const
        {
            TMVAssert(
                (!isunit || x == T2(1)) && 
                "Trying to assign to the diagonal of a UnitDiag TriMatrix.");
        }
        TMV_INLINE_ND void check() const
        {
            TMVAssert(
                !isunit && 
                "Trying to assign to the diagonal of a UnitDiag TriMatrix.");
        }
        TMV_INLINE reference ref() { check(); return itsref; }
        template <class T2>
        TMV_INLINE void assign(T2 x) { check(x); if (!isunit) itsref = x; }
        TMV_INLINE T val() const { return isunit ? T(1) : T(itsref); }

        const bool isunit;
        reference itsref;
    };

    // Overload some functions to work with TriRef
    template <class T, bool C>
    TMV_INLINE T TMV_CONJ(const TriRef<T,C>& x) 
    { return TMV_CONJ(T(x)); }
    template <class T, bool C>
    TMV_INLINE typename Traits<T>::real_type TMV_NORM(const TriRef<T,C>& x) 
    { return TMV_NORM(T(x)); }
    template <class T, bool C>
    TMV_INLINE typename Traits<T>::real_type TMV_ABS(const TriRef<T,C>& x) 
    { return TMV_ABS(T(x)); }
    template <class T, bool C>
    TMV_INLINE T TMV_SQR(const TriRef<T,C>& x) 
    { return TMV_SQR(T(x)); }
    template <class T, bool C>
    TMV_INLINE T TMV_SQRT(const TriRef<T,C>& x) 
    { return TMV_SQRT(T(x)); }
    template <class T, bool C>
    TMV_INLINE typename Traits<T>::real_type TMV_REAL(const TriRef<T,C>& x) 
    { return TMV_REAL(T(x)); }
    template <class T, bool C>
    TMV_INLINE typename Traits<T>::real_type TMV_IMAG(const TriRef<T,C>& x) 
    { return TMV_IMAG(T(x)); }

    template <class T, bool C1, bool C2>
    TMV_INLINE void TMV_SWAP(TriRef<T,C1> x1, TriRef<T,C2> x2)
    { TMV_SWAP(x1.getRef(),x2.getRef()); }
    template <class T, bool C2>
    TMV_INLINE void TMV_SWAP(T& x1, TriRef<T,C2> x2)
    { TMV_SWAP(x1,x2.getRef()); }
    template <class T, bool C1>
    TMV_INLINE void TMV_SWAP(TriRef<T,C1> x1, T& x2)
    { TMV_SWAP(x1.getRef(),x2); }

    template <class T, bool C, bool nonunit>
    struct TriRefHelper // D = NonUnitDiag
    {
        typedef T& reference;
        static TMV_INLINE reference makeRef(bool, T& r) { return r; }
    };
    template <class T, bool C>
    struct TriRefHelper<T,C,false>
    {
        typedef TriRef<T,C> reference;
        static TMV_INLINE reference makeRef(bool isunit, T& r)
        { return reference(isunit,r); }
    };

    // A Helper class to make a TriMatrix from a BaseMatrix
    // If the BaseMatrix is assignable to the Tri then we do the 
    // assign.  But if not, then we allow the construction - we
    // just make the TriMatrix from the correct portion of the matrix.
    template <bool assignable_to_tri>
    struct TriCopy // assignable
    {
        template <class M1, class M2>
        static TMV_INLINE void copy(
            const BaseMatrix<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
        { m2.noAlias() = m1; }
    };
    template <>
    struct TriCopy<false>
    {
        template <class M1, class M2>
        static void copy(
            const BaseMatrix<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
        {
            const bool up = M2::_upper;
            if (!m2.isunit()) {
                typename M2::nonunitdiag_type::noalias_type m2nu = 
                    m2.viewAsNonUnitDiag().noAlias();
                Maybe<up>::uppertri(m1.calc()).assignTo(m2nu);
            } else if (m2.size() > 1) {
                typename M2::offdiag_type::noalias_type m2o = 
                    m2.offDiag().noAlias();
                Maybe<up>::uppertri(m1.calc()).offDiag().assignTo(m2o);
            }
        }
    };

    template <class T, int A0, int A1, int A2>
    struct Traits<UpperTriMatrix<T,A0,A1,A2> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A012 = A0 | A1 | A2 };
        enum { A = (A012 & ~NoDivider & ~CheckAlias) | (
                ( Attrib<A012>::rowmajor ? 0 : ColMajor ) |
                ( Attrib<A012>::unitdiag ? 0 : NonUnitDiag ) )};
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor) &&
                (Attrib<A>::rowmajor != int(Attrib<A>::colmajor)) &&
                !Attrib<A>::diagmajor &&
                (Attrib<A>::unitdiag || Attrib<A>::nonunitdiag) &&
                (Attrib<A>::unitdiag != int(Attrib<A>::nonunitdiag)) &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::checkalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef UpperTriMatrix<T,A0,A1,A2> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef UpperTriMatrix<T,A012> copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _nlo = 0 };
        enum { _nhi = TMV_UNKNOWN };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = false };
        enum { _dt = A & AllDiagType };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { baseA = (_fort ? FortranStyle : 0) };
        enum { colA = baseA | (_colmajor ? Unit : 0) };
        enum { rowA = baseA | (_rowmajor ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { ndA = (A & ~AllDiagType) };
        enum { ndnmA = (ndA & ~AllStorageType) };
        enum { colpairA = _colmajor ? cmA : int(nmA) };
        enum { rowpairA = _rowmajor ? rmA : int(nmA) };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonunitA = ndA | NonUnitDiag };
        enum { unitA = ndA | UnitDiag };
        enum { nonconjA = A };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { An = (A & ~NoAlias) };

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,baseA> const_diag_type;
        typedef ConstVectorView<T,baseA> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,A> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmA> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmA> const_submatrix_step_type;
        typedef ConstVectorView<T,baseA> const_subvector_type;

        typedef ConstUpperTriMatrixView<T,A> const_view_type;
        typedef ConstUpperTriMatrixView<T,cstyleA> const_cview_type;
        typedef ConstUpperTriMatrixView<T,fstyleA> const_fview_type;
        typedef ConstUpperTriMatrixView<T> const_xview_type;
        typedef ConstUpperTriMatrixView<T,cmA> const_cmview_type;
        typedef ConstUpperTriMatrixView<T,rmA> const_rmview_type;
        typedef ConstUpperTriMatrixView<T,conjA> const_conjugate_type;
        typedef ConstLowerTriMatrixView<T,trA> const_transpose_type;
        typedef ConstLowerTriMatrixView<T,adjA> const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,nonunitA> const_offdiag_type;
        typedef ConstUpperTriMatrixView<T,unitA> const_unitdiag_type;
        typedef ConstUpperTriMatrixView<T,nonunitA> const_nonunitdiag_type;
        typedef ConstUpperTriMatrixView<T,ndA> const_unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallUpperTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                ConstUpperTriMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstUpperTriMatrixView<T,nonconjA> const_nonconj_type;
        typedef UpperTriMatrixView<T,A> nonconst_type;

        typedef typename TriRefHelper<T,false,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,baseA> diag_type;
        typedef VectorView<T,baseA> diag_sub_type;

        typedef UpperTriMatrixView<T,A> subtrimatrix_type;
        typedef UpperTriMatrixView<T,nmA> subtrimatrix_step_type;
        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,ndnmA> submatrix_step_type;
        typedef VectorView<T,baseA> subvector_type;

        typedef UpperTriMatrixView<T,A> view_type;
        typedef UpperTriMatrixView<T,cstyleA> cview_type;
        typedef UpperTriMatrixView<T,fstyleA> fview_type;
        typedef UpperTriMatrixView<T> xview_type;
        typedef UpperTriMatrixView<T,cmA> cmview_type;
        typedef UpperTriMatrixView<T,rmA> rmview_type;
        typedef UpperTriMatrixView<T,conjA> conjugate_type;
        typedef LowerTriMatrixView<T,trA> transpose_type;
        typedef LowerTriMatrixView<T,adjA> adjoint_type;

        typedef UpperTriMatrixView<T,nonunitA> offdiag_type;
        typedef UpperTriMatrixView<T,unitA> unitdiag_type;
        typedef UpperTriMatrixView<T,nonunitA> nonunitdiag_type;
        typedef UpperTriMatrixView<T,ndA> unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                SmallUpperTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                UpperTriMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef UpperTriMatrixView<T,nonconjA> nonconj_type;
        typedef UpperTriMatrixView<T,An|NoAlias> noalias_type;
        typedef UpperTriMatrixView<T,An> alias_type;
    };

    template <class T, int A0, int A1, int A2>
    class UpperTriMatrix : 
        public BaseMatrix_Tri_Mutable<UpperTriMatrix<T,A0,A1,A2> >
    {
    public:
        enum { A = A0 | A1 | A2 };

        typedef UpperTriMatrix<T,A0,A1,A2> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename Traits<type>::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _dt = Traits<type>::_dt };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        explicit UpperTriMatrix(size_t n=0) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(Traits<type>::okA);
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::constr_value());
#endif
        }

        UpperTriMatrix(size_t n, T x) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(Traits<type>::okA);
            Maybe<_unit>::offdiag(*this).setAllTo(x);
        }

        UpperTriMatrix(const type& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            this->noAlias() = m2;
        }

        template <class M2>
        UpperTriMatrix(const BaseMatrix<M2>& m2) :
            itss(m2.rowsize()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert(M2::_calc || assignable);
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        UpperTriMatrix(const BaseMatrix_Tri<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M2::_upper);
            typename Traits<type>::noalias_type na = this->noAlias();
            Maybe<_unit && !M2::_unit>::unitview(m2).assignTo(na);
        }

        template <class M2>
        UpperTriMatrix(const BaseMatrix_Diag<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            this->setZero();
            this->diag().noAlias() = m2.calc().diag();
        }

        ~UpperTriMatrix()
        {
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::destr_value());
#endif
        }


        //
        // Op=
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
        { base_mut::operator=(x); return *this; }


        // 
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i>j) ? T(0) :
                itsm[_rowmajor ? i*stepi() + j : i + j*stepj()]);
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,false,_nonunit>::makeRef(
                isunit() && i==j,
                itsm[_rowmajor ? i*stepi() + j : i + j*stepj()] ); 
        }

        void swapWith(type& m2)
        {
            TMVAssert(m2.size() == size());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }

        void resize(const size_t s)
        {
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::destr_value());
#endif
            itss = s;
            itsm.resize(s*s);
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::constr_value());
#endif
        }

        TMV_INLINE size_t size() const { return itss; }
        int nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE int stepi() const { return _rowmajor ? itss : 1; }
        TMV_INLINE int stepj() const { return _rowmajor ? 1 : itss; }
        TMV_INLINE DiagType dt() const { return static_cast<DiagType>(_dt); }
        TMV_INLINE bool isunit() const { return _unit; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    protected :

        size_t itss;
        AlignedArray<T> itsm;

    }; // UpperTriMatrix

    template <class T, int A0>
    struct Traits<ConstUpperTriMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~CheckAlias) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::checkalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstUpperTriMatrixView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _nlo = 0 };
        enum { _nhi = TMV_UNKNOWN };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? TMV_UNKNOWN : A & AllDiagType };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { baseA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) )};
        enum { colA = baseA | (_colmajor ? Unit : 0) };
        enum { rowA = baseA | (_rowmajor ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { ndA = (A & ~AllDiagType) };
        enum { ndnmA = (ndA & ~AllStorageType) };
        enum { colpairA = _colmajor ? cmA : int(nmA) };
        enum { rowpairA = _rowmajor ? rmA : int(nmA) };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonunitA = ndA | NonUnitDiag };
        enum { unitA = ndA | UnitDiag };
        enum { nonconjA = A & ~Conj };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                (_unit ? UnitDiag : NonUnitDiag) )};
        typedef UpperTriMatrix<T,copyA> copy_type;

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,baseA> const_diag_type;
        typedef ConstVectorView<T,baseA> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,A> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmA> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmA> const_submatrix_step_type;
        typedef ConstVectorView<T,baseA> const_subvector_type;

        typedef ConstUpperTriMatrixView<T,A> const_view_type;
        typedef ConstUpperTriMatrixView<T,cstyleA> const_cview_type;
        typedef ConstUpperTriMatrixView<T,fstyleA> const_fview_type;
        typedef ConstUpperTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstUpperTriMatrixView<T,cmA> const_cmview_type;
        typedef ConstUpperTriMatrixView<T,rmA> const_rmview_type;
        typedef ConstUpperTriMatrixView<T,conjA> const_conjugate_type;
        typedef ConstLowerTriMatrixView<T,trA> const_transpose_type;
        typedef ConstLowerTriMatrixView<T,adjA> const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,nonunitA> const_offdiag_type;
        typedef ConstUpperTriMatrixView<T,unitA> const_unitdiag_type;
        typedef ConstUpperTriMatrixView<T,nonunitA> const_nonunitdiag_type;
        typedef ConstUpperTriMatrixView<T,ndA> const_unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallUpperTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                ConstUpperTriMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstUpperTriMatrixView<T,nonconjA> const_nonconj_type;
        typedef UpperTriMatrixView<T,A> nonconst_type;

        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
    };

    template <class T, int A> 
    class ConstUpperTriMatrixView :
        public BaseMatrix_Tri<ConstUpperTriMatrixView<T,A> >
    {
    public:
        typedef ConstUpperTriMatrixView<T,A> type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _dt = Traits<type>::_dt };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        ConstUpperTriMatrixView(
            const T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstUpperTriMatrixView(const T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        ConstUpperTriMatrixView(const T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        ConstUpperTriMatrixView(const T* m, size_t s) :
            itsm(m), itss(s), itssi(_stepi), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != TMV_UNKNOWN);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        ConstUpperTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstUpperTriMatrixView(
            const ConstUpperTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstUpperTriMatrixView(
            const UpperTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstUpperTriMatrixView(
            const ConstSmallUpperTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstUpperTriMatrixView(
            const SmallUpperTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstUpperTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }

    private :
        void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }

        T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i>j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        int nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
        TMV_INLINE bool isunit() const { return dt() == UnitDiag; }
        TMV_INLINE bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        TMV_INLINE bool iscm() const
        {
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

        const T* itsm;
        const size_t itss;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // ConstUpperTriMatrixView

    template <class T, int A0>
    struct Traits<UpperTriMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~CheckAlias) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::checkalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef UpperTriMatrixView<T,A0> type;
        typedef ConstUpperTriMatrixView<T,A> calc_type;
        typedef const type& eval_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _nlo = 0 };
        enum { _nhi = TMV_UNKNOWN };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? TMV_UNKNOWN : A & AllDiagType };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { baseA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) )};
        enum { colA = baseA | (_colmajor ? Unit : 0) };
        enum { rowA = baseA | (_rowmajor ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { ndA = (A & ~AllDiagType) };
        enum { ndnmA = (ndA & ~AllStorageType) };
        enum { colpairA = _colmajor ? cmA : int(nmA) };
        enum { rowpairA = _rowmajor ? rmA : int(nmA) };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonunitA = ndA | NonUnitDiag };
        enum { unitA = ndA | UnitDiag };
        enum { nonconjA = A & ~Conj };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { An = (A & ~NoAlias) };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                (_unit ? UnitDiag : NonUnitDiag) )};
        typedef UpperTriMatrix<T,copyA> copy_type;

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,baseA> const_diag_type;
        typedef ConstVectorView<T,baseA> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,A> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmA> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmA> const_submatrix_step_type;
        typedef ConstVectorView<T,baseA> const_subvector_type;

        typedef ConstUpperTriMatrixView<T,A> const_view_type;
        typedef ConstUpperTriMatrixView<T,cstyleA> const_cview_type;
        typedef ConstUpperTriMatrixView<T,fstyleA> const_fview_type;
        typedef ConstUpperTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstUpperTriMatrixView<T,cmA> const_cmview_type;
        typedef ConstUpperTriMatrixView<T,rmA> const_rmview_type;
        typedef ConstUpperTriMatrixView<T,conjA> const_conjugate_type;
        typedef ConstLowerTriMatrixView<T,trA> const_transpose_type;
        typedef ConstLowerTriMatrixView<T,adjA> const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,nonunitA> const_offdiag_type;
        typedef ConstUpperTriMatrixView<T,unitA> const_unitdiag_type;
        typedef ConstUpperTriMatrixView<T,nonunitA> const_nonunitdiag_type;
        typedef ConstUpperTriMatrixView<T,ndA> const_unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallUpperTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                ConstUpperTriMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstUpperTriMatrixView<T,nonconjA> const_nonconj_type;
        typedef UpperTriMatrixView<T,A> nonconst_type;

        typedef typename TriRefHelper<T,_conj,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,baseA> diag_type;
        typedef VectorView<T,baseA> diag_sub_type;

        typedef UpperTriMatrixView<T,A> subtrimatrix_type;
        typedef UpperTriMatrixView<T,nmA> subtrimatrix_step_type;
        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,ndnmA> submatrix_step_type;
        typedef VectorView<T,baseA> subvector_type;

        typedef UpperTriMatrixView<T,A> view_type;
        typedef UpperTriMatrixView<T,cstyleA> cview_type;
        typedef UpperTriMatrixView<T,fstyleA> fview_type;
        typedef UpperTriMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef UpperTriMatrixView<T,cmA> cmview_type;
        typedef UpperTriMatrixView<T,rmA> rmview_type;
        typedef UpperTriMatrixView<T,conjA> conjugate_type;
        typedef LowerTriMatrixView<T,trA> transpose_type;
        typedef LowerTriMatrixView<T,adjA> adjoint_type;

        typedef UpperTriMatrixView<T,nonunitA> offdiag_type;
        typedef UpperTriMatrixView<T,unitA> unitdiag_type;
        typedef UpperTriMatrixView<T,nonunitA> nonunitdiag_type;
        typedef UpperTriMatrixView<T,ndA> unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                SmallUpperTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                UpperTriMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef UpperTriMatrixView<T,nonconjA> nonconj_type;
        typedef UpperTriMatrixView<T,An|NoAlias> noalias_type;
        typedef UpperTriMatrixView<T,An> alias_type;
    };

    template <class T, int A>
    class UpperTriMatrixView :
        public BaseMatrix_Tri_Mutable<UpperTriMatrixView<T,A> >
    {
    public:
        typedef UpperTriMatrixView<T,A> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename Traits<type>::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _dt = Traits<type>::_dt };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        UpperTriMatrixView(T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        UpperTriMatrixView(T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        UpperTriMatrixView(T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        UpperTriMatrixView(T* m, size_t s) :
            itsm(m), itss(s), itssi(_stepi), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != TMV_UNKNOWN);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        UpperTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        UpperTriMatrixView(UpperTriMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int Si2, int Sj2, int A2>
        UpperTriMatrixView(SmallUpperTriMatrixView<T,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~UpperTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }


        //
        // Op = 
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i>j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,_conj,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi()+j*stepj()]); 
        }

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        int nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
        TMV_INLINE bool isunit() const { return itsdt == UnitDiag; }
        TMV_INLINE bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        TMV_INLINE bool iscm() const
        {
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

        T* itsm;
        const size_t itss;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // UpperTriMatrixView

    template <class T, int A0, int A1, int A2>
    struct Traits<LowerTriMatrix<T,A0,A1,A2> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A012 = A0 | A1 | A2 };
        enum { A = (A012 & ~NoDivider & ~CheckAlias) | (
                ( Attrib<A012>::rowmajor ? 0 : ColMajor ) |
                ( Attrib<A012>::unitdiag ? 0 : NonUnitDiag ) )};
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor) &&
                (Attrib<A>::rowmajor != int(Attrib<A>::colmajor)) &&
                !Attrib<A>::diagmajor &&
                (Attrib<A>::unitdiag || Attrib<A>::nonunitdiag) &&
                (Attrib<A>::unitdiag != int(Attrib<A>::nonunitdiag)) &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::checkalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef LowerTriMatrix<T,A0,A1,A2> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef LowerTriMatrix<T,A012> copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = 0 };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = false };
        enum { _dt = A & AllDiagType };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { baseA = (_fort ? FortranStyle : 0) };
        enum { colA = baseA | (_colmajor ? Unit : 0) };
        enum { rowA = baseA | (_rowmajor ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { ndA = (A & ~AllDiagType) };
        enum { ndnmA = (ndA & ~AllStorageType) };
        enum { colpairA = _colmajor ? cmA : int(nmA) };
        enum { rowpairA = _rowmajor ? rmA : int(nmA) };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonunitA = ndA | NonUnitDiag };
        enum { unitA = ndA | UnitDiag };
        enum { nonconjA = A };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { An = (A & ~NoAlias) };

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,baseA> const_diag_type;
        typedef ConstVectorView<T,baseA> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,A> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmA> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmA> const_submatrix_step_type;
        typedef ConstVectorView<T,baseA> const_subvector_type;

        typedef ConstLowerTriMatrixView<T,A> const_view_type;
        typedef ConstLowerTriMatrixView<T,cstyleA> const_cview_type;
        typedef ConstLowerTriMatrixView<T,fstyleA> const_fview_type;
        typedef ConstLowerTriMatrixView<T> const_xview_type;
        typedef ConstLowerTriMatrixView<T,cmA> const_cmview_type;
        typedef ConstLowerTriMatrixView<T,rmA> const_rmview_type;
        typedef ConstLowerTriMatrixView<T,conjA> const_conjugate_type;
        typedef ConstUpperTriMatrixView<T,trA> const_transpose_type;
        typedef ConstUpperTriMatrixView<T,adjA> const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,nonunitA> const_offdiag_type;
        typedef ConstLowerTriMatrixView<T,unitA> const_unitdiag_type;
        typedef ConstLowerTriMatrixView<T,nonunitA> const_nonunitdiag_type;
        typedef ConstLowerTriMatrixView<T,ndA> const_unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallLowerTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                ConstLowerTriMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstLowerTriMatrixView<T,nonconjA> const_nonconj_type;
        typedef LowerTriMatrixView<T,A> nonconst_type;

        typedef typename TriRefHelper<T,false,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,baseA> diag_type;
        typedef VectorView<T,baseA> diag_sub_type;

        typedef LowerTriMatrixView<T,A> subtrimatrix_type;
        typedef LowerTriMatrixView<T,nmA> subtrimatrix_step_type;
        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,ndnmA> submatrix_step_type;
        typedef VectorView<T,baseA> subvector_type;

        typedef LowerTriMatrixView<T,A> view_type;
        typedef LowerTriMatrixView<T,cstyleA> cview_type;
        typedef LowerTriMatrixView<T,fstyleA> fview_type;
        typedef LowerTriMatrixView<T> xview_type;
        typedef LowerTriMatrixView<T,cmA> cmview_type;
        typedef LowerTriMatrixView<T,rmA> rmview_type;
        typedef LowerTriMatrixView<T,conjA> conjugate_type;
        typedef UpperTriMatrixView<T,trA> transpose_type;
        typedef UpperTriMatrixView<T,adjA> adjoint_type;

        typedef LowerTriMatrixView<T,nonunitA> offdiag_type;
        typedef LowerTriMatrixView<T,unitA> unitdiag_type;
        typedef LowerTriMatrixView<T,nonunitA> nonunitdiag_type;
        typedef LowerTriMatrixView<T,ndA> unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                SmallLowerTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                LowerTriMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef LowerTriMatrixView<T,nonconjA> nonconj_type;
        typedef LowerTriMatrixView<T,An|NoAlias> noalias_type;
        typedef LowerTriMatrixView<T,An> alias_type;
    };

    template <class T, int A0, int A1, int A2>
    class LowerTriMatrix : 
        public BaseMatrix_Tri_Mutable<LowerTriMatrix<T,A0,A1,A2> >
    {
    public:
        typedef LowerTriMatrix<T,A0,A1,A2> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename Traits<type>::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _dt = Traits<type>::_dt };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        explicit LowerTriMatrix(size_t n=0) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(Traits<type>::okA);
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::constr_value());
#endif
        }

        LowerTriMatrix(size_t n, T x) : itss(n), itsm(n*n)
        {
            TMVStaticAssert(Traits<type>::okA);
            Maybe<_unit>::offdiag(*this).setAllTo(x);
        }

        LowerTriMatrix(const type& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            this->noAlias() = m2;
        }

        template <class M2>
        LowerTriMatrix(const BaseMatrix<M2>& m2) :
            itss(m2.colsize()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert(M2::_calc || assignable);
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        LowerTriMatrix(const BaseMatrix_Tri<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M2::_lower);
            typename Traits<type>::noalias_type na = this->noAlias();
            Maybe<_unit && !M2::_unit>::unitview(m2).assignTo(na);
        }

        template <class M2>
        LowerTriMatrix(const BaseMatrix_Diag<M2>& m2) :
            itss(m2.size()), itsm(itss*itss)
        {
            TMVStaticAssert(Traits<type>::okA);
            this->setZero();
            this->diag().noAlias() = m2.calc().diag();
        }

        ~LowerTriMatrix()
        {
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::destr_value());
#endif
        }


        //
        // Op=
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
        { base_mut::operator=(x); return *this; }


        // 
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i<j) ? T(0) :
                itsm[_rowmajor ? i*stepi() + j : i + j*stepj()]);
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,false,_nonunit>::makeRef(
                isunit() && i==j,
                itsm[_rowmajor ? i*stepi() + j : i + j*stepj()] ); 
        }

        void swapWith(type& m2)
        {
            TMVAssert(m2.size() == size());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }

        void resize(const size_t s)
        {
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::destr_value());
#endif
            itss = s;
            itsm.resize(s*s);
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::constr_value());
#endif
        }

        TMV_INLINE size_t size() const { return itss; }
        int nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE int stepi() const { return _rowmajor ? itss : 1; }
        TMV_INLINE int stepj() const { return _rowmajor ? 1 : itss; }
        TMV_INLINE DiagType dt() const { return static_cast<DiagType>(_dt); }
        TMV_INLINE bool isunit() const { return _unit; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    protected :

        size_t itss;
        AlignedArray<T> itsm;

    }; // LowerTriMatrix

    template <class T, int A0>
    struct Traits<ConstLowerTriMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~CheckAlias) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::checkalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstLowerTriMatrixView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = 0 };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? TMV_UNKNOWN : A & AllDiagType };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { baseA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) )};
        enum { colA = baseA | (_colmajor ? Unit : 0) };
        enum { rowA = baseA | (_rowmajor ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { ndA = (A & ~AllDiagType) };
        enum { ndnmA = (ndA & ~AllStorageType) };
        enum { colpairA = _colmajor ? cmA : int(nmA) };
        enum { rowpairA = _rowmajor ? rmA : int(nmA) };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonunitA = ndA | NonUnitDiag };
        enum { unitA = ndA | UnitDiag };
        enum { nonconjA = A & ~Conj };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                (_unit ? UnitDiag : NonUnitDiag) )};
        typedef LowerTriMatrix<T,copyA> copy_type;

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,baseA> const_diag_type;
        typedef ConstVectorView<T,baseA> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,A> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmA> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmA> const_submatrix_step_type;
        typedef ConstVectorView<T,baseA> const_subvector_type;

        typedef ConstLowerTriMatrixView<T,A> const_view_type;
        typedef ConstLowerTriMatrixView<T,cstyleA> const_cview_type;
        typedef ConstLowerTriMatrixView<T,fstyleA> const_fview_type;
        typedef ConstLowerTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstLowerTriMatrixView<T,cmA> const_cmview_type;
        typedef ConstLowerTriMatrixView<T,rmA> const_rmview_type;
        typedef ConstLowerTriMatrixView<T,conjA> const_conjugate_type;
        typedef ConstUpperTriMatrixView<T,trA> const_transpose_type;
        typedef ConstUpperTriMatrixView<T,adjA> const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,nonunitA> const_offdiag_type;
        typedef ConstLowerTriMatrixView<T,unitA> const_unitdiag_type;
        typedef ConstLowerTriMatrixView<T,nonunitA> const_nonunitdiag_type;
        typedef ConstLowerTriMatrixView<T,ndA> const_unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallLowerTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                ConstLowerTriMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstLowerTriMatrixView<T,nonconjA> const_nonconj_type;
        typedef LowerTriMatrixView<T,A> nonconst_type;

        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
    };

    template <class T, int A>
    class ConstLowerTriMatrixView :
        public BaseMatrix_Tri<ConstLowerTriMatrixView<T,A> >
    {
    public:
        typedef ConstLowerTriMatrixView<T,A> type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _dt = Traits<type>::_dt };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //
        
        ConstLowerTriMatrixView(
            const T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstLowerTriMatrixView(const T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        ConstLowerTriMatrixView(const T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        ConstLowerTriMatrixView(const T* m, size_t s) :
            itsm(m), itss(s), itssi(_stepi), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != TMV_UNKNOWN);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        ConstLowerTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstLowerTriMatrixView(
            const ConstLowerTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstLowerTriMatrixView(
            const LowerTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstLowerTriMatrixView(
            const ConstSmallLowerTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstLowerTriMatrixView(
            const SmallLowerTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstLowerTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }

    private :
        void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }

        T cref(int i, int j) const 
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i<j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        int nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
        TMV_INLINE bool isunit() const { return itsdt == UnitDiag; }
        TMV_INLINE bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        TMV_INLINE bool iscm() const
        {
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

        const T* itsm;
        const size_t itss;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // ConstLowerTriMatrixView

    template <class T, int A0>
    struct Traits<LowerTriMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~CheckAlias) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::checkalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef LowerTriMatrixView<T,A0> type;
        typedef ConstLowerTriMatrixView<T,A> calc_type;
        typedef const type& eval_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = 0 };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? TMV_UNKNOWN : A & AllDiagType };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { baseA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) )};
        enum { colA = baseA | (_colmajor ? Unit : 0) };
        enum { rowA = baseA | (_rowmajor ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { ndA = (A & ~AllDiagType) };
        enum { ndnmA = (ndA & ~AllStorageType) };
        enum { colpairA = _colmajor ? cmA : int(nmA) };
        enum { rowpairA = _rowmajor ? rmA : int(nmA) };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonunitA = ndA | NonUnitDiag };
        enum { unitA = ndA | UnitDiag };
        enum { nonconjA = A & ~Conj };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { An = (A & ~NoAlias) };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                (_unit ? UnitDiag : NonUnitDiag) )};
        typedef LowerTriMatrix<T,copyA> copy_type;

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,baseA> const_diag_type;
        typedef ConstVectorView<T,baseA> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,A> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmA> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmA> const_submatrix_step_type;
        typedef ConstVectorView<T,baseA> const_subvector_type;

        typedef ConstLowerTriMatrixView<T,A> const_view_type;
        typedef ConstLowerTriMatrixView<T,cstyleA> const_cview_type;
        typedef ConstLowerTriMatrixView<T,fstyleA> const_fview_type;
        typedef ConstLowerTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstLowerTriMatrixView<T,cmA> const_cmview_type;
        typedef ConstLowerTriMatrixView<T,rmA> const_rmview_type;
        typedef ConstLowerTriMatrixView<T,conjA> const_conjugate_type;
        typedef ConstUpperTriMatrixView<T,trA> const_transpose_type;
        typedef ConstUpperTriMatrixView<T,adjA> const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,nonunitA> const_offdiag_type;
        typedef ConstLowerTriMatrixView<T,unitA> const_unitdiag_type;
        typedef ConstLowerTriMatrixView<T,nonunitA> const_nonunitdiag_type;
        typedef ConstLowerTriMatrixView<T,ndA> const_unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallLowerTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                ConstLowerTriMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstLowerTriMatrixView<T,nonconjA> const_nonconj_type;
        typedef LowerTriMatrixView<T,A> nonconst_type;

        typedef typename TriRefHelper<T,_conj,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,baseA> diag_type;
        typedef VectorView<T,baseA> diag_sub_type;

        typedef LowerTriMatrixView<T,A> subtrimatrix_type;
        typedef LowerTriMatrixView<T,nmA> subtrimatrix_step_type;
        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,ndnmA> submatrix_step_type;
        typedef VectorView<T,baseA> subvector_type;

        typedef LowerTriMatrixView<T,A> view_type;
        typedef LowerTriMatrixView<T,cstyleA> cview_type;
        typedef LowerTriMatrixView<T,fstyleA> fview_type;
        typedef LowerTriMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef LowerTriMatrixView<T,cmA> cmview_type;
        typedef LowerTriMatrixView<T,rmA> rmview_type;
        typedef LowerTriMatrixView<T,conjA> conjugate_type;
        typedef UpperTriMatrixView<T,trA> transpose_type;
        typedef UpperTriMatrixView<T,adjA> adjoint_type;

        typedef LowerTriMatrixView<T,nonunitA> offdiag_type;
        typedef LowerTriMatrixView<T,unitA> unitdiag_type;
        typedef LowerTriMatrixView<T,nonunitA> nonunitdiag_type;
        typedef LowerTriMatrixView<T,ndA> unknowndiag_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                SmallLowerTriMatrixView<real_type,TMV_UNKNOWN,twoSi,twoSj,twosAsm> ,
                LowerTriMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef LowerTriMatrixView<T,nonconjA> nonconj_type;
        typedef LowerTriMatrixView<T,An|NoAlias> noalias_type;
        typedef LowerTriMatrixView<T,An> alias_type;
    };

    template <class T, int A>
    class LowerTriMatrixView :
        public BaseMatrix_Tri_Mutable<LowerTriMatrixView<T,A> >
    {
    public:
        typedef LowerTriMatrixView<T,A> type;
        typedef BaseMatrix_Tri_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _dt = Traits<type>::_dt };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        LowerTriMatrixView(T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        LowerTriMatrixView(T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        LowerTriMatrixView(T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        LowerTriMatrixView(T* m, size_t s) :
            itsm(m), itss(s), itssi(_stepi), itssj(_stepj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != TMV_UNKNOWN);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
            TMVStaticAssert(_dt != TMV_UNKNOWN); 
        }

        LowerTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        LowerTriMatrixView(LowerTriMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int Si2, int Sj2, int A2>
        LowerTriMatrixView(SmallLowerTriMatrixView<T,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~LowerTriMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0;
#endif
        }


        //
        // Op = 
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        {
            return (
                (isunit() && i==j) ? T(1) :
                (i<j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,_conj,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi()+j*stepj()]); 
        }

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        int nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
        TMV_INLINE bool isunit() const { return itsdt == UnitDiag; }
        TMV_INLINE bool isrm() const
        {
            return Traits<type>::_rowmajor || 
                (!Traits<type>::_colmajor && stepj() == 1); 
        }
        TMV_INLINE bool iscm() const
        {
            return Traits<type>::_colmajor ||
                (!Traits<type>::_rowmajor && stepi() == 1); 
        }

    private :

        T* itsm;
        const size_t itss;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // LowerTriMatrixView



    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   UpperTriMatrixViewOf(T* m, n, S)
    //   UpperTriMatrixViewOf(T* m, n, si, sj)
    //   UpperTriMatrixViewOf(T* m, n, S, dt)
    //   UpperTriMatrixViewOf(T* m, n, si, sj, dt)
    //   UnitUpperTriMatrixViewOf(T* m, n, S)
    //   UnitUpperTriMatrixViewOf(T* m, n, si, sj)
    //   LowerTriMatrixViewOf(T* m, n, S)
    //   LowerTriMatrixViewOf(T* m, n, si, sj)
    //   LowerTriMatrixViewOf(T* m, n, S, dt)
    //   LowerTriMatrixViewOf(T* m, n, si, sj, dt)
    //   UnitLowerTriMatrixViewOf(T* m, n, S)
    //   UnitLowerTriMatrixViewOf(T* m, n, si, sj)
    //

    template <class T>
    inline UpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return UpperTriMatrixView<T,NonUnitDiag>(m,size,size,1);
        else
            return UpperTriMatrixView<T,NonUnitDiag>(m,size,1,size);
    }

    template <class T>
    inline ConstUpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,size,1);
        else
            return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,1,size);
    }

    template <class T>
    TMV_INLINE UpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return UpperTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

    template <class T>
    TMV_INLINE ConstUpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

    template <class T>
    inline UpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return UpperTriMatrixView<T,UnitDiag>(m,size,size,1);
        else
            return UpperTriMatrixView<T,UnitDiag>(m,size,1,size);
    }

    template <class T>
    inline ConstUpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T,UnitDiag>(m,size,size,1);
        else
            return ConstUpperTriMatrixView<T,UnitDiag>(m,size,1,size);
    }

    template <class T>
    TMV_INLINE UpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return UpperTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

    template <class T>
    TMV_INLINE ConstUpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstUpperTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

    template <class T>
    inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
        T* m, size_t size, StorageType stor, DiagType dt)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        if (stor == RowMajor) 
            return UpperTriMatrixView<T>(m,size,size,1,dt);
        else 
            return UpperTriMatrixView<T>(m,size,1,size,dt);
    }

    template <class T>
    inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
        const T* m, size_t size, StorageType stor, DiagType dt)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        if (stor == RowMajor)
            return ConstUpperTriMatrixView<T>(m,size,size,1,dt);
        else
            return ConstUpperTriMatrixView<T>(m,size,1,size,dt);
    }

    template <class T>
    TMV_INLINE_ND UpperTriMatrixView<T> UpperTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj, DiagType dt)
    {
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        return UpperTriMatrixView<T>(m,size,stepi,stepj,dt); 
    }

    template <class T>
    TMV_INLINE_ND ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj, DiagType dt)
    {
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        return ConstUpperTriMatrixView<T>(m,size,stepi,stepj,dt); 
    }
 
    template <class T>
    inline LowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return LowerTriMatrixView<T,NonUnitDiag>(m,size,size,1);
        else
            return LowerTriMatrixView<T,NonUnitDiag>(m,size,1,size);
    }

    template <class T>
    inline ConstLowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,size,1);
        else
            return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,1,size);
    }

    template <class T>
    TMV_INLINE LowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return LowerTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

    template <class T>
    TMV_INLINE ConstLowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

    template <class T>
    inline LowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return LowerTriMatrixView<T,UnitDiag>(m,size,size,1);
        else
            return LowerTriMatrixView<T,UnitDiag>(m,size,1,size);
    }

    template <class T>
    inline ConstLowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        const T* m, size_t size, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T,UnitDiag>(m,size,size,1);
        else
            return ConstLowerTriMatrixView<T,UnitDiag>(m,size,1,size);
    }

    template <class T>
    TMV_INLINE LowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj)
    { return LowerTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

    template <class T>
    TMV_INLINE ConstLowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj)
    { return ConstLowerTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

    template <class T>
    inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
        T* m, size_t size, StorageType stor, DiagType dt)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        if (stor == RowMajor)
            return LowerTriMatrixView<T>(m,size,size,1,dt);
        else
            return LowerTriMatrixView<T>(m,size,1,size,dt);
    }

    template <class T>
    inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
        const T* m, size_t size, StorageType stor, DiagType dt)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        if (stor == RowMajor)
            return ConstLowerTriMatrixView<T>(m,size,size,1,dt);
        else
            return ConstLowerTriMatrixView<T>(m,size,1,size,dt);
    }

    template <class T>
    TMV_INLINE_ND LowerTriMatrixView<T> LowerTriMatrixViewOf(
        T* m, size_t size, int stepi, int stepj, DiagType dt)
    {
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        return LowerTriMatrixView<T>(m,size,stepi,stepj,dt); 
    }

    template <class T>
    TMV_INLINE_ND ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
        const T* m, size_t size, int stepi, int stepj, DiagType dt)
    {
        TMVAssert(dt == UnitDiag || dt == NonUnitDiag);
        return ConstLowerTriMatrixView<T>(m,size,stepi,stepj,dt); 
    }

 


    //
    // Swap
    //

    template <class T, int A0, int A1, int A2>
    TMV_INLINE void Swap(
        UpperTriMatrix<T,A0,A1,A2>& m1, UpperTriMatrix<T,A0,A1,A2>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, int A>
    TMV_INLINE void Swap(
        BaseMatrix_Tri<M>& m1, UpperTriMatrixView<T,A> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, int A>
    TMV_INLINE void Swap(
        UpperTriMatrixView<T,A> m1, BaseMatrix_Tri<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int A1, int A2>
    TMV_INLINE void Swap(
        UpperTriMatrixView<T,A1> m1, UpperTriMatrixView<T,A2> m2)
    { DoSwap(m1,m2); }

    template <class T, int A0, int A1, int A2>
    TMV_INLINE void Swap(
        LowerTriMatrix<T,A0,A1,A2>& m1, LowerTriMatrix<T,A0,A1,A2>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, int A>
    TMV_INLINE void Swap(
        BaseMatrix_Tri<M>& m1, LowerTriMatrixView<T,A> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, int A>
    TMV_INLINE void Swap(
        LowerTriMatrixView<T,A> m1, BaseMatrix_Tri<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int A1, int A2>
    TMV_INLINE void Swap(
        LowerTriMatrixView<T,A1> m1, LowerTriMatrixView<T,A2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int A0, int A1, int A2>
    TMV_INLINE typename UpperTriMatrix<T,A0,A1,A2>::conjugate_type Conjugate(
        UpperTriMatrix<T,A0,A1,A2>& m)
    { return m.conjugate(); }
    template <class T, int A>
    TMV_INLINE typename UpperTriMatrixView<T,A>::conjugate_type Conjugate(
        UpperTriMatrixView<T,A> m)
    { return m.conjugate(); }

    template <class T, int A0, int A1, int A2>
    TMV_INLINE typename UpperTriMatrix<T,A0,A1,A2>::transpose_type Transpose(
        UpperTriMatrix<T,A0,A1,A2>& m)
    { return m.transpose(); }
    template <class T, int A>
    TMV_INLINE typename UpperTriMatrixView<T,A>::transpose_type Transpose(
        UpperTriMatrixView<T,A> m)
    { return m.transpose(); }

    template <class T, int A0, int A1, int A2>
    TMV_INLINE typename UpperTriMatrix<T,A0,A1,A2>::adjoint_type Adjoint(
        UpperTriMatrix<T,A0,A1,A2>& m)
    { return m.adjoint(); }
    template <class T, int A>
    TMV_INLINE typename UpperTriMatrixView<T,A>::adjoint_type Adjoint(
        UpperTriMatrixView<T,A> m)
    { return m.adjoint(); }

    template <class T, int A0, int A1, int A2>
    TMV_INLINE typename LowerTriMatrix<T,A0,A1,A2>::conjugate_type Conjugate(
        LowerTriMatrix<T,A0,A1,A2>& m)
    { return m.conjugate(); }
    template <class T, int A>
    TMV_INLINE typename LowerTriMatrixView<T,A>::conjugate_type Conjugate(
        LowerTriMatrixView<T,A> m)
    { return m.conjugate(); }

    template <class T, int A0, int A1, int A2>
    TMV_INLINE typename LowerTriMatrix<T,A0,A1,A2>::transpose_type Transpose(
        LowerTriMatrix<T,A0,A1,A2>& m)
    { return m.transpose(); }
    template <class T, int A>
    TMV_INLINE typename LowerTriMatrixView<T,A>::transpose_type Transpose(
        LowerTriMatrixView<T,A> m)
    { return m.transpose(); }

    template <class T, int A0, int A1, int A2>
    TMV_INLINE typename LowerTriMatrix<T,A0,A1,A2>::adjoint_type Adjoint(
        LowerTriMatrix<T,A0,A1,A2>& m)
    { return m.adjoint(); }
    template <class T, int A>
    TMV_INLINE typename LowerTriMatrixView<T,A>::adjoint_type Adjoint(
        LowerTriMatrixView<T,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class T, int A0, int A1, int A2>
    inline std::string TMV_Text(const UpperTriMatrix<T,A0,A1,A2>& m)
    {
        const int A = A0 | A1 | A2;
        std::ostringstream s;
        s << "UpperTriMatrix<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstUpperTriMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "ConstUpperTriMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int A>
    inline std::string TMV_Text(const UpperTriMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "UpperTriMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int A0, int A1, int A2>
    inline std::string TMV_Text(const LowerTriMatrix<T,A0,A1,A2>& m)
    {
        const int A = A0 | A1 | A2;
        std::ostringstream s;
        s << "LowerTriMatrix<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstLowerTriMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "ConstLowerTriMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int A>
    inline std::string TMV_Text(const LowerTriMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "LowerTriMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }
#endif


} // namespace tmv

#endif
