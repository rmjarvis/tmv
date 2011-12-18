

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
// the matrix.  The valid options are:
// ColMajor or RowMajor or DiagMajor
// CStyle or FortranStyle
// WithDivider or NoDivider
//
// The options are treated as a bit field, so you | them together to 
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
// size_t BandStorageLength(stor, colsize, rowsize, nlo, nhi);
//  
// Constructors:
//   
//    BandMatrix<T,A>(size_t colsize, size_t rowsize, int nlo, int nhi)
//        Makes a BandMatrix with column size = colsize and row size = rowsize 
//        with nhi non-zero superdiagonals and nlo non-zero subdiagonals
//        with _uninitialized_ values
//   
//    BandMatrix<T,A>(size_t colsize, size_t rowsize, int nlo, int nhi, T x)
//        The same as above, but all values are initialized to x
//    
//    BandMatrix<T,A>(const BaseMatrix<M>& m, int nlo, int nhi)
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
//            size_t colsize, size_t rowsize, int nlo, int nhi, 
//            StorageType stor)
//    BandMatrixView BandMatrixViewOf(const T* m, 
//            size_t colsize, size_t rowsize, int nlo, int nhi, 
//            StorageType stor)
//    ConstBandMatrixView BandMatrixViewOf(const T* m,
//            size_t colsize, size_t rowsize, int nlo, int nhi, 
//            int stepi, int stepj)
//    BandMatrixView BandMatrixViewOf(const T* m, 
//            size_t colsize, size_t rowsize, int nlo, int nhi, 
//            int stepi, int stepj)
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
//    value_type operator()(int i, int j) const
//    reference operator()(int i, int j)
//    value_type cref(int i, int j) const
//    reference refint i, int j)
//        Return the (i,j) element of the BandMatrix.
//        The first two respect the index-style of the underlying matrix.
//        The next two always use CStyle indexing and do not do any
//        checking of the validity of i,j.
//
//    const_row_sub_type row(int i, int j1, int j2) const
//    row_sub_type row(int i, int j1, int j2)
//        Return a subset of the ith row of the BandMatrix
//        The range (i,j1)..(i,j2-1) must be entirely within the band.
//
//    const_col_sub_type col(int j, int i1, int i2) const
//    col_sub_type col(int j, int i1, int i2)
//        Return a subset of the jth column of the BandMatrix
//        The range (i1,j)..(i2-1,j) must be entirely within the band.
//
//    const_diag_type diag() const
//    diag_type diag()
//        Return the main diagonal of the BandMatrix
//
//    const_diag_sub_type diag(int i) const
//    diag_sub_type diag(int i)
//    const_diag_sub_type diag(int i, int j1, int j2) const
//    diag_sub_type diag(int i, int j1, int j2)
//        Return the super or sub-diagonal i. 
//        If i > 0 return the super diagonal starting at (0,i)
//        If i < 0 return the sub diagonal starting at (|i|,0)
//        If j1,j2 are given, it returns the diagonal subvector 
//        either from (j1,i+j1) to (j2,i+j2) (for i>0) 
//        or from (|i|+j1,j1) to (|i|+j2,j2) (for i<0)
//
//
// Modifying Functions
//
//    type& setZero()
//    type& setAllTo(T x)
//    type& addToAll(T x)
//    type& clip(T x)
//    type& conjugateSelf()
//    type& transposeSelf() 
//        Must be square, and have nhi=nlo for this function
//    type& setToIdentity(T x = 1)
//    void Swap(BandMatrix& m1, BandMatrix& m2)
//        Must be the same size and have the same band structure (nlo,nhi)
//
//
// Functions of BandMatrices:
//        (These are all both member functions and functions of a BandMatrix,
//         so Norm(m) and m.norm() for example are equivalent.)
//
//    value_type m.det() const                          or Det(m)
//    float_type m.logDet(value_type* sign=0) const     or LogDet(m)
//    value_type m.trace()                              or Trace(m)
//    float_type m.norm() or m.normF()                  or Norm(m), NormF(m)
//    value_type m.sumElements() or SumElements(m)
//    float_type m.sumAbsElements()                     or SumAbsElements(m)
//    real_type m.sumAbs2Elements()                     or SumAbs2Elements(m)
//    real_type m.normSq()                              or NormSq(m)
//    float_type m.normSq(float_type scale)             or NormSq(m)
//    float_type m.norm1()                              or Norm1(m)
//    float_type m.norm2()                              or Norm2(m)
//    float_type m.normInf()                            or NormInf(m)
//    float_type m.condition()
//    float_type m.maxAbsElement()                      or MaxAbsElements(m)
//    real_type m.maxAbs2Element()                      or MaxAbs2Elements(m)
//
//    inverse_type m.inverse()                          or Inverse(m)
//    void m.makeInverse(minv)
//    void m.makeInverseATA(invata)
//
//
// Views of a BandMatrix:
//
//    subbandmatrix_type subBandMatrix(int i1, int i2, int j1, int j2, 
//            int nlo, int nhi)
//    subbandmatrix_step_type subBandMatrix(int i1, int i2, int j1, int j2, 
//            int nlo, int nhi, int istep, int jstep)
//        Returns a BandMatrixView with (i1,j1) in the upper left corner,
//        (i2,j2) in the lower right corder, and with nlo and nhi 
//        sub- and super-diagonals respectively.
//        All members of the new subBandMatrix must be within the 
//        original band.
//
//    subbandmatrix_type subBandMatrix(int i1, int i2, int j1, int j2)
//        This is normally equivalenet to 
//        b.subBandMatrix(i1,i2,j1,j2,b.nlo(),b.nhi())
//        However, it lowers the new nlo, nhi as necessary for small 
//        i2-i1 or j2-j1 if the full band-width doesn't fit into the 
//        new matrix size.
//
//    submatrix_type subMatrix(int i1, int i2, int j1, int j2)
//    submatrix_step_type subMatrix(int i1, int i2, int j1, int j2,
//            int istep, int jstep)
//    subvector_type subVector(int i1, int i2, int istep, int jstep,
//            int size)
//        Just like the regular Matrix version, but all entries in the 
//        submatrix or subvector must be within the band.
//
//    rowrange_type rowRange(int i1, int i2)
//        A bit different from the Matrix version, since the rowsize will
//        not be the full rowsize of the source BandMatrix.
//        Instead, this returns (as a BandMatrix) the nontrivial portion
//        of the BandMatrix in rows i1..i2 (not including i2)
//    colrange_type colRange(int j1, int j2)
//        As with rowRange, this returns the nontrivial portion of the 
//        BandMatrix in columns j1..j2 (not including j2)
//    diagrange_type diagRange(int i1, int i2)
//        Returns a thinner BandMatrix including the diagonals from i1..i2
//        (not including i2)
//
//    upperband_type upperBand()
//    lowerband_type lowerBand()
//        Returns a BandMatrix of only the upper or lower portion of the
//        matrix.
//
//    upperbandoff_type upperBandOff(int noff)
//    lowerbandoff_type lowerBandOff(int noff)
//        Returns a BandMatrix of only the strictly upper or lower portion
//        of the matrix.  (ie. not including the diagonal.)
//        Conceptually equivalent to upperBand().offDiag(noff).
//        Except that a BandMatrix doesn't have an offDiag(noff) function.
//
//    view_type m.view()
//    transpose_type m.transpose() or Transpose(m)
//    conjugate_type m.conjugate() or Conjugate(m)
//    adjoint_type m.adjoint() or Adjoint(m)
//    realpart_type m.realPart()
//    imagpart_type m.imagPart()
//    nonconj_type nonConj()
//    noalias_type noAlias()
//    alias_type alias()
//    nonconst_type nonConst()
//    const_view_type constView()
//        Just like the regular Matrix versions
//
//    linearview_type m.linearView()
//        Returns a VectorView with all the elements of the BandMatrix.
//        Extra care is required with this than with the corresponding
//        function for a regular matrix, since the view also includes the 
//        few elements in memory that fall outside the actual matrix.
//        These elements are generally uninitialized memory.  
//        So this function should not be used for anything that 
//        does a calculation on the elements, like maxElement.
//        It is mostly useful for implementing things like setZero
//        where it is ok to operate on the extra elements too.
//
//
//
// Operators:
//    Band matrices can do all the same arithmetic operations as normal
//    matrices whenever the meaning makes sense.
//
//    So for example, there is no implementation of b = v1 ^ v2, 
//    since the vector outer product is necessarily a full rectangular
//    matrix.
//
//    The return types may be assignable to a band matrix if that makes 
//    sense.  For example b1 * b2 is banded, so that can be assigned
//    to a new band matrix, as can b * U (where U is an upper triangular
//    matrix).  But b * m is not banded, so that needs to be assigned to 
//    a regular full matrix object.
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

#include "TMV_BaseMatrix_Band.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseMatrix_Diag.h"
#include "TMV_Divider.h"
#include "TMV_Array.h"
#include "TMV_VIt.h"
#include "TMV_MIt.h"

namespace tmv {


    //
    // Define the DivHelper class for rectangular matrices:
    //

    template <class M> class BandLUD;
    template <class T> class InstBandLUD;
    template <class M> class BandQRD;
    template <class T> class InstBandQRD;
    template <class M> class BandSVD;
    template <class T> class InstBandSVD;

    template <class T>
    class BandMatrixDivHelper2 :
        public DivHelper<T>
    {
    public:
        // There is no inline version of this.  These functions are defined
        // in TMV_BandMatrix.cpp for the instantiated types.

        typedef const InstBandLUD<T>& lud_type;
        typedef const InstBandQRD<T>& qrd_type;
        typedef const InstBandSVD<T>& svd_type;

        BandMatrixDivHelper2();
        ~BandMatrixDivHelper2();

        void setDiv() const;
        Matrix<T> getM() const;
        bool mIsSquare() const;

#if 0
        bandlud_type lud() const;
        bandqrd_type qrd() const;
        bandsvd_type svd() const;
#endif

        virtual ConstBandMatrixView<T> getConstView() const=0;

    }; // MatrixDivHelper2

    template <class M, bool hasdivider>
    class BandMatrixDivHelper1;

    template <class M>
    class BandMatrixDivHelper1<M,true> :
        public BandMatrixDivHelper2<typename Traits<M>::value_type>
    {
    public:

        typedef typename Traits<M>::value_type T;

        TMV_INLINE ConstBandMatrixView<T> getConstView() const
        { return mat2().constView(); }

        // Use mat2() rather than mat() to avoid having to disambiguate
        // it with the normal BaseMatrix mat() function.
        TMV_INLINE const M& mat2() const
        { return static_cast<const M&>(*this); }
    };

    template <class M>
    class BandMatrixDivHelper1<M,false>
    {
    public:
        typedef typename Traits<M>::lud_type lud_type;
        typedef typename Traits<M>::qrd_type qrd_type;
        typedef typename Traits<M>::svd_type svd_type;

        TMV_INLINE void resetDivType() const {} 
#if 0
        TMV_INLINE lud_type lud() const { return lud_type(mat2(),false); }
        TMV_INLINE qrd_type qrd() const { return qrd_type(mat2(),false); }
        TMV_INLINE svd_type svd() const { return svd_type(mat2(),false); }
#endif

        TMV_INLINE const M& mat2() const
        { return static_cast<const M&>(*this); }
    };

    template <class M>
    class BandMatrixDivHelper :
        public BandMatrixDivHelper1<M,Traits<M>::_hasdivider>
    {};

    static size_t BandStorageLength(
        int stor, size_t cs, size_t rs, int lo, int hi)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        if (cs == 0 || rs == 0) return 0;
        else if (cs == rs) return (cs-1)*(lo+hi)+cs;
        else {
            // correct cs, rs to be actual end of data
            if (cs > rs+lo) cs = rs+lo;
            if (rs > cs+hi) rs = cs+hi;

            if (stor == RowMajor)
                // si = lo+hi, sj = 1, size = (cs-1)*si + (rs-1)*sj + 1
                return (cs-1)*(lo+hi) + rs;
            else if (stor == ColMajor)
                // si = 1, sj = lo+hi, size = (cs-1)*si + (rs-1)*sj + 1
                return (rs-1)*(lo+hi) + cs;
            else if (cs > rs)
                // si = -rs, sj = 1+rs, size = (rs-lo-hi-1)*si + (rs-1)*sj + 1
                // size = (rs-lo-hi-1)(-rs) + (rs-1)(1+rs) + 1
                //      = rs(lo+hi+1)
                return rs*(lo+hi+1);
            else
                // si = 1-cs, sj = cs, size = (rs-lo-hi-1)*si + (rs-1)*sj + 1
                // size = (cs-lo-hi-1+rs-cs)(1-cs) + (rs-1)(cs) + 1
                //      = cs(lo+hi+1-rs+rs-1) + rs-lo-hi-1 + 1
                //      = (cs-1)(lo+hi) + rs
                return (cs-1)*(lo+hi) + rs;
        }
    }


    template <class T, int A0, int A1>
    struct Traits<BandMatrix<T,A0,A1> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A01 = A0 | A1 };
        enum { A = (A01 & ~CheckAlias) | (
                ((Attrib<A01>::rowmajor || Attrib<A01>::diagmajor) ? 
                 0 : ColMajor) |
                ( !Traits<real_type>::isinst ? NoDivider :
                  Traits<real_type>::isinteger ? NoDivider :
                  Attrib<A01>::nodivider ? 0 : WithDivider ) )};
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor || 
                 Attrib<A>::diagmajor) &&
                !(Attrib<A>::rowmajor && Attrib<A>::colmajor) &&
                !(Attrib<A>::rowmajor && Attrib<A>::diagmajor) &&
                !(Attrib<A>::colmajor && Attrib<A>::diagmajor) &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                !Attrib<A>::checkalias &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider) ) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger) ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef BandMatrix<T,A0,A1> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef BandMatrix<T,A01> copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = TMV_UNKNOWN };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = (_diagmajor ? 1 : TMV_UNKNOWN) };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = true };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { _hasdivider = Attrib<A>::withdivider };

#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , SVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? 0 : NoAlias) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA };
        enum { twosA = isreal ? int(ndA) : (ndA & ~AllStorageType) };
        enum { Asm = _checkalias ? (ndA | CheckAlias) : (ndA & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~AllStorageType) };
        enum { An = ndA & ~NoAlias };

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,diagA> const_diag_type;
        typedef ConstVectorView<T,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,nmA> const_submatrix_step_type;
        typedef ConstVectorView<T,vecA> const_subvector_type;
        typedef ConstBandMatrixView<T,ndA> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmA> const_subbandmatrix_step_type;

        typedef ConstBandMatrixView<T,ndA> const_colrange_type;
        typedef ConstBandMatrixView<T,ndA> const_rowrange_type;
        typedef ConstBandMatrixView<T,ndA> const_diagrange_type;

        typedef ConstBandMatrixView<T,ndA> const_view_type;
        typedef ConstBandMatrixView<T,cstyleA> const_cview_type;
        typedef ConstBandMatrixView<T,fstyleA> const_fview_type;
        typedef ConstBandMatrixView<T> const_xview_type;
        typedef ConstBandMatrixView<T,cmA> const_cmview_type;
        typedef ConstBandMatrixView<T,rmA> const_rmview_type;
        typedef ConstBandMatrixView<T,dmA> const_dmview_type;
        typedef ConstBandMatrixView<T,conjA> const_conjugate_type;
        typedef ConstBandMatrixView<T,trA> const_transpose_type;
        typedef ConstBandMatrixView<T,adjA> const_adjoint_type;

        typedef ConstBandMatrixView<T,ndA> const_upperband_type;
        typedef ConstBandMatrixView<T,ndA> const_lowerband_type;
        typedef ConstBandMatrixView<T,ndA> const_upperbandoff_type;
        typedef ConstBandMatrixView<T,ndA> const_lowerbandoff_type;

        enum { xx = TMV_UNKNOWN }; // Just for brevity.
        typedef typename TypeSelect< iscomplex ,
                ConstSmallBandMatrixView<real_type,xx,xx,xx,xx,twoSi,twoSj,twosAsm> ,
                ConstBandMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecA|Unit)> const_linearview_type;
        typedef ConstBandMatrixView<T,nonconjA> const_nonconj_type;
        typedef BandMatrixView<T,ndA> nonconst_type;

        typedef T& reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef CDMIt<type> const_diagmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef DMIt<type> diagmajor_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,diagA> diag_type;
        typedef VectorView<T,diagA> diag_sub_type;

        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,nmA> submatrix_step_type;
        typedef VectorView<T,vecA> subvector_type;
        typedef BandMatrixView<T,ndA> subbandmatrix_type;
        typedef BandMatrixView<T,nmA> subbandmatrix_step_type;

        typedef BandMatrixView<T,ndA> colrange_type;
        typedef BandMatrixView<T,ndA> rowrange_type;
        typedef BandMatrixView<T,ndA> diagrange_type;

        typedef BandMatrixView<T,ndA> view_type;
        typedef BandMatrixView<T,cstyleA> cview_type;
        typedef BandMatrixView<T,fstyleA> fview_type;
        typedef BandMatrixView<T> xview_type;
        typedef BandMatrixView<T,cmA> cmview_type;
        typedef BandMatrixView<T,rmA> rmview_type;
        typedef BandMatrixView<T,dmA> dmview_type;
        typedef BandMatrixView<T,conjA> conjugate_type;
        typedef BandMatrixView<T,trA> transpose_type;
        typedef BandMatrixView<T,adjA> adjoint_type;

        typedef BandMatrixView<T,ndA> upperband_type;
        typedef BandMatrixView<T,ndA> lowerband_type;
        typedef BandMatrixView<T,ndA> upperbandoff_type;
        typedef BandMatrixView<T,ndA> lowerbandoff_type;

        typedef typename TypeSelect< iscomplex ,
                SmallBandMatrixView<real_type,xx,xx,xx,xx,twoSi,twoSj,twosAsm> ,
                BandMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,(vecA|Unit)> linearview_type;
        typedef BandMatrixView<T,nonconjA> nonconj_type;
        typedef BandMatrixView<T,An|NoAlias> noalias_type;
        typedef BandMatrixView<T,An> alias_type;
    };

    template <class T, int A0, int A1>
    class BandMatrix : 
        public BaseMatrix_Band_Mutable<BandMatrix<T,A0,A1> >,
        public BandMatrixDivHelper<BandMatrix<T,A0,A1> >
    {
    public:

        typedef BandMatrix<T,A0,A1> type;
        typedef BaseMatrix_Band_Mutable<type> base_mut;
        typedef BandMatrixDivHelper<type> divhelper;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };
        enum { _stor = Attrib<_attrib>::stor };

        typedef typename Traits<type>::linearview_type linearview_type;
        typedef typename Traits<type>::const_linearview_type 
            const_linearview_type;

        //
        // Constructors
        //

        BandMatrix() : 
            itscs(0), itsrs(0), itsnlo(0), itsnhi(0), linsize(0),
            itssi(_colmajor ? 1 : 0), itssj(_rowmajor ? 1 : 0),
            itsm1(0), itsm(0)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        BandMatrix(size_t cs, size_t rs, int lo, int hi) :
            itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? lo+hi : _colmajor ? 1 : 
                  rs >= cs ? 1-int(cs) : -int(rs)),
            itssj(_rowmajor ? 1 : _colmajor ? lo+hi : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? lo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(lo >= 0 && lo < int(cs));
            TMVAssert(hi >= 0 && hi < int(rs));
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        BandMatrix(size_t cs, size_t rs, int lo, int hi, T x) :
            itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? lo+hi : _colmajor ? 1 : 
                  rs >= cs ? 1-int(cs) : -int(rs)),
            itssj(_rowmajor ? 1 : _colmajor ? lo+hi : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? lo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(lo >= 0 && lo < int(cs));
            TMVAssert(hi >= 0 && hi < int(rs));
            this->setAllTo(x);
        }

        BandMatrix(const type& m2) :
            itscs(m2.itscs), itsrs(m2.itsrs),
            itsnlo(m2.itsnlo), itsnhi(m2.itsnhi),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? itsnlo+itsnhi : _colmajor ? 1 : 
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? itsnlo+itsnhi : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? itsnlo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            typename Traits<type>::noalias_type na = this->noAlias();
            m2.assignTo(na);
        }

        template <class M2>
        BandMatrix(const BaseMatrix<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), 
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? itsnlo+itsnhi : _colmajor ? 1 : 
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? itsnlo+itsnhi : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? itsnlo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            typename Traits<type>::noalias_type na = this->noAlias();
            m2.assignTo(na);
        }

        template <class M2>
        BandMatrix(const BaseMatrix_Band<M2>& m2, int lo, int hi) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), itsnlo(lo), itsnhi(hi),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? lo+hi : _colmajor ? 1 : 
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? lo+hi : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? lo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(lo >= 0 && lo <= m2.nlo());
            TMVAssert(hi >= 0 && hi <= m2.nhi());
            typename Traits<type>::noalias_type na = this->noAlias();
            m2.cDiagRange(-lo,hi+1).assignTo(na);
        }

        template <class M2>
        BandMatrix(const BaseMatrix_Rec<M2>& m2, int lo, int hi) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), itsnlo(lo), itsnhi(hi),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? lo+hi : _colmajor ? 1 : 
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? lo+hi : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? lo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(lo >= 0 && lo < int(m2.colsize()));
            TMVAssert(hi >= 0 && hi < int(m2.rowsize()));
            typedef typename M2::value_type T2;
            ConstBandMatrixView<T2> m2b(
                m2.cptr(),m2.colsize(),m2.rowsize(),lo,hi,
                m2.stepi(),m2.stepj());
            typename Traits<type>::noalias_type na = this->noAlias();
            m2b.assignTo(na);
        }

        template <class M2>
        BandMatrix(const BaseMatrix_Tri<M2>& m2, int lohi) :
            itscs(m2.size()), itsrs(m2.size()), 
            itsnlo(Maybe<M2::_upper>::select(0,lohi)),
            itsnhi(Maybe<M2::_upper>::select(lohi,0)),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? lohi : _colmajor ? 1 : 1-int(itscs)),
            itssj(_rowmajor ? 1 : _colmajor ? lohi : int(itscs)),
            itsm1(linsize), itsm(itsm1.get() - (_diagmajor ? itsnlo*itssi : 0))
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(lohi >= 0 && lohi < int(m2.size()));
            typedef typename M2::value_type T2;
            ConstBandMatrixView<T2> m2b(
                m2.cptr(),m2.size(),m2.size(),itsnlo,itsnhi,
                m2.stepi(),m2.stepj());
            typename Traits<type>::noalias_type na = this->noAlias();
            m2b.assignTo(na);
        }

        template <class M2>
        BandMatrix(const BaseMatrix_Diag<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), itsnlo(0), itsnhi(0),
            linsize(BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi)),
            itssi(_rowmajor ? itsnlo+itsnhi : _colmajor ? 1 : 1-int(itscs)),
            itssj(_rowmajor ? 1 : _colmajor ? itsnlo+itsnhi : int(itscs)),
            itsm1(linsize), itsm(itsm1.get())
        {
            TMVStaticAssert(Traits<type>::okA);
            typename Traits<type>::noalias_type diagna = this->diag().noAlias();
            m2.diag().assignTo(diagna);
        }

        ~BandMatrix() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::destr_value());
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
        TMV_INLINE type& operator=(const BaseMatrix_Band<M2>& m2)
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
        // LinearView is only allowed for an actual BandMatrix
        // Views cannot be linearized.
        //
        
        TMV_INLINE_ND linearview_type linearView() 
        {
#ifdef TMV_USE_VALGRIND
            AlignedArray<int> ok(ls());
            for (size_t i=0;i<ls();++i) ok[i] = 0;
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) {
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = 1;
                }
            }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return linearview_type(start_mem(),ls(),1);
        }

        TMV_INLINE_ND const_linearview_type linearView() const
        {
#ifdef TMV_USE_VALGRIND
            AlignedArray<int> ok(ls());
            for (size_t i=0;i<ls();++i) ok[i] = 0;
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) {
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = 1;
                }
            }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return const_linearview_type(start_mem(),ls(),1);
        }

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }
        TMV_INLINE const T* start_mem() const { return itsm1; }
        TMV_INLINE T* start_mem() { return itsm1; }

        T cref(int i, int j) const
        { 
            if (InBand(i,j,itsnlo,itsnhi)) return itsm[i*itssi + j*itssj]; 
            else return T(0);
        }

        T& ref(int i, int j)
        { return itsm[i*itssi + j*itssj]; }

        void swapWith(type& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(m2.nlo() == nlo());
            TMVAssert(m2.nhi() == nhi());
            if (itsm1.get() == m2.itsm1.get()) return;
            itsm1.swapWith(m2.itsm1);
            TMV_SWAP(itsm,m2.itsm);
        }

        void resize(size_t cs, size_t rs, int lo, int hi)
        {
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::destr_value());
#endif
            itscs = cs;
            itsrs = rs;
            itsnlo = lo;
            itsnhi = hi;
            divhelper::resetDivType();
            linsize = BandStorageLength(_stor,itscs,itsrs,itsnlo,itsnhi);
            itsm1.resize(linsize);
            itsm = itsm1.get() - (_diagmajor ? lo*itssi : 0);
            itssi = _rowmajor ? lo+hi : _colmajor ? 1 : 
                  rs >= cs ? 1-int(cs) : -int(rs);
            itssj = _rowmajor ? 1 : _colmajor ? lo+hi : -itssi+1;
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        TMV_INLINE size_t ls() const { return linsize; }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nlo() const { return itsnlo; }
        TMV_INLINE int nhi() const { return itsnhi; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE int diagstep() const 
        { return _diagmajor ? 1 : stepi()+stepj(); }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }
        TMV_INLINE bool isdm() const { return _diagmajor; }

    private:

        size_t itscs;
        size_t itsrs;
        int itsnlo;
        int itsnhi;
        size_t linsize;
        CheckedInt<_stepi> itssi;
        CheckedInt<_stepj> itssj;
        AlignedArray<T> itsm1;
        T* itsm;

    }; // BandMatrix

    template <class T, int A0>
    struct Traits<ConstBandMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~CheckAlias) };
        enum { okA = (
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                !Attrib<A>::checkalias &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger ) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstBandMatrixView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = TMV_UNKNOWN };
        enum { _unit = false };
        enum { _nonunit = true };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = (_diagmajor ? 1 : TMV_UNKNOWN) };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = false };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { copyA = (
                (_rowmajor ? RowMajor : _colmajor ? ColMajor : DiagMajor) |
                (_fort ? FortranStyle : CStyle) |
                NoDivider ) };
        typedef BandMatrix<T,copyA> copy_type;
        enum { _hasdivider = Attrib<A>::withdivider };

#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , BandLUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , BandQRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , BandSVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? 0 : NoAlias) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (ndA | CheckAlias) : (ndA & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,diagA> const_diag_type;
        typedef ConstVectorView<T,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,nmA> const_submatrix_step_type;
        typedef ConstVectorView<T,vecA> const_subvector_type;
        typedef ConstBandMatrixView<T,ndA> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmA> const_subbandmatrix_step_type;

        typedef ConstBandMatrixView<T,ndA> const_colrange_type;
        typedef ConstBandMatrixView<T,ndA> const_rowrange_type;
        typedef ConstBandMatrixView<T,ndA> const_diagrange_type;

        typedef ConstBandMatrixView<T,ndA> const_view_type;
        typedef ConstBandMatrixView<T,cstyleA> const_cview_type;
        typedef ConstBandMatrixView<T,fstyleA> const_fview_type;
        typedef ConstBandMatrixView<T,(_conj?Conj:NonConj)> const_xview_type;
        typedef ConstBandMatrixView<T,cmA> const_cmview_type;
        typedef ConstBandMatrixView<T,rmA> const_rmview_type;
        typedef ConstBandMatrixView<T,dmA> const_dmview_type;
        typedef ConstBandMatrixView<T,conjA> const_conjugate_type;
        typedef ConstBandMatrixView<T,trA> const_transpose_type;
        typedef ConstBandMatrixView<T,adjA> const_adjoint_type;

        typedef ConstBandMatrixView<T,ndA> const_upperband_type;
        typedef ConstBandMatrixView<T,ndA> const_lowerband_type;
        typedef ConstBandMatrixView<T,ndA> const_upperbandoff_type;
        typedef ConstBandMatrixView<T,ndA> const_lowerbandoff_type;

        enum { xx = TMV_UNKNOWN }; // Just for brevity.
        typedef typename TypeSelect< iscomplex && (_colmajor||_rowmajor||_diagmajor) ,
                ConstSmallBandMatrixView<real_type,xx,xx,xx,xx,twoSi,twoSj,twosAsm> ,
                ConstBandMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstBandMatrixView<T,nonconjA> const_nonconj_type;
        typedef BandMatrixView<T,ndA> nonconst_type;

        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef CDMIt<type> const_diagmajor_iterator;
    };

    template <class T, int A>
    class ConstBandMatrixView : 
        public BaseMatrix_Band<ConstBandMatrixView<T,A> >,
        public BandMatrixDivHelper<ConstBandMatrixView<T,A> >
    {
    public:
        typedef ConstBandMatrixView<T,A> type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE ConstBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo, int hi, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(sj) 
        { 
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE ConstBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo, int hi, int si) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(_stepj) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo, int hi) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(_stepi), itssj(_stepj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != TMV_UNKNOWN);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstBandMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE ConstBandMatrixView(const ConstBandMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        TMV_INLINE ConstBandMatrixView(const BandMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int LO2, int HI2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstBandMatrixView(
            const ConstSmallBandMatrixView<T,M2,N2,LO2,HI2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int LO2, int HI2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstBandMatrixView(
            const SmallBandMatrixView<T,M2,N2,LO2,HI2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~ConstBandMatrixView() {
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
            if (InBand(i,j,itsnlo,itsnhi)) 
                return DoConj<_conj>(itsm[i*stepi() + j*stepj()]); 
            else return T(0);
        }

        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nlo() const { return itsnlo; }
        TMV_INLINE int nhi() const { return itsnhi; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE int diagstep() const 
        { return _diagmajor ? 1 : stepi()+stepj(); }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && !_diagmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && !_diagmajor && stepi() == 1); }
        TMV_INLINE bool isdm() const 
        { return _diagmajor || (!_rowmajor && !_colmajor && diagstep() == 1); }

    private :

        const T* itsm;
        const size_t itscs;
        const size_t itsrs;
        const int itsnlo;
        const int itsnhi;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;

    }; // ConstBandMatrixView

    template <class T, int A0>
    struct Traits<BandMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~CheckAlias) };
        enum { okA = (
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                !Attrib<A>::checkalias &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger ) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef BandMatrixView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = TMV_UNKNOWN };
        enum { _unit = false };
        enum { _nonunit = true };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = (_colmajor ? 1 : TMV_UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : TMV_UNKNOWN) };
        enum { _diagstep = TMV_UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = false };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                NoDivider ) };
        typedef BandMatrix<T,copyA> copy_type;
        enum { _hasdivider = Attrib<A>::withdivider };

#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , BandLUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , BandQRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , BandSVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? 0 : NoAlias) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { Asm = _checkalias ? (ndA | CheckAlias) : (ndA & ~NoAlias) };
        enum { twosAsm = isreal ? int(Asm) : (Asm & ~Conj & ~AllStorageType) };
        enum { An = ndA & ~NoAlias };

        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,diagA> const_diag_type;
        typedef ConstVectorView<T,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,nmA> const_submatrix_step_type;
        typedef ConstVectorView<T,vecA> const_subvector_type;
        typedef ConstBandMatrixView<T,ndA> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmA> const_subbandmatrix_step_type;

        typedef ConstBandMatrixView<T,ndA> const_colrange_type;
        typedef ConstBandMatrixView<T,ndA> const_rowrange_type;
        typedef ConstBandMatrixView<T,ndA> const_diagrange_type;

        typedef ConstBandMatrixView<T,ndA> const_view_type;
        typedef ConstBandMatrixView<T,cstyleA> const_cview_type;
        typedef ConstBandMatrixView<T,fstyleA> const_fview_type;
        typedef ConstBandMatrixView<T,(_conj?Conj:NonConj)> const_xview_type;
        typedef ConstBandMatrixView<T,cmA> const_cmview_type;
        typedef ConstBandMatrixView<T,rmA> const_rmview_type;
        typedef ConstBandMatrixView<T,dmA> const_dmview_type;
        typedef ConstBandMatrixView<T,conjA> const_conjugate_type;
        typedef ConstBandMatrixView<T,trA> const_transpose_type;
        typedef ConstBandMatrixView<T,adjA> const_adjoint_type;

        typedef ConstBandMatrixView<T,ndA> const_upperband_type;
        typedef ConstBandMatrixView<T,ndA> const_lowerband_type;
        typedef ConstBandMatrixView<T,ndA> const_upperbandoff_type;
        typedef ConstBandMatrixView<T,ndA> const_lowerbandoff_type;

        enum { xx = TMV_UNKNOWN }; // Just for brevity.
        typedef typename TypeSelect< iscomplex && (_colmajor||_rowmajor||_diagmajor) ,
                ConstSmallBandMatrixView<real_type,xx,xx,xx,xx,twoSi,twoSj,twosAsm> ,
                ConstBandMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstBandMatrixView<T,nonconjA> const_nonconj_type;
        typedef BandMatrixView<T,ndA> nonconst_type;

        typedef T& reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef CDMIt<type> const_diagmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef DMIt<type> diagmajor_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,diagA> diag_type;
        typedef VectorView<T,diagA> diag_sub_type;

        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,nmA> submatrix_step_type;
        typedef VectorView<T,vecA> subvector_type;
        typedef BandMatrixView<T,ndA> subbandmatrix_type;
        typedef BandMatrixView<T,nmA> subbandmatrix_step_type;

        typedef BandMatrixView<T,ndA> colrange_type;
        typedef BandMatrixView<T,ndA> rowrange_type;
        typedef BandMatrixView<T,ndA> diagrange_type;

        typedef BandMatrixView<T,ndA> view_type;
        typedef BandMatrixView<T,cstyleA> cview_type;
        typedef BandMatrixView<T,fstyleA> fview_type;
        typedef BandMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef BandMatrixView<T,cmA> cmview_type;
        typedef BandMatrixView<T,rmA> rmview_type;
        typedef BandMatrixView<T,dmA> dmview_type;
        typedef BandMatrixView<T,conjA> conjugate_type;
        typedef BandMatrixView<T,trA> transpose_type;
        typedef BandMatrixView<T,adjA> adjoint_type;

        typedef BandMatrixView<T,ndA> upperband_type;
        typedef BandMatrixView<T,ndA> lowerband_type;
        typedef BandMatrixView<T,ndA> upperbandoff_type;
        typedef BandMatrixView<T,ndA> lowerbandoff_type;

        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                SmallBandMatrixView<real_type,xx,xx,xx,xx,twoSi,twoSj,twosAsm> ,
                BandMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,(vecA|Unit)> linearview_type;
        typedef BandMatrixView<T,nonconjA> nonconj_type;
        typedef BandMatrixView<T,An|NoAlias> noalias_type;
        typedef BandMatrixView<T,An> alias_type;
    };

    template <class T, int A>
    class BandMatrixView : 
        public BaseMatrix_Band_Mutable<BandMatrixView<T,A> >,
        public BandMatrixDivHelper<BandMatrixView<T,A> >
    {
    public:

        typedef BandMatrixView<T,A> type;
        typedef BaseMatrix_Band_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE BandMatrixView(
            T* m, size_t cs, size_t rs, int lo, int hi, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE BandMatrixView(
            T* m, size_t cs, size_t rs, int lo, int hi, int si) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(_stepj) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
        }

        TMV_INLINE BandMatrixView(T* m, size_t cs, size_t rs, int lo, int hi) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(_stepi), itssj(_stepj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != TMV_UNKNOWN);
            TMVStaticAssert(_stepj != TMV_UNKNOWN); 
        }

        TMV_INLINE BandMatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE BandMatrixView(BandMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int LO2, int HI2, int Si2, int Sj2, int A2>
        TMV_INLINE BandMatrixView(
            SmallBandMatrixView<T,M2,N2,LO2,HI2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~BandMatrixView() {
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
        TMV_INLINE type& operator=(const BaseMatrix_Band<M2>& m2)
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
            if (InBand(i,j,itsnlo,itsnhi)) 
                return DoConj<_conj>(itsm[i*stepi() + j*stepj()]); 
            else return T(0);
        }

        reference ref(int i, int j) 
        { return reference(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nlo() const { return itsnlo; }
        TMV_INLINE int nhi() const { return itsnhi; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE int diagstep() const
        { return _diagmajor ? 1 : stepi()+stepj(); }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && !_diagmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && !_diagmajor && stepi() == 1); }
        TMV_INLINE bool isdm() const 
        { return _diagmajor || (!_rowmajor && !_colmajor && diagstep() == 1); }

    private :

        T* itsm;
        const size_t itscs;
        const size_t itsrs;
        const int itsnlo;
        const int itsnhi;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;

    }; // BandMatrixView


    //-------------------------------------------------------------------------

    //
    // Special Creators:
    //
    // View raw memory as a BandMatrix.
    //
    //   BandMatrixViewOf(T* m, cs, rs, lo, hi, stor)
    //   BandMatrixViewOf(T* m, cs, rs, lo, hi, stepi, stepj)
    //

    template <class T>
    inline BandMatrixView<T> BandMatrixViewOf(
        T* m, size_t cs, size_t rs, int lo, int hi, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(lo < int(cs));
        TMVAssert(hi < int(rs));
        const int stepi = (
            stor == RowMajor ? lo+hi :
            stor == ColMajor ? 1 :
            rs >= cs ? -int(cs)+1 : -int(rs) );
        const int stepj = (
            stor == RowMajor ? 1 :
            stor == ColMajor ? lo+hi :
            rs >= cs ? int(cs) : int(rs)+1 );
        T* m0 = (stor == DiagMajor) ? m - lo*stepi : m;
        return BandMatrixView<T>(m0,cs,rs,lo,hi,stepi,stepj);
    }

    template <class T>
    inline ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, size_t cs, size_t rs, int lo, int hi, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor || stor == DiagMajor);
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(lo < int(cs));
        TMVAssert(hi < int(rs));
        const int stepi = (
            stor == RowMajor ? lo+hi :
            stor == ColMajor ? 1 :
            rs >= cs ? -int(cs)+1 : -int(rs) );
        const int stepj = (
            stor == RowMajor ? 1 :
            stor == ColMajor ? lo+hi :
            rs >= cs ? int(cs) : int(rs)+1 );
        const T* m0 = (stor == DiagMajor) ? m - lo*stepi : m;
        return ConstBandMatrixView<T>(m0,cs,rs,lo,hi,stepi,stepj);
    }

    template <class T>
    TMV_INLINE BandMatrixView<T> BandMatrixViewOf(
        T* m, size_t cs, size_t rs, int lo, int hi, int stepi, int stepj)
    {
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(lo < int(cs));
        TMVAssert(hi < int(rs));
        return BandMatrixView<T>(m,cs,rs,lo,hi,stepi,stepj); 
    }

    template <class T>
    TMV_INLINE ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, size_t cs, size_t rs, int lo, int hi, int stepi, int stepj)
    {
        TMVAssert(cs > 0);
        TMVAssert(rs > 0);
        TMVAssert(lo < int(cs));
        TMVAssert(hi < int(rs));
        return ConstBandMatrixView<T>(m,cs,rs,lo,hi,stepi,stepj); 
    }


    //
    // Swap
    //

    template <class T, int A0, int A1>
    TMV_INLINE void Swap(BandMatrix<T,A0,A1>& m1, BandMatrix<T,A0,A1>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, int A>
    TMV_INLINE void Swap(
        BaseMatrix_Band_Mutable<M>& m1, BandMatrixView<T,A> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, int A>
    TMV_INLINE void Swap(
        BandMatrixView<T,A> m1, BaseMatrix_Band_Mutable<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int A1, int A2>
    TMV_INLINE void Swap(BandMatrixView<T,A1> m1, BandMatrixView<T,A2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int A0, int A1>
    TMV_INLINE typename BandMatrix<T,A0,A1>::conjugate_type Conjugate(
        BandMatrix<T,A0,A1>& m)
    { return m.conjugate(); }
    template <class T, int A>
    TMV_INLINE typename BandMatrixView<T,A>::conjugate_type Conjugate(
        BandMatrixView<T,A> m)
    { return m.conjugate(); }

    template <class T, int A0, int A1>
    TMV_INLINE typename BandMatrix<T,A0,A1>::transpose_type Transpose(
        BandMatrix<T,A0,A1>& m)
    { return m.transpose(); }
    template <class T, int A>
    TMV_INLINE typename BandMatrixView<T,A>::transpose_type Transpose(
        BandMatrixView<T,A> m)
    { return m.transpose(); }

    template <class T, int A0, int A1>
    TMV_INLINE typename BandMatrix<T,A0,A1>::adjoint_type Adjoint(
        BandMatrix<T,A0,A1>& m)
    { return m.adjoint(); }
    template <class T, int A>
    TMV_INLINE typename BandMatrixView<T,A>::adjoint_type Adjoint(
        BandMatrixView<T,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class T, int A0, int A1>
    inline std::string TMV_Text(const BandMatrix<T,A0,A1>& m)
    {
        const int A = A0 | A1;
        std::ostringstream s;
        s << "BandMatrix<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int A>
    inline std::string TMV_Text(const ConstBandMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "ConstBandMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int A>
    inline std::string TMV_Text(const BandMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "BandMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
