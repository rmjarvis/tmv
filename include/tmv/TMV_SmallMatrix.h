

//---------------------------------------------------------------------------
//
// This file defines the TMV SmallMatrix class.
//
// A SmallMatrix is a matrix whose size is known at compile time.
// This allows for a lot of optimizations by the compiler in implementing
// the various operations on it. 
// 
// Since the instantiation needs to know the size of the matrix, all 
// operations on the SmallMatrix are done inline, rather than precompiled.
// This gives another speed boost in most cases, since the compiler can
// avoid implementing a function call in most cases.
//
// Finally, another advantage is that we allocate the data on the stack
// rather than the heap, so we avoid new and delete calls as well.
// However, stack sizes are usually limited to hundreds of KB, 
// (mine is 8 MB), so we set a maximum size of 1KB for each SmallMatrix.
// (For double, this means up to 11x11 will be allocated on the stack.)
// Any bigger than that, and the performance drop from using the
// heap is pretty irrelevant.  (See TMV_Array.h to change this.)
// 
// One drawback of using a SmallMatrix is that it does not do any
// alias checking in the aritmetic statements.  So a statement like
// m = m2 * m will not produce a correct answer. 
// Normally this is a feature, since the alias checks can be a 
// significant fraction of the calculation time for small matrices.
// 
// You can workaround this when necessary by explicitly making a copy.
// The easiest way is with the .copy() method.  e.g. m = m2 * m.copy().
// 
//
// Constructors:
//
//    SmallMatrix<T,M,N,stor,I>()
//        Makes a SmallMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SmallMatrix<T,M,N,stor,I>(T x)
//        Makes a SmallMatrix of size n with all values = x
//
//    SmallMatrix<T,M,N,stor,I>(const GenMatrix<T>& m)
//        Make a SmallMatrix which copies the elements of m.
//
// Special Creators:
//
//    SmallMatrixView RowVectorViewOf(Vector& v)
//    SmallMatrixView RowVectorViewOf(VectorView v)
//    ConstSmallMatrixView RowVectorViewOf(const Vector& v)
//    ConstSmallMatrixView RowVectorViewOf(const ConstVectorView v)
//        Returns a 1xn MatrixView with v in the (only) row.
//        Unlike creating a Matrix with RowVector(v), this uses the same
//        data storage as the actual Vector v.
//        For a const Vector or a ConstVectorView, this returns a 
//        ConstMatrixView.
//
//    SmallMatrixView ColVectorViewOf(Vector& v)
//    SmallMatrixView ColVectorViewOf(VectorView v)
//    ConstSmallMatrixView ColVectorViewOf(const Vector& v)
//    ConstSmallMatrixView ColVectorViewOf(const ConstVectorView& v)
//        Returns an nx1 MatrixView with v in the (only) column.


#ifndef TMV_SmallMatrix_H
#define TMV_SmallMatrix_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"
#include "TMV_Matrix.h"

namespace tmv {

    template <class M> class LUD;
    template <class T> class InstLUD;
    template <class M> class QRD;
    template <class T> class InstQRD;
    template <class M> class QRPD;
    template <class T> class InstQRPD;
    template <class M> class SVD;
    template <class T> class InstSVD;

    //
    // SmallMatrix
    //

    template <class T, int M, int N, int A0>
    struct Traits<SmallMatrix<T,M,N,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                (Attrib<A0>::rowmajor ? 0 : ColMajor) |
                ( Attrib<A0>::withdivider ? 0 : NoDivider ) ) };
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor) &&
                (Attrib<A>::rowmajor != int(Attrib<A>::colmajor)) &&
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                !Attrib<A>::noalias &&
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

        typedef SmallMatrix<T,M,N,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallMatrix<T,M,N,A0> copy_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _nlo = IntTraits2<M-1,0>::max };
        enum { _nhi = IntTraits2<N-1,0>::max };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = _rowmajor ? N : 1 };
        enum { _stepj = _rowmajor ? 1 : M };
        enum { _diagstep = _stepi + _stepj };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = true };
        enum { _linsize = IntTraits2<M,N>::prod };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { minMN = IntTraits2<M,N>::min };

        enum { _hasdivider = Attrib<A>::withdivider };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<T>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<T>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRPD<T>& , QRPD<copy_type> >::type qrpd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstSVD<T>& , SVD<copy_type> >::type svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { colpairA = ndA & ~RowMajor };
        enum { rowpairA = ndA & ~ColMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA };
        enum { twosA = isreal ? int(nonconjA) : (nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        enum { xx = Unknown }; // for brevity
        typedef ConstSmallVectorView<T,M,_stepi,colA> const_col_type;
        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,rowA> const_row_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstVectorView<T,diagAr> const_diag_sub_type;

        typedef ConstMatrixView<T,Ar> const_submatrix_type;
        typedef ConstMatrixView<T,nmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,xx,colpairA>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,xx,_stepj,rowpairA>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,xx,_stepi,_stepj,ndA>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,xx,N,_stepi,_stepj,ndA>
            const_rowrange_type;

        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,ndA>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,cmA> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,rmA> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_lowertri_type;

        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecA|Unit)> 
            const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,ndA> nonconst_type;

        typedef T& reference;
        typedef CVIt<T,1,NonConj> const_linear_iterator;
        typedef VIt<T,1,NonConj> linear_iterator;
        typedef typename TypeSelect< _rowmajor , const_linear_iterator ,
                CRMIt<type> >::type const_rowmajor_iterator;
        typedef typename TypeSelect< _colmajor , const_linear_iterator ,
                CCMIt<type> >::type const_colmajor_iterator;
        typedef typename TypeSelect< _rowmajor , linear_iterator ,
                RMIt<type> >::type rowmajor_iterator;
        typedef typename TypeSelect< _colmajor , linear_iterator ,
                CMIt<type> >::type colmajor_iterator;

        typedef SmallVectorView<T,M,_stepi,colA> col_type;
        typedef VectorView<T,colAr> col_sub_type;
        typedef SmallVectorView<T,N,_stepj,rowA> row_type;
        typedef VectorView<T,rowAr> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,diagA> diag_type;
        typedef VectorView<T,diagAr> diag_sub_type;

        typedef MatrixView<T,Ar> submatrix_type;
        typedef MatrixView<T,nmAr> submatrix_step_type;
        typedef VectorView<T,vecAr> subvector_type;
        typedef SmallMatrixView<T,M,2,_stepi,xx,colpairA> colpair_type;
        typedef SmallMatrixView<T,2,N,xx,_stepj,rowpairA> rowpair_type;
        typedef SmallMatrixView<T,M,xx,_stepi,_stepj,ndA> colrange_type;
        typedef SmallMatrixView<T,xx,N,_stepi,_stepj,ndA> rowrange_type;

        typedef SmallMatrixView<T,M,N,_stepi,_stepj,ndA> view_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,fstyleA> fview_type;
        typedef MatrixView<T> xview_type;
        typedef SmallMatrixView<T,M,N,1,_stepj,cmA> cmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,1,rmA> rmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,trA> transpose_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,adjA> adjoint_type;

        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            unit_uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            unit_lowertri_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_lowertri_type;

        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,_linsize,1,(vecA|Unit)> linearview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,nonconjA> nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,An> noalias_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,An|CheckAlias> alias_type;
    };

    template <class T, int M, int N, int A>
    class SmallMatrix : 
        public BaseMatrix_Rec_Mutable<SmallMatrix<T,M,N,A> >,
        public MatrixDivHelper<SmallMatrix<T,M,N,A> >
    {
    public:

        typedef SmallMatrix<T,M,N,A> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        //typedef typename Traits<type>::lud_type lud_type;
        //typedef typename Traits<type>::qrd_type qrd_type;
        //typedef typename Traits<type>::qrpd_type qrpd_type;
        //typedef typename Traits<type>::svd_type svd_type;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _linsize = Traits<type>::_linsize };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        SmallMatrix() 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
#ifdef TMV_EXTRA_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        explicit SmallMatrix(T x) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            this->setAllTo(x);
        }

        SmallMatrix(const type& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            this->noAlias() = m2;
        }

        template <class M2>
        SmallMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            this->noAlias() = m2;
        }

        ~SmallMatrix()
        {
#ifdef TMV_EXTRA_DEBUG
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

        TMV_INLINE type& operator=(T x) 
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return itsm[i*stepi()+j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[i*stepi()+j*stepj()]; }

        TMV_INLINE int ls() const { return _linsize; }
        TMV_INLINE int colsize() const { return M; }
        TMV_INLINE int rowsize() const { return N; }
        TMV_INLINE int nElements() const { return _linsize; }
        TMV_INLINE int stepi() const { return _stepi; }
        TMV_INLINE int stepj() const { return _stepj; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    private:

        StackArray<T,_linsize> itsm;

    }; // SmallMatrix

    template <class T, int M, int N, int Si, int Sj, int A0>
    struct Traits<ConstSmallMatrixView<T,M,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) |
                ( Attrib<A0>::withdivider ? 0 : NoDivider ) ) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                !Attrib<A>::noalias &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider) ) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst &&
                    !Traits<real_type>::isinteger) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != Unknown && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != Unknown && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _nlo = IntTraits2<IntTraits<M>::Sm1,0>::max };
        enum { _nhi = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = (
                (Si == 1 && M != Unknown && Sj == M) || 
                (Sj == 1 && N != Unknown && Si == N) )};
        enum { _linsize = IntTraits2<M,N>::prod };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { minMN = IntTraits2<M,N>::min };

        enum { allknown = (M != Unknown && N != Unknown) };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename TypeSelect<allknown, 
                SmallMatrix<T,M,N,copyA|NoAlias>,
                Matrix<T,copyA|NoDivider> >::type copy_type;

        enum { _hasdivider = Attrib<A>::withdivider };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<T>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<T>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRPD<T>& , QRPD<copy_type> >::type qrpd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstSVD<T>& , SVD<copy_type> >::type svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { colpairA = ndA & ~RowMajor };
        enum { rowpairA = ndA & ~ColMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { nmAr = (Ar & ~AllStorageType) };

        enum { xx = Unknown }; // for brevity
        typedef ConstSmallVectorView<T,M,_stepi,colA> const_col_type;
        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,rowA> const_row_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstVectorView<T,diagAr> const_diag_sub_type;

        typedef ConstMatrixView<T,Ar> const_submatrix_type;
        typedef ConstMatrixView<T,nmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,xx,colpairA>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,xx,_stepj,rowpairA>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,xx,_stepi,_stepj,ndA>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,xx,N,_stepi,_stepj,ndA>
            const_rowrange_type;

        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,ndA>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,cmA> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,rmA> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_lowertri_type;

        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecA|Unit)> 
            const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,ndA> nonconst_type;

        typedef CVIt<T,1,_conj?Conj:NonConj> const_linear_iterator;
        typedef typename TypeSelect< _canlin && _rowmajor ,
                const_linear_iterator ,
                CRMIt<type> >::type const_rowmajor_iterator;
        typedef typename TypeSelect< _canlin && _colmajor ,
                const_linear_iterator ,
                CCMIt<type> >::type const_colmajor_iterator;
    };

    template <class T, int M, int N, int Si, int Sj, int A>
    class ConstSmallMatrixView :
        public BaseMatrix_Rec<ConstSmallMatrixView<T,M,N,Si,Sj,A> >,
        public MatrixDivHelper<ConstSmallMatrixView<T,M,N,Si,Sj,A> >
    {
    public:

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,A> type;
        //typedef typename Traits<type>::lud_type lud_type;
        //typedef typename Traits<type>::qrd_type qrd_type;
        //typedef typename Traits<type>::qrpd_type qrpd_type;
        //typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _linsize = Traits<type>::_linsize };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE ConstSmallMatrixView(
            const T* m, int cs, int rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE ConstSmallMatrixView(
            const T* m, int cs, int rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE ConstSmallMatrixView(const T* m, int cs, int rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { 
            TMVStaticAssert(Si != Unknown);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE ConstSmallMatrixView(const T* m, int cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != Unknown);
            TMVStaticAssert(Si != Unknown);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE ConstSmallMatrixView(const T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != Unknown);
            TMVStaticAssert(N != Unknown);
            TMVStaticAssert(Si != Unknown);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE ConstSmallMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE ConstSmallMatrixView(const ConstMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        TMV_INLINE ConstSmallMatrixView(const MatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstSmallMatrixView(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstSmallMatrixView(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~ConstSmallMatrixView() 
        {
#ifdef TMV_EXTRA_DEBUG
            itsm = 0; 
#endif
        }

    private:
        void operator=(const type& v2);
    public:


        // 
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE int ls() const { return itscs*itsrs; }
        TMV_INLINE int colsize() const { return itscs; }
        TMV_INLINE int rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        const T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstSmallMatrixView

    template <class T, int M, int N, int Si, int Sj, int A0>
    struct Traits<SmallMatrixView<T,M,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) |
                ( Attrib<A0>::withdivider ? 0 : NoDivider ) ) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                !Attrib<A>::noalias &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider) ) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst &&
                    !Traits<real_type>::isinteger) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != Unknown && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != Unknown && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallMatrixView<T,M,N,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _nlo = IntTraits2<IntTraits<M>::Sm1,0>::max };
        enum { _nhi = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = (
                (Si == 1 && M != Unknown && Sj == M) || 
                (Sj == 1 && N != Unknown && Si == N) )};
        enum { _linsize = IntTraits2<M,N>::prod };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { minMN = IntTraits2<M,N>::min };

        enum { allknown = (M != Unknown && N != Unknown) };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename TypeSelect<allknown, 
                SmallMatrix<T,M,N,copyA|NoAlias>,
                Matrix<T,copyA|NoDivider> >::type copy_type;

        enum { _hasdivider = Attrib<A>::withdivider };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<T>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<T>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRPD<T>& , QRPD<copy_type> >::type qrpd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstSVD<T>& , SVD<copy_type> >::type svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { colpairA = ndA & ~RowMajor };
        enum { rowpairA = ndA & ~ColMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        enum { xx = Unknown }; // for brevity
        typedef ConstSmallVectorView<T,M,_stepi,colA> const_col_type;
        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,rowA> const_row_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstVectorView<T,diagAr> const_diag_sub_type;

        typedef ConstMatrixView<T,Ar> const_submatrix_type;
        typedef ConstMatrixView<T,nmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,xx,colpairA>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,xx,_stepj,rowpairA>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,xx,_stepi,_stepj,ndA>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,xx,N,_stepi,_stepj,ndA>
            const_rowrange_type;

        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,ndA>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,cmA> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,rmA> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_lowertri_type;

        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecA|Unit)> 
            const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,ndA> nonconst_type;

        typedef typename AuxRef<T,_conj>::reference reference;
        typedef CVIt<T,1,_conj?Conj:NonConj> const_linear_iterator;
        typedef VIt<T,1,_conj?Conj:NonConj> linear_iterator;
        typedef typename TypeSelect< _canlin && _rowmajor ,
                const_linear_iterator ,
                CRMIt<type> >::type const_rowmajor_iterator;
        typedef typename TypeSelect< _canlin && _colmajor ,
                const_linear_iterator ,
                CCMIt<type> >::type const_colmajor_iterator;
        typedef typename TypeSelect< _canlin && _rowmajor , 
                linear_iterator ,
                RMIt<type> >::type rowmajor_iterator;
        typedef typename TypeSelect< _canlin && _colmajor ,
                linear_iterator ,
                CMIt<type> >::type colmajor_iterator;


        typedef SmallVectorView<T,M,_stepi,colA> col_type;
        typedef VectorView<T,colAr> col_sub_type;
        typedef SmallVectorView<T,N,_stepj,rowA> row_type;
        typedef VectorView<T,rowAr> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,diagA> diag_type;
        typedef VectorView<T,diagAr> diag_sub_type;

        typedef MatrixView<T,Ar> submatrix_type;
        typedef MatrixView<T,nmAr> submatrix_step_type;
        typedef VectorView<T,vecAr> subvector_type;
        typedef SmallMatrixView<T,M,2,_stepi,xx,colpairA> colpair_type;
        typedef SmallMatrixView<T,2,N,xx,_stepj,rowpairA> rowpair_type;
        typedef SmallMatrixView<T,M,xx,_stepi,_stepj,ndA> colrange_type;
        typedef SmallMatrixView<T,xx,N,_stepi,_stepj,ndA> rowrange_type;

        typedef SmallMatrixView<T,M,N,_stepi,_stepj,ndA> view_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,fstyleA> fview_type;
        typedef MatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef SmallMatrixView<T,M,N,1,_stepj,cmA> cmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,1,rmA> rmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,trA> transpose_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,adjA> adjoint_type;

        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            unit_uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            unit_lowertri_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_lowertri_type;

        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,_linsize,1,(vecA|Unit)> linearview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,nonconjA> nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,An> noalias_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,An|CheckAlias> alias_type;
    };

    template <class T, int M, int N, int Si, int Sj, int A>
    class SmallMatrixView :
        public BaseMatrix_Rec_Mutable<SmallMatrixView<T,M,N,Si,Sj,A> >,
        public MatrixDivHelper<SmallMatrixView<T,M,N,Si,Sj,A> >
    {
    public:

        typedef SmallMatrixView<T,M,N,Si,Sj,A> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;
        //typedef typename Traits<type>::lud_type lud_type;
        //typedef typename Traits<type>::qrd_type qrd_type;
        //typedef typename Traits<type>::qrpd_type qrpd_type;
        //typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _linsize = Traits<type>::_linsize };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE SmallMatrixView(T* m, int cs, int rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE SmallMatrixView(T* m, int cs, int rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE SmallMatrixView(T* m, int cs, int rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != Unknown);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE SmallMatrixView(T* m, int cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != Unknown);
            TMVStaticAssert(Si != Unknown); 
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE SmallMatrixView(T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != Unknown);
            TMVStaticAssert(N != Unknown);
            TMVStaticAssert(Si != Unknown);
            TMVStaticAssert(Sj != Unknown); 
        }

        TMV_INLINE SmallMatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE SmallMatrixView(MatrixView<T,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        TMV_INLINE SmallMatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~SmallMatrixView() 
        {
#ifdef TMV_EXTRA_DEBUG
            itsm = 0; 
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

        TMV_INLINE type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        reference ref(int i, int j)
        { return reference(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE int ls() const { return itscs*itsrs; }
        TMV_INLINE int colsize() const { return itscs; }
        TMV_INLINE int rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // SmallMatrixView


    // Special Creators: 
    //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
    //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage

    // A quick helper class to get the return types correct for
    // RowVectorView and ColVectorView.
    template <class V>
    struct VVO
    {
        typedef typename V::value_type T;
        enum { N = V::_size };
        enum { S = V::_step };
        enum { vecA = (
                ( V::_conj ? Conj : NonConj ) |
                ( V::_fort ? FortranStyle : CStyle ) |
                ( V::_checkalias ? CheckAlias : NoAlias ) ) };
        enum { colA = vecA | (S == 1 ? ColMajor : 0 ) };
        enum { rowA = vecA | (S == 1 ? RowMajor : 0 ) };
        typedef ConstSmallMatrixView<T,1,N,N,S,rowA> crv;
        typedef SmallMatrixView<T,1,N,N,S,rowA> rv;
        typedef ConstSmallMatrixView<T,N,1,S,N,colA> ccv;
        typedef SmallMatrixView<T,N,1,S,N,colA> cv;
    };

    template <class V>
    TMV_INLINE typename VVO<V>::crv RowVectorViewOf(const BaseVector_Calc<V>& v)
    { return typename VVO<V>::crv(v.cptr(),1,v.size(),v.size(),v.step()); }

    template <class V>
    TMV_INLINE typename VVO<V>::rv RowVectorViewOf(BaseVector_Mutable<V>& v)
    { return typename VVO<V>::rv(v.ptr(),1,v.size(),v.size(),v.step()); }

    template <class V>
    TMV_INLINE typename VVO<V>::ccv ColVectorViewOf(const BaseVector_Calc<V>& v)
    { return typename VVO<V>::ccv(v.cptr(),v.size(),1,v.step(),v.size()); }

    template <class V>
    TMV_INLINE typename VVO<V>::cv ColVectorViewOf(BaseVector_Mutable<V>& v)
    { return typename VVO<V>::cv(v.ptr(),v.size(),1,v.step(),v.size()); }

    // Repeat for VectorView and SmallVectorView so we can have them 
    // work without requiring the non-const reference.
    template <class T, int A>
    TMV_INLINE typename VVO<VectorView<T,A> >::rv RowVectorViewOf(
        VectorView<T,A> v)
    {
        return typename VVO<VectorView<T,A> >::rv(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int N, int S, int A>
    TMV_INLINE typename VVO<SmallVectorView<T,N,S,A> >::rv RowVectorViewOf(
        SmallVectorView<T,N,S,A> v)
    {
        return typename VVO<SmallVectorView<T,N,S,A> >::rv(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int S, int A>
    TMV_INLINE typename VVO<VectorView<T,A> >::cv ColVectorViewOf(
        VectorView<T,A> v)
    {
        return typename VVO<VectorView<T,A> >::cv(
            v.ptr(),v.size(),1,v.step(),v.size());
    }

    template <class T, int N, int S, int A>
    TMV_INLINE typename VVO<SmallVectorView<T,N,S,A> >::cv ColVectorViewOf(
        SmallVectorView<T,N,S,A> v)
    {
        return typename VVO<SmallVectorView<T,N,S,A> >::cv(
            v.ptr(),v.size(),1,v.step(),v.size());
    }


    //
    // Swap
    //

    template <class T, int M, int N, int Si, int Sj, int A, class MM>
    TMV_INLINE void Swap(
        BaseMatrix_Rec_Mutable<MM>& m1, SmallMatrixView<T,M,N,Si,Sj,A> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si, int Sj, int A, class MM>
    TMV_INLINE void Swap(
        SmallMatrixView<T,M,N,Si,Sj,A> m1, BaseMatrix_Rec_Mutable<MM>& m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, int A1, int Si2, int Sj2, int A2>
    TMV_INLINE void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,A1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, int A1, int A2>
    TMV_INLINE void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,A1> m1, MatrixView<T,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int A1, int Si2, int Sj2, int A2>
    TMV_INLINE void Swap(
        MatrixView<T,A1> m1, SmallMatrixView<T,M,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int M, int N, int A>
    TMV_INLINE typename SmallMatrix<T,M,N,A>::conjugate_type Conjugate(
        SmallMatrix<T,M,N,A>& m)
    { return m.conjugate(); }
    template <class T, int M, int N, int Si, int Sj, int A>
    TMV_INLINE typename SmallMatrixView<T,M,N,Si,Sj,A>::conjugate_type Conjugate(
        SmallMatrixView<T,M,N,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, int M, int N, int A>
    TMV_INLINE typename SmallMatrix<T,M,N,A>::transpose_type Transpose(
        SmallMatrix<T,M,N,A>& m)
    { return m.transpose(); }
    template <class T, int M, int N, int Si, int Sj, int A>
    TMV_INLINE typename SmallMatrixView<T,M,N,Si,Sj,A>::transpose_type Transpose(
        SmallMatrixView<T,M,N,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, int M, int N, int A>
    TMV_INLINE typename SmallMatrix<T,M,N,A>::adjoint_type Adjoint(
        SmallMatrix<T,M,N,A>& m)
    { return m.adjoint(); }
    template <class T, int M, int N, int Si, int Sj, int A>
    TMV_INLINE typename SmallMatrixView<T,M,N,Si,Sj,A>::adjoint_type Adjoint(
        SmallMatrixView<T,M,N,Si,Sj,A> m)
    { return m.adjoint(); }


    // 
    // TMV_Text
    //

    template <class T, int M, int N, int A>
    inline std::string TMV_Text(const SmallMatrix<T,M,N,A>& m)
    {
        std::ostringstream s;
        s << "SmallMatrix<"<<TMV_Text(T());
        s << ','<<M<<','<<N;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    inline std::string TMV_Text(
        const ConstSmallMatrixView<T,M,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "ConstSmallMatrixView<"<<TMV_Text(T());
        s << ','<<IntTraits<M>::text()<<','<<IntTraits<N>::text();
        s << ','<<IntTraits<Si>::text()<<','<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    inline std::string TMV_Text(
        const SmallMatrixView<T,M,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "SmallMatrixView<"<<TMV_Text(T());
        s << ','<<IntTraits<M>::text()<<','<<IntTraits<N>::text();
        s << ','<<IntTraits<Si>::text()<<','<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

} // namespace tmv

#endif
