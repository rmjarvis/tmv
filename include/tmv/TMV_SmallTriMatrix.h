

//---------------------------------------------------------------------------
//
// This file defines the TMV SmallTriMatrix class.
//
// Constructors:
//
//    As with the comments in TMV_TriMatrix.h, I omit the Upper or Lower
//    before TriMatrix in the below constructors, since they are the same
//    form for each.
//
//    SmallTriMatrix<T,n,A>()
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SmallTriMatrix<T,n,A>(T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    SmallTriMatrix<T,n,A>(const Matrix<T>& m)
//    SmallTriMatrix<T,n,A>(const SmallTriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//

#ifndef TMV_SmallTriMatrix_H
#define TMV_SmallTriMatrix_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_VIt.h"
#include "TMV_MIt.h"
#include "TMV_Array.h"

namespace tmv {

    template <class T, ptrdiff_t N, int A0>
    struct Traits<SmallUpperTriMatrix<T,N,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::rowmajor ? 0 : ColMajor ) |
                ( Attrib<A0>::unitdiag ? 0 : NonUnitDiag ) )};
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
                !Attrib<A>::noalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallUpperTriMatrix<T,N,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallUpperTriMatrix<T,N,A0> copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = 0 };
        enum { _nhi = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = _rowmajor ? N : 1 };
        enum { _stepj = _rowmajor ? 1 : N };
        enum { _diagstep = IntTraits<N>::Sp1 };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = false };
        enum { _dt = A & AllDiagType };
        enum { twoSi = isreal ? ptrdiff_t(_stepi) : ptrdiff_t(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? ptrdiff_t(_stepj) : ptrdiff_t(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
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
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { ndAr = (Ar & ~AllDiagType) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { ndnmAr = (ndAr & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAr> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,Ar> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmAr> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAr> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;

        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstUpperTriMatrixView<T> const_xview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,1,_stepj,cmA>
            const_cmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,1,rmA>
            const_rmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallUpperTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            const_offdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,unitA>
            const_unitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            const_nonunitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,ndA>
            const_unknowndiag_type;
        typedef ConstSmallUpperTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,A> nonconst_type;

        typedef typename TriRefHelper<T,false,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colAr> col_sub_type;
        typedef VectorView<T,rowAr> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAr> diag_sub_type;

        typedef UpperTriMatrixView<T,Ar> subtrimatrix_type;
        typedef UpperTriMatrixView<T,nmAr> subtrimatrix_step_type;
        typedef MatrixView<T,ndAr> submatrix_type;
        typedef MatrixView<T,ndnmAr> submatrix_step_type;
        typedef VectorView<T,vecAr> subvector_type;

        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,A> view_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,fstyleA> fview_type;
        typedef UpperTriMatrixView<T> xview_type;
        typedef SmallUpperTriMatrixView<T,N,1,_stepj,cmA> cmview_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,1,rmA> rmview_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallLowerTriMatrixView<T,N,_stepj,_stepi,trA> transpose_type;
        typedef SmallLowerTriMatrixView<T,N,_stepj,_stepi,adjA> adjoint_type;

        typedef SmallUpperTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            offdiag_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,unitA>
            unitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            nonunitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,ndA>
            unknowndiag_type;
        typedef SmallUpperTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,An> noalias_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,An|CheckAlias> alias_type;
    };

    template <class T, ptrdiff_t N, int A>
    class SmallUpperTriMatrix : 
        public BaseMatrix_Tri_Mutable<SmallUpperTriMatrix<T,N,A> >
    {
    public:

        typedef SmallUpperTriMatrix<T,N,A> type;
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
        enum { _dt = Traits<type>::_dt };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE_ND SmallUpperTriMatrix()
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
#ifdef TMV_EXTRA_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::constr_value());
#endif
        }

        explicit SmallUpperTriMatrix(T x)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            this->setAllTo(x);
        }

        SmallUpperTriMatrix(const type& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            this->noAlias() = m2;
        }

        template <class M2>
        SmallUpperTriMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert(M2::_calc || assignable);
            TMVAssert(m2.colsize() == N);
            TMVAssert(m2.rowsize() == N);
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        SmallUpperTriMatrix(const BaseMatrix_Tri<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(M2::_upper);
            TMVAssert(m2.size() == N);
            typename Traits<type>::noalias_type na = this->noAlias();
            Maybe<_unit && !M2::_unit>::unitview(m2).assignTo(na);
        }

        template <class M2>
        SmallUpperTriMatrix(const BaseMatrix_Diag<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVAssert(m2.size() == N);
            this->setZero();
            this->diag().noAlias() = m2.calc().diag();
        }

        TMV_INLINE_ND ~SmallUpperTriMatrix()
        {
#ifdef TMV_EXTRA_DEBUG
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

        T cref(ptrdiff_t i, ptrdiff_t j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                itsm[i*stepi() + j*stepj()]);
        }

        reference ref(ptrdiff_t i, ptrdiff_t j)
        {
            return TriRefHelper<T,false,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }

        TMV_INLINE ptrdiff_t size() const { return N; }
        ptrdiff_t nElements() const { return N*(N+1)/2; }
        TMV_INLINE ptrdiff_t stepi() const { return _stepi; }
        TMV_INLINE ptrdiff_t stepj() const { return _stepj; }
        TMV_INLINE DiagType dt() const { return static_cast<DiagType>(_dt); }
        TMV_INLINE bool isunit() const { return _unit; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    protected :

        StackArray<T,N*N> itsm;

    }; // SmallUpperTriMatrix

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A0>
    struct Traits<ConstSmallUpperTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::noalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != Unknown && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != Unknown && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallUpperTriMatrixView<T,N,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = 0 };
        enum { _nhi = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? Unknown : A & AllDiagType };
        enum { twoSi = isreal ? ptrdiff_t(_stepi) : ptrdiff_t(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? ptrdiff_t(_stepj) : ptrdiff_t(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        enum { known = N != Unknown };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_unit ? UnitDiag : _nonunit ? NonUnitDiag : UnknownDiag) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename TypeSelect<known, 
                SmallUpperTriMatrix<T,N,copyA|NoAlias>,
                UpperTriMatrix<T,copyA|NoDivider> >::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
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
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { ndAr = (Ar & ~AllDiagType) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { ndnmAr = (ndAr & ~AllStorageType) };

        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAr> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,Ar> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmAr> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAr> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;

        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstUpperTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,1,_stepj,cmA>
            const_cmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,1,rmA>
            const_rmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallUpperTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            const_offdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,unitA>
            const_unitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            const_nonunitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,ndA>
            const_unknowndiag_type;
        typedef ConstSmallUpperTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,A> nonconst_type;

        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
    };

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    class ConstSmallUpperTriMatrixView :
        public BaseMatrix_Tri<ConstSmallUpperTriMatrixView<T,N,Si,Sj,A> >
    {
    public:

        typedef ConstSmallUpperTriMatrixView<T,N,Si,Sj,A> type;

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
        enum { _dt = Traits<type>::_dt };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //
        
        ConstSmallUpperTriMatrixView(
            const T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallUpperTriMatrixView(const T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallUpperTriMatrixView(const T* m, ptrdiff_t s, ptrdiff_t si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallUpperTriMatrixView(const T* m, ptrdiff_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != Unknown); TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallUpperTriMatrixView(const T* m) :
            itsm(m), itss(N), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != Unknown);
            TMVStaticAssert(Si != Unknown); TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallUpperTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstSmallUpperTriMatrixView(
            const ConstUpperTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstSmallUpperTriMatrixView(
            const UpperTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
        ConstSmallUpperTriMatrixView(
            const ConstSmallUpperTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
        ConstSmallUpperTriMatrixView(
            const SmallUpperTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstSmallUpperTriMatrixView() 
        {
#ifdef TMV_EXTRA_DEBUG
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

        T cref(ptrdiff_t i, ptrdiff_t j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        TMV_INLINE ptrdiff_t colsize() const { return itss; }
        TMV_INLINE ptrdiff_t rowsize() const { return itss; }
        TMV_INLINE ptrdiff_t size() const { return itss; }
        ptrdiff_t nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE ptrdiff_t stepi() const { return itssi; }
        TMV_INLINE ptrdiff_t stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<ptrdiff_t>(itsdt)); }
        TMV_INLINE bool isunit() const { return dt() == UnitDiag; }
        TMV_INLINE bool isrm() const
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        const T* itsm;
        const CheckedInt<N> itss;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // ConstSmallUpperTriMatrixView

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A0>
    struct Traits<SmallUpperTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::noalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != Unknown && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != Unknown && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallUpperTriMatrixView<T,N,Si,Sj,A0> type;
        typedef ConstSmallUpperTriMatrixView<T,N,Si,Sj,A> calc_type;
        typedef const type& eval_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = 0 };
        enum { _nhi = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? Unknown : A & AllDiagType };
        enum { twoSi = isreal ? ptrdiff_t(_stepi) : ptrdiff_t(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? ptrdiff_t(_stepj) : ptrdiff_t(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        enum { known = N != Unknown };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_unit ? UnitDiag : _nonunit ? NonUnitDiag : UnknownDiag) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename TypeSelect<known, 
                SmallUpperTriMatrix<T,N,copyA|NoAlias>,
                UpperTriMatrix<T,copyA|NoDivider> >::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
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
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { ndAr = (Ar & ~AllDiagType) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { ndnmAr = (ndAr & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAr> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,Ar> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmAr> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAr> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;

        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstUpperTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,1,_stepj,cmA>
            const_cmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,1,rmA>
            const_rmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallUpperTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            const_offdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,unitA>
            const_unitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            const_nonunitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,ndA>
            const_unknowndiag_type;
        typedef ConstSmallUpperTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,A> nonconst_type;

        typedef typename TriRefHelper<T,_conj,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colAr> col_sub_type;
        typedef VectorView<T,rowAr> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAr> diag_sub_type;

        typedef UpperTriMatrixView<T,Ar> subtrimatrix_type;
        typedef UpperTriMatrixView<T,nmAr> subtrimatrix_step_type;
        typedef MatrixView<T,ndAr> submatrix_type;
        typedef MatrixView<T,ndnmAr> submatrix_step_type;
        typedef VectorView<T,vecAr> subvector_type;

        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,A> view_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,fstyleA> fview_type;
        typedef UpperTriMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef SmallUpperTriMatrixView<T,N,1,_stepj,cmA> cmview_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,1,rmA> rmview_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallLowerTriMatrixView<T,N,_stepj,_stepi,trA> transpose_type;
        typedef SmallLowerTriMatrixView<T,N,_stepj,_stepi,adjA> adjoint_type;

        typedef SmallUpperTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            offdiag_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,unitA>
            unitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            nonunitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,ndA>
            unknowndiag_type;
        typedef SmallUpperTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,An> noalias_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,An|CheckAlias> alias_type;
    };

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    class SmallUpperTriMatrixView :
        public BaseMatrix_Tri_Mutable<SmallUpperTriMatrixView<T,N,Si,Sj,A> >
    {
    public:

        typedef SmallUpperTriMatrixView<T,N,Si,Sj,A> type;
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
        enum { _dt = Traits<type>::_dt };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        SmallUpperTriMatrixView(T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallUpperTriMatrixView(T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != Unknown); 
        }

        SmallUpperTriMatrixView(T* m, ptrdiff_t s, ptrdiff_t si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        SmallUpperTriMatrixView(T* m, ptrdiff_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != Unknown); TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        SmallUpperTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        SmallUpperTriMatrixView(
            UpperTriMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
        SmallUpperTriMatrixView(
            SmallUpperTriMatrixView<T,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~SmallUpperTriMatrixView() 
        {
#ifdef TMV_EXTRA_DEBUG
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

        T cref(ptrdiff_t i, ptrdiff_t j) const
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        reference ref(ptrdiff_t i, ptrdiff_t j)
        {
            return TriRefHelper<T,_conj,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }


        TMV_INLINE ptrdiff_t colsize() const { return itss; }
        TMV_INLINE ptrdiff_t rowsize() const { return itss; }
        TMV_INLINE ptrdiff_t size() const { return itss; }
        ptrdiff_t nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE ptrdiff_t stepi() const { return itssi; }
        TMV_INLINE ptrdiff_t stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<ptrdiff_t>(itsdt)); }
        TMV_INLINE bool isunit() const { return dt() == UnitDiag; }
        TMV_INLINE bool isrm() const
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        T* itsm;
        const CheckedInt<N> itss;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // SmallUpperTriMatrixView

    template <class T, ptrdiff_t N, int A0>
    struct Traits<SmallLowerTriMatrix<T,N,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::rowmajor ? 0 : ColMajor ) |
                ( Attrib<A0>::unitdiag ? 0 : NonUnitDiag ) )};
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
                !Attrib<A>::noalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallLowerTriMatrix<T,N,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallLowerTriMatrix<T,N,A0> copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _nhi = 0 };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = _rowmajor ? N : 1 };
        enum { _stepj = _rowmajor ? 1 : N };
        enum { _diagstep = IntTraits<N>::Sp1 };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = false };
        enum { _dt = A & AllDiagType };
        enum { twoSi = isreal ? ptrdiff_t(_stepi) : ptrdiff_t(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? ptrdiff_t(_stepj) : ptrdiff_t(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
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
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { ndAr = (Ar & ~AllDiagType) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { ndnmAr = (ndAr & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAr> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,Ar> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmAr> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAr> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;

        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstLowerTriMatrixView<T> const_xview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,1,_stepj,cmA>
            const_cmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,1,rmA>
            const_rmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallLowerTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            const_offdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,unitA>
            const_unitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            const_nonunitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,ndA>
            const_unknowndiag_type;
        typedef ConstSmallLowerTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,A> nonconst_type;

        typedef typename TriRefHelper<T,false,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colAr> col_sub_type;
        typedef VectorView<T,rowAr> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAr> diag_sub_type;

        typedef LowerTriMatrixView<T,Ar> subtrimatrix_type;
        typedef LowerTriMatrixView<T,nmAr> subtrimatrix_step_type;
        typedef MatrixView<T,ndAr> submatrix_type;
        typedef MatrixView<T,ndnmAr> submatrix_step_type;
        typedef VectorView<T,vecAr> subvector_type;

        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,A> view_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,fstyleA> fview_type;
        typedef LowerTriMatrixView<T> xview_type;
        typedef SmallLowerTriMatrixView<T,N,1,_stepj,cmA> cmview_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,1,rmA> rmview_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallUpperTriMatrixView<T,N,_stepj,_stepi,trA> transpose_type;
        typedef SmallUpperTriMatrixView<T,N,_stepj,_stepi,adjA> adjoint_type;

        typedef SmallLowerTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            offdiag_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,unitA>
            unitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            nonunitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,ndA>
            unknowndiag_type;
        typedef SmallLowerTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,An> noalias_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,An|CheckAlias> alias_type;
    };

    template <class T, ptrdiff_t N, int A>
    class SmallLowerTriMatrix : 
        public BaseMatrix_Tri_Mutable<SmallLowerTriMatrix<T,N,A> >
    {
    public:

        typedef SmallLowerTriMatrix<T,N,A> type;
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
        enum { _dt = Traits<type>::_dt };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE_ND SmallLowerTriMatrix()
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
#ifdef TMV_EXTRA_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(
                Traits<T>::constr_value());
#endif
        }

        explicit SmallLowerTriMatrix(T x)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            this->setAllTo(x);
        }

        SmallLowerTriMatrix(const type& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            this->noAlias() = m2;
        }

        template <class M2>
        SmallLowerTriMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert(M2::_calc || assignable);
            TMVAssert(m2.colsize() == N);
            TMVAssert(m2.rowsize() == N);
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        SmallLowerTriMatrix(const BaseMatrix_Tri<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(M2::_lower);
            TMVAssert(m2.size() == N);
            typename Traits<type>::noalias_type na = this->noAlias();
            Maybe<_unit && !M2::_unit>::unitview(m2).assignTo(na);
        }

        template <class M2>
        SmallLowerTriMatrix(const BaseMatrix_Diag<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVAssert(m2.size() == N);
            this->setZero();
            this->diag().noAlias() = m2.calc().diag();
        }

        TMV_INLINE_ND ~SmallLowerTriMatrix()
        {
#ifdef TMV_EXTRA_DEBUG
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

        T cref(ptrdiff_t i, ptrdiff_t j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                itsm[i*stepi() + j*stepj()]);
        }

        reference ref(ptrdiff_t i, ptrdiff_t j)
        {
            return TriRefHelper<T,false,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }

        TMV_INLINE ptrdiff_t size() const { return N; }
        ptrdiff_t nElements() const { return N*(N+1)/2; }
        TMV_INLINE ptrdiff_t stepi() const { return _stepi; }
        TMV_INLINE ptrdiff_t stepj() const { return _stepj; }
        TMV_INLINE DiagType dt() const { return static_cast<DiagType>(_dt); }
        TMV_INLINE bool isunit() const { return _unit; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    protected :

        StackArray<T,N*N> itsm;

    }; // SmallLowerTriMatrix

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A0>
    struct Traits<ConstSmallLowerTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::noalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != Unknown && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != Unknown && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallLowerTriMatrixView<T,N,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _nhi = 0 };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? Unknown : A & AllDiagType };
        enum { twoSi = isreal ? ptrdiff_t(_stepi) : ptrdiff_t(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? ptrdiff_t(_stepj) : ptrdiff_t(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        enum { known = N != Unknown };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_unit ? UnitDiag : _nonunit ? NonUnitDiag : UnknownDiag) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename TypeSelect<known, 
                SmallLowerTriMatrix<T,N,copyA|NoAlias>,
                LowerTriMatrix<T,copyA|NoDivider> >::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
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
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { ndAr = (Ar & ~AllDiagType) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { ndnmAr = (ndAr & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAr> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,Ar> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmAr> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAr> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;

        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstLowerTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,1,_stepj,cmA>
            const_cmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,1,rmA>
            const_rmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallLowerTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            const_offdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,unitA>
            const_unitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            const_nonunitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,ndA>
            const_unknowndiag_type;
        typedef ConstSmallLowerTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,A> nonconst_type;

        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
    };

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    class ConstSmallLowerTriMatrixView :
        public BaseMatrix_Tri<ConstSmallLowerTriMatrixView<T,N,Si,Sj,A> >
    {
    public:

        typedef ConstSmallLowerTriMatrixView<T,N,Si,Sj,A> type;

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
        enum { _dt = Traits<type>::_dt };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //
        
        ConstSmallLowerTriMatrixView(
            const T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallLowerTriMatrixView(const T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallLowerTriMatrixView(const T* m, ptrdiff_t s, ptrdiff_t si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallLowerTriMatrixView(const T* m, ptrdiff_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != Unknown); TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        ConstSmallLowerTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstSmallLowerTriMatrixView(
            const ConstLowerTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstSmallLowerTriMatrixView(
            const LowerTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
        ConstSmallLowerTriMatrixView(
            const ConstSmallLowerTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
        ConstSmallLowerTriMatrixView(
            const SmallLowerTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstSmallLowerTriMatrixView() 
        {
#ifdef TMV_EXTRA_DEBUG
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

        T cref(ptrdiff_t i, ptrdiff_t j) const 
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        TMV_INLINE ptrdiff_t colsize() const { return itss; }
        TMV_INLINE ptrdiff_t rowsize() const { return itss; }
        TMV_INLINE ptrdiff_t size() const { return itss; }
        ptrdiff_t nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE ptrdiff_t stepi() const { return itssi; }
        TMV_INLINE ptrdiff_t stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<ptrdiff_t>(itsdt)); }
        TMV_INLINE bool isunit() const { return dt() == UnitDiag; }
        TMV_INLINE bool isrm() const
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        const T* itsm;
        const CheckedInt<N> itss;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // ConstSmallLowerTriMatrixView

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A0>
    struct Traits<SmallLowerTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::noalias &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != Unknown && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != Unknown && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallLowerTriMatrixView<T,N,Si,Sj,A0> type;
        typedef ConstSmallLowerTriMatrixView<T,N,Si,Sj,A> calc_type;
        typedef const type& eval_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = IntTraits2<IntTraits<N>::Sm1,0>::max };
        enum { _nhi = 0 };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? Unknown : A & AllDiagType };
        enum { twoSi = isreal ? ptrdiff_t(_stepi) : ptrdiff_t(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? ptrdiff_t(_stepj) : ptrdiff_t(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        enum { known = N != Unknown };
        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_unit ? UnitDiag : _nonunit ? NonUnitDiag : UnknownDiag) |
                (_fort ? FortranStyle : CStyle) ) };
        typedef typename TypeSelect<known, 
                SmallLowerTriMatrix<T,N,copyA|NoAlias>,
                LowerTriMatrix<T,copyA|NoDivider> >::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
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
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAr = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAr = vecAr | (_colmajor ? Unit : 0) };
        enum { rowAr = vecAr | (_rowmajor ? Unit : 0) };
        enum { diagAr = vecAr | (_diagstep == 1 ? Unit : 0) };
        enum { ndAr = (Ar & ~AllDiagType) };
        enum { nmAr = (Ar & ~AllStorageType) };
        enum { ndnmAr = (ndAr & ~AllStorageType) };
        enum { An = ndA & ~CheckAlias };

        typedef ConstVectorView<T,colAr> const_col_sub_type;
        typedef ConstVectorView<T,rowAr> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAr> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,Ar> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmAr> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAr> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAr> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAr> const_subvector_type;

        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstLowerTriMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,1,_stepj,cmA>
            const_cmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,1,rmA>
            const_rmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepj,_stepi,adjA>
            const_adjoint_type;

        typedef ConstSmallLowerTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            const_offdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,unitA>
            const_unitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            const_nonunitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,ndA>
            const_unknowndiag_type;
        typedef ConstSmallLowerTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,A> nonconst_type;

        typedef typename TriRefHelper<T,_conj,_nonunit>::reference reference;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;

        typedef VectorView<T,colAr> col_sub_type;
        typedef VectorView<T,rowAr> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAr> diag_sub_type;

        typedef LowerTriMatrixView<T,Ar> subtrimatrix_type;
        typedef LowerTriMatrixView<T,nmAr> subtrimatrix_step_type;
        typedef MatrixView<T,ndAr> submatrix_type;
        typedef MatrixView<T,ndnmAr> submatrix_step_type;
        typedef VectorView<T,vecAr> subvector_type;

        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,A> view_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,fstyleA> fview_type;
        typedef LowerTriMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef SmallLowerTriMatrixView<T,N,1,_stepj,cmA> cmview_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,1,rmA> rmview_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallUpperTriMatrixView<T,N,_stepj,_stepi,trA> transpose_type;
        typedef SmallUpperTriMatrixView<T,N,_stepj,_stepi,adjA> adjoint_type;

        typedef SmallLowerTriMatrixView<T,Nm1,_stepi,_stepj,nonunitA>
            offdiag_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,unitA>
            unitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,nonunitA>
            nonunitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,ndA>
            unknowndiag_type;
        typedef SmallLowerTriMatrixView<real_type,N,twoSi,twoSj,twosA>
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,nonconjA>
            nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,An> noalias_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,An|CheckAlias> alias_type;
    };

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    class SmallLowerTriMatrixView :
        public BaseMatrix_Tri_Mutable<SmallLowerTriMatrixView<T,N,Si,Sj,A> >
    {
    public:

        typedef SmallLowerTriMatrixView<T,N,Si,Sj,A> type;
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
        enum { _dt = Traits<type>::_dt };
        enum { _unit = Traits<type>::_unit };
        enum { _nonunit = Traits<type>::_nonunit };
        enum { _unknowndiag = Traits<type>::_unknowndiag };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        SmallLowerTriMatrixView(T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallLowerTriMatrixView(T* m, ptrdiff_t s, ptrdiff_t si, ptrdiff_t sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != Unknown); 
        }

        SmallLowerTriMatrixView(T* m, ptrdiff_t s, ptrdiff_t si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        SmallLowerTriMatrixView(T* m, ptrdiff_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != Unknown); TMVStaticAssert(Sj != Unknown); 
            TMVStaticAssert(_dt != Unknown); 
        }

        SmallLowerTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), 
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        SmallLowerTriMatrixView(
            LowerTriMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
        SmallLowerTriMatrixView(
            SmallLowerTriMatrixView<T,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~SmallLowerTriMatrixView() 
        {
#ifdef TMV_EXTRA_DEBUG
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

        T cref(ptrdiff_t i, ptrdiff_t j) const
        {
            return (
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        reference ref(ptrdiff_t i, ptrdiff_t j)
        {
            return TriRefHelper<T,_conj,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }
 

        TMV_INLINE ptrdiff_t colsize() const { return itss; }
        TMV_INLINE ptrdiff_t rowsize() const { return itss; }
        TMV_INLINE ptrdiff_t size() const { return itss; }
        ptrdiff_t nElements() const { return itss*(itss+1)/2; }
        TMV_INLINE ptrdiff_t stepi() const { return itssi; }
        TMV_INLINE ptrdiff_t stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<ptrdiff_t>(itsdt)); }
        TMV_INLINE bool isunit() const { return dt() == UnitDiag; }
        TMV_INLINE bool isrm() const
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        T* itsm;
        const CheckedInt<N> itss;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;
        const CheckedInt<_dt> itsdt;

    }; // SmallLowerTriMatrixView


    //
    // Swap
    //

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, class M>
    TMV_INLINE void Swap(
        BaseMatrix_Tri_Mutable<M>& m1,
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m2)
    { DoSwap(m1,m2); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, class M>
    TMV_INLINE void Swap(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m1,
        BaseMatrix_Tri_Mutable<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    TMV_INLINE void Swap(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> m1,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, int A2>
    TMV_INLINE void Swap(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> m1,
        UpperTriMatrixView<T,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, ptrdiff_t N, int A1, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    TMV_INLINE void Swap(
        UpperTriMatrixView<T,A1> m1,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, class M>
    TMV_INLINE void Swap(
        BaseMatrix_Tri_Mutable<M>& m1,
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A, class M>
    TMV_INLINE void Swap(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m1,
        BaseMatrix_Tri_Mutable<M>& m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    TMV_INLINE void Swap(
        SmallLowerTriMatrixView<T,N,Si1,Sj1,A1> m1,
        SmallLowerTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, ptrdiff_t N, ptrdiff_t Si1, ptrdiff_t Sj1, int A1, int A2>
    TMV_INLINE void Swap(
        SmallLowerTriMatrixView<T,N,Si1,Sj1,A1> m1,
        LowerTriMatrixView<T,A2> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, ptrdiff_t N, int A1, ptrdiff_t Si2, ptrdiff_t Sj2, int A2>
    TMV_INLINE void Swap(
        LowerTriMatrixView<T,A1> m1,
        SmallLowerTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { Swap(m1.transpose(),m2.transpose()); }


    //
    // Conjugate, Transpose, Adjoint
    //
    
    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallUpperTriMatrix<T,N,A>::conjugate_type Conjugate(
        SmallUpperTriMatrix<T,N,A>& m)
    { return m.conjugate(); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE typename SmallUpperTriMatrixView<T,N,Si,Sj,A>::conjugate_type Conjugate(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallUpperTriMatrix<T,N,A>::transpose_type Transpose(
        SmallUpperTriMatrix<T,N,A>& m)
    { return m.transpose(); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE typename SmallUpperTriMatrixView<T,N,Si,Sj,A>::transpose_type Transpose(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallUpperTriMatrix<T,N,A>::adjoint_type Adjoint(
        SmallUpperTriMatrix<T,N,A>& m)
    { return m.adjoint(); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE typename SmallUpperTriMatrixView<T,N,Si,Sj,A>::adjoint_type Adjoint(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    { return m.adjoint(); }

    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallLowerTriMatrix<T,N,A>::conjugate_type Conjugate(
        SmallLowerTriMatrix<T,N,A>& m)
    { return m.conjugate(); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE typename SmallLowerTriMatrixView<T,N,Si,Sj,A>::conjugate_type Conjugate(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallLowerTriMatrix<T,N,A>::transpose_type Transpose(
        SmallLowerTriMatrix<T,N,A>& m)
    { return m.transpose(); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE typename SmallLowerTriMatrixView<T,N,Si,Sj,A>::transpose_type Transpose(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallLowerTriMatrix<T,N,A>::adjoint_type Adjoint(
        SmallLowerTriMatrix<T,N,A>& m)
    { return m.adjoint(); }
    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE typename SmallLowerTriMatrixView<T,N,Si,Sj,A>::adjoint_type Adjoint(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

    template <class T, ptrdiff_t N, int A>
    inline std::string TMV_Text(const SmallUpperTriMatrix<T,N,A>& m)
    {
        std::ostringstream s;
        s << "SmallUpperTriMatrix<"<<TMV_Text(T());
        s << ","<<N;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    inline std::string TMV_Text(
        const ConstSmallUpperTriMatrixView<T,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "ConstSmallUpperTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<Si>::text()<<","<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    inline std::string TMV_Text(
        const SmallUpperTriMatrixView<T,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "SmallUpperTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<Si>::text()<<","<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, int A>
    inline std::string TMV_Text(const SmallLowerTriMatrix<T,N,A>& m)
    {
        std::ostringstream s;
        s << "SmallLowerTriMatrix<"<<TMV_Text(T());
        s << ","<<N;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    inline std::string TMV_Text(
        const ConstSmallLowerTriMatrixView<T,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "ConstSmallLowerTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<Si>::text()<<","<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    inline std::string TMV_Text(
        const SmallLowerTriMatrixView<T,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "SmallLowerTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<Si>::text()<<","<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }


} // namespace tmv

#endif
