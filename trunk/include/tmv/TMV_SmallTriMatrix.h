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
//    SmallTriMatrix<T,n,A>(T* vv)
//    SmallTriMatrix<T,n,A>(const std::vector<T>& vv)
//        Makes a Triangular Matrix with column size = row size = n
//        which copies the values from vv.
//
//    SmallTriMatrix<T,n,A>(const Matrix<T>& m)
//    SmallTriMatrix<T,n,A>(const SmallTriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//

#ifndef TMV_SmallTriMatrix_H
#define TMV_SmallTriMatrix_H

#include <vector>
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"

namespace tmv {

    template <class T, int N, int A0, int A1, int A2>
    struct Traits<SmallUpperTriMatrix<T,N,A0,A1,A2> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A012 = A0 | A1 | A2 };
        enum { A = (A012 & ~NoDivider) | (
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
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallUpperTriMatrix<T,N,A0,A1,A2> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallUpperTriMatrix<T,N,A012> copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = _rowmajor ? N : 1 };
        enum { _stepj = _rowmajor ? 1 : N };
        enum { _diagstep = IntTraits<N>::Sp1 };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = false };
        enum { _dt = A & AllDiagType };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
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
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { ndAn = (An & ~AllDiagType) };
        enum { nmAn = (An & ~AllStorageType) };
        enum { ndnmAn = (ndAn & ~AllStorageType) };

        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAn> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,An> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmAn> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAn> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;

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

        typedef VectorView<T,colAn> col_sub_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAn> diag_sub_type;

        typedef UpperTriMatrixView<T,An> subtrimatrix_type;
        typedef UpperTriMatrixView<T,nmAn> subtrimatrix_step_type;
        typedef MatrixView<T,ndAn> submatrix_type;
        typedef MatrixView<T,ndnmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;

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
    };

    template <class T, int N, int A0, int A1, int A2>
    class SmallUpperTriMatrix : 
        public BaseMatrix_Tri_Mutable<SmallUpperTriMatrix<T,N,A0,A1,A2> >
    {
    public:

        typedef SmallUpperTriMatrix<T,N,A0,A1,A2> type;
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
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
        }

        explicit SmallUpperTriMatrix(T x)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            this->setAllTo(x);
        }

        explicit SmallUpperTriMatrix(const T* vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            SmallVectorView<T,N*N,1> lv(itsm);
            ConstSmallVectorView<T,N*N,1>(vv).newAssignTo(lv);
        }

        explicit SmallUpperTriMatrix(const std::vector<T>& vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVAssert(vv.size() == N*N);
            SmallVectorView<T,N*N,1> lv(itsm);
            ConstSmallVectorView<T,N*N,1>(&vv[0]).newAssignTo(lv);
        }

        SmallUpperTriMatrix(const type& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            m2.newAssignTo(*this);
        }

        template <class M2>
        SmallUpperTriMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert((
                (M2::_calc && ShapeTraits<M2::_shape>::upper) || assignable));
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
            Maybe<_unit && !M2::_unit>::unitview(m2).newAssignTo(*this);
        }

        template <class M2>
        SmallUpperTriMatrix(const BaseMatrix_Diag<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVAssert(m2.size() == N);
            typename type::diag_type d = this->diag();
            this->setZero();
            m2.calc().diag().newAssignTo(d);
        }

        TMV_INLINE_ND ~SmallUpperTriMatrix()
        {
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(999));
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
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                itsm[i*stepi() + j*stepj()]);
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,false,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }

        TMV_INLINE size_t size() const { return N; }
        TMV_INLINE int nElements() const { return N*(N+1)/2; }
        TMV_INLINE int stepi() const { return _stepi; }
        TMV_INLINE int stepj() const { return _stepj; }
        TMV_INLINE DiagType dt() const { return static_cast<DiagType>(_dt); }
        TMV_INLINE bool isunit() const { return _unit; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    protected :

        StackArray<T,N*N> itsm;

    }; // SmallUpperTriMatrix

    template <class T, int N, int Si, int Sj, int A0>
    struct Traits<ConstSmallUpperTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != UNKNOWN && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != UNKNOWN && Sj != 1 && Attrib<A>::rowmajor ) )};
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
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? UNKNOWN : A & AllDiagType };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        typedef typename MCopyHelper<T,_shape,N,N,_rowmajor,_fort>::type 
            copy_type;

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
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { ndAn = (An & ~AllDiagType) };
        enum { nmAn = (An & ~AllStorageType) };
        enum { ndnmAn = (ndAn & ~AllStorageType) };

        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAn> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,An> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmAn> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAn> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;

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
    };

    template <class T, int N, int Si, int Sj, int A>
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
            const T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallUpperTriMatrixView(const T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        ConstSmallUpperTriMatrixView(const T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        ConstSmallUpperTriMatrixView(const T* m, size_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        ConstSmallUpperTriMatrixView(const T* m) :
            itsm(m), itss(N), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
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
        }

        template <int A2>
        ConstSmallUpperTriMatrixView(
            const UpperTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstSmallUpperTriMatrixView(
            const ConstSmallUpperTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstSmallUpperTriMatrixView(
            const SmallUpperTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ~ConstSmallUpperTriMatrixView() {
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
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        TMV_INLINE int nElements() const { return int(itss)*int(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
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

    template <class T, int N, int Si, int Sj, int A0>
    struct Traits<SmallUpperTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != UNKNOWN && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != UNKNOWN && Sj != 1 && Attrib<A>::rowmajor ) )};
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
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? UNKNOWN : A & AllDiagType };
        enum { _shape = _unit ? UnitUpperTri : UpperTri };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        typedef typename MCopyHelper<T,_shape,N,N,_rowmajor,_fort>::type 
            copy_type;

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
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { ndAn = (An & ~AllDiagType) };
        enum { nmAn = (An & ~AllStorageType) };
        enum { ndnmAn = (ndAn & ~AllStorageType) };

        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAn> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,An> const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,nmAn> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAn> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;

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

        typedef VectorView<T,colAn> col_sub_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAn> diag_sub_type;

        typedef UpperTriMatrixView<T,An> subtrimatrix_type;
        typedef UpperTriMatrixView<T,nmAn> subtrimatrix_step_type;
        typedef MatrixView<T,ndAn> submatrix_type;
        typedef MatrixView<T,ndnmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;

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
    };

    template <class T, int N, int Si, int Sj, int A>
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

        SmallUpperTriMatrixView(T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallUpperTriMatrixView(T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        SmallUpperTriMatrixView(T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        SmallUpperTriMatrixView(T* m, size_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
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
        }

        template <int N2, int Si2, int Sj2, int A2>
        SmallUpperTriMatrixView(
            SmallUpperTriMatrixView<T,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ~SmallUpperTriMatrixView() {
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
                (isunit() && i==j ) ? T(1) :
                (i>j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,_conj,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }


        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        TMV_INLINE int nElements() const { return int(itss)*int(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
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

    template <class T, int N, int A0, int A1, int A2>
    struct Traits<SmallLowerTriMatrix<T,N,A0,A1,A2> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A012 = A0 | A1 | A2 };
        enum { A = (A012 & ~NoDivider) | (
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
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallLowerTriMatrix<T,N,A0,A1,A2> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallLowerTriMatrix<T,N,A012> copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = _rowmajor ? N : 1 };
        enum { _stepj = _rowmajor ? 1 : N };
        enum { _diagstep = IntTraits<N>::Sp1 };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = false };
        enum { _dt = A & AllDiagType };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
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
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { ndAn = (An & ~AllDiagType) };
        enum { nmAn = (An & ~AllStorageType) };
        enum { ndnmAn = (ndAn & ~AllStorageType) };

        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAn> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,An> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmAn> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAn> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;

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

        typedef VectorView<T,colAn> col_sub_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAn> diag_sub_type;

        typedef LowerTriMatrixView<T,An> subtrimatrix_type;
        typedef LowerTriMatrixView<T,nmAn> subtrimatrix_step_type;
        typedef MatrixView<T,ndAn> submatrix_type;
        typedef MatrixView<T,ndnmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;

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
    };

    template <class T, int N, int A0, int A1, int A2>
    class SmallLowerTriMatrix : 
        public BaseMatrix_Tri_Mutable<SmallLowerTriMatrix<T,N,A0,A1,A2> >
    {
    public:

        typedef SmallLowerTriMatrix<T,N,A0,A1,A2> type;
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
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
        }

        explicit SmallLowerTriMatrix(T x)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            this->setAllTo(x);
        }

        explicit SmallLowerTriMatrix(const T* vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            SmallVectorView<T,N*N,1> lv(ptr());
            ConstSmallVectorView<T,N*N,1>(vv).newAssignTo(lv);
        }

        explicit SmallLowerTriMatrix(const std::vector<T>& vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVAssert(vv.size() == N*N);
            SmallVectorView<T,N*N,1> lv(ptr());
            ConstSmallVectorView<T,N*N,1>(&vv[0]).newAssignTo(lv);
        }

        SmallLowerTriMatrix(const type& m2) 
        {
            TMVStaticAssert(N>=0);
            m2.newAssignTo(*this);
        }

        template <class M2>
        SmallLowerTriMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert((
                (M2::_calc && ShapeTraits<M2::_shape>::lower) || assignable));
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
            Maybe<_unit && !M2::_unit>::unitview(m2).newAssignTo(*this);
        }

        template <class M2>
        SmallLowerTriMatrix(const BaseMatrix_Diag<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N>=0);
            TMVAssert(m2.size() == N);
            typename type::diag_type d = this->diag();
            this->setZero();
            m2.calc().diag().newAssignTo(d);
        }

        TMV_INLINE_ND ~SmallLowerTriMatrix()
        {
#ifdef TMV_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
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
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                itsm[i*stepi() + j*stepj()]);
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,false,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }

        TMV_INLINE size_t size() const { return N; }
        TMV_INLINE int nElements() const { return N*(N+1)/2; }
        TMV_INLINE int stepi() const { return _stepi; }
        TMV_INLINE int stepj() const { return _stepj; }
        TMV_INLINE DiagType dt() const { return static_cast<DiagType>(_dt); }
        TMV_INLINE bool isunit() const { return _unit; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    protected :

        StackArray<T,N*N> itsm;

    }; // SmallLowerTriMatrix

    template <class T, int N, int Si, int Sj, int A0>
    struct Traits<ConstSmallLowerTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != UNKNOWN && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != UNKNOWN && Sj != 1 && Attrib<A>::rowmajor ) )};
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
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? UNKNOWN : A & AllDiagType };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        typedef typename MCopyHelper<T,_shape,N,N,_rowmajor,_fort>::type 
            copy_type;

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
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { ndAn = (An & ~AllDiagType) };
        enum { nmAn = (An & ~AllStorageType) };
        enum { ndnmAn = (ndAn & ~AllStorageType) };

        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAn> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,An> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmAn> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAn> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;

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
    };

    template <class T, int N, int Si, int Sj, int A>
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
            const T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallLowerTriMatrixView(const T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        ConstSmallLowerTriMatrixView(const T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        ConstSmallLowerTriMatrixView(const T* m, size_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
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
        }

        template <int A2>
        ConstSmallLowerTriMatrixView(
            const LowerTriMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstSmallLowerTriMatrixView(
            const ConstSmallLowerTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int Si2, int Sj2, int A2>
        ConstSmallLowerTriMatrixView(
            const SmallLowerTriMatrixView<T,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ~ConstSmallLowerTriMatrixView() {
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
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        TMV_INLINE int nElements() const { return int(itss)*int(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
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

    template <class T, int N, int Si, int Sj, int A0>
    struct Traits<SmallLowerTriMatrixView<T,N,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::nodivider &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != UNKNOWN && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != UNKNOWN && Sj != 1 && Attrib<A>::rowmajor ) )};
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
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unitdiag };
        enum { _nonunit = Attrib<A>::nonunitdiag };
        enum { _unknowndiag = !(_unit || _nonunit) };
        enum { _dt = _unknowndiag ? UNKNOWN : A & AllDiagType };
        enum { _shape = _unit ? UnitLowerTri : LowerTri };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { Nm1 = IntTraits<N>::Sm1 };

        typedef typename MCopyHelper<T,_shape,N,N,_rowmajor,_fort>::type 
            copy_type;

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
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { ndAn = (An & ~AllDiagType) };
        enum { nmAn = (An & ~AllStorageType) };
        enum { ndnmAn = (ndAn & ~AllStorageType) };

        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,N,_diagstep,vecA> const_diag_type;
        typedef ConstVectorView<T,vecAn> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,An> const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,nmAn> const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,ndAn> const_submatrix_type;
        typedef ConstMatrixView<T,ndnmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;

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

        typedef VectorView<T,colAn> col_sub_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,N,_diagstep,vecA> diag_type;
        typedef VectorView<T,vecAn> diag_sub_type;

        typedef LowerTriMatrixView<T,An> subtrimatrix_type;
        typedef LowerTriMatrixView<T,nmAn> subtrimatrix_step_type;
        typedef MatrixView<T,ndAn> submatrix_type;
        typedef MatrixView<T,ndnmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;

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
    };

    template <class T, int N, int Si, int Sj, int A>
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

        SmallLowerTriMatrixView(T* m, size_t s, int si, int sj, DiagType dt) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(dt)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallLowerTriMatrixView(T* m, size_t s, int si, int sj) :
            itsm(m), itss(s), itssi(si), itssj(sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        SmallLowerTriMatrixView(T* m, size_t s, int si) :
            itsm(m), itss(s), itssi(si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
        }

        SmallLowerTriMatrixView(T* m, size_t s) :
            itsm(m), itss(s), itssi(Si), itssj(Sj), itsdt(_dt)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
            TMVStaticAssert(_dt != UNKNOWN); 
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
        }

        template <int N2, int Si2, int Sj2, int A2>
        SmallLowerTriMatrixView(
            SmallLowerTriMatrixView<T,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itssi(m2.stepi()), itssj(m2.stepj()), itsdt(m2.dt())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ~SmallLowerTriMatrixView() {
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
                (isunit() && i==j ) ? T(1) :
                (i<j) ? T(0) :
                DoConj<_conj>(itsm[i*stepi() + j*stepj()]));
        }

        reference ref(int i, int j)
        {
            return TriRefHelper<T,_conj,_nonunit>::makeRef(
                isunit() && i==j, itsm[i*stepi() + j*stepj()] ); 
        }
 

        TMV_INLINE size_t colsize() const { return itss; }
        TMV_INLINE size_t rowsize() const { return itss; }
        TMV_INLINE size_t size() const { return itss; }
        TMV_INLINE int nElements() const { return int(itss)*int(itss+1)/2; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE DiagType dt() const 
        { return static_cast<DiagType>(static_cast<int>(itsdt)); }
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

    template <class T, int N, int Si, int Sj, int A, class M>
    static TMV_INLINE void Swap(
        BaseMatrix_Tri_Mutable<M>& m1,
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m2)
    { DoSwap(m1,m2); }
    template <class T, int N, int Si, int Sj, int A, class M>
    static TMV_INLINE void Swap(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m1,
        BaseMatrix_Tri_Mutable<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int N, int Si1, int Sj1, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> m1,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int N, int Si1, int Sj1, int A1, int A2>
    static TMV_INLINE void Swap(
        SmallUpperTriMatrixView<T,N,Si1,Sj1,A1> m1,
        UpperTriMatrixView<T,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int N, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        UpperTriMatrixView<T,A1> m1,
        SmallUpperTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }

    template <class T, int N, int Si, int Sj, int A, class M>
    static TMV_INLINE void Swap(
        BaseMatrix_Tri_Mutable<M>& m1,
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, int Si, int Sj, int A, class M>
    static TMV_INLINE void Swap(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m1,
        BaseMatrix_Tri_Mutable<M>& m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, int Si1, int Sj1, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        SmallLowerTriMatrixView<T,N,Si1,Sj1,A1> m1,
        SmallLowerTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, int Si1, int Sj1, int A1, int A2>
    static TMV_INLINE void Swap(
        SmallLowerTriMatrixView<T,N,Si1,Sj1,A1> m1,
        LowerTriMatrixView<T,A2> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        LowerTriMatrixView<T,A1> m1,
        SmallLowerTriMatrixView<T,N,Si2,Sj2,A2> m2)
    { Swap(m1.transpose(),m2.transpose()); }


    //
    // Conjugate, Transpose, Adjoint
    //
    
    template <class T, int N, int A0, int A1, int A2>
    static TMV_INLINE typename SmallUpperTriMatrix<T,N,A0,A1,A2>::conjugate_type Conjugate(
        SmallUpperTriMatrix<T,N,A0,A1,A2>& m)
    { return m.conjugate(); }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallUpperTriMatrixView<T,N,Si,Sj,A>::conjugate_type Conjugate(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, int N, int A0, int A1, int A2>
    static TMV_INLINE typename SmallUpperTriMatrix<T,N,A0,A1,A2>::transpose_type Transpose(
        SmallUpperTriMatrix<T,N,A0,A1,A2>& m)
    { return m.transpose(); }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallUpperTriMatrixView<T,N,Si,Sj,A>::transpose_type Transpose(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, int N, int A0, int A1, int A2>
    static TMV_INLINE typename SmallUpperTriMatrix<T,N,A0,A1,A2>::adjoint_type Adjoint(
        SmallUpperTriMatrix<T,N,A0,A1,A2>& m)
    { return m.adjoint(); }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallUpperTriMatrixView<T,N,Si,Sj,A>::adjoint_type Adjoint(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    { return m.adjoint(); }

    template <class T, int N, int A0, int A1, int A2>
    static TMV_INLINE typename SmallLowerTriMatrix<T,N,A0,A1,A2>::conjugate_type Conjugate(
        SmallLowerTriMatrix<T,N,A0,A1,A2>& m)
    { return m.conjugate(); }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallLowerTriMatrixView<T,N,Si,Sj,A>::conjugate_type Conjugate(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, int N, int A0, int A1, int A2>
    static TMV_INLINE typename SmallLowerTriMatrix<T,N,A0,A1,A2>::transpose_type Transpose(
        SmallLowerTriMatrix<T,N,A0,A1,A2>& m)
    { return m.transpose(); }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallLowerTriMatrixView<T,N,Si,Sj,A>::transpose_type Transpose(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, int N, int A0, int A1, int A2>
    static TMV_INLINE typename SmallLowerTriMatrix<T,N,A0,A1,A2>::adjoint_type Adjoint(
        SmallLowerTriMatrix<T,N,A0,A1,A2>& m)
    { return m.adjoint(); }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallLowerTriMatrixView<T,N,Si,Sj,A>::adjoint_type Adjoint(
        SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_DEBUG
    template <class T, int N, int A0, int A1, int A2>
    static inline std::string TMV_Text(
        const SmallUpperTriMatrix<T,N,A0,A1,A2>& m)
    {
        const int A = A0 | A1 | A2;
        std::ostringstream s;
        s << "SmallUpperTriMatrix<"<<TMV_Text(T());
        s << ","<<N;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int N, int Si, int Sj, int A>
    static inline std::string TMV_Text(
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

    template <class T, int N, int Si, int Sj, int A>
    static inline std::string TMV_Text(
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

    template <class T, int N, int A0, int A1, int A2>
    static inline std::string TMV_Text(
        const SmallLowerTriMatrix<T,N,A0,A1,A2>& m)
    {
        const int A = A0 | A1 | A2;
        std::ostringstream s;
        s << "SmallLowerTriMatrix<"<<TMV_Text(T());
        s << ","<<N;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.size()<<","<<m.stepi()<<","<<m.stepj();
        s << ","<<(m.isunit()?"Unit":"NonUnit")<<")";
        return s.str();
    }

    template <class T, int N, int Si, int Sj, int A>
    static inline std::string TMV_Text(
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

    template <class T, int N, int Si, int Sj, int A>
    static inline std::string TMV_Text(
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
#endif


} // namespace tmv

#endif
