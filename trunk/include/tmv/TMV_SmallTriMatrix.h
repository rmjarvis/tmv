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
//    SmallTriMatrix<T,n,dt,stor,I>()
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SmallTriMatrix<T,n,dt,stor,I>(T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    SmallTriMatrix<T,n,dt,stor,I>(T* vv)
//    SmallTriMatrix<T,n,dt,stor,I>(const std::vector<T>& vv)
//        Makes a Triangular Matrix with column size = row size = n
//        which copies the values from vv.
//
//    SmallTriMatrix<T,n,dt,stor,I>(const Matrix<T>& m)
//    SmallTriMatrix<T,n,dt,stor,I>(const SmallTriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//

#ifndef TMV_SmallTriMatrix_H
#define TMV_SmallTriMatrix_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_TriMatrix.h"
#include <vector>

namespace tmv {

    template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
    struct Traits<SmallUpperTriMatrix<T,N,D,S,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallUpperTriMatrix<T,N,D,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==RowMajor ? N : 1) };
        enum { _stepj = (S==RowMajor ? 1 : N) };
        enum { _diagstep = N+1 };
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
        typedef ConstSmallVectorView<T,N,_diagstep,false,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_sub_type;

        typedef ConstUpperTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_subtrimatrix_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_view_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,CStyle> 
            const_cview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,
                false,FortranStyle> const_fview_type;
        typedef ConstUpperTriMatrixView<T,D> const_xview_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag> const_xdview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,1,_stepj,false,I> 
            const_cmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,1,false,I> 
            const_rmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepj,_stepi,false,I> 
            const_transpose_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_offdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,_stepi,_stepj,
                false,I> const_unitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,
                false,I> const_nonunitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,
                false,I> const_unknowndiag_type;
        typedef ConstSmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            nonconst_type;

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

        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> view_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,CStyle> 
            cview_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,FortranStyle> 
            fview_type;
        typedef UpperTriMatrixView<T,D> xview_type;
        typedef UpperTriMatrixView<T,UnknownDiag> xdview_type;
        typedef SmallUpperTriMatrixView<T,N,D,1,_stepj,false,I> cmview_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,1,false,I> rmview_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            conjugate_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepj,_stepi,false,I> 
            transpose_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            adjoint_type;

        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            offdiag_type;
        typedef SmallUpperTriMatrixView<T,N,UnitDiag,_stepi,_stepj,false,I> 
            unitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,false,I> 
            nonunitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,false,I> 
            unknowndiag_type;
        typedef SmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            nonconj_type;
    };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
    class SmallUpperTriMatrix : 
        public BaseMatrix_Tri_Mutable<SmallUpperTriMatrix<T,N,D,S,I> >
    {
    public:

        typedef SmallUpperTriMatrix<T,N,D,S,I> type;
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

        inline SmallUpperTriMatrix(size_t n=N)
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
        }

        explicit inline SmallUpperTriMatrix(size_t n, T x)
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
            this->setAllTo(x);
        }

        explicit inline SmallUpperTriMatrix(size_t n, const T* vv) 
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            SmallVectorView<T,N*N,1> lv(itsm);
            ConstSmallVectorView<T,N*N,1>(vv).newAssignTo(lv);
        }

        explicit inline SmallUpperTriMatrix(size_t n, const std::vector<T>& vv) 
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
            TMVAssert(vv.size() == N*N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            SmallVectorView<T,N*N,1> lv(itsm);
            ConstSmallVectorView<T,N*N,1>(&vv[0]).newAssignTo(lv);
        }

        inline SmallUpperTriMatrix(const type& m2) 
        { 
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            m2.newAssignTo(*this);
        }

        template <class M2> 
        inline SmallUpperTriMatrix(const BaseMatrix<M2>& m2) 
        { 
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert((
                (M2::_calc && ShapeTraits<M2::_shape>::upper) || assignable));
            TMVAssert(m2.colsize() == N);
            TMVAssert(m2.rowsize() == N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        inline SmallUpperTriMatrix(const BaseMatrix_Tri<M2>& m2) 
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(M2::_upper);
            TMVAssert(m2.size() == N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            Maybe<_unit && !M2::_unit>::unitview(m2).newAssignTo(*this);
        }

        template <class M2> 
        inline SmallUpperTriMatrix(const BaseMatrix_Diag<M2>& m2) 
        { 
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(m2.size() == N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            typename type::diag_type d = this->diag();
            this->setZero();
            m2.calc().diag().newAssignTo(d);
        }

        inline ~SmallUpperTriMatrix()
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

        inline size_t size() const { return N; }
        inline int stepi() const { return _stepi; }
        inline int stepj() const { return _stepj; }
        inline DiagType dt() const { return D; }
        inline bool isunit() const { return D == UnitDiag; }
        inline bool isconj() const { return false; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline StorageType stor() const { return S; }

    protected :

        StackArray<T,N*N> itsm;

    }; // SmallUpperTriMatrix

    template <class T, int N, DiagType D, StorageType S>
    class SmallUpperTriMatrixF : 
        public SmallUpperTriMatrix<T,N,D,S,FortranStyle>
    {
    public:

        typedef SmallUpperTriMatrixF<T,N,D,S> type;
        typedef SmallUpperTriMatrix<T,N,D,S,FortranStyle> mtype;

        inline SmallUpperTriMatrixF() : mtype() {}
        explicit inline SmallUpperTriMatrixF(T x) : mtype(x) {}
        explicit inline SmallUpperTriMatrixF(const T* vv) : mtype(vv) {}
        explicit inline SmallUpperTriMatrixF(const std::vector<T>& vv) : 
            mtype(vv) {}
        template <class M2> 
        inline SmallUpperTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        inline SmallUpperTriMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        inline SmallUpperTriMatrixF(const BaseMatrix_Tri<M2>& m2) :
            mtype(m2) {}
        template <class M2>
        inline SmallUpperTriMatrixF(const BaseMatrix_Rec<M2>& m2) :
            mtype(m2) {}
        inline ~SmallUpperTriMatrixF() {}

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

    }; // SmallUpperTriMatrixF

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
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

        // In case N == UNKNOWN
        typedef typename MCopyHelper<T,UpperTri,N,N,_rowmajor,_fort>::type 
            copy_type;

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

        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,I> 
            const_view_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,
                FortranStyle> const_fview_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,1,_stepj,C,I> 
            const_cmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,1,C,I> 
            const_rmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_offdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,_stepi,_stepj,C,I> 
            const_unitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,C,I> 
            const_nonunitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknowndiag_type;
        typedef ConstSmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,I> nonconst_type;
    };

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class ConstSmallUpperTriMatrixView :
        public BaseMatrix_Tri<ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
    public:

        typedef ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;

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
        inline ConstSmallUpperTriMatrixView(
            const T* m, size_t s, bool u, int si, int sj
        ) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline ConstSmallUpperTriMatrixView(
            const T* m, size_t s, bool u, int si
        ) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstSmallUpperTriMatrixView(const T* m, size_t s) :
            itsm(m), itss(s), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline ConstSmallUpperTriMatrixView(const T* m) :
            itsm(m), itss(N), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline ConstSmallUpperTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixView(
            const ConstUpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixView(
            const UpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixView(
            const ConstSmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixView(
            const SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~ConstSmallUpperTriMatrixView() {
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
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        inline bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

#ifdef TMV_DEBUG
        const T* itsm;
#else
        const T*const itsm;
#endif
        const CheckedInt<N> itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstSmallUpperTriMatrixView

    template <class T, int N, DiagType D, int Si, int Sj, bool C>
    class ConstSmallUpperTriMatrixViewF :
        public ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef ConstSmallUpperTriMatrixViewF<T,N,D,Si,Sj,C> type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

        inline ConstSmallUpperTriMatrixViewF(
            const T* m, size_t s, bool u, int si, int sj
        ) : 
            mtype(m,s,si,sj) {}
        inline ConstSmallUpperTriMatrixViewF(const T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline ConstSmallUpperTriMatrixViewF(const T* m, size_t s) :
            mtype(m,s) {}
        inline ConstSmallUpperTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixViewF(
            const ConstUpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixViewF(
            const UpperTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixViewF(
            const ConstSmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallUpperTriMatrixViewF(
            const SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        inline ~ConstSmallUpperTriMatrixViewF() {}

    private :
        inline void operator=(const type& m2);

    }; // ConstSmallUpperTriMatrixViewF

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
        typedef const ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
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

        typedef typename MCopyHelper<T,UpperTri,N,N,_rowmajor,_fort>::type 
            copy_type;

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

        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,I> 
            const_view_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,
                FortranStyle> const_fview_type;
        typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,1,_stepj,C,I> 
            const_cmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,1,C,I> 
            const_rmview_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,1,C,I> 
            const_offdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,_stepi,1,C,I> 
            const_unitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,1,C,I> 
            const_nonunitdiag_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,1,C,I> 
            const_unknowndiag_type;
        typedef ConstSmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,I> nonconst_type;

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

        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,I> view_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,CStyle> 
            cview_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,C,FortranStyle> 
            fview_type;
        typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
        typedef UpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            xdview_type;
        typedef SmallUpperTriMatrixView<T,N,D,1,_stepj,C,I> cmview_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,1,C,I> rmview_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            conjugate_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepj,_stepi,C,I> 
            transpose_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            adjoint_type;

        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            offdiag_type;
        typedef SmallUpperTriMatrixView<T,N,UnitDiag,_stepi,_stepj,C,I> 
            unitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,C,I> 
            nonunitdiag_type;
        typedef SmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,C,I> 
            unknowndiag_type;
        typedef SmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            nonconj_type;
    };

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class SmallUpperTriMatrixView :
        public BaseMatrix_Tri_Mutable<SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
    public:

        typedef SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
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

        inline SmallUpperTriMatrixView(T* m, size_t s, bool u, int si, int sj) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline SmallUpperTriMatrixView(T* m, size_t s, bool u, int si) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline SmallUpperTriMatrixView(T* m, size_t s) :
            itsm(m), itss(s), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline SmallUpperTriMatrixView(T* m) :
            itsm(m), itss(N), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline SmallUpperTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), 
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallUpperTriMatrixView(
            UpperTriMatrixView<T,D2,Si2,Sj2,C,I2> m2
        ) :
            itsm(m2.ptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallUpperTriMatrixView(
            SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~SmallUpperTriMatrixView() {
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
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        inline bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

#ifdef TMV_DEBUG
        T* itsm;
#else
        T*const itsm;
#endif
        const CheckedInt<N> itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // SmallUpperTriMatrixView

    template <class T, int N, DiagType D, int Si, int Sj, bool C>
    class SmallUpperTriMatrixViewF :
        public SmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef SmallUpperTriMatrixViewF<T,N,D,Si,Sj,C> type;
        typedef SmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

        inline SmallUpperTriMatrixViewF(T* m, size_t s, bool u, int si, int sj) :
            mtype(m,s,si,sj) {}
        inline SmallUpperTriMatrixViewF(T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline SmallUpperTriMatrixViewF(T* m, size_t s) : mtype(m,s) {}
        inline SmallUpperTriMatrixViewF(T* m) : mtype(m) {}
        inline SmallUpperTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallUpperTriMatrixViewF(
            UpperTriMatrixView<T,D2,Si2,Sj2,C,I2> m2
        ) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallUpperTriMatrixViewF(
            SmallUpperTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2) : mtype(m2) {}
        inline ~SmallUpperTriMatrixViewF() {}

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

    }; // SmallUpperTriMatrixViewF


    template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
    struct Traits<SmallLowerTriMatrix<T,N,D,S,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallLowerTriMatrix<T,N,D,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==RowMajor ? N : 1) };
        enum { _stepj = (S==RowMajor ? 1 : N) };
        enum { _diagstep = N+1 };
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
        typedef ConstSmallVectorView<T,N,_diagstep,false,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_sub_type;

        typedef ConstLowerTriMatrixView<T,D,_stepi,_stepj,false,I> 
            const_subtrimatrix_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> 
            const_subtrimatrix_step_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_view_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,CStyle> 
            const_cview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,
                FortranStyle> const_fview_type;
        typedef ConstLowerTriMatrixView<T,D> const_xview_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag> const_xdview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,1,_stepj,false,I> 
            const_cmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,1,false,I> 
            const_rmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepj,_stepi,false,I> 
            const_transpose_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_offdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnitDiag,_stepi,_stepj,
                false,I> const_unitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,
                false,I> const_nonunitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,
                false,I> const_unknowndiag_type;
        typedef ConstSmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            nonconst_type;

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

        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> view_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,CStyle> 
            cview_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,
                FortranStyle> fview_type;
        typedef LowerTriMatrixView<T,D> xview_type;
        typedef LowerTriMatrixView<T,UnknownDiag> xdview_type;
        typedef SmallLowerTriMatrixView<T,N,D,1,_stepj,false,I> cmview_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,1,false,I> rmview_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            conjugate_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepj,_stepi,false,I> 
            transpose_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            adjoint_type;

        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            offdiag_type;
        typedef SmallLowerTriMatrixView<T,N,UnitDiag,_stepi,_stepj,false,I> 
            unitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,false,I> 
            nonunitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,false,I> 
            unknowndiag_type;
        typedef SmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            nonconj_type;
    };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
    class SmallLowerTriMatrix : 
        public BaseMatrix_Tri_Mutable<SmallLowerTriMatrix<T,N,D,S,I> >
    {
    public:

        typedef SmallLowerTriMatrix<T,N,D,S,I> type;
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

        inline SmallLowerTriMatrix(size_t n=N)
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
        }

        explicit inline SmallLowerTriMatrix(size_t n, T x)
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
            this->setAllTo(x);
        }

        explicit inline SmallLowerTriMatrix(size_t n, const T* vv) 
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            SmallVectorView<T,N*N,1> lv(ptr());
            ConstSmallVectorView<T,N*N,1>(vv).newAssignTo(lv);
        }

        explicit inline SmallLowerTriMatrix(size_t n, const std::vector<T>& vv) 
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(n==N);
            TMVAssert(vv.size() == N*N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            SmallVectorView<T,N*N,1> lv(ptr());
            ConstSmallVectorView<T,N*N,1>(&vv[0]).newAssignTo(lv);
        }

        inline SmallLowerTriMatrix(const type& m2) 
        { 
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            m2.newAssignTo(*this);
        }

        template <class M2> 
        inline SmallLowerTriMatrix(const BaseMatrix<M2>& m2) 
        { 
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            const bool assignable = 
                ShapeTraits2<M2::_shape,_shape>::assignable;
            TMVStaticAssert((
                (M2::_calc && ShapeTraits<M2::_shape>::lower) || assignable));
            TMVAssert(m2.colsize() == N);
            TMVAssert(m2.rowsize() == N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            TriCopy<assignable>::copy(m2,*this);
        }

        template <class M2>
        inline SmallLowerTriMatrix(const BaseMatrix_Tri<M2>& m2) 
        {
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(M2::_lower);
            TMVAssert(m2.size() == N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            Maybe<_unit && !M2::_unit>::unitview(m2).newAssignTo(*this);
        }

        template <class M2> 
        inline SmallLowerTriMatrix(const BaseMatrix_Diag<M2>& m2) 
        { 
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor); 
            TMVStaticAssert(D != UnknownDiag);
            TMVAssert(m2.size() == N);
#ifdef XTEST_DEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
#endif
            typename type::diag_type d = this->diag();
            this->setZero();
            m2.calc().diag().newAssignTo(d);
        }

        inline ~SmallLowerTriMatrix()
        {
#ifdef TMVDEBUG
            Maybe<_unit>::offdiag(*this).setAllTo(T(888));
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

        inline size_t size() const { return N; }
        inline int stepi() const { return _stepi; }
        inline int stepj() const { return _stepj; }
        inline DiagType dt() const { return D; }
        inline bool isunit() const { return D == UnitDiag; }
        inline bool isconj() const { return false; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline StorageType stor() const { return S; }

    protected :

        StackArray<T,N*N> itsm;

    }; // SmallLowerTriMatrix

    template <class T, int N, DiagType D, StorageType S>
    class SmallLowerTriMatrixF : 
        public SmallLowerTriMatrix<T,N,D,S,FortranStyle>
    {
    public:

        typedef SmallLowerTriMatrixF<T,N,D,S> type;
        typedef SmallLowerTriMatrix<T,N,D,S,FortranStyle> mtype;

        inline SmallLowerTriMatrixF() : mtype() {}
        explicit inline SmallLowerTriMatrixF(T x) : mtype(x) {}
        explicit inline SmallLowerTriMatrixF(const T* vv) : mtype(vv) {}
        explicit inline SmallLowerTriMatrixF(const std::vector<T>& vv) :
            mtype(vv) {}
        template <class M2> 
        inline SmallLowerTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        inline SmallLowerTriMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        inline SmallLowerTriMatrixF(const BaseMatrix_Tri<M2>& m2) :
            mtype(m2) {}
        template <class M2>
        inline SmallLowerTriMatrixF(const BaseMatrix_Rec<M2>& m2) : 
            mtype(m2) {}
        inline ~SmallLowerTriMatrixF() {}

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

    }; // SmallLowerTriMatrixF

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
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

        // In case N == UNKNOWN
        typedef typename MCopyHelper<T,LowerTri,N,N,_rowmajor,_fort>::type 
            copy_type;

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

        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,I> 
            const_view_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,
                FortranStyle> const_fview_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,1,_stepj,C,I> 
            const_cmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,1,C,I> 
            const_rmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_offdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnitDiag,_stepi,_stepj,C,I> 
            const_unitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,C,I> 
            const_nonunitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknowndiag_type;
        typedef ConstSmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,I> nonconst_type;
    };

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class ConstSmallLowerTriMatrixView :
        public BaseMatrix_Tri<ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
    public:

        typedef ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;

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
        inline ConstSmallLowerTriMatrixView(
            const T* m, size_t s, bool u, int si, int sj
        ) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline ConstSmallLowerTriMatrixView(
            const T* m, size_t s, bool u, int si
        ) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstSmallLowerTriMatrixView(const T* m, size_t s) :
            itsm(m), itss(s), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline ConstSmallLowerTriMatrixView(const T* m) :
            itsm(m), itss(N), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline ConstSmallLowerTriMatrixView(const type& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixView(
            const ConstLowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixView(
            const LowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixView(
            const ConstSmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixView(
            const SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~ConstSmallLowerTriMatrixView() {
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
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        inline bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

#ifdef TMV_DEBUG
        const T* itsm;
#else
        const T*const itsm;
#endif
        const CheckedInt<N> itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstSmallLowerTriMatrixView

    template <class T, int N, DiagType D, int Si, int Sj, bool C>
    class ConstSmallLowerTriMatrixViewF :
        public ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef ConstSmallLowerTriMatrixViewF<T,N,D,Si,Sj,C> type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

        inline ConstSmallLowerTriMatrixViewF(
            const T* m, size_t s, bool u, int si, int sj
        ) :
            mtype(m,s,si,sj) {}
        inline ConstSmallLowerTriMatrixViewF(const T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline ConstSmallLowerTriMatrixViewF(const T* m, size_t s) :
            mtype(m,s) {}
        inline ConstSmallLowerTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixViewF(
            const ConstLowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixViewF(
            const LowerTriMatrixView<T,D2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixViewF(
            const ConstSmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallLowerTriMatrixViewF(
            const SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        inline ~ConstSmallLowerTriMatrixViewF() {}

    private :
        inline void operator=(const type& m2);

    }; // ConstSmallLowerTriMatrixViewF

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
        typedef const ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
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

        typedef typename MCopyHelper<T,LowerTri,N,N,_rowmajor,_fort>::type 
            copy_type;

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

        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,I> 
            const_view_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,CStyle> 
            const_cview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,
                FortranStyle> const_fview_type;
        typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> 
            const_xview_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            const_xdview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,1,_stepj,C,I> 
            const_cmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,1,C,I> 
            const_rmview_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepj,_stepi,C,I> 
            const_transpose_type;
        typedef ConstSmallUpperTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            const_adjoint_type;

        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,1,C,I> 
            const_offdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnitDiag,_stepi,1,C,I> 
            const_unitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,NonUnitDiag,_stepi,1,C,I> 
            const_nonunitdiag_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,1,C,I> 
            const_unknowndiag_type;
        typedef ConstSmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            const_nonconj_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,I> nonconst_type;

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

        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,I> view_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,CStyle> 
            cview_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,C,FortranStyle> 
            fview_type;
        typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
        typedef LowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> 
            xdview_type;
        typedef SmallLowerTriMatrixView<T,N,D,1,_stepj,C,I> cmview_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,1,C,I> rmview_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,notC,I> 
            conjugate_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepj,_stepi,C,I> 
            transpose_type;
        typedef SmallUpperTriMatrixView<T,N,D,_stepj,_stepi,notC,I> 
            adjoint_type;

        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            offdiag_type;
        typedef SmallLowerTriMatrixView<T,N,UnitDiag,_stepi,_stepj,C,I> 
            unitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,C,I> 
            nonunitdiag_type;
        typedef SmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,C,I> 
            unknowndiag_type;
        typedef SmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallLowerTriMatrixView<T,N,D,_stepi,_stepj,false,I> 
            nonconj_type;
    };

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    class SmallLowerTriMatrixView :
        public BaseMatrix_Tri_Mutable<SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
    {
    public:

        typedef SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
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

        inline SmallLowerTriMatrixView(T* m, size_t s, bool u, int si, int sj) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(sj) {}

        inline SmallLowerTriMatrixView(T* m, size_t s, bool u, int si) :
            itsm(m), itss(s), itsu(u), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline SmallLowerTriMatrixView(T* m, size_t s) :
            itsm(m), itss(s), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline SmallLowerTriMatrixView(T* m) :
            itsm(m), itss(N), itsu(D==UnitDiag), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(D != UnknownDiag);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline SmallLowerTriMatrixView(const type& m2) :
            itsm(m2.itsm), itss(m2.size()), 
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallLowerTriMatrixView(
            LowerTriMatrixView<T,D2,Si2,Sj2,C,I2> m2
        ) :
            itsm(m2.ptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallLowerTriMatrixView(
            SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itss(m2.size()),
            itsu(m2.isunit()), itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~SmallLowerTriMatrixView() {
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
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        inline bool iscm() const
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

#ifdef TMV_DEBUG
        T* itsm;
#else
        T*const itsm;
#endif
        const CheckedInt<N> itss;
        const CheckedInt<D==UnknownDiag ? UNKNOWN : D==UnitDiag> itsu;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // SmallLowerTriMatrixView

    template <class T, int N, DiagType D, int Si, int Sj, bool C>
    class SmallLowerTriMatrixViewF :
        public SmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
    {
    public:

        typedef SmallLowerTriMatrixViewF<T,N,D,Si,Sj,C> type;
        typedef SmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

        inline SmallLowerTriMatrixViewF(T* m, size_t s, bool u, int si, int sj) :
            mtype(m,s,si,sj) {}
        inline SmallLowerTriMatrixViewF(T* m, size_t s, bool u, int si) :
            mtype(m,s,si) {}
        inline SmallLowerTriMatrixViewF(T* m, size_t s) : mtype(m,s) {}
        inline SmallLowerTriMatrixViewF(T* m) : mtype(m) {}
        inline SmallLowerTriMatrixViewF(const type& m2) : mtype(m2) {}
        template <DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallLowerTriMatrixViewF(
            LowerTriMatrixView<T,D2,Si2,Sj2,C,I2> m2
        ) :
            mtype(m2) {}
        template <int N2, DiagType D2, int Si2, int Sj2, IndexStyle I2>
        inline SmallLowerTriMatrixViewF(
            SmallLowerTriMatrixView<T,N2,D2,Si2,Sj2,C,I2> m2) : mtype(m2) {}
        inline ~SmallLowerTriMatrixViewF() {}

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

    }; // SmallLowerTriMatrixViewF


    //
    // Swap
    //

    template <class T, int N, DiagType D, int Si, int Sj, bool C, 
              IndexStyle I, class M>
    inline void Swap(
        BaseMatrix_Tri_Mutable<M>& m1,
        SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> m2)
    { DoSwap(m1,m2); }
    template <class T, int N, DiagType D, int Si, int Sj, bool C, 
              IndexStyle I, class M>
    inline void Swap(
        SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> m1,
        BaseMatrix_Tri_Mutable<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallUpperTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
        SmallUpperTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }
    template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallUpperTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
        UpperTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }
    template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        UpperTriMatrixView<T,D,Si1,Sj1,C1,I1> m1,
        SmallUpperTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }

    template <class T, int N, DiagType D, int Si, int Sj, bool C, 
              IndexStyle I, class M>
    inline void Swap(
        BaseMatrix_Tri_Mutable<M>& m1,
        SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, DiagType D, int Si, int Sj, bool C, 
              IndexStyle I, class M>
    inline void Swap(
        SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> m1,
        BaseMatrix_Tri_Mutable<M>& m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallLowerTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
        SmallLowerTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallLowerTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
        LowerTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
    { Swap(m1.transpose(),m2.transpose()); }
    template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        LowerTriMatrixView<T,D,Si1,Sj1,C1,I1> m1,
        SmallLowerTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
    { Swap(m1.transpose(),m2.transpose()); }


    //
    // TMV_Text 
    //

    template <class T, int N, DiagType D, StorageType S, IndexStyle I>
    inline std::string TMV_Text(const SmallUpperTriMatrix<T,N,D,S,I>& )
    {
        std::ostringstream s;
        s << "SmallUpperTriMatrix<"<<TMV_Text(T())<<","<<N<<","<<TMV_Text(D);
        s << ","<<TMV_Text(S)<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
        const ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstSmallUpperTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.stepi()<<")";
        s <<","<<TMV_Text(D);
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
        const SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "SmallUpperTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.stepi()<<")";
        s <<","<<TMV_Text(D);
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, DiagType D, StorageType S, IndexStyle I>
    inline std::string TMV_Text(const SmallLowerTriMatrix<T,N,D,S,I>& )
    {
        std::ostringstream s;
        s << "SmallLowerTriMatrix<"<<TMV_Text(T())<<","<<N<<","<<TMV_Text(D);
        s << ","<<TMV_Text(S)<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
        const ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstSmallLowerTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.stepi()<<")";
        s <<","<<TMV_Text(D);
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
        const SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "SmallLowerTriMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.stepi()<<")";
        s <<","<<TMV_Text(D);
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }


} // namespace tmv

#endif
