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
//    SmallMatrix<T,M,N,stor,I>(const vector<vector<T> >& m)
//        Makes a SmallMatrix with a_ij = m[i][j]
//
//    SmallMatrix<T,M,N,stor,I>(const T* m)
//    SmallMatrix<T,M,N,stor,I>(const vector<T>& m)
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

#include "TMV_SmallVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"

namespace tmv {

    //
    // SmallMatrix
    //

    template <class T, int M, int N, StorageType S, IndexStyle I>
    struct Traits<SmallMatrix<T,M,N,S,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallMatrix<T,M,N,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _shape = Rec };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==RowMajor ? N : 1) };
        enum { _stepj = (S==RowMajor ? 1 : M) };
        enum { _diagstep = _stepi + _stepj };
        enum { _conj = false };
        enum { _canlin = true };

        enum { _hasdivider = false };
        typedef LUD<copy_type> lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;

        enum { minMN = M > N ? N : M };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { notC = iscomplex };

        typedef ConstSmallVectorView<T,M,_stepi,false,I> const_col_type;
        typedef ConstVectorView<T,_stepi,false,I> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,false,I> const_row_type;
        typedef ConstVectorView<T,_stepj,false,I> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,false,I>
            const_diag_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_sub_type;

        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I>
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,UNKNOWN,false,I>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,_stepj,false,I>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,false,I>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,false,I>
            const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,false,I>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,false,CStyle>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,false,FortranStyle>
            const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,false,I> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,false,I> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,notC,I>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,false,I>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,notC,I>
            const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,
                false,I> const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,_stepi,_stepj,
                false,I> const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,
                false,I> const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,NonUnitDiag,_stepi,_stepj,
                false,I> const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,UnitDiag,_stepi,_stepj,
                false,I> const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,
                false,I> const_unknown_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,false,I>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,M*N,1,false,I> const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,false,I>
            const_nonconj_type;

        typedef SmallMatrixView<T,M,N,_stepi,_stepj,false,I> nonconst_type;

        typedef T& reference;

        typedef SmallVectorView<T,M,_stepi,false,I> col_type;
        typedef VectorView<T,_stepi,false,I> col_sub_type;
        typedef SmallVectorView<T,N,_stepj,false,I> row_type;
        typedef VectorView<T,_stepj,false,I> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,false,I> diag_type;
        typedef VectorView<T,_diagstep,false,I> diag_sub_type;

        typedef MatrixView<T,_stepi,_stepj,false,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_type;
        typedef SmallMatrixView<T,M,2,_stepi,UNKNOWN,false,I> colpair_type;
        typedef SmallMatrixView<T,2,N,UNKNOWN,_stepj,false,I> rowpair_type;
        typedef SmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,false,I> colrange_type;
        typedef SmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,false,I> rowrange_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,false,I> view_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,false,CStyle> cview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,false,FortranStyle>
            fview_type;
        typedef MatrixView<T> xview_type;
        typedef SmallMatrixView<T,M,N,1,_stepj,false,I> cmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,1,false,I> rmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,notC,I> conjugate_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,false,I> transpose_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,notC,I> adjoint_type;
        typedef SmallUpperTriMatrixView<T,N,NonUnitDiag,_stepi,_stepj,false,I>
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,UnitDiag,_stepi,_stepj,false,I>
            unit_uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,false,I>
            unknown_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,NonUnitDiag,_stepi,_stepj,false,I>
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,UnitDiag,_stepi,_stepj,false,I>
            unit_lowertri_type;
        typedef SmallLowerTriMatrixView<T,N,UnknownDiag,_stepi,_stepj,false,I>
            unknown_lowertri_type;
        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,false,I>
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,M*N,1,false,I> linearview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,false,I> nonconj_type;
    };

    template <class T, int M, int N, StorageType S, IndexStyle I>
    class SmallMatrix : 
        public BaseMatrix_Rec_Mutable<SmallMatrix<T,M,N,S,I> >,
        public MatrixDivHelper<SmallMatrix<T,M,N,S,I>,false>
    {
    public:

        typedef SmallMatrix<T,M,N,S,I> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename Traits<type>::lud_type lud_type;
        typedef typename Traits<type>::qrd_type qrd_type;
        typedef typename Traits<type>::qrpd_type qrpd_type;
        typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _canlin = Traits<type>::_canlin };

        //
        // Constructors
        //

        SmallMatrix() 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        explicit SmallMatrix(T x) 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            this->setAllTo(x);
        }

        explicit SmallMatrix(const T* vv) 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            typename type::linearview_type lv = this->linearView();
            ConstSmallVectorView<T,M*N,1>(vv).newAssignTo(lv);
        }

        explicit SmallMatrix(const std::vector<T>& vv) 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == M*N);
            typename type::linearview_type lv = this->linearView();
            ConstSmallVectorView<T,M*N,1>(&vv[0]).newAssignTo(lv);
        }

        explicit SmallMatrix(const std::vector<std::vector<T> >& vv) 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == M);
            for(int i=0;i<M;++i) {
                TMVAssert(vv[i].size() == N);
                typename type::row_type mi = this->row(i);
                ConstSmallVectorView<T,N,1>(&vv[i][0]).newAssignTo(mi);
            }
        }

        SmallMatrix(const type& m2) 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            m2.newAssignTo(*this);
        }

        template <class M2>
        SmallMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            m2.newAssignTo(*this);
        }

        ~SmallMatrix()
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
#endif
        }


        //
        // Op=
        //

        type& operator=(const type& m2)
        {
            if (&m2 != this) base_mut::operator=(m2);
            return *this;
        }

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this;
        }

        type& operator=(T x) 
        {
            base_mut::operator=(x); 
            return *this;
        }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        size_t ls() const { return M*N; }
        size_t colsize() const { return M; }
        size_t rowsize() const { return N; }
        int nElements() const { return M*N; }
        int stepi() const { return _stepi; }
        int stepj() const { return _stepj; }
        bool isconj() const { return false; }
        bool isrm() const { return S == RowMajor; }
        bool iscm() const { return S == ColMajor; }
        StorageType stor() const { return S; }

    private:

        StackArray<T,M*N> itsm;

    }; // SmallMatrix

    template <class T, int M, int N, StorageType S>
    class SmallMatrixF : 
        public SmallMatrix<T,M,N,S,FortranStyle>
    {
    public:
        typedef SmallMatrixF<T,M,N,S> type;
        typedef SmallMatrix<T,M,N,S,FortranStyle> mtype;

        SmallMatrixF() : mtype() {}
        explicit SmallMatrixF(T x) : mtype(x) {}
        explicit SmallMatrixF(const T* vv) : mtype(vv) {}
        explicit SmallMatrixF(const std::vector<T>& vv) : mtype(vv) {}
        explicit SmallMatrixF(
            const std::vector<std::vector<T> >& vv) : mtype(vv) {}
        SmallMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        SmallMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        ~SmallMatrixF() {}

        type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(T x) 
        { mtype::operator=(x); return *this; }

    }; // SmallMatrixF

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstSmallMatrixView<T,M,N,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _shape = Rec };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _canlin = (Si == 1 && Sj == M) || (Sj == 1 && Si == N) };

        // In case M,N == UNKNOWN
        typedef typename MCopyHelper<T,Rec,M,N,_rowmajor,_fort>::type 
            copy_type;

        enum { _hasdivider = false };
        typedef LUD<copy_type> lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;

        enum { minMN = M > N ? N : M };
        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };
        enum { prodMN = IntTraits2<M,N>::prod };

        typedef ConstSmallVectorView<T,M,Si,C,I> const_col_type;
        typedef ConstVectorView<T,Si,C,I> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,Sj,C,I> const_row_type;
        typedef ConstVectorView<T,Sj,C,I> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstMatrixView<T,Si,Sj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I>
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,Si,UNKNOWN,C,I> const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,Sj,C,I> const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,Si,Sj,C,I> const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,Si,Sj,C,I> const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,I> const_view_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle> const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>
            const_fview_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,Sj,C,I> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,Si,1,C,I> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,notC,I> const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,Sj,Si,C,I> const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,Sj,Si,notC,I> const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,Si,Sj,C,I>
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,Si,Sj,C,I>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnknownDiag,Si,Sj,C,I>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,NonUnitDiag,Si,Sj,C,I>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,UnitDiag,Si,Sj,C,I>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnknownDiag,Si,Sj,C,I>
            const_unknown_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,false,I>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,prodMN,1,C,I> const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,false,I>
            const_nonconj_type;

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> nonconst_type;
    };

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    class ConstSmallMatrixView :
        public BaseMatrix_Rec<ConstSmallMatrixView<T,M,N,Si,Sj,C,I> >,
        public MatrixDivHelper<ConstSmallMatrixView<T,M,N,Si,Sj,C,I>,false>
    {
    public:

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef typename Traits<type>::lud_type lud_type;
        typedef typename Traits<type>::qrd_type qrd_type;
        typedef typename Traits<type>::qrpd_type qrpd_type;
        typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _canlin = Traits<type>::_canlin };

        //
        // Constructors
        //

        ConstSmallMatrixView(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

        ConstSmallMatrixView(const T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        ConstSmallMatrixView(const T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        ConstSmallMatrixView(const T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        ConstSmallMatrixView(const T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(M != UNKNOWN); TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        ConstSmallMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixView(
            const ConstMatrixView<T,Si2,Sj2,C,I2>& m2
        ) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixView(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixView(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixView(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        ~ConstSmallMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

    private:
        void operator=(const type& v2);
    public:


        // 
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }

        T cref(int i, int j) const
        { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

        size_t ls() const { return int(itscs)*int(itsrs); }
        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return int(itscs)*int(itsrs); }
        int stepi() const { return itssi; }
        int stepj() const { return itssj; }
        bool isconj() const { return C; }
        bool isrm() const 
        { return _rowmajor || (!_colmajor &&  stepj() == 1); }
        bool iscm() const 
        { return _colmajor || (!_rowmajor &&  stepi() == 1); }

    private :

        const T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstSmallMatrixView

    template <class T, int M, int N, int Si, int Sj, bool C>
    class ConstSmallMatrixViewF :
        public ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>
    {
    public:

        typedef ConstSmallMatrixViewF<T,M,N,Si,Sj,C> type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> mtype;

        ConstSmallMatrixViewF(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            mtype(m,cs,rs,si,sj) {}
        ConstSmallMatrixViewF(
            const T* m, size_t cs, size_t rs, int si) : mtype(m,cs,rs,si) {}
        ConstSmallMatrixViewF(const T* m, size_t cs, size_t rs) :
            mtype(m,cs,rs) {}
        ConstSmallMatrixViewF(const T* m, size_t cs) : mtype(m,cs) {}
        ConstSmallMatrixViewF(const T* m) : mtype(m) {}
        ConstSmallMatrixViewF(const type& m2) : mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixViewF(
            const ConstMatrixView<T,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixViewF(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixViewF(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstSmallMatrixViewF(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        ~ConstSmallMatrixViewF() {}

    private:
        void operator=(const type& v2);
    }; // ConstSmallMatrixViewF

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<SmallMatrixView<T,M,N,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef const ConstSmallMatrixView<T,M,N,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _shape = Rec };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _canlin = (Si == 1 && Sj == M) || (Sj == 1 && Si == N) };

        // In case M,N == UNKNOWN
        typedef typename MCopyHelper<T,Rec,M,N,_rowmajor,_fort>::type 
            copy_type;

        enum { minMN = M > N ? N : M };
        enum { _hasdivider = false };
        typedef LUD<copy_type> lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;

        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };
        enum { prodMN = IntTraits2<M,N>::prod };

        typedef ConstSmallVectorView<T,M,Si,C,I> const_col_type;
        typedef ConstVectorView<T,Si,C,I> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,Sj,C,I> const_row_type;
        typedef ConstVectorView<T,Sj,C,I> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstMatrixView<T,Si,Sj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I>
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,Si,UNKNOWN,C,I> const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,Sj,C,I> const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,Si,Sj,C,I> const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,Si,Sj,C,I> const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,I> const_view_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle> const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>
            const_fview_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,Sj,C,I> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,Si,1,C,I> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,notC,I> const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,Sj,Si,C,I> const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,Sj,Si,notC,I> const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,Si,Sj,C,I>
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,Si,Sj,C,I>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnknownDiag,Si,Sj,C,I>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,NonUnitDiag,Si,Sj,C,I>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,UnitDiag,Si,Sj,C,I>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,UnknownDiag,Si,Sj,C,I>
            const_unknown_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,false,I>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,prodMN,1,C,I> const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,false,I>
            const_nonconj_type;

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> nonconst_type;

        typedef typename AuxRef<T,C>::reference reference;

        typedef SmallVectorView<T,M,Si,C,I> col_type;
        typedef VectorView<T,Si,C,I> col_sub_type;
        typedef SmallVectorView<T,N,Sj,C,I> row_type;
        typedef VectorView<T,Sj,C,I> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,C,I> diag_type;
        typedef VectorView<T,_diagstep,C,I> diag_sub_type;

        typedef MatrixView<T,Si,Sj,C,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,C,I> subvector_type;
        typedef SmallMatrixView<T,M,2,Si,UNKNOWN,C,I> colpair_type;
        typedef SmallMatrixView<T,2,N,UNKNOWN,Sj,C,I> rowpair_type;
        typedef SmallMatrixView<T,M,UNKNOWN,Si,Sj,C,I> colrange_type;
        typedef SmallMatrixView<T,UNKNOWN,N,Si,Sj,C,I> rowrange_type;
        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> view_type;
        typedef SmallMatrixView<T,M,N,Si,Sj,C,CStyle> cview_type;
        typedef SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> fview_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,C> xview_type;
        typedef SmallMatrixView<T,M,N,1,Sj,C,I> cmview_type;
        typedef SmallMatrixView<T,M,N,Si,1,C,I> rmview_type;
        typedef SmallMatrixView<T,M,N,Si,Sj,notC,I> conjugate_type;
        typedef SmallMatrixView<T,N,M,Sj,Si,C,I> transpose_type;
        typedef SmallMatrixView<T,N,M,Sj,Si,notC,I> adjoint_type;
        typedef SmallUpperTriMatrixView<T,N,NonUnitDiag,Si,Sj,C,I>
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,UnitDiag,Si,Sj,C,I>
            unit_uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,UnknownDiag,Si,Sj,C,I>
            unknown_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,NonUnitDiag,Si,Sj,C,I>
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,UnitDiag,Si,Sj,C,I>
            unit_lowertri_type;
        typedef SmallLowerTriMatrixView<T,N,UnknownDiag,Si,Sj,C,I>
            unknown_lowertri_type;
        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,false,I>
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,prodMN,1,C,I> linearview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,false,I> nonconj_type;
    };

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    class SmallMatrixView :
        public BaseMatrix_Rec_Mutable<SmallMatrixView<T,M,N,Si,Sj,C,I> >,
        public MatrixDivHelper<SmallMatrixView<T,M,N,Si,Sj,C,I>,false>
    {
    public:

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;
        typedef typename Traits<type>::lud_type lud_type;
        typedef typename Traits<type>::qrd_type qrd_type;
        typedef typename Traits<type>::qrpd_type qrpd_type;
        typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _canlin = Traits<type>::_canlin };

        //
        // Constructors
        //

        SmallMatrixView(T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

        SmallMatrixView(T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        SmallMatrixView(T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        SmallMatrixView(T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        SmallMatrixView(T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(M != UNKNOWN); TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        SmallMatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        SmallMatrixView(MatrixView<T,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        SmallMatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        ~SmallMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

        //
        // Op=
        //

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        type& operator=(const type& m2)
        { base_mut::operator=(m2); return *this; }

        type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

        reference ref(int i, int j)
        { return reference(itsm[i*stepi()+j*stepj()]); }

        size_t ls() const { return int(itscs)*int(itsrs); }
        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return int(itscs)*int(itsrs); }
        int stepi() const { return itssi; }
        int stepj() const { return itssj; }
        bool isconj() const { return C; }
        bool isrm() const 
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        bool iscm() const 
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // SmallMatrixView

    template <class T, int M, int N, int Si, int Sj, bool C>
    class SmallMatrixViewF :
        public SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>
    {
    public:
        typedef SmallMatrixViewF<T,M,N,Si,Sj,C> type;
        typedef SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> mtype;

        SmallMatrixViewF(T* m, size_t cs, size_t rs, int si, int sj) :
            mtype(m,cs,rs,si,sj) {}
        SmallMatrixViewF(T* m, size_t cs, size_t rs, int si) :
            mtype(m,cs,rs,si) {}
        SmallMatrixViewF(T* m, size_t cs, size_t rs) : mtype(m,cs,rs) {}
        SmallMatrixViewF(T* m, size_t cs) : mtype(m,cs) {}
        SmallMatrixViewF(T* m) : mtype(m) {}
        SmallMatrixViewF(const type& m2) : mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        SmallMatrixViewF(MatrixView<T,Si2,Sj2,C,I2> m2) : mtype(m2) {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        SmallMatrixViewF(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
            mtype(m2) {}
        ~SmallMatrixViewF() {}

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2);        return *this; }
        type& operator=(const type& m2)
        { mtype::operator=(m2);        return *this; }
        type& operator=(const T x)
        { mtype::operator=(x);        return *this; }

    }; // SmallMatrixViewF


    // Special Creators: 
    //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
    //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage

#define VT typename V::value_type
#define VN V::_size
#define VS V::_step
#define VC V::_conj
#define VI (V::_fort ? FortranStyle : CStyle)

    template <class V>
    static inline ConstSmallMatrixView<VT,1,VN,VN,VS,VC,VI> RowVectorViewOf(
        const BaseVector_Calc<V>& v)
    {
        return ConstSmallMatrixView<VT,1,VN,VN,VS,VC,VI>(
            v.cptr(),1,v.size(),v.size(),v.step());
    }

    template <class V>
    static inline SmallMatrixView<VT,1,VN,VN,VS,VC,VI> RowVectorViewOf(
        BaseVector_Mutable<V>& v)
    {
        return SmallMatrixView<VT,1,VN,VN,VS,VC,VI>(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class V>
    static inline ConstSmallMatrixView<VT,VN,1,VS,VN,VC,VI> ColVectorViewOf(
        const BaseVector_Calc<V>& v)
    {
        return ConstSmallMatrixView<VT,VN,1,VS,VN,VC,VI>(
            v.cptr(),v.size(),1,v.step(),v.size());
    }

    template <class V>
    static inline SmallMatrixView<VT,VN,1,VS,VN,VC,VI> ColVectorViewOf(
        BaseVector_Mutable<V>& v)
    {
        return SmallMatrixView<VT,VN,1,VS,VN,VC,VI>(
            v.ptr(),v.size(),1,v.step(),v.size());
    }

#undef VT
#undef VN
#undef VS
#undef VC
#undef VI


    // Repeat for VectorView and SmallVectorView so we can have them 
    // work without requiring the non-const reference.
    template <class T, int S, bool C, IndexStyle I>
    static inline SmallMatrixView<T,1,UNKNOWN,UNKNOWN,S,C,I> RowVectorViewOf(
        VectorView<T,S,C,I> v)
    {
        return SmallMatrixView<T,1,UNKNOWN,UNKNOWN,S,C,I>(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    static inline SmallMatrixView<T,1,N,N,S,C,I> RowVectorViewOf(
        SmallVectorView<T,N,S,C,I> v)
    {
        return SmallMatrixView<T,1,N,N,S,C,I>(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int S, bool C, IndexStyle I>
    static inline SmallMatrixView<T,UNKNOWN,1,S,UNKNOWN,C,I> ColVectorViewOf(
        VectorView<T,S,C,I> v)
    {
        return SmallMatrixView<T,UNKNOWN,1,S,UNKNOWN,C,I>(
            v.ptr(),v.size(),1,v.step(),v.size());
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    static inline SmallMatrixView<T,N,1,S,N,C,I> ColVectorViewOf(
        SmallVectorView<T,N,S,C,I> v)
    {
        return SmallMatrixView<T,N,1,S,N,C,I>(
            v.ptr(),v.size(),1,v.step(),v.size());
    }


    //
    // Swap
    //

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I, class MM>
    static inline void Swap(
        BaseMatrix_Rec_Mutable<MM>& m1, SmallMatrixView<T,M,N,Si,Sj,C,I> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I, class MM>
    static inline void Swap(
        SmallMatrixView<T,M,N,Si,Sj,C,I> m1, BaseMatrix_Rec_Mutable<MM>& m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    static inline void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,C1,I1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    static inline void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,C1,I1> m1, 
        MatrixView<T,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    static inline void Swap(
        MatrixView<T,Si1,Sj1,C1,I1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //
    
    template <class T, int M, int N, StorageType S, IndexStyle I>
    static inline typename SmallMatrix<T,M,N,S,I>::conjugate_type Conjugate(
        SmallMatrix<T,M,N,S,I>& m)
    { return m.conjugate(); }
    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    static inline typename SmallMatrixView<T,M,N,Si,Sj,C,I>::conjugate_type Conjugate(
        SmallMatrixView<T,M,N,Si,Sj,C,I> m)
    { return m.conjugate(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    static inline typename SmallMatrix<T,M,N,S,I>::transpose_type Transpose(
        SmallMatrix<T,M,N,S,I>& m)
    { return m.transpose(); }
    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    static inline typename SmallMatrixView<T,M,N,Si,Sj,C,I>::transpose_type Transpose(
        SmallMatrixView<T,M,N,Si,Sj,C,I> m)
    { return m.transpose(); }

    template <class T, int M, int N, StorageType S, IndexStyle I>
    static inline typename SmallMatrix<T,M,N,S,I>::adjoint_type Adjoint(
        SmallMatrix<T,M,N,S,I>& m)
    { return m.adjoint(); }
    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    static inline typename SmallMatrixView<T,M,N,Si,Sj,C,I>::adjoint_type Adjoint(
        SmallMatrixView<T,M,N,Si,Sj,C,I> m)
    { return m.adjoint(); }


    // 
    // TMV_Text
    //

    template <class T, int M, int N, StorageType S, IndexStyle I>
    static inline std::string TMV_Text(const SmallMatrix<T,M,N,S,I>& )
    {
        std::ostringstream s;
        s << "SmallMatrix<"<<TMV_Text(T())<<','<<M<<','<<N;
        s <<','<<TMV_Text(S)<<','<<TMV_Text(I)<<">"; 
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    static inline std::string TMV_Text(
        const SmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "SmallMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<M>::text();
        if (M == UNKNOWN) s << "("<<m.colsize()<<")";
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.rowsize()<<")";
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    static inline std::string TMV_Text(
        const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstSmallMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<M>::text();
        if (M == UNKNOWN) s << "("<<m.colsize()<<")";
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.rowsize()<<")";
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

} // namespace tmv

#endif
