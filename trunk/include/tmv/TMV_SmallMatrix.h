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
// heap is pretty irrelevant.  (See TMV_StackArray.h to change this.)
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
        enum { misreal = Traits<T>::isreal };
        enum { miscomplex = Traits<T>::iscomplex };

        typedef SmallMatrix<T,M,N,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;
        typedef InvalidType inverse_type;

        enum { mcolsize = M };
        enum { mrowsize = N };
        enum { mshape = Rec };
        enum { mfort = (I == FortranStyle) };
        enum { mcalc = true };
        enum { mrowmajor = (S == RowMajor) };
        enum { mcolmajor = (S == ColMajor) };
        enum { mstor = S };
        enum { mstepi = (S==RowMajor ? N : 1) };
        enum { mstepj = (S==RowMajor ? 1 : M) };
        enum { mdiagstep = mstepi + mstepj };
        enum { mconj = false };
        enum { mcanlin = true };

        enum { minMN = M > N ? N : M };
        enum { twoSi = misreal ? int(mstepi) : int(IntTraits<mstepi>::twoS) };
        enum { twoSj = misreal ? int(mstepj) : int(IntTraits<mstepj>::twoS) };
        enum { notC = miscomplex };

        typedef ConstSmallVectorView<T,M,mstepi,false,I> const_col_type;
        typedef ConstVectorView<T,mstepi,false,I> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,mstepj,false,I> const_row_type;
        typedef ConstVectorView<T,mstepj,false,I> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,mdiagstep,false,I> 
            const_diag_type;
        typedef ConstVectorView<T,mdiagstep,false,I> const_diag_sub_type;

        typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,mstepi,UNKNOWN,false,I> 
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,mstepj,false,I> 
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,mstepi,mstepj,false,I> 
            const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,mstepi,mstepj,false,I> 
            const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,false,I> 
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,false,CStyle> 
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,false,FortranStyle> 
            const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,mstepj,false,I> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,1,false,I> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,notC,I> 
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,mstepj,mstepi,false,I> 
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,mstepj,mstepi,notC,I> 
            const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,NonUnitDiag,mstepi,mstepj,
                false,I> const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,mstepi,mstepj,
                false,I> const_unit_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,NonUnitDiag,mstepi,mstepj,
                false,I> const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,UnitDiag,mstepi,mstepj,
                false,I> const_unit_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,M*N,1,false,I> const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,false,I> 
            const_nonconj_type;

        typedef SmallMatrixView<T,M,N,mstepi,mstepj,false,I> nonconst_type;

        typedef T& reference;

        typedef SmallVectorView<T,M,mstepi,false,I> col_type;
        typedef VectorView<T,mstepi,false,I> col_sub_type;
        typedef SmallVectorView<T,N,mstepj,false,I> row_type;
        typedef VectorView<T,mstepj,false,I> row_sub_type;
        typedef SmallVectorView<T,minMN,mdiagstep,false,I> diag_type;
        typedef VectorView<T,mdiagstep,false,I> diag_sub_type;

        typedef MatrixView<T,mstepi,mstepj,false,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_type;
        typedef SmallMatrixView<T,M,2,mstepi,UNKNOWN,false,I> colpair_type;
        typedef SmallMatrixView<T,2,N,UNKNOWN,mstepj,false,I> rowpair_type;
        typedef SmallMatrixView<T,M,UNKNOWN,mstepi,mstepj,false,I> colrange_type;
        typedef SmallMatrixView<T,UNKNOWN,N,mstepi,mstepj,false,I> rowrange_type;
        typedef SmallMatrixView<T,M,N,mstepi,mstepj,false,I> view_type;
        typedef SmallMatrixView<T,M,N,mstepi,mstepj,false,CStyle> cview_type;
        typedef SmallMatrixView<T,M,N,mstepi,mstepj,false,FortranStyle> 
            fview_type;
        typedef MatrixView<T> xview_type;
        typedef SmallMatrixView<T,M,N,1,mstepj,false,I> cmview_type;
        typedef SmallMatrixView<T,M,N,mstepi,1,false,I> rmview_type;
        typedef SmallMatrixView<T,M,N,mstepi,mstepj,notC,I> conjugate_type;
        typedef SmallMatrixView<T,N,M,mstepj,mstepi,false,I> transpose_type;
        typedef SmallMatrixView<T,N,M,mstepj,mstepi,notC,I> adjoint_type;
        typedef SmallUpperTriMatrixView<T,N,NonUnitDiag,mstepi,mstepj,false,I> 
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,UnitDiag,mstepi,mstepj,false,I> 
            unit_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,NonUnitDiag,mstepi,mstepj,false,I> 
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,UnitDiag,mstepi,mstepj,false,I> 
            unit_lowertri_type;
        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,M*N,1,false,I> linearview_type;
        typedef SmallMatrixView<T,M,N,mstepi,mstepj,false,I> nonconj_type;
    };

#ifdef XTEST
#ifdef TMV_DEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    class SmallMatrix : public BaseMatrix_Rec_Mutable<SmallMatrix<T,M,N,S,I> >
    {
    public:

        typedef SmallMatrix<T,M,N,S,I> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;

        enum { mcolsize = Traits<type>::mcolsize };
        enum { mrowsize = Traits<type>::mrowsize };
        enum { mshape = Traits<type>::mshape };
        enum { mfort = Traits<type>::mfort };
        enum { mrowmajor = Traits<type>::mrowmajor };
        enum { mcolmajor = Traits<type>::mcolmajor };
        enum { mstor = Traits<type>::mstor };
        enum { mcalc = Traits<type>::mcalc };
        enum { mstepi = Traits<type>::mstepi };
        enum { mstepj = Traits<type>::mstepj };
        enum { mdiagstep = Traits<type>::mdiagstep };
        enum { mconj = Traits<type>::mconj };
        enum { mcanlin = Traits<type>::mcanlin };

        //
        // Constructors
        //

        inline SmallMatrix() 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        inline SmallMatrix(size_t cs, size_t rs) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(cs==M);
            TMVAssert(rs==N);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        explicit inline SmallMatrix(size_t cs, size_t rs, T x) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(cs==M);
            TMVAssert(rs==N);
            this->setAllTo(x);
        }

        explicit inline SmallMatrix(size_t cs, size_t rs, const T* vv) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(cs==M);
            TMVAssert(rs==N);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            typename type::linearview_type lv = this->linearView();
            ConstSmallVectorView<T,M*N,1>(vv).newAssignTo(lv);
        }

        explicit inline SmallMatrix(
            size_t cs, size_t rs, const std::vector<T>& vv) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(cs==M);
            TMVAssert(rs==N);
            TMVAssert(vv.size() == M*N);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            typename type::linearview_type lv = this->linearView();
            ConstSmallVectorView<T,M*N,1>(&vv[0]).newAssignTo(lv);
        }

        explicit inline SmallMatrix(const std::vector<std::vector<T> >& vv) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == M);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            for(int i=0;i<M;++i) {
                TMVAssert(vv[i].size() == N);
                typename type::row_type mi = this->row(i);
                ConstSmallVectorView<T,N,1>(&vv[i][0]).newAssignTo(mi);
            }
        }

        inline SmallMatrix(const type& m2) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            m2.newAssignTo(*this);
        }

        template <class M2>
        inline SmallMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(M>0);
            TMVStaticAssert(N>0);
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            m2.newAssignTo(*this);
        }

        inline ~SmallMatrix()
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
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
        { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        inline T& ref(int i, int j)
        { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        inline size_t ls() const { return M*N; }
        inline size_t colsize() const { return M; }
        inline size_t rowsize() const { return N; }
        inline int stepi() const { return mstepi; }
        inline int stepj() const { return mstepj; }
        inline bool isconj() const { return false; }
        inline bool isrm() const { return S == RowMajor; }
        inline bool iscm() const { return S == ColMajor; }
        inline StorageType stor() const { return S; }

    private:

        StackArray<T,M*N> itsm;

    }; // SmallMatrix

    template <class T, int M, int N, StorageType S>
    class SmallMatrixF : public SmallMatrix<T,M,N,S,FortranStyle>
    {
    public:
        typedef SmallMatrixF<T,M,N,S> type;
        typedef SmallMatrix<T,M,N,S,FortranStyle> mtype;

        inline SmallMatrixF() : mtype() {}
        explicit inline SmallMatrixF(T x) : mtype(x) {}
        explicit inline SmallMatrixF(const T* vv) : mtype(vv) {}
        explicit inline SmallMatrixF(const std::vector<T>& vv) : mtype(vv) {}
        explicit inline SmallMatrixF(
            const std::vector<std::vector<T> >& vv) : mtype(vv) {}
        inline SmallMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        inline SmallMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        inline ~SmallMatrixF() {}

        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(T x) 
        { mtype::operator=(x); return *this; }

    }; // SmallMatrixF

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstSmallMatrixView<T,M,N,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { misreal = Traits<T>::isreal };
        enum { miscomplex = Traits<T>::iscomplex };

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef InvalidType inverse_type;

        enum { mcolsize = M };
        enum { mrowsize = N };
        enum { mshape = Rec };
        enum { mfort = (I == FortranStyle) };
        enum { mcalc = true };
        enum { mrowmajor = (Sj == 1) };
        enum { mcolmajor = (Si == 1) };
        enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
        enum { mstepi = Si };
        enum { mstepj = Sj };
        enum { mdiagstep = IntTraits2<Si,Sj>::sum };
        enum { mconj = C };
        enum { mcanlin = (Si == 1 && Sj == M) || (Sj == 1 && Si == N) };

        // In case M,N == UNKNOWN
        typedef typename MCopyHelper<T,Rec,M,N,mrowmajor,mfort>::type 
            copy_type;

        enum { minMN = M > N ? N : M };
        enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && miscomplex };
        enum { prodMN = IntTraits2<M,N>::prod };

        typedef ConstSmallVectorView<T,M,Si,C,I> const_col_type;
        typedef ConstVectorView<T,Si,C,I> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,Sj,C,I> const_row_type;
        typedef ConstVectorView<T,Sj,C,I> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,mdiagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,mdiagstep,C,I> const_diag_sub_type;

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
        typedef ConstSmallLowerTriMatrixView<T,M,NonUnitDiag,Si,Sj,C,I> 
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,UnitDiag,Si,Sj,C,I> 
            const_unit_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,prodMN,1,C,I> const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,false,I> 
            const_nonconj_type;

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> nonconst_type;
    };

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    class ConstSmallMatrixView :
        public BaseMatrix_Rec<ConstSmallMatrixView<T,M,N,Si,Sj,C,I> >
    {
    public:

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,C,I> type;

        enum { mcolsize = Traits<type>::mcolsize };
        enum { mrowsize = Traits<type>::mrowsize };
        enum { mshape = Traits<type>::mshape };
        enum { mfort = Traits<type>::mfort };
        enum { mrowmajor = Traits<type>::mrowmajor };
        enum { mcolmajor = Traits<type>::mcolmajor };
        enum { mstor = Traits<type>::mstor };
        enum { mcalc = Traits<type>::mcalc };
        enum { mstepi = Traits<type>::mstepi };
        enum { mstepj = Traits<type>::mstepj };
        enum { mdiagstep = Traits<type>::mdiagstep };
        enum { mconj = Traits<type>::mconj };
        enum { mcanlin = Traits<type>::mcanlin };

        //
        // Constructors
        //

        inline ConstSmallMatrixView(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

        inline ConstSmallMatrixView(const T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstSmallMatrixView(const T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        inline ConstSmallMatrixView(const T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline ConstSmallMatrixView(const T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(M != UNKNOWN); TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline ConstSmallMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixView(
            const ConstMatrixView<T,Si2,Sj2,C,I2>& m2
        ) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixView(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixView(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixView(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~ConstSmallMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

    private:
        inline void operator=(const type& v2);
    public:


        // 
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }

        inline T cref(int i, int j) const
        { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

        inline size_t ls() const { return int(itscs)*int(itsrs); }
        inline size_t colsize() const { return itscs; }
        inline size_t rowsize() const { return itsrs; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline bool isconj() const { return C; }
        inline bool isrm() const 
        { return mrowmajor || (!mcolmajor &&  stepj() == 1); }
        inline bool iscm() const 
        { return mcolmajor || (!mrowmajor &&  stepi() == 1); }

    private :

#ifdef TMV_DEBUG
        const T* itsm;
#else
        const T*const itsm;
#endif
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

        inline ConstSmallMatrixViewF(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            mtype(m,cs,rs,si,sj) {}
        inline ConstSmallMatrixViewF(
            const T* m, size_t cs, size_t rs, int si
        ) :
            mtype(m,cs,rs,si) {}
        inline ConstSmallMatrixViewF(const T* m, size_t cs, size_t rs) :
            mtype(m,cs,rs) {}
        inline ConstSmallMatrixViewF(const T* m, size_t cs) : mtype(m,cs) {}
        inline ConstSmallMatrixViewF(const T* m) : mtype(m) {}
        inline ConstSmallMatrixViewF(const type& m2) : mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixViewF(
            const ConstMatrixView<T,Si2,Sj2,C,I2>& m2
        ) :
            mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixViewF(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
            mtype(m2) {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixViewF(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2
        ) : 
            mtype(m2) {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        inline ConstSmallMatrixViewF(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
        inline ~ConstSmallMatrixViewF() {}

    private:
        inline void operator=(const type& v2);
    }; // ConstSmallMatrixViewF

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<SmallMatrixView<T,M,N,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { misreal = Traits<T>::isreal };
        enum { miscomplex = Traits<T>::iscomplex };

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef const ConstSmallMatrixView<T,M,N,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;
        typedef InvalidType inverse_type;

        enum { mcolsize = M };
        enum { mrowsize = N };
        enum { mshape = Rec };
        enum { mfort = (I == FortranStyle) };
        enum { mcalc = true };
        enum { mrowmajor = (Sj == 1) };
        enum { mcolmajor = (Si == 1) };
        enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
        enum { mstepi = Si };
        enum { mstepj = Sj };
        enum { mdiagstep = IntTraits2<Si,Sj>::sum };
        enum { mconj = C };
        enum { mcanlin = (Si == 1 && Sj == M) || (Sj == 1 && Si == N) };

        // In case M,N == UNKNOWN
        typedef typename MCopyHelper<T,Rec,M,N,mrowmajor,mfort>::type 
            copy_type;

        enum { minMN = M > N ? N : M };
        enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && miscomplex };
        enum { prodMN = IntTraits2<M,N>::prod };

        typedef ConstSmallVectorView<T,M,Si,C,I> const_col_type;
        typedef ConstVectorView<T,Si,C,I> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,Sj,C,I> const_row_type;
        typedef ConstVectorView<T,Sj,C,I> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,mdiagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,mdiagstep,C,I> const_diag_sub_type;

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
        typedef ConstSmallLowerTriMatrixView<T,M,NonUnitDiag,Si,Sj,C,I> 
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,UnitDiag,Si,Sj,C,I> 
            const_unit_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,prodMN,1,C,I> const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,mstepi,mstepj,false,I> 
            const_nonconj_type;

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> nonconst_type;

        typedef typename AuxRef<T,C>::reference reference;

        typedef SmallVectorView<T,M,Si,C,I> col_type;
        typedef VectorView<T,Si,C,I> col_sub_type;
        typedef SmallVectorView<T,N,Sj,C,I> row_type;
        typedef VectorView<T,Sj,C,I> row_sub_type;
        typedef SmallVectorView<T,minMN,mdiagstep,C,I> diag_type;
        typedef VectorView<T,mdiagstep,C,I> diag_sub_type;

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
        typedef SmallLowerTriMatrixView<T,M,NonUnitDiag,Si,Sj,C,I> 
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,UnitDiag,Si,Sj,C,I> 
            unit_lowertri_type;
        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,false,I> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,prodMN,1,C,I> linearview_type;
        typedef SmallMatrixView<T,M,N,mstepi,mstepj,false,I> nonconj_type;
    };

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    class SmallMatrixView :
        public BaseMatrix_Rec_Mutable<SmallMatrixView<T,M,N,Si,Sj,C,I> >
    {
    public:

        typedef SmallMatrixView<T,M,N,Si,Sj,C,I> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { mcolsize = Traits<type>::mcolsize };
        enum { mrowsize = Traits<type>::mrowsize };
        enum { mshape = Traits<type>::mshape };
        enum { mfort = Traits<type>::mfort };
        enum { mrowmajor = Traits<type>::mrowmajor };
        enum { mcolmajor = Traits<type>::mcolmajor };
        enum { mstor = Traits<type>::mstor };
        enum { mcalc = Traits<type>::mcalc };
        enum { mstepi = Traits<type>::mstepi };
        enum { mstepj = Traits<type>::mstepj };
        enum { mdiagstep = Traits<type>::mdiagstep };
        enum { mconj = Traits<type>::mconj };
        enum { mcanlin = Traits<type>::mcanlin };

        //
        // Constructors
        //

        inline SmallMatrixView(T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

        inline SmallMatrixView(T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        { TMVStaticAssert(Sj != UNKNOWN); }

        inline SmallMatrixView(T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        inline SmallMatrixView(T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline SmallMatrixView(T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(M != UNKNOWN); TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        inline SmallMatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        inline SmallMatrixView(MatrixView<T,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        inline SmallMatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        inline ~SmallMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

        //
        // Op=
        //

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        inline type& operator=(const type& m2)
        { base_mut::operator=(m2); return *this; }

        inline type& operator=(const T x)
        { base_mut::operator=(x); return *this; }

        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const
        { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

        inline reference ref(int i, int j)
        { return reference(itsm[i*stepi()+j*stepj()]); }

        inline size_t ls() const { return int(itscs)*int(itsrs); }
        inline size_t colsize() const { return itscs; }
        inline size_t rowsize() const { return itsrs; }
        inline int stepi() const { return itssi; }
        inline int stepj() const { return itssj; }
        inline bool isconj() const { return C; }
        inline bool isrm() const 
        { return mrowmajor || (!mcolmajor && stepj() == 1); }
        inline bool iscm() const 
        { return mcolmajor || (!mrowmajor && stepi() == 1); }

    private :

#ifdef TMV_DEBUG
        T* itsm;
#else
        T*const itsm;
#endif
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

        inline SmallMatrixViewF(T* m, size_t cs, size_t rs, int si, int sj) :
            mtype(m,cs,rs,si,sj) {}
        inline SmallMatrixViewF(T* m, size_t cs, size_t rs, int si) :
            mtype(m,cs,rs,si) {}
        inline SmallMatrixViewF(T* m, size_t cs, size_t rs) : mtype(m,cs,rs) {}
        inline SmallMatrixViewF(T* m, size_t cs) : mtype(m,cs) {}
        inline SmallMatrixViewF(T* m) : mtype(m) {}
        inline SmallMatrixViewF(const type& m2) : mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        inline SmallMatrixViewF(MatrixView<T,Si2,Sj2,C,I2> m2) : mtype(m2) {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        inline SmallMatrixViewF(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
            mtype(m2) {}
        inline ~SmallMatrixViewF() {}

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2);        return *this; }
        inline type& operator=(const type& m2)
        { mtype::operator=(m2);        return *this; }
        inline type& operator=(const T x)
        { mtype::operator=(x);        return *this; }

    }; // SmallMatrixViewF


    // Special Creators: 
    //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
    //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage

#define VT typename V::value_type
#define VN V::vsize
#define VS V::vstep
#define VC V::visconj
#define VI (V::vfort ? FortranStyle : CStyle)

    template <class V> 
    inline ConstSmallMatrixView<VT,1,VN,VN,VS,VC,VI> RowVectorViewOf(
        const BaseVector_Calc<V>& v)
    {
        return ConstSmallMatrixView<VT,1,VN,VN,VS,VC,VI>(
            v.cptr(),1,v.size(),v.size(),v.step());
    }

    template <class V> 
    inline SmallMatrixView<VT,1,VN,VN,VS,VC,VI> RowVectorViewOf(
        BaseVector_Mutable<V>& v)
    {
        return SmallMatrixView<VT,1,VN,VN,VS,VC,VI>(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class V> 
    inline ConstSmallMatrixView<VT,VN,1,VS,VN,VC,VI> ColVectorViewOf(
        const BaseVector_Calc<V>& v)
    {
        return ConstSmallMatrixView<VT,VN,1,VS,VN,VC,VI>(
            v.cptr(),v.size(),1,v.step(),v.size());
    }

    template <class V> 
    inline SmallMatrixView<VT,VN,1,VS,VN,VC,VI> ColVectorViewOf(
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
    inline SmallMatrixView<T,1,UNKNOWN,UNKNOWN,S,C,I> RowVectorViewOf(
        VectorView<T,S,C,I> v)
    {
        return SmallMatrixView<T,1,UNKNOWN,UNKNOWN,S,C,I>(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int N, int S, bool C, IndexStyle I> 
    inline SmallMatrixView<T,1,N,N,S,C,I> RowVectorViewOf(
        SmallVectorView<T,N,S,C,I> v)
    {
        return SmallMatrixView<T,1,N,N,S,C,I>(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int S, bool C, IndexStyle I> 
    inline SmallMatrixView<T,UNKNOWN,1,S,UNKNOWN,C,I> ColVectorViewOf(
        VectorView<T,S,C,I> v)
    {
        return SmallMatrixView<T,UNKNOWN,1,S,UNKNOWN,C,I>(
            v.ptr(),v.size(),1,v.step(),v.size());
    }

    template <class T, int N, int S, bool C, IndexStyle I> 
    inline SmallMatrixView<T,N,1,S,N,C,I> ColVectorViewOf(
        SmallVectorView<T,N,S,C,I> v)
    {
        return SmallMatrixView<T,N,1,S,N,C,I>(
            v.ptr(),v.size(),1,v.step(),v.size());
    }


    //
    // Swap
    //

    template <class T, int M, int N, int Si, int Sj, bool C, 
              IndexStyle I, class MM>
    inline void Swap(
        BaseMatrix_Rec_Mutable<MM>& m1, SmallMatrixView<T,M,N,Si,Sj,C,I> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si, int Sj, bool C, 
              IndexStyle I, class MM>
    inline void Swap(
        SmallMatrixView<T,M,N,Si,Sj,C,I> m1, BaseMatrix_Rec_Mutable<MM>& m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,C1,I1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,C1,I1> m1, 
        MatrixView<T,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, bool C1, 
              IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    inline void Swap(
        MatrixView<T,Si1,Sj1,C1,I1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }


    // 
    // TMV_Text
    //

    template <class T, int M, int N, StorageType S, IndexStyle I>
    inline std::string TMV_Text(const SmallMatrix<T,M,N,S,I>& )
    {
        std::ostringstream s;
        s << "SmallMatrix<"<<TMV_Text(T())<<','<<M<<','<<N;
        s <<','<<TMV_Text(S)<<','<<TMV_Text(I)<<">"; 
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, bool C, IndexStyle I>
    inline std::string TMV_Text(
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
    inline std::string TMV_Text(
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
