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
// This file defines the TMV SimpleMatrix class.
//
// It is used when I want some of the simpler functionality of 
// the Matrix class without some of the overhead that the regular 
// Matrix class has.
//
// So this class is completely header-only.  It doesn't have any 
// any of the arithmetic operations defined, so everything has to 
// be done on the elements directly.  Similary, it doesn't have
// all of the normal functions and methods.  Just a few that are easy
// to implement inline.
//
// Constructors:
//
//    SimpleMatrix<T,stor>(int cs, int rs)
//        Makes a SimpleMatrix with column size = cs and row size = rs
//        with _uninitialized_ values
//
//    SimpleMatrix<T,stor>(int cs, int rs, T x)
//        Makes a SimpleMatrix with column size = cs and row size = rs
//        with values initially set to x.
//
//    SimpleMatrix<T,stor>(const BaseMatrix<M2>& m)
//        Make a SimpleMatrix which copies the elements of m.

#ifndef TMV_SimpleMatrix_H
#define TMV_SimpleMatrix_H

#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_BaseMatrix_Tri.h"
#include "tmv/TMV_Array.h"

namespace tmv {

    template <class T, StorageType S=ColMajor>
    class SimpleMatrix;

    template <class T, StorageType S>
    struct Traits<SimpleMatrix<T,S> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SimpleMatrix<T,S> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        typedef InvalidType inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = false };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==ColMajor ? 1 : UNKNOWN) };
        enum { _stepj = (S==RowMajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = false };
        enum { _canlin = true };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { notC = iscomplex };

        enum { _hasdivider = false };

        typedef ConstVectorView<T,_stepi,false> const_col_type;
        typedef ConstVectorView<T,_stepi,false> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,false> const_row_type;
        typedef ConstVectorView<T,_stepj,false> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,false> const_diag_type;
        typedef ConstVectorView<T,_diagstep,false> const_diag_sub_type;

        typedef ConstMatrixView<T,_stepi,_stepj,false> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false>
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,false>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,false>
            const_rowpair_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false> const_colrange_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false> const_rowrange_type;

        typedef ConstMatrixView<T,_stepi,_stepj,false> const_view_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,CStyle> const_cview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,FortranStyle>
            const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstMatrixView<T,1,_stepj,false> const_cmview_type;
        typedef ConstMatrixView<T,_stepi,1,false> const_rmview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,notC> const_conjugate_type;
        typedef ConstMatrixView<T,_stepj,_stepi,false> const_transpose_type;
        typedef ConstMatrixView<T,_stepj,_stepi,notC> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false>
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,_stepj,false>
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,false>
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false>
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,_stepj,false>
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,false>
            const_unknown_lowertri_type;
        typedef ConstMatrixView<real_type,twoSi,twoSj,false>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,1,false> const_linearview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false> const_nonconj_type;
        typedef MatrixView<T,_stepi,_stepj,false> nonconst_type;

        typedef T& reference;

        typedef VectorView<T,_stepi,false> col_type;
        typedef VectorView<T,_stepi,false> col_sub_type;
        typedef VectorView<T,_stepj,false> row_type;
        typedef VectorView<T,_stepj,false> row_sub_type;
        typedef VectorView<T,_diagstep,false> diag_type;
        typedef VectorView<T,_diagstep,false> diag_sub_type;

        typedef MatrixView<T,_stepi,_stepj,false> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,false> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,false> subvector_type;
        typedef SmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,false>
            colpair_type;
        typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,false>
            rowpair_type;
        typedef MatrixView<T,_stepi,_stepj,false> colrange_type;
        typedef MatrixView<T,_stepi,_stepj,false> rowrange_type;

        typedef MatrixView<T,_stepi,_stepj,false> view_type;
        typedef MatrixView<T,_stepi,_stepj,false,CStyle> cview_type;
        typedef MatrixView<T,_stepi,_stepj,false,FortranStyle> fview_type;
        typedef MatrixView<T> xview_type;
        typedef MatrixView<T,1,_stepj,false> cmview_type;
        typedef MatrixView<T,_stepi,1,false> rmview_type;
        typedef MatrixView<T,_stepi,_stepj,notC> conjugate_type;
        typedef MatrixView<T,_stepj,_stepi,false> transpose_type;
        typedef MatrixView<T,_stepj,_stepi,notC> adjoint_type;
        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false>
            uppertri_type;
        typedef UpperTriMatrixView<T,UnitDiag,_stepi,_stepj,false>
            unit_uppertri_type;
        typedef UpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,false>
            unknown_uppertri_type;
        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false>
            lowertri_type;
        typedef LowerTriMatrixView<T,UnitDiag,_stepi,_stepj,false>
            unit_lowertri_type;
        typedef LowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,false>
            unknown_lowertri_type;
        typedef MatrixView<real_type,twoSi,twoSj,false> realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,1,false> linearview_type;
        typedef MatrixView<T,_stepi,_stepj,false> nonconj_type;
    };

    template <class T, StorageType S>
    class SimpleMatrix :
        public BaseMatrix<SimpleMatrix<T,S> >
    {
        typedef SimpleMatrix<T,S> type;
        typedef BaseMatrix<type> base;

    public:

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

        typedef typename Traits<type>::value_type value_type;
        typedef typename Traits<type>::calc_type calc_type;
        typedef typename Traits<type>::eval_type eval_type;
        typedef typename Traits<type>::copy_type copy_type;

        typedef typename Traits<type>::const_row_type const_row_type;
        typedef typename Traits<type>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<type>::const_col_type const_col_type;
        typedef typename Traits<type>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<type>::const_diag_type const_diag_type;
        typedef typename Traits<type>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<type>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<type>::const_submatrix_step_type
            const_submatrix_step_type;
        typedef typename Traits<type>::const_subvector_type const_subvector_type;
        typedef typename Traits<type>::const_colpair_type const_colpair_type;
        typedef typename Traits<type>::const_rowpair_type const_rowpair_type;
        typedef typename Traits<type>::const_colrange_type const_colrange_type;
        typedef typename Traits<type>::const_rowrange_type const_rowrange_type;

        typedef typename Traits<type>::const_view_type const_view_type;
        typedef typename Traits<type>::const_cview_type const_cview_type;
        typedef typename Traits<type>::const_fview_type const_fview_type;
        typedef typename Traits<type>::const_xview_type const_xview_type;
        typedef typename Traits<type>::const_cmview_type const_cmview_type;
        typedef typename Traits<type>::const_rmview_type const_rmview_type;
        typedef typename Traits<type>::const_transpose_type const_transpose_type;
        typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<type>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<type>::const_uppertri_type const_uppertri_type;
        typedef typename Traits<type>::const_unit_uppertri_type
            const_unit_uppertri_type;
        typedef typename Traits<type>::const_unknown_uppertri_type
            const_unknown_uppertri_type;
        typedef typename Traits<type>::const_lowertri_type const_lowertri_type;
        typedef typename Traits<type>::const_unit_lowertri_type
            const_unit_lowertri_type;
        typedef typename Traits<type>::const_unknown_lowertri_type
            const_unknown_lowertri_type;
        typedef typename Traits<type>::const_linearview_type
            const_linearview_type;
        typedef typename Traits<type>::const_realpart_type const_realpart_type;
        typedef typename Traits<type>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<type>::const_nonconj_type const_nonconj_type;

        typedef typename Traits<type>::row_type row_type;
        typedef typename Traits<type>::row_sub_type row_sub_type;
        typedef typename Traits<type>::col_type col_type;
        typedef typename Traits<type>::col_sub_type col_sub_type;
        typedef typename Traits<type>::diag_type diag_type;
        typedef typename Traits<type>::diag_sub_type diag_sub_type;

        typedef typename Traits<type>::submatrix_type submatrix_type;
        typedef typename Traits<type>::submatrix_step_type submatrix_step_type;
        typedef typename Traits<type>::subvector_type subvector_type;
        typedef typename Traits<type>::colpair_type colpair_type;
        typedef typename Traits<type>::rowpair_type rowpair_type;
        typedef typename Traits<type>::colrange_type colrange_type;
        typedef typename Traits<type>::rowrange_type rowrange_type;

        typedef typename Traits<type>::view_type view_type;
        typedef typename Traits<type>::cview_type cview_type;
        typedef typename Traits<type>::fview_type fview_type;
        typedef typename Traits<type>::xview_type xview_type;
        typedef typename Traits<type>::cmview_type cmview_type;
        typedef typename Traits<type>::rmview_type rmview_type;
        typedef typename Traits<type>::transpose_type transpose_type;
        typedef typename Traits<type>::conjugate_type conjugate_type;
        typedef typename Traits<type>::adjoint_type adjoint_type;
        typedef typename Traits<type>::uppertri_type uppertri_type;
        typedef typename Traits<type>::unit_uppertri_type unit_uppertri_type;
        typedef typename Traits<type>::unknown_uppertri_type
            unknown_uppertri_type;
        typedef typename Traits<type>::lowertri_type lowertri_type;
        typedef typename Traits<type>::unit_lowertri_type unit_lowertri_type;
        typedef typename Traits<type>::unknown_lowertri_type
            unknown_lowertri_type;
        typedef typename Traits<type>::linearview_type linearview_type;
        typedef typename Traits<type>::realpart_type realpart_type;
        typedef typename Traits<type>::imagpart_type imagpart_type;
        typedef typename Traits<type>::nonconj_type nonconj_type;

        typedef typename Traits<type>::reference reference;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;
        typedef typename Traits<value_type>::complex_type complex_type;


        //
        // Constructors
        //

        SimpleMatrix() : itscs(0), itsrs(0), linsize(0), itsm(0) 
        { TMVAssert(S==RowMajor || S==ColMajor); }

        SimpleMatrix(int cs, int rs) :
            itscs(cs), itsrs(rs),
            linsize((cs)*(rs)), itsm(linsize)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_DEBUG
            setAllTo(T(888));
#endif
        }

        SimpleMatrix(int cs, int rs, T x) :
            itscs(cs), itsrs(rs),
            linsize((cs)*(rs)), itsm(linsize)
        { setAllTo(x); }

        SimpleMatrix(const type& m2) :
            itscs(m2.itscs), itsrs(m2.itsrs),
            linsize(m2.linsize), itsm(linsize)
        { for(int i=0;i<linsize;++i) itsm[i] = m2.itsm[i]; }

        template <class M2>
        SimpleMatrix(const BaseMatrix<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()),
            linsize(itscs * itsrs), itsm(linsize)
        {
            typename M2::calc_type m2c = m2.calc();
            for(int i=0;i<itscs;++i) for(int j=0;j<itsrs;++j) 
                ref(i,j) = m2c.cref(i,j);
        }

        ~SimpleMatrix()
        {
#ifdef TMV_DEBUG
            setAllTo(T(999));
#endif
        }


        //
        // Access
        //

        T operator()(int i,int j) const
        {
            TMVAssert(i>=0 && i<itscs);
            TMVAssert(j>=0 && j<itsrs);
            return cref(i,j); 
        }

        reference operator()(int i,int j) 
        {
            TMVAssert(i>=0 && i<itscs);
            TMVAssert(j>=0 && j<itsrs);
            return ref(i,j);
        }

        const_row_type row(int i) const 
        { return row_type(ptr()+i*stepi(),rowsize(),stepj()); }

        const_col_type col(int j) const 
        { return col_type(ptr()+j*stepj(),colsize(),stepi()); }

        // No need for a get_ routine for diag()
        const_diag_type diag() const 
        { return diag_type(ptr(),TMV_MIN(colsize(),rowsize()),diagstep()); }

        const_row_sub_type row(int i, int j1, int j2) const 
        { return row_sub_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

        const_col_sub_type col(int j, int i1, int i2) const 
        { return col_sub_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

        const_diag_sub_type diag(int i) const 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ?
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i) ),
                diagstep());
        }

        const_diag_sub_type diag(int i, int j1, int j2) const 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }

        const_row_type operator[](int i) const 
        { return row(i); }

        row_type row(int i)
        { return row_type(ptr()+i*stepi(),rowsize(),stepj()); }

        col_type col(int j)
        { return col_type(ptr()+j*stepj(),colsize(),stepi()); }

        // No need for a get_ routine for diag()
        diag_type diag()
        { return diag_type(ptr(),TMV_MIN(colsize(),rowsize()),diagstep()); }

        row_sub_type row(int i, int j1, int j2)
        { return row_sub_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

        col_sub_type col(int j, int i1, int i2)
        { return col_sub_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

        diag_sub_type diag(int i)
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ?
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i) ),
                diagstep());
        }

        diag_sub_type diag(int i, int j1, int j2)
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }

        row_type operator[](int i)
        { return row(i); }


        //
        // Op=
        //

        type& operator=(const type& m2)
        {
            if (&m2 != this) {
                for(int i=0;i<linsize;++i) itsm[i] = m2.itsm[i];
            }
            return *this; 
        }

        template <class M2>
        type& operator=(const M2& m2)
        {
            for(int i=0;i<itscs;++i) for(int j=0;j<itsrs;++j) 
                ref(i,j) = m2.cref(i,j);
            return *this; 
        }

        type& operator=(const T& x) 
        { return setToIdentity(x); }

        typedef typename linearview_type::iterator lin_iter;
        ListAssigner<value_type,lin_iter> operator<<(value_type x)
        {
            linearview_type v = linearView();
            return ListAssigner<value_type,lin_iter>(v.begin(),v.size(),x);
        }


        //
        // Functions of Matrix
        //

        T trace() const
        {
            TMVAssert(itscs == itsrs);
            T sum(0);
            for(int i=0; i<itsrs; ++i) sum += cref(i,i);
            return sum;
        }

        T sumElements() const
        {
            T sum(0);
            for(int i=0;i<linsize;++i) sum += itsm[i];
            return sum;
        }

        real_type sumAbsElements() const
        {
            real_type sum(0);
            for(int i=0;i<linsize;++i) sum += TMV_ABS(itsm[i]);
            return sum;
        }

        real_type sumAbs2Elements() const
        {
            real_type sum(0);
            for(int i=0;i<linsize;++i) sum += TMV_ABS2(itsm[i]);
            return sum;
        }

        real_type norm() const 
        { return normF(); }

        real_type normF() const
        {
            TMVAssert(!std::numeric_limits<T>::is_integer);
            return TMV_SQRT(normSq()); 
        }

        real_type normSq(real_type scale=real_type(1)) const
        {
            real_type sum(0);
            if (scale == real_type(1))
                for(int i=0;i<linsize;++i) sum += TMV_NORM(itsm[i]);
            else
                for(int i=0;i<linsize;++i) sum += TMV_NORM(itsm[i]*scale);
            return sum;
        }

        real_type norm1() const
        {
            real_type max(0);
            for(int j=0;j<itsrs;++j) {
                real_type temp(0);
                for(int i=0;i<itscs;++i) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // inf-norm = max_i (sum_j |a_ij|)
        real_type normInf() const
        {
            real_type max(0);
            for(int i=0;i<itscs;++i) {
                real_type temp(0);
                for(int j=0;j<itsrs;++j) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|a_ij|)
        real_type maxAbsElement() const
        {
            real_type max(0);
            for(int i=0;i<linsize;++i) {
                real_type temp = TMV_ABS(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|real(a_ij)|+|imag(a_ij)|)
        real_type maxAbs2Element() const
        {
            real_type max(0);
            for(int i=0;i<linsize;++i) {
                real_type temp = TMV_ABS2(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        //
        // Modifying Functions
        //

        type& setZero() 
        {
            for(int i=0;i<linsize;++i) itsm[i] = T(0);
            return *this;
        }

        type& clip(real_type thresh)
        {
            for(int i=0;i<linsize;++i) 
                if (TMV_ABS(itsm[i]) < thresh) itsm[i] = T(0);
            return *this;
        }

        type& setAllTo(const T& x) 
        {
            for(int i=0;i<linsize;++i) itsm[i] = x;
            return *this;
        }

        type& addToAll(const T& x) 
        {
            for(int i=0;i<linsize;++i) itsm[i] += x;
            return *this;
        }

        type& transposeSelf() 
        {
            TMVAssert(itscs == itsrs);
            for(int i=1; i<itscs; ++i) 
                for(int j=0; j<i; ++j) TMV_SWAP(ref(i,j),ref(j,i));
            return *this;
        }

        type& conjugateSelf() 
        {
            if (isComplex(T())) {
                for(int i=0;i<linsize;++i) itsm[i] = TMV_CONJ(itsm[i]);
            }
            return *this;
        }

        type& setToIdentity(const T& x=T(1)) 
        {
            TMVAssert(itscs == itsrs);
            setZero();
            for(int i=0; i<itsrs; ++i) ref(i,i) = T(1);
            return *this;
        }

        type& swapRows(int i1, int i2)
        {
            TMVAssert(i1 >= 0 && i1 < itscs);
            TMVAssert(i2 >= 0 && i2 < itscs);
            if (i1 != i2)
                for(int j=0; j<itsrs; ++j) TMV_SWAP(ref(i1,j),ref(i2,j));
            return *this;
        }

        type& swapCols(int j1, int j2)
        {
            TMVAssert(j1 >= 0 && j1 < itsrs);
            TMVAssert(j2 >= 0 && j2 < itsrs);
            if (j1 != j2)
                for(int i=0; i<itscs; ++i) TMV_SWAP(ref(i,j1),ref(i,j2));
            return *this;
        }


        //
        // subMatrix, etc.
        //

        const_submatrix_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            return const_submatrix_type(
                cptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj());
        }

        const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            return const_submatrix_step_type(
                cptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        const_subvector_type subVector( 
            int i, int j, int istep, int jstep, int s) const
        {
            return const_subvector_type(
                cptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        const_colpair_type colPair(int j1, int j2) const
        {
            return const_colpair_type(
                cptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
        }

        const_rowpair_type rowPair(int i1, int i2) const
        {
            return const_rowpair_type(
                cptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
        }

        const_colrange_type colRange(int j1, int j2) const
        {
            return const_colrange_type(
                cptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
        }

        const_rowrange_type rowRange(int i1, int i2) const
        {
            return const_rowrange_type(
                cptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
        }

        submatrix_type subMatrix(int i1, int i2, int j1, int j2)
        {
            return submatrix_type(
                ptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj());
        }

        submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep)
        {
            return submatrix_step_type(
                ptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        subvector_type subVector(int i, int j, int istep, int jstep, int s)
        {
            return subvector_type(
                ptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        colpair_type colPair(int j1, int j2)
        {
            return colpair_type(
                ptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
        }

        rowpair_type rowPair(int i1, int i2)
        {
            return rowpair_type(
                ptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
        }

        colrange_type colRange(int j1, int j2)
        {
            return colrange_type(
                ptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
        }

        rowrange_type rowRange(int i1, int i2)
        {
            return rowrange_type(
                ptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
        }


        //
        // Views
        //

        const_view_type view() const
        { return const_view_type(cptr(),colsize(),rowsize(),stepi(),stepj()); }

        const_cview_type cView() const
        { return view(); }

        const_fview_type fView() const
        { return view(); }

        const_xview_type xView() const
        { return view(); }

        const_cmview_type cmView() const
        { return view(); }

        const_rmview_type rmView() const
        { return view(); }

        const_view_type constView() const
        { return view(); }

        const_transpose_type transpose() const
        {
            return const_transpose_type(
                cptr(),rowsize(),colsize(),stepj(),stepi()); 
        }

        const_conjugate_type conjugate() const
        {
            return const_conjugate_type(
                cptr(),colsize(),rowsize(),stepi(),stepj()); 
        }

        const_adjoint_type adjoint() const
        {
            return const_adjoint_type(
                cptr(),rowsize(),colsize(),stepj(),stepi()); 
        }

        const_uppertri_type upperTri() const
        { return const_uppertri_type(cptr(),rowsize(),false,stepi(),stepj()); }

        const_unit_uppertri_type unitUpperTri() const
        {
            return const_unit_uppertri_type(
                cptr(),rowsize(),true,stepi(),stepj()); 
        }

        const_unknown_uppertri_type upperTri(DiagType dt) const
        {
            return const_unknown_uppertri_type(
                cptr(),rowsize(),dt==UnitDiag,stepi(),stepj());
        }

        const_lowertri_type lowerTri() const
        { return const_lowertri_type(cptr(),colsize(),false,stepi(),stepj()); }

        const_unit_lowertri_type unitLowerTri() const
        {
            return const_unit_lowertri_type(
                cptr(),colsize(),true,stepi(),stepj()); 
        }

        const_unknown_lowertri_type lowerTri(DiagType dt) const
        {
            return const_unknown_lowertri_type(
                cptr(),rowsize(),dt==UnitDiag,stepi(),stepj());
        }

        const_linearview_type linearView() const
        { return const_linearview_type(cptr(),ls(),1); }

        const_realpart_type realPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            return const_realpart_type(
                reinterpret_cast<real_type*>(cptr()), colsize(), rowsize(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        const_imagpart_type imagPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<real_type*>(cptr())+1, colsize(), rowsize(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        const_nonconj_type nonConj() const
        {
            return const_nonconj_type(
                cptr(),colsize(),rowsize(),stepi(),stepj()); 
        }

        view_type view()
        { return view_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }

        cview_type cView()
        { return view(); }

        fview_type fView()
        { return view(); }

        xview_type xView()
        { return view(); }

        cmview_type cmView()
        { return view(); }

        rmview_type rmView()
        { return view(); }

        transpose_type transpose()
        { return transpose_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

        conjugate_type conjugate()
        { return conjugate_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }

        adjoint_type adjoint()
        { return adjoint_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

        uppertri_type upperTri()
        { return uppertri_type(ptr(),rowsize(),false,stepi(),stepj()); }

        unit_uppertri_type unitUpperTri()
        { return unit_uppertri_type(ptr(),rowsize(),true,stepi(),stepj()); }

        unknown_uppertri_type upperTri(DiagType dt)
        {
            return unknown_uppertri_type(
                ptr(),rowsize(),dt==UnitDiag,stepi(),stepj());
        }

        lowertri_type lowerTri()
        { return lowertri_type(ptr(),colsize(),false,stepi(),stepj()); }

        unit_lowertri_type unitLowerTri()
        { return unit_lowertri_type(ptr(),colsize(),true,stepi(),stepj()); }

        unknown_lowertri_type lowerTri(DiagType dt)
        {
            return unknown_lowertri_type(
                ptr(),rowsize(),dt==UnitDiag,stepi(),stepj());
        }

        linearview_type linearView()
        { return linearview_type(ptr(),ls(),1); }

        realpart_type realPart()
        {
            const bool isreal = Traits<value_type>::isreal;
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()), colsize(), rowsize(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        imagpart_type imagPart()
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1, colsize(), rowsize(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        nonconj_type nonConj()
        { return nonconj_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }


        //
        // I/O
        //

        void write(std::ostream& os) const
        {
            os << itscs <<"  "<< itsrs <<std::endl;
            for(int i=0;i<itscs;++i) {
                os << "( ";
                for(int j=0;j<itsrs;++j) os << ' '<<cref(i,j)<<' ';
                os << " )\n";
            }
        }

        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        bool isSquare() const { return itscs == itsrs; }
        int stepi() const { return S == RowMajor ? itsrs : 1; }
        int stepj() const { return S == RowMajor ? 1 : itscs; }
        int diagstep() const { return S == RowMajor ? itsrs+1 : itscs+1; }
        bool isrm() const { return S == RowMajor; }
        bool iscm() const { return S == ColMajor; }
        bool isconj() const { return false; }
        StorageType stor() const { return S; }
        size_t ls() const { return linsize; }
        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return S == RowMajor ? itsm[i*itsrs+j] : itsm[j*itscs+i]; }

        reference ref(int i, int j)
        { return S == RowMajor ? itsm[i*itsrs+j] : itsm[j*itscs+i]; }

        void resize(size_t cs, size_t rs)
        {
            linsize = cs*rs;
            itsm.resize(linsize);
            itscs = cs;
            itsrs = rs;
#ifdef TMV_DEBUG
            setAllTo(T(888));
#endif
        }

    protected :

        int itscs;
        int itsrs;
        int linsize;
        AlignedArray<T> itsm;

    }; // SimpleMatrix

    template <class T, StorageType S>
    static inline std::ostream& operator<<(
        std::ostream& os, const SimpleMatrix<T,S>& m)
    { m.write(os); return os; }

    template <class T, StorageType S>
    static inline std::string TMV_Text(const SimpleMatrix<T,S>& )
    {
        std::ostringstream s;
        s << std::string("SimpleMatrix<")<<TMV_Text(T());
        s << ','<<TMV_Text(S)<<">"; 
        return s.str();
    }

} // namespace tmv

#endif
