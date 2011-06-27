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

//-----------------------------------------------------------------------------
//
// This file defines the BaseMatrix_Rec and BaseMatrix_Rec_Mutable
// classes.  These are the base classes for the dense rectangular
// matrices, Matrix and SmallMatrix.
//
// BaseMatrix_Rec is the base class for all dense rectangular matrices.
// It adds the rest of the views described in TMV_Matrix.h like
// subMatrix, subVector, rowRange, colRange, etc. that aren't defined for 
// a general BaseMatrix.
//
// BaseMatrix_Rec_Mutable is the base class for mutable dense 
// rectangular matrices and adds the non-const views along with the
// other modifying functions like transposeSelf(), setToIdentity(), etc.

#ifndef TMV_BaseMatrix_Rec_H
#define TMV_BaseMatrix_Rec_H

#include "TMV_BaseMatrix.h"

namespace tmv {

    // BaseMatrix_Rec is derived from BaseMatrix_Calc, and is used
    // for dense rectangular matrices.
    template <class M>
    class BaseMatrix_Rec;

    // BaseMatrix_Rec adds the following requirements to Traits<M>:
    //
    //  _stepi = the step size along column if known (else UNKNOWN)
    //  _stepj = the step size along row if known (else UNKNOWN)
    //  _diagstep = the step size along row if known (else UNKNOWN)
    //
    //  const_row_type = return type from row(i) const
    //  const_col_type = return type from col(j) const
    //  const_diag_type = return type from diag() const
    //  const_row_sub_type = return type from row(i,j1,j2) const
    //  const_col_sub_type = return type from col(j,i1,i2) const
    //  const_diag_sub_type = return type from diag(i) and diag(i,j1,j2) const
    //
    //  const_submatrix_type = return type from subMatrix(i1,i2,j1,j2) const
    //  const_submatrix_step_type = return type from 
    //      subMatrix(i1,i2,j1,j2,istep,jstep) const
    //  const_subvector_type = return type from 
    //      subVector(i,j,istep,jstep,size) const
    //  const_colpair_type = return type from colPair(j1,j2) const
    //  const_rowpair_type = return type from rowPair(i1,i2) const
    //  const_colrange_type = return type from colRange(j1,j2) const
    //  const_rowrange_type = return type from rowRange(i1,i2) const
    //
    //  const_uppertri_type = return type from upperTri() const
    //  const_unit_uppertri_type = return type from unitUpperTri() const
    //  const_unknown_uppertri_type = return type from upperTri(dt) const
    //  const_lowertri_type = return type from lowerTri() const
    //  const_unit_lowertri_type = return type from unitLowerTri() const
    //  const_unknown_lowertri_type = return type from lowerTri(dt) const
    //  const_linearview_type = return type from linearView() const
    //  


    // BaseMatrix_Rec_Mutable is derived from both BaseMatrix_Rec and
    // BaseMatrix_Mutable, and is used for mutable, dense, rectangular matrices.
    template <class M>
    class BaseMatrix_Rec_Mutable;

    // BaseMatrix_Rec_Mutable adds:
    //
    //  col_type = the return type of col(i) 
    //  row_type = the return type of row(j) 
    //  diag_type = return type from diag() 
    //  row_sub_type = return type from row(i,j1,j2) 
    //  col_sub_type = return type from col(j,i1,i2) 
    //  diag_sub_type = return type from diag(i) and diag(i,j1,j2) 
    //
    //  submatrix_type = return type from subMatrix(i1,i2,j1,j2) 
    //  submatrix_step_type = return type from 
    //      subMatrix(i1,i2,j1,j2,istep,jstep) 
    //  subvector_type = return type from subVector(i,j,istep,jstep,size)
    //  colpair_type = return type from colPair(j1,j2) 
    //  rowpair_type = return type from rowPair(i1,i2) 
    //  colrange_type = return type from colRange(j1,j2) 
    //  rowrange_type = return type from rowRange(i1,i2) 
    //
    //  uppertri_type = return type from upperTri() 
    //  unit_uppertri_type = return type from unitUpperTri() 
    //  unknown_uppertri_type = return type from upperTri(dt) 
    //  lowertri_type = return type from lowerTri() 
    //  unit_lowertri_type = return type from unitLowerTri() 
    //  unknown_lowertri_type = return type from lowerTri(dt) 
    //  realpart_type = return type from realPart()
    //  imagpart_type = return type from imagPart() 
    //  nonconj_type = return type from nonConj() 
    //  linearview_type = return type from linearView() 
    //  

    // The following all derive from BaseMatrix_Rec or BaseMatrix_Rec_Mutable.
    // See TMV_Matrix.h and TMV_SmallMatrix.h for their definitions:
    template <class T, int A=0, int A2=0>
    class Matrix;
    template <class T, int A=0>
    class ConstMatrixView;
    template <class T, int A=0>
    class MatrixView;
    template <class T, int M, int N, int A=0, int A2=0> 
    class SmallMatrix;
    template <class T, int M, int N, int Si=UNKNOWN, int Sj=UNKNOWN, int A=0>
    class ConstSmallMatrixView;
    template <class T, int M, int N, int Si=UNKNOWN, int Sj=UNKNOWN, int A=0>
    class SmallMatrixView;

    // In TMV_Norm.h
    template <class M>
    static inline typename M::float_type DoNormF(const BaseMatrix_Rec<M>& m);

    //
    // Helper functions and values:
    //

    // Specify ExactSameStorage for rectangular matrices:
    template <class M1, class M2>
    static TMV_INLINE bool ExactSameStorage(
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()); 
    }

    // These helper functions check the validity of indices according
    // to whether the matrix uses CStyle or FortranStyle indexing.
    // They also update the indices to be consistent with CStyle.
    template <bool _fort>
    static TMV_INLINE_ND void CheckDiagIndex(int& i, int m, int n)
    { // CStyle or FortranStyle
        TMVAssert(i >= -m && "negative diag index must be <= nrows");
        TMVAssert(i <= n && "positive diag index must be <= ncols");
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckDiagIndex(int& i, int& j1, int& j2, int m, int n)
    { // CStyle
        TMVAssert(i >= -m && "negative diag index must be <= nrows");
        TMVAssert(i <= n && "positive diag index must be <= ncols");
        TMVAssert(j1 >= 0 && "first element must be in matrix");
        TMVAssert(j2 >= j1 && 
                  "range must have a non-negative number of elements");
        // if i >= 0, i + j2 <= n  and j2 <= m
        TMVAssert((i < 0 || i+j2 <= n) && "last element must be in matrix");
        TMVAssert((i < 0 || j2 <= m) && "last element must be in matrix");
        // if i <= 0, -i + j2 <= m  and j2 <= n
        TMVAssert((i > 0 || -i+j2 <= m) && "last element must be in matrix");
        TMVAssert((i > 0 || j2 <= n) && "last element must be in matrix");
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckRowRange(int& i1, int i2, int m)
    { // CStyle
        TMVAssert(i1 >= 0 && "first row must be in matrix");
        TMVAssert(i2 <= m && "last row must be in matrix");
        TMVAssert(i2 >= i1 && "range must have a non-negative number of rows");
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckColRange(int& j1, int j2, int n)
    { // CStyle
        TMVAssert(j1 >= 0 && "first column must be in matrix");
        TMVAssert(j2 <= n && "last column must be in matrix");
        TMVAssert(j2 >= j1 && "range must have a non-negative number of columns");
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckRowRange(int& i1, int& i2, int istep, int m)
    { // CStyle
        TMVAssert(istep != 0 && "istep cannot be 0");
        TMVAssert(((i1 >= 0 && i1 < m) || i1==i2) && 
                  "first row must be in matrix");
        TMVAssert(((i2-istep >= 0 && i2-istep < m) || i1==i2) && 
                  "last row must be in matrix");
        TMVAssert((i2-i1) % istep == 0 &&
                  "row range must be an integral multiple of istep");
        TMVAssert((i2-i1) / istep >= 0 &&
                  "must have a non-negative number of rows");
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckColRange(int& j1, int& j2, int jstep, int n)
    { // CStyle
        TMVAssert(jstep != 0 && "jstep cannot be 0");
        TMVAssert(((j1 >= 0 && j1 < n) || j1==j2) && 
                  "first column must be in matrix");
        TMVAssert(((j2-jstep >= 0 && j2-jstep < n) || j1==j2) &&
                  "last column must be in matrix");
        TMVAssert((j2-j1) % jstep == 0 &&
                  "column range must be an integral multiple of jstep");
        TMVAssert((j2-j1) / jstep >= 0 &&
                  "must have a non-negative number of columns");
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckMatSubVector(
        int& i, int& j, int istep, int jstep, int size, int m, int n)
    { // CStyle
        TMVAssert(!(istep == 0 && jstep == 0) && 
                  "istep and jstep cannot both be 0");
        TMVAssert(((i >= 0 && i < m) || size==0) && 
                  "first element must be in matrix");
        TMVAssert(((j >= 0 && j < n) || size==0) && 
                  "first element must be in matrix");
        TMVAssert(( (i+istep*(size-1) >= 0 && i+istep*(size-1) < m) || 
                    size==0 ) && 
                  "last element must be in matrix");
        TMVAssert(( (j+jstep*(size-1) >= 0 && j+jstep*(size-1) < n) || 
                    size==0 ) && 
                  "last element must be in matrix");
    }
    template <>
    TMV_INLINE_ND void CheckDiagIndex<true>(
        int& i, int& j1, int& j2, int m, int n)
    { // FortranStyle
        TMVAssert(i >= -m && "negative diag index must be <= nrows");
        TMVAssert(i <= n && "positive diag index must be <= ncols");
        TMVAssert(j1 >= 1 && "first element must be in matrix");
        TMVAssert(j2 >= j1 && "range must have a positive number of elements");
        TMVAssert((i < 0 || i+j2 <= n) && "last element must be in matrix");
        TMVAssert((i < 0 || j2 <= m) && "last element must be in matrix");
        TMVAssert((i > 0 || -i+j2 <= m) && "last element must be in matrix");
        TMVAssert((i > 0 || j2 <= n) && "last element must be in matrix");
        --j1;
    }
    template <>
    TMV_INLINE_ND void CheckRowRange<true>(int& i1, int i2, int m)
    { // FortranStyle
        TMVAssert(i1 >= 1 && "first row must be in matrix");
        TMVAssert(i2 <= m && "last row must be in matrix");
        TMVAssert(i2 >= i1 && "range must have a positive number of rows");
        --i1;
    }
    template <>
    TMV_INLINE_ND void CheckColRange<true>(int& j1, int j2, int n)
    { // FortranStyle
        TMVAssert(j1 >= 1 && "first column must be in matrix");
        TMVAssert(j2 <= n && "last column must be in matrix");
        TMVAssert(j2 >= j1 && "range must have a positive number of columns");
        --j1;
    }
    template <>
    TMV_INLINE_ND void CheckRowRange<true>(int& i1, int& i2, int istep, int m)
    { // FortranStyle
        TMVAssert(istep != 0 && "istep cannot be 0");
        TMVAssert(i1 >= 1 && i1 <= m && "first row must be in matrix");
        TMVAssert(i2 >= 1 && i2 <= m && "last row must be in matrix");
        TMVAssert((i2-i1) % istep == 0 &&
                  "row range must be an integral multiple of istep");
        TMVAssert((i2-i1) / istep >= 0 &&
                  "must have a positive number of rows");
        --i1; i2 += istep-1;
    }
    template <>
    TMV_INLINE_ND void CheckColRange<true>(int& j1, int& j2, int jstep, int n)
    { // FortranStyle
        TMVAssert(jstep != 0 && "jstep cannot be 0");
        TMVAssert(j1 >= 1 && j1 <= n && "first column must be in matrix");
        TMVAssert(j2 >= 1 && j2 <= n && "last column must be in matrix");
        TMVAssert((j2-j1) % jstep == 0 &&
                  "column range must be an integral multiple of jstep");
        TMVAssert((j2-j1) / jstep >= 0 &&
                  "must have a positive number of columns");
        --j1; j2 += jstep-1;
    }
    template <>
    TMV_INLINE_ND void CheckMatSubVector<true>(
        int& i, int& j, int istep, int jstep, int size, int m, int n)
    { // FortranStyle
        TMVAssert(!(istep == 0 && jstep == 0) && 
                  "istep and jstep cannot both be 0");
        TMVAssert(i >= 1 && i <= m && "first element must be in matrix");
        TMVAssert(j >= 1 && j <= n && "first element must be in matrix");
        TMVAssert(i+istep*(size-1) >= 1 && i+istep*(size-1) <= m && 
                  "last element must be in matrix");
        TMVAssert(j+jstep*(size-1) >= 1 && j+jstep*(size-1) <= n && 
                  "last element must be in matrix");
        --i; --j;
    }

    template <class T, int cs, int rs, bool rm, bool fort>
    struct MCopyHelper<T,Rec,cs,rs,rm,fort>
    {
        enum { A2 = (
                (rm ? RowMajor : ColMajor) |
                (fort ? FortranStyle : CStyle) |
                NoDivider | NoAlias ) };
        typedef SmallMatrix<T,cs,rs,A2> type;
    };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,Rec,UNKNOWN,rs,rm,fort>
    {
        enum { A2 = (
                (rm ? RowMajor : ColMajor) |
                (fort ? FortranStyle : CStyle) |
                NoDivider | NoAlias ) };
        typedef Matrix<T,A2> type; 
    };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,Rec,cs,UNKNOWN,rm,fort>
    { 
        enum { A2 = (
                (rm ? RowMajor : ColMajor) |
                (fort ? FortranStyle : CStyle) |
                NoDivider | NoAlias ) };
        typedef Matrix<T,A2> type; 
    };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,Rec,UNKNOWN,UNKNOWN,rm,fort>
    { 
        enum { A2 = (
                (rm ? RowMajor : ColMajor) |
                (fort ? FortranStyle : CStyle) |
                NoDivider | NoAlias ) };
        typedef Matrix<T,A2> type; 
    };


    template <class T, int cs, int rs, int si, int sj, int c>
    struct MViewHelper<T,Rec,cs,rs,si,sj,c>
    { 
        typedef SmallMatrixView<T,cs,rs,si,sj,c> type; 
        typedef ConstSmallMatrixView<T,cs,rs,si,sj,c> ctype; 
    };
    template <class T, int si, int sj, int c>
    struct MViewHelper<T,Rec,UNKNOWN,UNKNOWN,si,sj,c>
    {
        enum { A2 = c | (si == 1 ? ColMajor : sj == 1 ? RowMajor : NonMajor) };
        typedef MatrixView<T,A2> type; 
        typedef ConstMatrixView<T,A2> ctype; 
    };

    // A quick helper class to get the return types correct for
    // RowVectorView and ColVectorView.
    // These are in SmallMatrix.h, but since they are used elsewhere
    // we put the VVO helper here.
    template <class V>
    struct VVO
    {
        typedef typename V::value_type T;
        enum { N = V::_size };
        enum { S = V::_step };
        enum { vecA = (
                ( V::_conj ? Conj : NonConj ) |
                ( V::_fort ? FortranStyle : CStyle ) |
                ( Traits<V>::A & AllAliasStatus ) )};
        enum { colA = vecA | (S == 1 ? ColMajor : 0 ) };
        enum { rowA = vecA | (S == 1 ? RowMajor : 0 ) };
        typedef ConstSmallMatrixView<T,1,N,N,S,rowA> crv;
        typedef SmallMatrixView<T,1,N,N,S,rowA> rv;
        typedef ConstSmallMatrixView<T,N,1,S,N,colA> ccv;
        typedef SmallMatrixView<T,N,1,S,N,colA> cv;
    };


    // A quick auxilliary function for canLinearize.
    // (It only accesses the steps that are unknown at compile time.)
    template <bool canlin, int Si, int Sj, class M>
    struct AuxCanLinearize;

    template <int Si, int Sj, class M>
    struct AuxCanLinearize<true,Si,Sj,M>
    {
        static TMV_INLINE bool ok(const M& m) 
        { return true; }
    };
    template <int Si, int Sj, class M>
    struct AuxCanLinearize<false,Si,Sj,M>
    {
        static inline bool ok(const M& m) 
        {
            return (
                (m.stepi() == 1 && m.stepj() == int(m.colsize())) ||
                (m.stepj() == 1 && m.stepi() == int(m.rowsize())) );
        }
    };
    template <int Sj, class M>
    struct AuxCanLinearize<false,1,Sj,M>
    {
        static inline bool ok(const M& m) 
        { return m.stepj() == int(m.colsize()); } 
    };
    template <int Si, class M>
    struct AuxCanLinearize<false,Si,1,M>
    {
        static inline bool ok(const M& m) 
        { return m.stepi() == int(m.rowsize()); } 
    };
    template <class M>
    struct AuxCanLinearize<false,1,1,M>
    { static TMV_INLINE bool ok(const M& m) { return true; } };

    // Defined in TMV_CopyM.h
    template <class M1, class M2>
    static inline void Copy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2>
    static inline void NoAliasCopy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);

    // Defined in TMV_SwapM.h
    template <class M1, class M2>
    static inline void Swap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2>
    static inline void NoAliasSwap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M>
    static inline void TransposeSelf(BaseMatrix_Rec_Mutable<M>& m);

    // Defined in TMV_Permute.h
    template <class M>
    static inline void PermuteRows(
        BaseMatrix_Rec_Mutable<M>& m, const int* p, int i1, int i2);
    template <class M>
    static inline void ReversePermuteRows(
        BaseMatrix_Rec_Mutable<M>& m, const int* p, int i1, int i2);

    // Defined in TMV_NormM.h
    template <class M>
    static inline typename M::real_type DoNormSq(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::float_type DoNormSq(
        const BaseMatrix_Rec<M>& m, const typename M::float_type scale);
    template <class M>
    static inline typename M::float_type DoNorm1(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::float_type DoNormInf(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::value_type DoSumElements(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::float_type DoSumAbsElements(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::real_type DoSumAbs2Elements(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::float_type DoMaxAbsElement(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::real_type DoMaxAbs2Element(const BaseMatrix_Rec<M>& m);

    // Defined in TMV_MatrixIO.h
    template <class M>
    static inline void Read(std::istream& is, BaseMatrix_Rec_Mutable<M>& m);

    // Defined below:
    template <class M>
    static inline void SetZero(BaseMatrix_Rec_Mutable<M>& m);
    template <class M, class RT>
    static inline void Clip(BaseMatrix_Rec_Mutable<M>& m, const RT& thresh);
    template <class M, class T>
    static inline void SetAllTo(BaseMatrix_Rec_Mutable<M>& m, const T& val);
    template <class M, class T>
    static inline void AddToAll(BaseMatrix_Rec_Mutable<M>& m, const T& val);
    template <class M>
    static inline void ConjugateSelf(BaseMatrix_Rec_Mutable<M>& m);
    template <class M, class F>
    static inline void ApplyToAll(BaseMatrix_Rec_Mutable<M>& m, const F& f);


#if 0
    // Defined in TMV_MatrixDiv.h
    template <class M>
    static inline typename M::float_type DoNorm2(const BaseMatrix_Rec<M>& m);
    template <class M>
    static inline typename M::float_type DoCondition(const BaseMatrix_Rec<M>& m);
#endif

    // A helper class for returning views without necessarily
    // making a new object.
    template <bool ref, class type, class view_type>
    struct MakeRecView_Helper;

    template <class type, class view_type>
    struct MakeRecView_Helper<true,type,view_type>
    {
        typedef type& ret_type;
        typedef const type& const_ret_type;
        static TMV_INLINE ret_type call(type& m) { return m; }
        static TMV_INLINE const_ret_type call(const type& m) { return m; }
    };

    template <class type, class view_type>
    struct MakeRecView_Helper<false,type,view_type>
    {
        typedef view_type ret_type;
        typedef view_type const_ret_type;
        static TMV_INLINE ret_type call(type& m) 
        {
            return view_type(
                m.ptr(),m.colsize(),m.rowsize(),m.stepi(),m.stepj()); 
        }
        static TMV_INLINE const_ret_type call(const type& m) 
        {
            return view_type(
                m.cptr(),m.colsize(),m.rowsize(),m.stepi(),m.stepj()); 
        }
    };

    template <class type, class view_type>
    struct MakeRecView
    {
        enum { ref = Traits2<type,view_type>::sametype };
        typedef MakeRecView_Helper<ref,type,view_type> helper;

        static TMV_INLINE typename helper::ret_type call(type& m)
        { return helper::call(m); }
        static TMV_INLINE typename helper::const_ret_type call(const type& m)
        { return helper::call(m); }
    };

 
    template <class M>
    class BaseMatrix_Rec : 
        public BaseMatrix_Calc<M>
    {
    public:
        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _calc = Traits<M>::_calc };
        enum { _stepi = Traits<M>::_stepi };
        enum { _stepj = Traits<M>::_stepj };
        enum { _diagstep = Traits<M>::_diagstep };
        enum { _conj = Traits<M>::_conj };
        enum { _canlin = Traits<M>::_canlin };

        typedef M type;
        typedef BaseMatrix_Calc<M> base;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;
        typedef typename base::inverse_type inverse_type;
        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename base::const_view_type const_view_type;
        typedef typename base::const_cview_type const_cview_type;
        typedef typename base::const_fview_type const_fview_type;
        typedef typename base::const_xview_type const_xview_type;
        typedef typename base::const_transpose_type const_transpose_type;
        typedef typename base::const_conjugate_type const_conjugate_type;
        typedef typename base::const_adjoint_type const_adjoint_type;
        typedef typename base::const_realpart_type const_realpart_type;
        typedef typename base::const_imagpart_type const_imagpart_type;
        typedef typename base::const_nonconj_type const_nonconj_type;
        typedef typename base::nonconst_type nonconst_type;

        typedef typename Traits<M>::const_row_type const_row_type;
        typedef typename Traits<M>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<M>::const_col_type const_col_type;
        typedef typename Traits<M>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<M>::const_diag_type const_diag_type;
        typedef typename Traits<M>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;

        typedef typename Traits<M>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<M>::const_submatrix_step_type 
            const_submatrix_step_type;

        typedef typename Traits<M>::const_subvector_type const_subvector_type;

        typedef typename Traits<M>::const_colpair_type const_colpair_type;
        typedef typename Traits<M>::const_rowpair_type const_rowpair_type;
        typedef typename Traits<M>::const_colrange_type const_colrange_type;
        typedef typename Traits<M>::const_rowrange_type const_rowrange_type;

        typedef typename Traits<M>::const_uppertri_type const_uppertri_type;
        typedef typename Traits<M>::const_unit_uppertri_type 
            const_unit_uppertri_type;
        typedef typename Traits<M>::const_unknown_uppertri_type 
            const_unknown_uppertri_type;
        typedef typename Traits<M>::const_lowertri_type const_lowertri_type;
        typedef typename Traits<M>::const_unit_lowertri_type 
            const_unit_lowertri_type;
        typedef typename Traits<M>::const_unknown_lowertri_type 
            const_unknown_lowertri_type;

        typedef typename Traits<M>::const_linearview_type 
            const_linearview_type;



        //
        // Constructor
        //

        TMV_INLINE BaseMatrix_Rec() {}
        TMV_INLINE BaseMatrix_Rec(const BaseMatrix_Rec<M>&) {}
        TMV_INLINE ~BaseMatrix_Rec() {}

    private:
        void operator=(const BaseMatrix_Rec<M>&);
    public:


        //
        // Access 
        //

        // The get_ routines always use CStyle indexing.
        const_row_type get_row(int i) const
        { return const_row_type(cptr()+i*stepi(),rowsize(),stepj()); }

        const_col_type get_col(int j) const
        { return const_col_type(cptr()+j*stepj(),colsize(),stepi()); }

        const_row_sub_type get_row(int i, int j1, int j2) const
        {
            return const_row_sub_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); 
        }

        const_col_sub_type get_col(int j, int i1, int i2) const
        {
            return const_col_sub_type(
                cptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); 
        }

        const_diag_sub_type get_diag(int i) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ? 
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i)),
                diagstep());
        }

        const_diag_sub_type get_diag(int i, int j1, int j2) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        TMV_INLINE const_row_type row(int i) const
        {
            CheckRowIndex<_fort>(i,colsize());
            return get_row(i);
        }

        TMV_INLINE const_row_type operator[](int i) const
        { return row(i); }

        TMV_INLINE const_col_type col(int j) const
        {
            CheckColIndex<_fort>(j,rowsize());
            return get_col(j);
        }

        TMV_INLINE_ND const_row_sub_type row(int i, int j1, int j2) const
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return get_row(i,j1,j2); 
        }

        TMV_INLINE_ND const_col_sub_type col(int j, int i1, int i2) const
        {
            CheckColIndex<_fort>(j,rowsize());
            CheckRowRange<_fort>(i1,i2,colsize());
            return get_col(j,i1,i2); 
        }

        // No need for a get_ routine for diag()
        const_diag_type diag() const
        {
            return const_diag_type(
                cptr(),TMV_MIN(colsize(),rowsize()),diagstep()); 
        }

        TMV_INLINE const_diag_sub_type diag(int i) const
        {
            CheckDiagIndex<_fort>(i,colsize(),rowsize());
            return get_diag(i);
        }

        TMV_INLINE const_diag_sub_type diag(int i, int j1, int j2) const
        {
            CheckDiagIndex<_fort>(i,j1,j2,colsize(),rowsize());
            return get_diag(i,j1,j2);
        }


        //
        // Functions
        //

        TMV_INLINE value_type sumElements() const
        { return tmv::DoSumElements(mat()); }

        TMV_INLINE float_type sumAbsElements() const
        { return tmv::DoSumAbsElements(mat()); }

        TMV_INLINE real_type sumAbs2Elements() const
        { return tmv::DoSumAbs2Elements(mat()); }

        TMV_INLINE float_type maxAbsElement() const
        { return tmv::DoMaxAbsElement(mat()); }

        TMV_INLINE real_type maxAbs2Element() const
        { return tmv::DoMaxAbs2Element(mat()); }

        TMV_INLINE real_type normSq() const
        { return tmv::DoNormSq(mat()); }

        TMV_INLINE float_type normSq(const float_type scale) const
        { return tmv::DoNormSq(mat(),scale); }

        TMV_INLINE float_type normF() const 
        { return tmv::DoNormF(mat()); }

        TMV_INLINE float_type norm() const
        { return normF(); }

        TMV_INLINE float_type norm1() const
        { return tmv::DoNorm1(mat()); }

        TMV_INLINE float_type normInf() const
        { return tmv::DoNormInf(mat()); }



        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        const_submatrix_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        {
            return const_submatrix_type(
                cptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
        }

        const_submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            return const_submatrix_step_type(
                cptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        const_subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        {
            return const_subvector_type(
                cptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        const_colpair_type cColPair(int j1, int j2) const
        {
            return const_colpair_type(
                cptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
        }

        const_rowpair_type cRowPair(int i1, int i2) const
        {
            return const_rowpair_type(
                cptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
        }

        const_colrange_type cColRange(int j1, int j2) const
        {
            return const_colrange_type(
                cptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
        }

        const_rowrange_type cRowRange(int i1, int i2) const
        {
            return const_rowrange_type(
                cptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
        }


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        const_submatrix_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return cSubMatrix(i1,i2,j1,j2);
        }

        const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            CheckRowRange<_fort>(i1,i2,istep,colsize());
            CheckColRange<_fort>(j1,j2,jstep,rowsize());
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        const_subvector_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            CheckMatSubVector<_fort>(i,j,istep,jstep,s,colsize(),rowsize());
            return cSubVector(i,j,istep,jstep,s); 
        }

        const_colpair_type colPair(int j1, int j2) const
        {
            CheckColIndex<_fort>(j1,rowsize());
            CheckColIndex<_fort>(j2,rowsize());
            return cColPair(j1,j2);
        }

        const_rowpair_type rowPair(int i1, int i2) const
        {
            CheckRowIndex<_fort>(i1,colsize());
            CheckRowIndex<_fort>(i2,colsize());
            return cRowPair(i1,i2);
        }

        const_colrange_type colRange(int j1, int j2) const
        {
            CheckColRange<_fort>(j1,j2,rowsize());
            return cColRange(j1,j2);
        }

        const_rowrange_type rowRange(int i1, int i2) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            return cRowRange(i1,i2);
        }


        //
        // Views
        //

        TMV_MAYBE_CREF(type,const_view_type) view() const
        { return MakeRecView<type,const_view_type>::call(mat()); }

        TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return MakeRecView<type,const_cview_type>::call(mat()); }

        TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return MakeRecView<type,const_fview_type>::call(mat()); }

        TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return MakeRecView<type,const_xview_type>::call(mat()); }

        TMV_MAYBE_CREF(type,const_cmview_type) cmView() const
        { return MakeRecView<type,const_cmview_type>::call(mat()); }

        TMV_MAYBE_CREF(type,const_rmview_type) rmView() const
        { return MakeRecView<type,const_rmview_type>::call(mat()); }

        TMV_MAYBE_CREF(type,const_view_type) constView() const
        { return MakeRecView<type,const_view_type>::call(mat()); }

        const_transpose_type transpose() const
        {
            return const_transpose_type(
                cptr(),rowsize(),colsize(),stepj(),stepi()); 
        }

        TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return MakeRecView<type,const_conjugate_type>::call(mat()); }

        const_adjoint_type adjoint() const
        {
            return const_adjoint_type(
                cptr(),rowsize(),colsize(),stepj(),stepi()); 
        }

        const_uppertri_type upperTri() const
        {
            return const_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),NonUnitDiag); 
        }

        const_unit_uppertri_type unitUpperTri() const
        {
            return const_unit_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),UnitDiag); 
        }

        const_unknown_uppertri_type upperTri(DiagType dt) const
        {
            TMVAssert(dt == NonUnitDiag || dt == UnitDiag);
            return const_unknown_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),dt); 
        }

        const_lowertri_type lowerTri() const
        {
            return const_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),NonUnitDiag); 
        }

        const_unit_lowertri_type unitLowerTri() const
        {
            return const_unit_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),UnitDiag); 
        }

        const_unknown_lowertri_type lowerTri(DiagType dt) const
        {
            TMVAssert(dt == NonUnitDiag || dt == UnitDiag);
            return const_unknown_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),dt); 
        }

        const_linearview_type linearView() const
        {
            TMVAssert(
                (stepi() == 1 && stepj() == int(colsize())) ||
                (stepj() == 1 && stepi() == int(rowsize())) );
            return const_linearview_type(cptr(),ls(),1);
        }

        const_realpart_type realPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()),
                colsize(), rowsize(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        const_imagpart_type imagPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1,
                colsize(), rowsize(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return MakeRecView<type,const_nonconj_type>::call(mat()); }

        nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<value_type*>(cptr()),
                colsize(),rowsize(),stepi(),stepj());
        }



        //
        // Auxilliary routines
        //

        template <class M2>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::Copy(mat(),m2.mat()); 
        }

        template <class M2>
        TMV_INLINE_ND void newAssignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::NoAliasCopy(mat(),m2.mat()); 
        }

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }

        TMV_INLINE int diagstep() const 
        { return _diagstep == UNKNOWN ? stepi() + stepj() : _diagstep; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool canLinearize() const 
        { return (AuxCanLinearize<_canlin,_colmajor,_rowmajor,M>::ok(mat())); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref from BaseMatrix.

        TMV_INLINE size_t colsize() const { return mat().colsize(); }
        TMV_INLINE size_t rowsize() const { return mat().rowsize(); }
        TMV_INLINE size_t ls() const { return mat().ls(); }
        TMV_INLINE int stepi() const { return mat().stepi(); }
        TMV_INLINE int stepj() const { return mat().stepj(); }
        TMV_INLINE bool isrm() const { return mat().isrm(); }
        TMV_INLINE bool iscm() const { return mat().iscm(); }

        TMV_INLINE const value_type* cptr() const { return mat().cptr(); }

    }; // BaseMatrix_Rec

    template <class M>
    class BaseMatrix_Rec_Mutable : 
        public BaseMatrix_Rec<M>,
        public BaseMatrix_Mutable<M>
    {
    public:
        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _stepi = Traits<M>::_stepi };
        enum { _stepj = Traits<M>::_stepj };
        enum { _diagstep = Traits<M>::_diagstep };
        enum { _conj = Traits<M>::_conj };
        enum { _canlin = Traits<M>::_canlin };

        typedef M type;
        typedef BaseMatrix_Rec<M> base;
        typedef BaseMatrix_Mutable<M> base_mut;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;
        typedef typename base::inverse_type inverse_type;
        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename base::const_view_type const_view_type;
        typedef typename base::const_cview_type const_cview_type;
        typedef typename base::const_fview_type const_fview_type;
        typedef typename base::const_xview_type const_xview_type;
        typedef typename base::const_transpose_type const_transpose_type;
        typedef typename base::const_conjugate_type const_conjugate_type;
        typedef typename base::const_adjoint_type const_adjoint_type;
        typedef typename base::const_realpart_type const_realpart_type;
        typedef typename base::const_imagpart_type const_imagpart_type;
        typedef typename base::const_nonconj_type const_nonconj_type;
        typedef typename base::nonconst_type nonconst_type;

        typedef typename base_mut::view_type view_type;
        typedef typename base_mut::cview_type cview_type;
        typedef typename base_mut::fview_type fview_type;
        typedef typename base_mut::xview_type xview_type;
        typedef typename base_mut::transpose_type transpose_type;
        typedef typename base_mut::conjugate_type conjugate_type;
        typedef typename base_mut::adjoint_type adjoint_type;
        typedef typename base_mut::realpart_type realpart_type;
        typedef typename base_mut::imagpart_type imagpart_type;
        typedef typename base_mut::nonconj_type nonconj_type;
        typedef typename base_mut::reference reference;

        typedef typename base::const_row_type const_row_type;
        typedef typename base::const_row_sub_type const_row_sub_type;
        typedef typename base::const_col_type const_col_type;
        typedef typename base::const_col_sub_type const_col_sub_type;
        typedef typename base::const_diag_type const_diag_type;
        typedef typename base::const_diag_sub_type const_diag_sub_type;
        typedef typename base::const_cmview_type const_cmview_type;
        typedef typename base::const_rmview_type const_rmview_type;
        typedef typename base::const_submatrix_type const_submatrix_type;
        typedef typename base::const_submatrix_step_type 
            const_submatrix_step_type;
        typedef typename base::const_subvector_type const_subvector_type;
        typedef typename base::const_colpair_type const_colpair_type;
        typedef typename base::const_rowpair_type const_rowpair_type;
        typedef typename base::const_colrange_type const_colrange_type;
        typedef typename base::const_rowrange_type const_rowrange_type;
        typedef typename base::const_uppertri_type const_uppertri_type;
        typedef typename base::const_unit_uppertri_type 
            const_unit_uppertri_type;
        typedef typename base::const_unknown_uppertri_type 
            const_unknown_uppertri_type;
        typedef typename base::const_lowertri_type const_lowertri_type;
        typedef typename base::const_unit_lowertri_type 
            const_unit_lowertri_type;
        typedef typename base::const_unknown_lowertri_type 
            const_unknown_lowertri_type;
        typedef typename base::const_linearview_type const_linearview_type;

        typedef typename Traits<M>::row_type row_type;
        typedef typename Traits<M>::row_sub_type row_sub_type;
        typedef typename Traits<M>::col_type col_type;
        typedef typename Traits<M>::col_sub_type col_sub_type;
        typedef typename Traits<M>::diag_type diag_type;
        typedef typename Traits<M>::diag_sub_type diag_sub_type;

        typedef typename Traits<M>::cmview_type cmview_type;
        typedef typename Traits<M>::rmview_type rmview_type;

        typedef typename Traits<M>::submatrix_type submatrix_type;
        typedef typename Traits<M>::submatrix_step_type submatrix_step_type;

        typedef typename Traits<M>::subvector_type subvector_type;

        typedef typename Traits<M>::colpair_type colpair_type;
        typedef typename Traits<M>::rowpair_type rowpair_type;
        typedef typename Traits<M>::colrange_type colrange_type;
        typedef typename Traits<M>::rowrange_type rowrange_type;

        typedef typename Traits<M>::uppertri_type uppertri_type;
        typedef typename Traits<M>::unit_uppertri_type unit_uppertri_type;
        typedef typename Traits<M>::unknown_uppertri_type unknown_uppertri_type;
        typedef typename Traits<M>::lowertri_type lowertri_type;
        typedef typename Traits<M>::unit_lowertri_type unit_lowertri_type;
        typedef typename Traits<M>::unknown_lowertri_type unknown_lowertri_type;

        typedef typename Traits<M>::linearview_type linearview_type;
        typedef typename Traits<M>::linear_iterator linear_iterator;

        //
        // Constructor
        //

        TMV_INLINE BaseMatrix_Rec_Mutable() {}
        TMV_INLINE BaseMatrix_Rec_Mutable(const BaseMatrix_Rec_Mutable<M>&) {}
        TMV_INLINE ~BaseMatrix_Rec_Mutable() {}


        //
        // Access 
        //

        TMV_INLINE reference operator()(int i, int j)
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColIndex<_fort>(j,rowsize());
            return ref(i,j);
        }

        // The get_ routines always use CStyle indexing.
        row_type get_row(int i) 
        { return row_type(ptr()+i*stepi(),rowsize(),stepj()); }

        col_type get_col(int j) 
        { return col_type(ptr()+j*stepj(),colsize(),stepi()); }

        // No need for a get_ routine for diag()
        diag_type diag() 
        { return diag_type(ptr(),TMV_MIN(colsize(),rowsize()),diagstep()); }

        row_sub_type get_row(int i, int j1, int j2) 
        { return row_sub_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

        col_sub_type get_col(int j, int i1, int i2) 
        { return col_sub_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

        diag_sub_type get_diag(int i) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ? 
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i) ),
                diagstep());
        }

        diag_sub_type get_diag(int i, int j1, int j2) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        TMV_INLINE row_type row(int i) 
        {
            CheckRowIndex<_fort>(i,colsize());
            return get_row(i);
        }

        TMV_INLINE row_type operator[](int i) 
        { return row(i); }

        TMV_INLINE col_type col(int j) 
        {
            CheckColIndex<_fort>(j,rowsize());
            return get_col(j);
        }

        TMV_INLINE row_sub_type row(int i, int j1, int j2) 
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return get_row(i,j1,j2); 
        }

        TMV_INLINE col_sub_type col(int j, int i1, int i2) 
        {
            CheckColIndex<_fort>(j,rowsize());
            CheckRowRange<_fort>(i1,i2,colsize());
            return get_col(j,i1,i2); 
        }

        TMV_INLINE diag_sub_type diag(int i) 
        {
            CheckDiagIndex<_fort>(i,colsize(),rowsize());
            return get_diag(i);
        }

        TMV_INLINE diag_sub_type diag(int i, int j1, int j2) 
        {
            CheckDiagIndex<_fort>(i,j1,j2,colsize(),rowsize());
            return get_diag(i,j1,j2);
        }


        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        TMV_INLINE value_type operator()(int i, int j) const
        { return base::operator()(i,j); }

        TMV_INLINE const_row_type get_row(int i) const
        { return base::get_row(i); }
        TMV_INLINE const_col_type get_col(int j) const
        { return base::get_col(j); }
        TMV_INLINE const_diag_type diag() const
        { return base::diag(); }
        TMV_INLINE const_row_sub_type get_row(int i, int j1, int j2) const
        { return base::get_row(i,j1,j2); }
        TMV_INLINE const_col_sub_type get_col(int j, int i1, int i2) const
        { return base::get_col(j,i1,i2); }
        TMV_INLINE const_diag_sub_type get_diag(int i) const
        { return base::get_diag(i); }
        TMV_INLINE const_diag_sub_type get_diag(int i, int j1, int j2) const
        { return base::get_diag(i,j1,j2); }

        TMV_INLINE const_row_type row(int i) const
        { return base::row(i); }
        TMV_INLINE const_row_type operator[](int i) const
        { return base::row(i); }
        TMV_INLINE const_col_type col(int j) const
        { return base::col(j); }
        TMV_INLINE const_row_sub_type row(int i, int j1, int j2) const
        { return base::row(i,j1,j2); }
        TMV_INLINE const_col_sub_type col(int j, int i1, int i2) const
        { return base::col(j,i1,i2); }
        TMV_INLINE const_diag_sub_type diag(int i) const
        { return base::diag(i); }
        TMV_INLINE const_diag_sub_type diag(int i, int j1, int j2) const
        { return base::diag(i,j1,j2); }


        //
        // Op =
        //

        TMV_INLINE_ND type& operator=(BaseMatrix_Rec_Mutable<M>& m2) 
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        TMV_INLINE_ND type& operator=(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_colsize,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_rowsize,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            m2.assignTo(mat());
            return mat(); 
        }

        TMV_INLINE_ND type& operator=(const value_type x)
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            setToIdentity(x);
            return mat();
        }

        ListAssigner<value_type,linear_iterator> operator<<(value_type x)
        {
            TMVAssert(this->canLinearize());
            linearview_type v = linearView();
            return ListAssigner<value_type,linear_iterator>(
                v.begin(),v.size(),x); 
        }


        //
        // Modifying Functions
        //

        // First the ones from BaseMatrix_Mutable:
        TMV_INLINE type& setZero()
        { tmv::SetZero(mat()); return mat(); }

        TMV_INLINE type& setAllTo(value_type val)
        { tmv::SetAllTo(mat(),val); return mat(); }

        TMV_INLINE type& addToAll(value_type val)
        { tmv::AddToAll(mat(),val); return mat(); }

        TMV_INLINE type& clip(float_type thresh)
        { tmv::Clip(mat(),thresh); return mat(); }

        template <class F>
        TMV_INLINE type& applyToAll(const F& f)
        { tmv::ApplyToAll(mat(),f); return mat(); }

        TMV_INLINE type& conjugateSelf()
        { tmv::ConjugateSelf(mat()); return mat(); }

        // Some more that are added for Rec shape:
        TMV_INLINE type& transposeSelf() 
        { tmv::TransposeSelf(mat()); return mat(); }

        type& setToIdentity(const value_type x=value_type(1))
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            this->setZero(); diag().setAllTo(x);
            return mat();
        }

        type& cSwapRows(int i1, int i2) 
        {
            if (i1 != i2) {
                row_type row1 = get_row(i1);
                row_type row2 = get_row(i2);
                tmv::NoAliasSwap(row1,row2);
            }
            return mat();
        }
        TMV_INLINE type& swapRows(int i1, int i2) 
        {
            CheckRowIndex<_fort>(i1,colsize());
            CheckRowIndex<_fort>(i2,colsize());
            return cSwapRows(i1,i2);
        }

        type& cSwapCols(int j1, int j2) 
        {
            if (j1 != j2) {
                col_type col1 = get_col(j1);
                col_type col2 = get_col(j2);
                tmv::NoAliasSwap(col1,col2);
            }
            return mat();
        }
        TMV_INLINE type& swapCols(int j1, int j2) 
        {
            CheckColIndex<_fort>(j1,rowsize());
            CheckColIndex<_fort>(j2,rowsize());
            return cSwapCols(j1,j2);
        }

        TMV_INLINE type& cPermuteRows(const int* p, int i1, int i2) 
        { tmv::PermuteRows(mat(),p,i1,i2); return mat(); }
        TMV_INLINE type& permuteRows(const int* p, int i1, int i2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            return cPermuteRows(p,i1,i2);
        }
        TMV_INLINE type& permuteRows(const int* p) 
        { return cPermuteRows(p,0,colsize()); }

        TMV_INLINE type& cReversePermuteRows(const int* p, int i1, int i2) 
        { tmv::ReversePermuteRows(mat(),p,i1,i2); return mat(); }
        TMV_INLINE type& reversePermuteRows(const int* p, int i1, int i2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            return cReversePermuteRows(p,i1,i2);
        }
        TMV_INLINE type& reversePermuteRows(const int* p) 
        { return cReversePermuteRows(p,0,colsize()); }

        TMV_INLINE type& permuteCols(const int* p, int j1, int j2) 
        {
            CheckColRange<_fort>(j1,j2,rowsize());
            transpose().cPermuteRows(p,j1,j2);
            return mat();
        }
        TMV_INLINE type& permuteCols(const int* p) 
        { transpose().cPermuteRows(p,0,rowsize()); return mat(); }

        TMV_INLINE type& reversePermuteCols(const int* p, int j1, int j2) 
        {
            CheckColRange<_fort>(j1,j2,rowsize());
            transpose().cReversePermuteRows(p,j1,j2);
            return mat();
        }
        TMV_INLINE type& reversePermuteCols(const int* p) 
        { transpose().cReversePermuteRows(p,0,rowsize()); return mat(); }


        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        submatrix_type cSubMatrix(int i1, int i2, int j1, int j2) 
        {
            return submatrix_type(
                ptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
        }

        submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            return submatrix_step_type(
                ptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) 
        {
            return subvector_type(
                ptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        colpair_type cColPair(int j1, int j2) 
        {
            return colpair_type(
                ptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
        }

        rowpair_type cRowPair(int i1, int i2) 
        {
            return rowpair_type(
                ptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
        }

        colrange_type cColRange(int j1, int j2) 
        {
            return colrange_type(
                ptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
        }

        rowrange_type cRowRange(int i1, int i2) 
        {
            return rowrange_type(
                ptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
        }


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return cSubMatrix(i1,i2,j1,j2);
        }

        submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            CheckRowRange<_fort>(i1,i2,istep,colsize());
            CheckColRange<_fort>(j1,j2,jstep,rowsize());
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        subvector_type subVector(
            int i, int j, int istep, int jstep, int s) 
        {
            CheckMatSubVector<_fort>(i,j,istep,jstep,s,colsize(),rowsize());
            return cSubVector(i,j,istep,jstep,s); 
        }

        colpair_type colPair(int j1, int j2) 
        {
            CheckColIndex<_fort>(j1,rowsize());
            CheckColIndex<_fort>(j2,rowsize());
            return cColPair(j1,j2);
        }

        rowpair_type rowPair(int i1, int i2) 
        {
            CheckRowIndex<_fort>(i1,colsize());
            CheckRowIndex<_fort>(i2,colsize());
            return cRowPair(i1,i2);
        }

        colrange_type colRange(int j1, int j2) 
        {
            CheckColRange<_fort>(j1,j2,rowsize());
            return cColRange(j1,j2);
        }

        rowrange_type rowRange(int i1, int i2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            return cRowRange(i1,i2);
        }


        // Repeat the const versions:
        TMV_INLINE const_submatrix_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        { return base::cSubMatrix(i1,i2,j1,j2); }
        TMV_INLINE const_submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return base::cSubMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_INLINE const_subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return base::cSubVector(i,j,istep,jstep,s); }
        TMV_INLINE const_colpair_type cColPair(int j1, int j2) const
        { return base::cColPair(j1,j2); }
        TMV_INLINE const_rowpair_type cRowPair(int i1, int i2) const
        { return base::cRowPair(i1,i2); }
        TMV_INLINE const_colrange_type cColRange(int j1, int j2) const
        { return base::cColRange(j1,j2); }
        TMV_INLINE const_rowrange_type cRowRange(int i1, int i2) const
        { return base::cRowRange(i1,i2); }

        TMV_INLINE const_submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) const
        { return base::subMatrix(i1,i2,j1,j2); }
        TMV_INLINE const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return base::subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_INLINE const_subvector_type subVector(
            int i, int j, int istep, int jstep, int s) const
        { return base::subVector(i,j,istep,jstep,s); }
        TMV_INLINE const_colpair_type colPair(int j1, int j2) const
        { return base::colPair(j1,j2); }
        TMV_INLINE const_rowpair_type rowPair(int i1, int i2) const
        { return base::rowPair(i1,i2); }
        TMV_INLINE const_colrange_type colRange(int j1, int j2) const
        { return base::colRange(j1,j2); }
        TMV_INLINE const_rowrange_type rowRange(int i1, int i2) const
        { return base::rowRange(i1,i2); }


        //
        // Views
        //

        TMV_MAYBE_REF(type,view_type) view() 
        { return MakeRecView<type,view_type>::call(mat()); }

        TMV_MAYBE_REF(type,cview_type) cView() 
        { return MakeRecView<type,cview_type>::call(mat()); }

        TMV_MAYBE_REF(type,fview_type) fView() 
        { return MakeRecView<type,fview_type>::call(mat()); }

        TMV_MAYBE_REF(type,xview_type) xView() 
        { return MakeRecView<type,xview_type>::call(mat()); }

        TMV_MAYBE_REF(type,cmview_type) cmView() 
        {
            TMVAssert(iscm() && "Called cmView on non-ColMajor matrix");
            return MakeRecView<type,cmview_type>::call(mat()); 
        }

        TMV_MAYBE_REF(type,rmview_type) rmView() 
        {
            TMVAssert(isrm() && "Called rmView on non-RowMajor matrix");
            return MakeRecView<type,rmview_type>::call(mat()); 
        }

        transpose_type transpose() 
        { return transpose_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

        TMV_MAYBE_REF(type,conjugate_type) conjugate() 
        { return MakeRecView<type,conjugate_type>::call(mat()); }

        adjoint_type adjoint() 
        { return adjoint_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

        uppertri_type upperTri() 
        { return uppertri_type(ptr(),rowsize(),stepi(),stepj(),NonUnitDiag); }

        unit_uppertri_type unitUpperTri() 
        {
            return unit_uppertri_type(
                ptr(),rowsize(),stepi(),stepj(),UnitDiag); 
        }

        unknown_uppertri_type upperTri(DiagType dt)
        {
            TMVAssert(dt == NonUnitDiag || dt == UnitDiag);
            return unknown_uppertri_type(
                ptr(),rowsize(),stepi(),stepj(),dt); 
        }

        lowertri_type lowerTri() 
        { return lowertri_type(ptr(),colsize(),stepi(),stepj(),NonUnitDiag); }

        unit_lowertri_type unitLowerTri() 
        {
            return unit_lowertri_type(
                ptr(),colsize(),stepi(),stepj(),UnitDiag); 
        }

        unknown_lowertri_type lowerTri(DiagType dt)
        {
            TMVAssert(dt == NonUnitDiag || dt == UnitDiag);
            return unknown_lowertri_type(
                ptr(),rowsize(),stepi(),stepj(),dt); 
        }

        linearview_type linearView() 
        {
            TMVAssert(this->canLinearize());
            return linearview_type(ptr(),ls(),1);
        }

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

        TMV_MAYBE_REF(type,nonconj_type) nonConj()
        { return MakeRecView<type,nonconj_type>::call(mat()); }


        // Repeat the const versions:
        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return base::view(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return base::cView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return base::fView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return base::xView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_cmview_type) cmView() const
        { return base::cmView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_rmview_type) rmView() const
        { return base::rmView(); }
        TMV_INLINE const_transpose_type transpose() const
        { return base::transpose(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return base::conjugate(); }
        TMV_INLINE const_adjoint_type adjoint() const
        { return base::adjoint(); }
        TMV_INLINE const_uppertri_type upperTri() const
        { return base::upperTri(); }
        TMV_INLINE const_unit_uppertri_type unitUpperTri() const
        { return base::unitUpperTri(); }
        TMV_INLINE const_unknown_uppertri_type upperTri(DiagType dt) const
        { return base::upperTri(dt); }
        TMV_INLINE const_lowertri_type lowerTri() const
        { return base::lowerTri(); }
        TMV_INLINE const_unit_lowertri_type unitLowerTri() const
        { return base::unitLowerTri(); }
        TMV_INLINE const_unknown_lowertri_type lowerTri(DiagType dt) const
        { return base::lowerTri(dt); }
        TMV_INLINE const_linearview_type linearView() const
        { return base::linearView(); }
        TMV_INLINE const_realpart_type realPart() const
        { return base::realPart(); }
        TMV_INLINE const_imagpart_type imagPart() const
        { return base::imagPart(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return base::nonConj(); }



        //
        // I/O
        //

        TMV_INLINE void read(std::istream& is)
        { tmv::Read(is,mat()); }

        //
        // Auxilliary routines
        //

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }
        TMV_INLINE type& mat()
        { return static_cast<type&>(*this); }

        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE int diagstep() const 
        { return _diagstep == UNKNOWN ? stepi() + stepj() : _diagstep; }
        TMV_INLINE bool canLinearize() const 
        { return (AuxCanLinearize<_canlin,_colmajor,_rowmajor,M>::ok(mat())); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref and cptr from above.

        TMV_INLINE size_t colsize() const { return mat().colsize(); }
        TMV_INLINE size_t rowsize() const { return mat().rowsize(); }
        TMV_INLINE size_t ls() const { return mat().ls(); }
        TMV_INLINE int stepi() const { return mat().stepi(); }
        TMV_INLINE int stepj() const { return mat().stepj(); }
        TMV_INLINE bool isrm() const { return mat().isrm(); }
        TMV_INLINE bool iscm() const { return mat().iscm(); }

        TMV_INLINE value_type* ptr() { return mat().ptr(); }
        TMV_INLINE reference ref(int i, int j) { return mat().ref(i,j); }

    }; // BaseMatrix_Rec_Mutable


    //
    // setZero
    //

    template <int algo, class M>
    struct SetZeroM_Helper;

    // algo 1: Linearize to vector version
    template <class M>
    struct SetZeroM_Helper<1,M>
    {
        static inline void call(M& m) 
        { m.linearView().setZero(); } 
    };

    // algo 2: RowMajor
    template <class M>
    struct SetZeroM_Helper<2,M>
    {
        static void call(M& m1) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).setZero();
        } 
    };

    // algo 3: ColMajor
    template <class M>
    struct SetZeroM_Helper<3,M>
    {
        static void call(M& m) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).setZero();
        } 
    };

    template <class M>
    static TMV_INLINE void SetZero(BaseMatrix_Rec_Mutable<M>& m)
    {
        const int algo = M::_canlin ? 1 : M::_colmajor ? 3 : 2;
        SetZeroM_Helper<algo,M>::call(m.mat());
    }

    //
    // setAllTo
    //

    template <int algo, class M, class T>
    struct SetAllToM_Helper;

    // algo 1: Linearize to vector version
    template <class M, class T>
    struct SetAllToM_Helper<1,M,T> // algo 1, linearize
    {
        static inline void call(M& m, const T& val) 
        { m.linearView().setAllTo(val); } 
    };

    // algo 2: RowMajor
    template <class M, class T>
    struct SetAllToM_Helper<2,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).setAllTo(val);
        } 
    };

    // algo 3: ColMajor
    template <class M, class T>
    struct SetAllToM_Helper<3,M,T>
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).setAllTo(val);
        } 
    };

    template <class M, class T>
    static TMV_INLINE void SetAllTo(BaseMatrix_Rec_Mutable<M>& m, const T& val)
    {
        const int algo = M::_canlin ? 1 : M::_colmajor ? 3 : 2;
        SetAllToM_Helper<algo,M,T>::call(m.mat(),val);
    }

    //
    // addToAll
    //

    template <int algo, class M, class T>
    struct AddToAllM_Helper;

    // algo 1: Linearize to vector version
    template <class M, class T>
    struct AddToAllM_Helper<1,M,T> // algo 1, linearize
    {
        static inline void call(M& m, const T& val) 
        { m.linearView().addToAll(val); } 
    };

    // algo 2: RowMajor
    template <class M, class T>
    struct AddToAllM_Helper<2,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).addToAll(val);
        } 
    };

    // algo 3: ColMajor
    template <class M, class T>
    struct AddToAllM_Helper<3,M,T>
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).addToAll(val);
        } 
    };

    template <class M, class T>
    static TMV_INLINE void AddToAll(BaseMatrix_Rec_Mutable<M>& m, const T& val)
    {
        const int algo = M::_canlin ? 1 : M::_colmajor ? 3 : 2;
        AddToAllM_Helper<algo,M,T>::call(m.mat(),val);
    }

    //
    // Clip
    //

    template <int algo, class M, class RT>
    struct ClipM_Helper;

    // algo 1: Linearize to vector version
    template <class M, class RT>
    struct ClipM_Helper<1,M,RT> // algo 1, linearize
    {
        static inline void call(M& m, const RT& thresh) 
        { m.linearView().clip(thresh); } 
    };

    // algo 2: RowMajor
    template <class M, class RT>
    struct ClipM_Helper<2,M,RT>
    {
        static void call(M& m1, const RT& thresh) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).clip(thresh);
        } 
    };

    // algo 3: ColMajor
    template <class M, class RT>
    struct ClipM_Helper<3,M,RT>
    {
        static void call(M& m, const RT& thresh) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).clip(thresh);
        } 
    };

    template <class M, class RT>
    static TMV_INLINE void Clip(BaseMatrix_Rec_Mutable<M>& m, const RT& thresh)
    {
        const int algo = M::_canlin ? 1 : M::_colmajor ? 3 : 2;
        ClipM_Helper<algo,M,RT>::call(m.mat(),thresh);
    }

    //
    // applyToAll
    //

    template <int algo, class M, class F>
    struct ApplyToAllM_Helper;

    // algo 1: Linearize to vector version
    template <class M, class F>
    struct ApplyToAllM_Helper<1,M,F> // algo 1, linearize
    {
        static inline void call(M& m, const F& f) 
        { m.linearView().applyToAll(f); } 
    };

    // algo 2: RowMajor
    template <class M, class F>
    struct ApplyToAllM_Helper<2,M,F>
    {
        static inline void call(M& m1, const F& f) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).applyToAll(f);
        } 
    };

    // algo 3: ColMajor
    template <class M, class F>
    struct ApplyToAllM_Helper<3,M,F>
    {
        static inline void call(M& m, const F& f) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).applyToAll(f);
        } 
    };

    template <class M, class F>
    static TMV_INLINE void ApplyToAll(BaseMatrix_Mutable<M>& m, const F& f)
    {
        const int algo = M::_canlin ? 1 : M::_colmajor ? 3 : 2;
        ApplyToAllM_Helper<algo,M,F>::call(m.mat(),f);
    }

    //
    // ConjugateSelf
    //

    template <int algo, class M>
    struct ConjugateM_Helper;

    // algo 0: Not complex, nothing to do
    template <class M>
    struct ConjugateM_Helper<0,M>
    { static inline void call(M& ) {} };

    // algo 1: Linearize to vector version
    template <class M>
    struct ConjugateM_Helper<1,M>
    {
        static inline void call(M& m)
        { m.linearView().conjugateSelf(); }
    };

    // In TMV_ScaleM.h
    template <int ix, class T, class M>
    static inline void Scale(const Scaling<ix,T>& x, BaseMatrix_Mutable<M>& m);

    // algo 2: m.imagPart() *= -1
    template <class M>
    struct ConjugateM_Helper<2,M>
    {
        static inline void call(M& m)
        {
            typedef typename M::real_type RT;
            typedef typename M::imagpart_type Mi;
            const Scaling<-1,RT> mone;
            Mi mi = m.imagPart();
            Scale(mone,mi);
        }
    };

    template <class M>
    static TMV_INLINE void ConjugateSelf(BaseMatrix_Rec_Mutable<M>& m)
    {
        const bool isreal = Traits<typename M::value_type>::isreal;
        const int algo = isreal ? 0 : M::_canlin ? 1 : 2;
        ConjugateM_Helper<algo,M>::call(m.mat());
    }



    // The BaseMatrix Trace call is efficient for composite types, since
    // it avoid calculating all the elements to do the sum.
    // But if we do have the elements calculated, this overloaded 
    // version will be faster:
    template <class M>
    static TMV_INLINE typename M::value_type DoTrace(const BaseMatrix_Rec<M>& m)
    { return m.diag().sumElements(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class M>
    static inline std::string TMV_Text(const BaseMatrix_Rec<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Rec< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    static inline std::string TMV_Text(const BaseMatrix_Rec_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Rec_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }
#endif

} // namespace tmv

#endif
