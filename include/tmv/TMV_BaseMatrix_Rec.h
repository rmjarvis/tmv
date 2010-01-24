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

#include <sstream>
#include "TMV_BaseMatrix.h"

namespace tmv {

    // BaseMatrix_Rec is derived from BaseMatrix_Calc, and is used
    // for dense rectangular matrices.
    template <class M>
    class BaseMatrix_Rec;

    // BaseMatrix_Rec adds the following requirements to Traits<M>:
    //
    //  mstepi = the step size along column if known (else UNKNOWN)
    //  mstepj = the step size along row if known (else UNKNOWN)
    //  mdiagstep = the step size along row if known (else UNKNOWN)
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
    //  const_lowertri_type = return type from lowerTri() const
    //  const_unit_lowertri_type = return type from unitLowerTri() const
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
    //  lowertri_type = return type from lowerTri() 
    //  unit_lowertri_type = return type from unitLowerTri() 
    //  realpart_type = return type from realPart()
    //  imagpart_type = return type from imagPart() 
    //  nonconj_type = return type from nonConj() 
    //  linearview_type = return type from linearView() 
    //  

    // The following all derive from BaseMatrix_Rec or BaseMatrix_Rec_Mutable.
    // See TMV_Matrix.h and TMV_SmallMatrix.h for their definitions:
    template <class T, StorageType S=ColMajor, IndexStyle I=CStyle>
    class Matrix;
    template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false,
              IndexStyle I=CStyle>
    class ConstMatrixView;
    template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false,
              IndexStyle I=CStyle>
    class MatrixView;
    template <class T, int M, int N, StorageType S=ColMajor,
              IndexStyle I=CStyle>
    class SmallMatrix;
    template <class T, int M, int N, int Si, int Sj, bool C=false,
              IndexStyle I=CStyle>
    class ConstSmallMatrixView;
    template <class T, int M, int N, int Si, int Sj, bool C=false,
              IndexStyle I=CStyle>
    class SmallMatrixView;

    // These are effectively aliases for I = FortranStyle
    // so you don' thave to write all the Si,Sj,C values if you are using
    // the defaults when you want to spectify FortranStyle.
    template <class T, StorageType S=ColMajor>
    class MatrixF;
    template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
    class ConstMatrixViewF;
    template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
    class MatrixViewF;
    template <class T, int M, int N, StorageType S=ColMajor>
    class SmallMatrixF;
    template <class T, int M, int N, int Si, int Sj, bool C=false>
    class ConstSmallMatrixViewF;
    template <class T, int M, int N, int Si, int Sj, bool C=false>
    class SmallMatrixViewF;


    //
    // Helper functions and values:
    //

    // Specify ExactSameStorage for rectangular matrices:
    template <class M1, class M2>
    inline bool ExactSameStorage(
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
    template <bool mfort>
    inline void CheckDiagIndex(int& i, int m, int n)
    { // CStyle or FortranStyle
        TMVAssert(i >= -m && "negative diag index must be <= nrows");
        TMVAssert(i <= n && "positive diag index must be <= ncols");
    }
    template <bool mfort>
    inline void CheckDiagIndex(int& i, int& j1, int& j2, int m, int n)
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
    template <bool mfort>
    inline void CheckRowRange(int& i1, int i2, int m)
    { // CStyle
        TMVAssert(i1 >= 0 && "first row must be in matrix");
        TMVAssert(i2 <= m && "last row must be in matrix");
        TMVAssert(i2 >= i1 && "range must have a non-negative number of rows");
    }
    template <bool mfort>
    inline void CheckColRange(int& j1, int j2, int n)
    { // CStyle
        TMVAssert(j1 >= 0 && "first column must be in matrix");
        TMVAssert(j2 <= n && "last column must be in matrix");
        TMVAssert(j2 >= j1 && 
                  "range must have a non-negative number of columns");
    }
    template <bool mfort>
    inline void CheckRowRange(int& i1, int& i2, int istep, int m)
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
    template <bool mfort>
    inline void CheckColRange(int& j1, int& j2, int jstep, int n)
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
    template <bool mfort>
    inline void CheckMatSubVector(int& i, int& j, int istep, int jstep, 
                                  int size, int m, int n)
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
    inline void CheckDiagIndex<true>(int& i, int& j1, int& j2, int m, int n)
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
    inline void CheckRowRange<true>(int& i1, int i2, int m)
    { // FortranStyle
        TMVAssert(i1 >= 1 && "first row must be in matrix");
        TMVAssert(i2 <= m && "last row must be in matrix");
        TMVAssert(i2 >= i1 && "range must have a positive number of rows");
        --i1;
    }
    template <>
    inline void CheckColRange<true>(int& j1, int j2, int n)
    { // FortranStyle
        TMVAssert(j1 >= 1 && "first column must be in matrix");
        TMVAssert(j2 <= n && "last column must be in matrix");
        TMVAssert(j2 >= j1 && "range must have a positive number of columns");
        --j1;
    }
    template <>
    inline void CheckRowRange<true>(int& i1, int& i2, int istep, int m)
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
    inline void CheckColRange<true>(int& j1, int& j2, int jstep, int n)
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
    inline void CheckMatSubVector<true>(int& i, int& j, int istep, int jstep, 
                                        int size, int m, int n)
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
        typedef SmallMatrix<T,cs,rs,rm?RowMajor:ColMajor,(
            fort?FortranStyle:CStyle)> type; 
    };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,Rec,UNKNOWN,rs,rm,fort>
    { typedef Matrix<T,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,Rec,cs,UNKNOWN,rm,fort>
    { typedef Matrix<T,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,Rec,UNKNOWN,UNKNOWN,rm,fort>
    { typedef Matrix<T,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };

    // A quick auxilliary function for canLinearize.
    // (It only accesses the steps that are unknown at compile time.)
    template <bool canlin, int Si, int Sj, class M>
    struct AuxCanLinearize;

    template <int Si, int Sj, class M>
    struct AuxCanLinearize<true,Si,Sj,M>
    {
        static inline bool ok(const M& m) 
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
    { static inline bool ok(const M& m) { return true; } };

    // Defined in TMV_CopyM.h
    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void NoAliasCopy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);

    // Defined in TMV_SwapM.h
    template <class M1, class M2>
    inline void Swap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M>
    inline void TransposeSelf(BaseMatrix_Rec_Mutable<M>& m);

    // Defined in TMV_Permute.h
    template <class M>
    inline void PermuteRows(BaseMatrix_Rec_Mutable<M>& m,
                            const int*const p, int i1, int i2);
    template <class M>
    inline void ReversePermuteRows(BaseMatrix_Rec_Mutable<M>& m,
                                   const int*const p, int i1, int i2);

    // Defined in TMV_NormM.h
    template <class M>
    inline typename M::real_type NormSq(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type NormSq(const BaseMatrix_Rec<M>& m, 
                                        const typename M::real_type scale);
    template <class M>
    inline typename M::real_type NormF(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type Norm1(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type NormInf(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::value_type SumElements(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type SumAbsElements(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type MaxAbsElement(const BaseMatrix_Rec<M>& m);

    // Defined in TMV_MatrixIO.h
    template <class M>
    inline void Read(std::istream& is, BaseMatrix_Rec_Mutable<M>& m);

    // Defined below:
    template <class M>
    inline void SetZero(BaseMatrix_Rec_Mutable<M>& m);
    template <class M, class RT>
    inline void Clip(BaseMatrix_Rec_Mutable<M>& m, const RT& thresh);
    template <class M, class T>
    inline void SetAllTo(BaseMatrix_Rec_Mutable<M>& m, const T& val);
    template <class M, class T>
    inline void AddToAll(BaseMatrix_Rec_Mutable<M>& m, const T& val);
    template <class M>
    inline void ConjugateSelf(BaseMatrix_Rec_Mutable<M>& m);
    template <class M, class F>
    inline void ApplyToAll(BaseMatrix_Rec_Mutable<M>& m, const F& f);


#if 0
    // Defined in TMV_MatrixDiv.h
    template <class M>
    inline typename M::value_type Det(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type LogDet(
        const BaseMatrix_Rec<M>& m, typename M::value_type* sign=0);
    template <class M, class M2>
    inline void DoInverse(
        const BaseMatrix_Rec<M>& m, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M, class M2>
    inline void DoInverseATA(
        const BaseMatrix_Rec<M>& m, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M>
    inline typename M::real_type Norm2(const BaseMatrix_Rec<M>& m);
    template <class M>
    inline typename M::real_type Condition(const BaseMatrix_Rec<M>& m);
#endif

    template <class M> 
    class BaseMatrix_Rec : public BaseMatrix_Calc<M>
    {
    public:
        enum { mcolsize = Traits<M>::mcolsize };
        enum { mrowsize = Traits<M>::mrowsize };
        enum { mshape = Traits<M>::mshape };
        enum { mfort = Traits<M>::mfort };
        enum { mrowmajor = Traits<M>::mrowmajor }; 
        enum { mcolmajor = Traits<M>::mcolmajor }; 
        enum { mstor = Traits<M>::mstor };
        enum { mcalc = Traits<M>::mcalc };
        enum { mstepi = Traits<M>::mstepi };
        enum { mstepj = Traits<M>::mstepj };
        enum { mdiagstep = Traits<M>::mdiagstep };
        enum { mconj = Traits<M>::mconj };
        enum { mcanlin = Traits<M>::mcanlin };

        typedef M type;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;

        typedef typename Traits<M>::const_row_type const_row_type;
        typedef typename Traits<M>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<M>::const_col_type const_col_type;
        typedef typename Traits<M>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<M>::const_diag_type const_diag_type;
        typedef typename Traits<M>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<M>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<M>::const_submatrix_step_type 
            const_submatrix_step_type;
        typedef typename Traits<M>::const_subvector_type const_subvector_type;
        typedef typename Traits<M>::const_colpair_type const_colpair_type;
        typedef typename Traits<M>::const_rowpair_type const_rowpair_type;
        typedef typename Traits<M>::const_colrange_type const_colrange_type;
        typedef typename Traits<M>::const_rowrange_type const_rowrange_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<M>::const_uppertri_type const_uppertri_type;
        typedef typename Traits<M>::const_unit_uppertri_type 
            const_unit_uppertri_type;
        typedef typename Traits<M>::const_lowertri_type const_lowertri_type;
        typedef typename Traits<M>::const_unit_lowertri_type 
            const_unit_lowertri_type;
        typedef typename Traits<M>::const_linearview_type 
            const_linearview_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;
        typedef typename Traits<M>::nonconst_type nonconst_type;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;



        //
        // Constructor
        //

        inline BaseMatrix_Rec() {}
        inline BaseMatrix_Rec(const BaseMatrix_Rec<M>&) {}
        inline ~BaseMatrix_Rec() {}

    private:
        void operator=(const BaseMatrix_Rec<M>&);
    public:


        //
        // Access 
        //

        // The get_ routines always use CStyle indexing.
        inline const_row_type get_row(int i) const
        { return const_row_type(cptr()+i*stepi(),rowsize(),stepj()); }

        inline const_col_type get_col(int j) const
        { return const_col_type(cptr()+j*stepj(),colsize(),stepi()); }

        inline const_row_sub_type get_row(int i, int j1, int j2) const
        {
            return const_row_sub_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); 
        }

        inline const_col_sub_type get_col(int j, int i1, int i2) const
        { 
            return const_col_sub_type(
                cptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); 
        }

        inline const_diag_sub_type get_diag(int i) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ? 
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i)),
                diagstep());
        }

        inline const_diag_sub_type get_diag(int i, int j1, int j2) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        inline const_row_type row(int i) const
        {
            CheckRowIndex<mfort>(i,colsize());
            return get_row(i);
        }

        inline const_row_type operator[](int i) const
        { return row(i); }

        inline const_col_type col(int j) const
        {
            CheckColIndex<mfort>(j,rowsize());
            return get_col(j);
        }

        inline const_row_sub_type row(int i, int j1, int j2) const
        {
            CheckRowIndex<mfort>(i,colsize());
            CheckColRange<mfort>(j1,j2,rowsize());
            return get_row(i,j1,j2); 
        }

        inline const_col_sub_type col(int j, int i1, int i2) const
        {
            CheckColIndex<mfort>(j,rowsize());
            CheckRowRange<mfort>(i1,i2,colsize());
            return get_col(j,i1,i2); 
        }

        // No need for a get_ routine for diag()
        inline const_diag_type diag() const
        { 
            return const_diag_type(
                cptr(),TMV_MIN(colsize(),rowsize()),diagstep()); 
        }

        inline const_diag_sub_type diag(int i) const
        {
            CheckDiagIndex<mfort>(i,colsize(),rowsize());
            return get_diag(i);
        }

        inline const_diag_sub_type diag(int i, int j1, int j2) const
        {
            CheckDiagIndex<mfort>(i,j1,j2,colsize(),rowsize());
            return get_diag(i,j1,j2);
        }


        //
        // Functions
        //

        inline value_type sumElements() const
        { return tmv::SumElements(cView()); }

        inline real_type sumAbsElements() const
        { return tmv::SumAbsElements(cView()); }

        inline real_type maxAbsElement() const
        { return tmv::MaxAbsElement(cView()); }

        inline real_type normSq() const
        { return tmv::NormSq(cView()); }

        inline real_type normSq(const real_type scale) const
        { return tmv::NormSq(cView(),scale); }

        inline real_type normF() const 
        { return tmv::NormF(cView()); }

        inline real_type norm() const
        { return normF(); }

        inline real_type norm1() const
        { return tmv::Norm1(cView()); }

        inline real_type normInf() const
        { return tmv::NormInf(cView()); }



        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        inline const_submatrix_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        {
            return const_submatrix_type(
                cptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
        }

        inline const_submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            return const_submatrix_step_type(
                cptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        inline const_subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        {
            return const_subvector_type(
                cptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        inline const_colpair_type cColPair(int j1, int j2) const
        {
            return const_colpair_type(
                cptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
        }

        inline const_rowpair_type cRowPair(int i1, int i2) const
        {
            return const_rowpair_type(
                cptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
        }

        inline const_colrange_type cColRange(int j1, int j2) const
        {
            return const_colrange_type(
                cptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
        }

        inline const_rowrange_type cRowRange(int i1, int i2) const
        {
            return const_rowrange_type(
                cptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
        }


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        inline const_submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            CheckRowRange<mfort>(i1,i2,colsize());
            CheckColRange<mfort>(j1,j2,rowsize());
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            CheckRowRange<mfort>(i1,i2,istep,colsize());
            CheckColRange<mfort>(j1,j2,jstep,rowsize());
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_subvector_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            CheckMatSubVector<mfort>(i,j,istep,jstep,s,colsize(),rowsize());
            return cSubVector(i,j,istep,jstep,s); 
        }

        inline const_colpair_type colPair(int j1, int j2) const
        {
            CheckColIndex<mfort>(j1,rowsize());
            CheckColIndex<mfort>(j2,rowsize());
            return cColPair(j1,j2);
        }

        inline const_rowpair_type rowPair(int i1, int i2) const
        {
            CheckRowIndex<mfort>(i1,colsize());
            CheckRowIndex<mfort>(i2,colsize());
            return cRowPair(i1,i2);
        }

        inline const_colrange_type colRange(int j1, int j2) const
        {
            CheckColRange<mfort>(j1,j2,rowsize());
            return cColRange(j1,j2);
        }

        inline const_rowrange_type rowRange(int i1, int i2) const
        {
            CheckRowRange<mfort>(i1,i2,colsize());
            return cRowRange(i1,i2);
        }


        //
        // Views
        //

        inline const_view_type view() const
        { return const_view_type(cptr(),colsize(),rowsize(),stepi(),stepj()); }

        inline const_cview_type cView() const
        { return view(); }

        inline const_fview_type fView() const
        { return view(); }

        inline const_xview_type xView() const
        { return view(); }

        inline const_cmview_type cmView() const
        { 
            TMVAssert(iscm() && "Called cmView on non-ColMajor matrix");
            return view(); 
        }

        inline const_rmview_type rmView() const
        {
            TMVAssert(isrm() && "Called rmView on non-RowMajor matrix");
            return view(); 
        }

        inline const_view_type constView() const
        { return view(); }

        inline const_transpose_type transpose() const
        { 
            return const_transpose_type(
                cptr(),rowsize(),colsize(),stepj(),stepi()); 
        }

        inline const_conjugate_type conjugate() const
        {
            return const_conjugate_type(
                cptr(),colsize(),rowsize(),stepi(),stepj()); 
        }

        inline const_adjoint_type adjoint() const
        { 
            return const_adjoint_type(
                cptr(),rowsize(),colsize(),stepj(),stepi()); 
        }

        inline const_uppertri_type upperTri() const
        { return const_uppertri_type(cptr(),rowsize(),false,stepi(),stepj()); }

        inline const_unit_uppertri_type unitUpperTri() const
        {
            return const_unit_uppertri_type(
                cptr(),rowsize(),true,stepi(),stepj()); 
        }

        inline const_lowertri_type lowerTri() const
        { return const_lowertri_type(cptr(),colsize(),false,stepi(),stepj()); }

        inline const_unit_lowertri_type unitLowerTri() const
        {
            return const_unit_lowertri_type(
                cptr(),colsize(),true,stepi(),stepj()); 
        }

        inline const_linearview_type linearView() const
        {
            TMVAssert(
                (stepi() == 1 && stepj() == colsize()) ||
                (stepj() == 1 && stepi() == rowsize()) );
            return const_linearview_type(cptr(),ls(),1);
        }

        inline const_realpart_type realPart() const
        {
            const bool misreal = Traits<value_type>::isreal;
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()),
                colsize(), rowsize(),
                misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
        }

        inline const_imagpart_type imagPart() const
        {
            const bool misreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1,
                colsize(), rowsize(),
                misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
        }

        inline const_nonconj_type nonConj() const
        { 
            return const_nonconj_type(
                cptr(),colsize(),rowsize(),stepi(),stepj()); 
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<value_type*>(cptr()),
                colsize(),rowsize(),stepi(),stepj());
        }



        //
        // Auxilliary routines
        //

        template <class M2>
        inline void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<mshape,M2::mshape>::assignable));
            tmv::Copy(mat(),m2.mat()); 
        }

        template <class M2>
        inline void newAssignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<mshape,M2::mshape>::assignable));
            tmv::NoAliasCopy(mat(),m2.mat()); 
        }

        inline const type& mat() const
        { return *static_cast<const type*>(this); }

        inline int diagstep() const 
        { return mdiagstep == UNKNOWN ? stepi() + stepj() : mdiagstep; }
        inline bool isconj() const { return mconj; }
        inline bool canLinearize() const 
        { return (AuxCanLinearize<mcanlin,mcolmajor,mrowmajor,M>::ok(mat())); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref from BaseMatrix.

        inline size_t colsize() const { return mat().colsize(); }
        inline size_t rowsize() const { return mat().rowsize(); }
        inline size_t ls() const { return mat().ls(); }
        inline int stepi() const { return mat().stepi(); }
        inline int stepj() const { return mat().stepj(); }
        inline bool isrm() const { return mat().isrm(); }
        inline bool iscm() const { return mat().iscm(); }

        inline const value_type* cptr() const { return mat().cptr(); }

    }; // BaseMatrix_Rec

    template <class M> 
    class BaseMatrix_Rec_Mutable : public BaseMatrix_Rec<M>,
                                   public BaseMatrix_Mutable<M>
    {
    public:
        enum { mcolsize = Traits<M>::mcolsize };
        enum { mrowsize = Traits<M>::mrowsize };
        enum { mshape = Traits<M>::mshape };
        enum { mfort = Traits<M>::mfort };
        enum { mcalc = Traits<M>::mcalc };
        enum { mrowmajor = Traits<M>::mrowmajor }; 
        enum { mcolmajor = Traits<M>::mcolmajor }; 
        enum { mstor = Traits<M>::mstor };
        enum { mstepi = Traits<M>::mstepi };
        enum { mstepj = Traits<M>::mstepj };
        enum { mdiagstep = Traits<M>::mdiagstep };
        enum { mconj = Traits<M>::mconj };
        enum { mcanlin = Traits<M>::mcanlin };

        typedef M type;
        typedef BaseMatrix_Rec<M> base_rec;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;

        typedef typename Traits<M>::const_row_type const_row_type;
        typedef typename Traits<M>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<M>::const_col_type const_col_type;
        typedef typename Traits<M>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<M>::const_diag_type const_diag_type;
        typedef typename Traits<M>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<M>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<M>::const_submatrix_step_type 
            const_submatrix_step_type;
        typedef typename Traits<M>::const_subvector_type const_subvector_type;
        typedef typename Traits<M>::const_colpair_type const_colpair_type;
        typedef typename Traits<M>::const_rowpair_type const_rowpair_type;
        typedef typename Traits<M>::const_colrange_type const_colrange_type;
        typedef typename Traits<M>::const_rowrange_type const_rowrange_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<M>::const_uppertri_type const_uppertri_type;
        typedef typename Traits<M>::const_unit_uppertri_type 
            const_unit_uppertri_type;
        typedef typename Traits<M>::const_lowertri_type const_lowertri_type;
        typedef typename Traits<M>::const_unit_lowertri_type 
            const_unit_lowertri_type;
        typedef typename Traits<M>::const_linearview_type 
            const_linearview_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;

        typedef typename Traits<M>::row_type row_type;
        typedef typename Traits<M>::row_sub_type row_sub_type;
        typedef typename Traits<M>::col_type col_type;
        typedef typename Traits<M>::col_sub_type col_sub_type;
        typedef typename Traits<M>::diag_type diag_type;
        typedef typename Traits<M>::diag_sub_type diag_sub_type;

        typedef typename Traits<M>::submatrix_type submatrix_type;
        typedef typename Traits<M>::submatrix_step_type submatrix_step_type;
        typedef typename Traits<M>::subvector_type subvector_type;
        typedef typename Traits<M>::colpair_type colpair_type;
        typedef typename Traits<M>::rowpair_type rowpair_type;
        typedef typename Traits<M>::colrange_type colrange_type;
        typedef typename Traits<M>::rowrange_type rowrange_type;

        typedef typename Traits<M>::view_type view_type;
        typedef typename Traits<M>::cview_type cview_type;
        typedef typename Traits<M>::fview_type fview_type;
        typedef typename Traits<M>::xview_type xview_type;
        typedef typename Traits<M>::cmview_type cmview_type;
        typedef typename Traits<M>::rmview_type rmview_type;
        typedef typename Traits<M>::transpose_type transpose_type;
        typedef typename Traits<M>::conjugate_type conjugate_type;
        typedef typename Traits<M>::adjoint_type adjoint_type;
        typedef typename Traits<M>::uppertri_type uppertri_type;
        typedef typename Traits<M>::unit_uppertri_type unit_uppertri_type;
        typedef typename Traits<M>::lowertri_type lowertri_type;
        typedef typename Traits<M>::unit_lowertri_type unit_lowertri_type;
        typedef typename Traits<M>::linearview_type linearview_type;
        typedef typename Traits<M>::realpart_type realpart_type;
        typedef typename Traits<M>::imagpart_type imagpart_type;
        typedef typename Traits<M>::nonconj_type nonconj_type;

        typedef typename Traits<M>::reference reference;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;


        //
        // Constructor
        //

        inline BaseMatrix_Rec_Mutable() {}
        inline BaseMatrix_Rec_Mutable(const BaseMatrix_Rec_Mutable<M>&) {}
        inline ~BaseMatrix_Rec_Mutable() {}


        //
        // Access 
        //

        inline reference operator()(int i, int j)
        {
            CheckRowIndex<mfort>(i,colsize());
            CheckColIndex<mfort>(j,rowsize());
            return ref(i,j);
        }

        // The get_ routines always use CStyle indexing.
        inline row_type get_row(int i) 
        { return row_type(ptr()+i*stepi(),rowsize(),stepj()); }

        inline col_type get_col(int j) 
        { return col_type(ptr()+j*stepj(),colsize(),stepi()); }

        // No need for a get_ routine for diag()
        inline diag_type diag() 
        {
            return diag_type(ptr(),TMV_MIN(colsize(),rowsize()),diagstep());
        }

        inline row_sub_type get_row(int i, int j1, int j2) 
        { return row_sub_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

        inline col_sub_type get_col(int j, int i1, int i2) 
        { return col_sub_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

        inline diag_sub_type get_diag(int i) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ? 
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i) ),
                diagstep());
        }

        inline diag_sub_type get_diag(int i, int j1, int j2) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        inline row_type row(int i) 
        {
            CheckRowIndex<mfort>(i,colsize());
            return get_row(i);
        }

        inline row_type operator[](int i) 
        { return row(i); }

        inline col_type col(int j) 
        {
            CheckColIndex<mfort>(j,rowsize());
            return get_col(j);
        }

        inline row_sub_type row(int i, int j1, int j2) 
        {
            CheckRowIndex<mfort>(i,colsize());
            CheckColRange<mfort>(j1,j2,rowsize());
            return get_row(i,j1,j2); 
        }

        inline col_sub_type col(int j, int i1, int i2) 
        {
            CheckColIndex<mfort>(j,rowsize());
            CheckRowRange<mfort>(i1,i2,colsize());
            return get_col(j,i1,i2); 
        }

        inline diag_sub_type diag(int i) 
        {
            CheckDiagIndex<mfort>(i,colsize(),rowsize());
            return get_diag(i);
        }

        inline diag_sub_type diag(int i, int j1, int j2) 
        {
            CheckDiagIndex<mfort>(i,j1,j2,colsize(),rowsize());
            return get_diag(i,j1,j2);
        }


        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        inline value_type operator()(int i, int j) const
        { return base_rec::operator()(i,j); }

        inline const_row_type get_row(int i) const
        { return base_rec::get_row(i); }
        inline const_col_type get_col(int j) const
        { return base_rec::get_col(j); }
        inline const_diag_type diag() const
        { return base_rec::diag(); }
        inline const_row_sub_type get_row(int i, int j1, int j2) const
        { return base_rec::get_row(i,j1,j2); }
        inline const_col_sub_type get_col(int j, int i1, int i2) const
        { return base_rec::get_col(j,i1,i2); }
        inline const_diag_sub_type get_diag(int i) const
        { return base_rec::get_diag(i); }
        inline const_diag_sub_type get_diag(int i, int j1, int j2) const
        { return base_rec::get_diag(i,j1,j2); }

        inline const_row_type row(int i) const
        { return base_rec::row(i); }
        inline const_row_type operator[](int i) const
        { return base_rec::row(i); }
        inline const_col_type col(int j) const
        { return base_rec::col(j); }
        inline const_row_sub_type row(int i, int j1, int j2) const
        { return base_rec::row(i,j1,j2); }
        inline const_col_sub_type col(int j, int i1, int i2) const
        { return base_rec::col(j,i1,i2); }
        inline const_diag_sub_type diag(int i) const
        { return base_rec::diag(i); }
        inline const_diag_sub_type diag(int i, int j1, int j2) const
        { return base_rec::diag(i,j1,j2); }


        //
        // Op =
        //

        inline type& operator=(BaseMatrix_Rec_Mutable<M>& m2) 
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert((Sizes<mcolsize,M2::mcolsize>::same));
            TMVStaticAssert((Sizes<mrowsize,M2::mrowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
            m2.assignTo(mat());
            return mat(); 
        }

        inline type& operator=(const value_type x)
        {
            TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
            TMVAssert(colsize() == rowsize());
            setToIdentity(x);
            return mat();
        }

        typedef typename linearview_type::iterator lin_iter;
        inline ListAssigner<value_type,lin_iter> operator<<(value_type x)
        {
            TMVAssert(this->canLinearize());
            linearview_type v = linearView();
            return ListAssigner<value_type,lin_iter>(v.begin(),v.size(),x); 
        }


        //
        // Modifying Functions
        //

        // First the ones from BaseMatrix_Mutable:
        inline type& setZero()
        { tmv::SetZero(mat()); return mat(); }

        inline type& setAllTo(value_type val)
        { tmv::SetAllTo(mat(),val); return mat(); }

        inline type& addToAll(value_type val)
        { tmv::AddToAll(mat(),val); return mat(); }

        inline type& clip(real_type thresh)
        { tmv::Clip(mat(),thresh); return mat(); }

        template <class F>
        inline type& applyToAll(const F& f)
        { tmv::ApplyToAll(mat(),f); return mat(); }

        inline type& conjugateSelf()
        { tmv::ConjugateSelf(mat()); return mat(); }

        // Some more that are added for Rec shape:
        inline type& transposeSelf() 
        { tmv::TransposeSelf(mat()); return mat(); }

        inline type& setToIdentity(const value_type x=value_type(1))
        {
            TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
            TMVAssert(colsize() == rowsize());
            this->setZero(); diag().setAllTo(x);
            return mat();
        }

        inline type& cSwapRows(int i1, int i2) 
        {
            if (i1 != i2) {
                tmv::Swap(get_row(i1),get_row(i2));
            }
            return mat();
        }
        inline type& swapRows(int i1, int i2) 
        {
            CheckRowIndex<mfort>(i1,colsize());
            CheckRowIndex<mfort>(i2,colsize());
            return cSwapRows(i1,i2);
        }

        inline type& cSwapCols(int j1, int j2) 
        {
            if (j1 != j2) {
                tmv::Swap(get_col(j1),get_col(j2));
            }
            return mat();
        }
        inline type& swapCols(int j1, int j2) 
        {
            CheckColIndex<mfort>(j1,rowsize());
            CheckColIndex<mfort>(j2,rowsize());
            return cSwapCols(j1,j2);
        }

        inline type& cPermuteRows(const int*const p, int i1, int i2) 
        { tmv::PermuteRows(mat(),p,i1,i2); return mat(); }
        inline type& permuteRows(const int*const p, int i1, int i2) 
        {
            CheckRowRange<mfort>(i1,i2,colsize());
            return cPermuteRows(p,i1,i2);
        }
        inline type& permuteRows(const int*const p) 
        { cPermuteRows(p,0,colsize()); return mat(); }

        inline type& cReversePermuteRows(const int*const p, int i1, int i2) 
        { tmv::ReversePermuteRows(mat(),p,i1,i2); }
        inline type& reversePermuteRows(const int*const p, int i1, int i2) 
        {
            CheckRowRange<mfort>(i1,i2,colsize());
            return cReversePermuteRows(p,i1,i2);
        }
        inline type& reversePermuteRows(const int*const p) 
        { cReversePermuteRows(p,0,colsize()); return mat(); }

        inline type& permuteCols(const int*const p, int j1, int j2) 
        {
            CheckColRange<mfort>(j1,j2,rowsize());
            transpose().cPermuteRows(p,j1,j2);
            return mat();
        }
        inline type& permuteCols(const int*const p) 
        { transpose().cPermuteRows(p,0,rowsize()); return mat(); }

        inline type& reversePermuteCols(const int*const p, int j1, int j2) 
        {
            CheckColRange<mfort>(j1,j2,rowsize());
            transpose().cReversePermuteRows(p,j1,j2);
            return mat();
        }
        inline type& reversePermuteCols(const int*const p) 
        { transpose().cReversePermuteRows(p,0,rowsize()); return mat(); }


        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        inline submatrix_type cSubMatrix(int i1, int i2, int j1, int j2) 
        {
            return submatrix_type(
                ptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
        }

        inline submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            return submatrix_step_type(
                ptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        inline subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) 
        {
            return subvector_type(
                ptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        inline colpair_type cColPair(int j1, int j2) 
        {
            return colpair_type(
                ptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
        }

        inline rowpair_type cRowPair(int i1, int i2) 
        {
            return rowpair_type(
                ptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
        }

        inline colrange_type cColRange(int j1, int j2) 
        {
            return colrange_type(
                ptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
        }

        inline rowrange_type cRowRange(int i1, int i2) 
        {
            return rowrange_type(
                ptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
        }


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        inline submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) 
        {
            CheckRowRange<mfort>(i1,i2,colsize());
            CheckColRange<mfort>(j1,j2,rowsize());
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            CheckRowRange<mfort>(i1,i2,istep,colsize());
            CheckColRange<mfort>(j1,j2,jstep,rowsize());
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline subvector_type subVector(
            int i, int j, int istep, int jstep, int s) 
        {
            CheckMatSubVector<mfort>(i,j,istep,jstep,s,colsize(),rowsize());
            return cSubVector(i,j,istep,jstep,s); 
        }

        inline colpair_type colPair(int j1, int j2) 
        {
            CheckColIndex<mfort>(j1,rowsize());
            CheckColIndex<mfort>(j2,rowsize());
            return cColPair(j1,j2);
        }

        inline rowpair_type rowPair(int i1, int i2) 
        {
            CheckRowIndex<mfort>(i1,colsize());
            CheckRowIndex<mfort>(i2,colsize());
            return cRowPair(i1,i2);
        }

        inline colrange_type colRange(int j1, int j2) 
        {
            CheckColRange<mfort>(j1,j2,rowsize());
            return cColRange(j1,j2);
        }

        inline rowrange_type rowRange(int i1, int i2) 
        {
            CheckRowRange<mfort>(i1,i2,colsize());
            return cRowRange(i1,i2);
        }


        // Repeat the const versions:
        inline const_submatrix_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        { return base_rec::cSubMatrix(i1,i2,j1,j2); }
        inline const_submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return base_rec::cSubMatrix(i1,i2,j1,j2,istep,jstep); }
        inline const_subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return base_rec::cSubVector(i,j,istep,jstep,s); }
        inline const_colpair_type cColPair(int j1, int j2) const
        { return base_rec::cColPair(j1,j2); }
        inline const_rowpair_type cRowPair(int i1, int i2) const
        { return base_rec::cRowPair(i1,i2); }
        inline const_colrange_type cColRange(int j1, int j2) const
        { return base_rec::cColRange(j1,j2); }
        inline const_rowrange_type cRowRange(int i1, int i2) const
        { return base_rec::cRowRange(i1,i2); }

        inline const_submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) const
        { return base_rec::subMatrix(i1,i2,j1,j2); }
        inline const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return base_rec::subMatrix(i1,i2,j1,j2,istep,jstep); }
        inline const_subvector_type subVector(
            int i, int j, int istep, int jstep, int s) const
        { return base_rec::subVector(i,j,istep,jstep,s); }
        inline const_colpair_type colPair(int j1, int j2) const
        { return base_rec::colPair(j1,j2); }
        inline const_rowpair_type rowPair(int i1, int i2) const
        { return base_rec::rowPair(i1,i2); }
        inline const_colrange_type colRange(int j1, int j2) const
        { return base_rec::colRange(j1,j2); }
        inline const_rowrange_type rowRange(int i1, int i2) const
        { return base_rec::rowRange(i1,i2); }


        //
        // Views
        //

        inline view_type view() 
        { return view_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }

        inline cview_type cView() 
        { return view(); }

        inline fview_type fView() 
        { return view(); }

        inline xview_type xView() 
        { return view(); }

        inline cmview_type cmView() 
        {
            TMVAssert(iscm() && "Called cmView on non-ColMajor matrix");
            return view(); 
        }

        inline rmview_type rmView() 
        {
            TMVAssert(isrm() && "Called rmView on non-RowMajor matrix");
            return view(); 
        }

        inline transpose_type transpose() 
        { return transpose_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

        inline conjugate_type conjugate() 
        { return conjugate_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }

        inline adjoint_type adjoint() 
        { return adjoint_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

        inline uppertri_type upperTri() 
        { return uppertri_type(ptr(),rowsize(),false,stepi(),stepj()); }

        inline unit_uppertri_type unitUpperTri() 
        { return unit_uppertri_type(ptr(),rowsize(),true,stepi(),stepj()); }

        inline lowertri_type lowerTri() 
        { return lowertri_type(ptr(),colsize(),false,stepi(),stepj()); }

        inline unit_lowertri_type unitLowerTri() 
        { return unit_lowertri_type(ptr(),colsize(),true,stepi(),stepj()); }

        inline linearview_type linearView() 
        {
            TMVAssert(this->canLinearize());
            return linearview_type(ptr(),ls(),1);
        }

        inline realpart_type realPart() 
        {
            const bool misreal = Traits<value_type>::isreal;
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()), colsize(), rowsize(),
                misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
        }

        inline imagpart_type imagPart() 
        {
            const bool misreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1, colsize(), rowsize(),
                misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
        }

        inline nonconj_type nonConj()
        { return nonconj_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }


        // Repeat the const versions:
        inline const_view_type view() const
        { return base_rec::view(); }
        inline const_cview_type cView() const
        { return base_rec::view(); }
        inline const_fview_type fView() const
        { return base_rec::view(); }
        inline const_xview_type xView() const
        { return base_rec::view(); }
        inline const_cmview_type cmView() const
        { return base_rec::view(); }
        inline const_rmview_type rmView() const
        { return base_rec::view(); }
        inline const_transpose_type transpose() const
        { return base_rec::transpose(); }
        inline const_conjugate_type conjugate() const
        { return base_rec::conjugate(); }
        inline const_adjoint_type adjoint() const
        { return base_rec::adjoint(); }
        inline const_uppertri_type upperTri() const
        { return base_rec::upperTri(); }
        inline const_unit_uppertri_type unitUpperTri() const
        { return base_rec::unitUpperTri(); }
        inline const_lowertri_type lowerTri() const
        { return base_rec::lowerTri(); }
        inline const_unit_lowertri_type unitLowerTri() const
        { return base_rec::unitLowerTri(); }
        inline const_linearview_type linearView() const
        { return base_rec::linearView(); }
        inline const_realpart_type realPart() const
        { return base_rec::realPart(); }
        inline const_imagpart_type imagPart() const
        { return base_rec::imagPart(); }
        inline const_nonconj_type nonConj() const
        { return base_rec::nonConj(); }



        //
        // I/O
        //

        inline void read(std::istream& is)
        {
            cview_type mcv = cView();
            tmv::Read(is,mcv); 
        }

        //
        // Auxilliary routines
        //

        inline const type& mat() const
        { return *static_cast<const type*>(this); }
        inline type& mat()
        { return *static_cast<type*>(this); }

        inline bool isconj() const { return mconj; }
        inline int diagstep() const 
        { return mdiagstep == UNKNOWN ? stepi() + stepj() : mdiagstep; }
        inline bool canLinearize() const 
        { return (AuxCanLinearize<mcanlin,mcolmajor,mrowmajor,M>::ok(mat())); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref and cptr from above.

        inline size_t colsize() const { return mat().colsize(); }
        inline size_t rowsize() const { return mat().rowsize(); }
        inline size_t ls() const { return mat().ls(); }
        inline int stepi() const { return mat().stepi(); }
        inline int stepj() const { return mat().stepj(); }
        inline bool isrm() const { return mat().isrm(); }
        inline bool iscm() const { return mat().iscm(); }

        inline value_type* ptr() { return mat().ptr(); }
        inline reference ref(int i, int j) { return mat().ref(i,j); }

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
        static inline void call(M& m1) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).setZero();
        } 
    };

    // algo 3: ColMajor
    template <class M> 
    struct SetZeroM_Helper<3,M> 
    {
        static inline void call(M& m) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).setZero();
        } 
    };

    // algo 4: Unknown sizes, determine which algorithm to use
    template <class M> 
    struct SetZeroM_Helper<4,M> 
    {
        static inline void call(M& m) 
        {
            const int algo2 = M::mcolmajor ? 3 : 2;
#if TMV_OPT >= 2
            if (m.canLinearize())
                SetZeroM_Helper<1,M>::call(m);
            else
#endif
                SetZeroM_Helper<algo2,M>::call(m);
        } 
    };

    template <class M>
    inline void SetZero(BaseMatrix_Rec_Mutable<M>& m)
    {
#if TMV_OPT == 0
        const int algo = M::mcolmajor ? 3 : 2;
#else
        const int algo = M::mcanlin ? 1 : 4;
#endif
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

    // algo 4: Unknown sizes, determine which algorithm to use
    template <class M, class T> 
    struct SetAllToM_Helper<4,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int algo2 = M::mcolmajor ? 3 : 2;
#if TMV_OPT >= 2
            if (m.canLinearize())
                SetAllToM_Helper<1,M,T>::call(m,val);
            else
#endif
                SetAllToM_Helper<algo2,M,T>::call(m,val);
        } 
    };

    template <class M, class T>
    inline void SetAllTo(BaseMatrix_Rec_Mutable<M>& m, const T& val)
    {
#if TMV_OPT == 0
        const int algo = M::mcolmajor ? 3 : 2;
#else
        const int algo = M::mcanlin ? 1 : 4;
#endif
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

    // algo 4: Unknown sizes, determine which algorithm to use
    template <class M, class T> 
    struct AddToAllM_Helper<4,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int algo2 = M::mcolmajor ? 3 : 2;
#if TMV_OPT >= 2
            if (m.canLinearize())
                AddToAllM_Helper<1,M,T>::call(m,val);
            else
#endif
                AddToAllM_Helper<algo2,M,T>::call(m,val);
        } 
    };

    template <class M, class T>
    inline void AddToAll(BaseMatrix_Rec_Mutable<M>& m, const T& val)
    {
#if TMV_OPT == 0
        const int algo = M::mcolmajor ? 3 : 2;
#else
        const int algo = M::mcanlin ? 1 : 4;
#endif
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

    // algo 4: Unknown sizes, determine which algorithm to use
    template <class M, class RT> 
    struct ClipM_Helper<4,M,RT> 
    {
        static void call(M& m, const RT& thresh) 
        {
            const int algo2 = M::mcolmajor ? 3 : 2;
#if TMV_OPT >= 2
            if (m.canLinearize())
                ClipM_Helper<1,M,RT>::call(m,thresh);
            else
#endif
                ClipM_Helper<algo2,M,RT>::call(m,thresh);
        } 
    };

    template <class M, class RT>
    inline void Clip(BaseMatrix_Rec_Mutable<M>& m, const RT& thresh)
    {
#if TMV_OPT == 0
        const int algo = M::mcolmajor ? 3 : 2;
#else
        const int algo = M::mcanlin ? 1 : 4;
#endif
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
        static void call(M& m1, const F& f) 
        {
            const int m = m1.colsize();
            for(int i=0;i<m;++i) m1.get_row(i).applyToAll(f);
        } 
    };

    // algo 3: ColMajor
    template <class M, class F> 
    struct ApplyToAllM_Helper<3,M,F> 
    {
        static void call(M& m, const F& f) 
        {
            const int n = m.rowsize();
            for(int j=0;j<n;++j) m.get_col(j).applyToAll(f);
        } 
    };

    // algo 4: Unknown sizes, determine which algorithm to use
    template <class M, class F> 
    struct ApplyToAllM_Helper<4,M,F>
    {
        static void call(M& m, const F& f) 
        {
            const int algo2 = M::mcolmajor ? 3 : 2;
#if TMV_OPT >= 2
            if (m.canLinearize())
                ApplyToAllM_Helper<1,M,F>::call(m,f);
            else
#endif
                ApplyToAllM_Helper<algo2,M,F>::call(m,f);
        } 
    };

    template <class M, class F>
    inline void ApplyToAll(BaseMatrix_Mutable<M>& m, const F& f)
    {
#if TMV_OPT == 0
        const int algo = M::mcolmajor ? 3 : 2;
#else
        const int algo = M::mcanlin ? 1 : 4;
#endif
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
    inline void Scale(const Scaling<ix,T>& x, BaseMatrix_Mutable<M>& m);

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

    // algo 4: Unknown sizes, determine which algorithm to use
    template <class M>
    struct ConjugateM_Helper<4,M>
    {
        static void call(M& m)
        {
#if TMV_OPT >= 2
            if (m.canLinearize())
                ConjugateM_Helper<1,M>::call(m);
            else
#endif
                ConjugateM_Helper<2,M>::call(m);
        }
    };

    template <class M>
    inline void ConjugateSelf(BaseMatrix_Rec_Mutable<M>& m)
    {
        const bool misreal = Traits<typename M::value_type>::isreal;
#if TMV_OPT == 0
        const int algo = misreal ? 0 : 2;
#else
        const int algo = misreal ? 0 : M::mcanlin ? 1 : 4;
#endif
        ConjugateM_Helper<algo,M>::call(m.mat());
    }


    //
    // TMV_Text 
    //

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Rec<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Rec< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Rec_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Rec_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
