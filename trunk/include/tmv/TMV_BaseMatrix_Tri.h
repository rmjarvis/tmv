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
// This file defines the BaseMatrix_Tri and BaseMatrix_Tri classes.
//
// See TMV_TriMatrix.h for the functions that are defined for these objects.
//

#ifndef TMV_BaseMatrix_Tri_H
#define TMV_BaseMatrix_Tri_H

#include "TMV_BaseMatrix.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Diag.h"

namespace tmv {

    template <class M>
    class BaseMatrix_Tri;
    template <class M>
    class BaseMatrix_Tri_Mutable;

    template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor,
              IndexStyle I=CStyle>
    class UpperTriMatrix;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN,
              bool C=false, IndexStyle I=CStyle>
    class ConstUpperTriMatrixView;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false, IndexStyle I=CStyle>
    class UpperTriMatrixView;
    template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor,
              IndexStyle I=CStyle>
    class SmallUpperTriMatrix;
    template <class T, int N, DiagType D, int Si, int Sj, 
              bool C=false, IndexStyle I=CStyle>
    class ConstSmallUpperTriMatrixView;
    template <class T, int N, DiagType D, int Si, int Sj, 
              bool C=false, IndexStyle I=CStyle>
    class SmallUpperTriMatrixView;

    template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor,
              IndexStyle I=CStyle>
    class LowerTriMatrix;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false, IndexStyle I=CStyle>
    class ConstLowerTriMatrixView;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false, IndexStyle I=CStyle>
    class LowerTriMatrixView;
    template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor,
              IndexStyle I=CStyle>
    class SmallLowerTriMatrix;
    template <class T, int N, DiagType D, int Si, int Sj, 
              bool C=false, IndexStyle I=CStyle>
    class ConstSmallLowerTriMatrixView;
    template <class T, int N, DiagType D, int Si, int Sj, 
              bool C=false, IndexStyle I=CStyle>
    class SmallLowerTriMatrixView;

    template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor>
    class UpperTriMatrixF;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false>
    class ConstUpperTriMatrixViewF;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false>
    class UpperTriMatrixViewF;
    template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor>
    class SmallUpperTriMatrixF;
    template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
    class ConstSmallUpperTriMatrixViewF;
    template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
    class SmallUpperTriMatrixViewF;

    template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor>
    class LowerTriMatrixF;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false>
    class ConstLowerTriMatrixViewF;
    template <class T, DiagType D=UnknownDiag, int Si=UNKNOWN, int Sj=UNKNOWN, 
              bool C=false>
    class LowerTriMatrixViewF;
    template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor>
    class SmallLowerTriMatrixF;
    template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
    class ConstSmallLowerTriMatrixViewF;
    template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
    class SmallLowerTriMatrixViewF;

    // Specify ExactSameStorage for triangle matrices:
    template <class M1, class M2>
    static bool ExactSameStorage(
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()); 
    }
    template <class M1, class M2>
    static bool OppositeStorage(
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepj() && m1.stepj() == m2.stepi()); 
    }

    template <class M1, class M2>
    static bool ExactSameStorage(
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()); 
    }
    template <class M1, class M2>
    static bool OppositeStorage(
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepj() && m1.stepj() == m2.stepi()); 
    }

    template <class M1, class M2>
    static bool ExactSameStorage(
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()); 
    }
    template <class M1, class M2>
    static bool OppositeStorage(
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepj() && m1.stepj() == m2.stepi()); 
    }

    template <int shape>
    static void CheckTri(const bool unit, int i, int j)
    { // UpperTri
        TMVAssert(i <= j && "i,j must be in upper triangle");
        TMVAssert((!unit || i < j) && "i,j must be in strict upper triangle");
    }
    template <>
    static void CheckTri<LowerTri>(const bool unit, int i, int j)
    { // LowerTri
        TMVAssert(i >= j && "i,j must be in lower triangle");
        TMVAssert((!unit || i > j) && "i,j must be in strict lower triangle");
    }
    template <>
    static void CheckTri<UnitUpperTri>(const bool , int i, int j)
    { // UnitUpperTri
        TMVAssert(i < j && "i,j must be in strict upper triangle");
    }
    template <>
    static void CheckTri<UnitLowerTri>(const bool , int i, int j)
    { // UnitLowerTri
        TMVAssert(i > j && "i,j must be in strict upper triangle");
    }
    template <int shape>
    static void CheckDiagTri(const bool unit, int i)
    { // UpperTri
        TMVAssert(i >= 0 && 
                  "sub-diagonal not allowed for upper triangle matrix");
        TMVAssert(!(unit && i == 0) && 
                  "main diagonal not allowed for UnitDiag triangle matrix");
    }
    template <>
    static void CheckDiagTri<LowerTri>(const bool unit, int i)
    { // LowerTri
        TMVAssert(i <= 0 && 
                  "super-diagonal not allowed for lower triangle matrix");
        TMVAssert(!(unit && i == 0) && 
                  "main diagonal not allowed for UnitDiag triangle matrix");
    }
    template <>
    static void CheckDiagTri<UnitUpperTri>(const bool , int i)
    { // UpperTri
        TMVAssert(i >= 0 && 
                  "sub-diagonal not allowed for upper triangle matrix");
        TMVAssert(i != 0 && 
                  "main diagonal not allowed for UnitDiag triangle matrix");
    }
    template <>
    static void CheckDiagTri<UnitLowerTri>(const bool , int i)
    { // LowerTri
        TMVAssert(i <= 0 && 
                  "super-diagonal not allowed for lower triangle matrix");
        TMVAssert(i != 0 && 
                  "main diagonal not allowed for UnitDiag triangle matrix");
    }
    static inline void CheckOffDiag(int i, int n) 
    { TMVAssert(i>0 && i<=n && "offDiag index is not valid"); } 
    static inline void CheckOffDiag(int n) 
    { TMVAssert(n>0 && "offDiag not possible for zero-sized TriMatrix"); } 

    template <class T, int cs, int rs, bool rm, bool fort>
    struct MCopyHelper<T,UpperTri,cs,rs,rm,fort>
    {
        typedef SmallUpperTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,UpperTri,UNKNOWN,rs,rm,fort>
    {
        typedef SmallUpperTriMatrix<T,rs,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,UpperTri,cs,UNKNOWN,rm,fort>
    { 
        typedef SmallUpperTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,UpperTri,UNKNOWN,UNKNOWN,rm,fort>
    { 
        typedef UpperTriMatrix<T,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };

    template <class T, int cs, int rs, bool rm, bool fort>
    struct MCopyHelper<T,LowerTri,cs,rs,rm,fort>
    { 
        typedef SmallLowerTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,LowerTri,UNKNOWN,rs,rm,fort>
    { 
        typedef SmallLowerTriMatrix<T,rs,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,LowerTri,cs,UNKNOWN,rm,fort>
    { 
        typedef SmallLowerTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,LowerTri,UNKNOWN,UNKNOWN,rm,fort>
    { 
        typedef LowerTriMatrix<T,NonUnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };

    template <class T, int cs, int rs, bool rm, bool fort>
    struct MCopyHelper<T,UnitUpperTri,cs,rs,rm,fort>
    { 
        typedef SmallUpperTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,UnitUpperTri,UNKNOWN,rs,rm,fort>
    { 
        typedef SmallUpperTriMatrix<T,rs,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,UnitUpperTri,cs,UNKNOWN,rm,fort>
    { 
        typedef SmallUpperTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,UnitUpperTri,UNKNOWN,UNKNOWN,rm,fort>
    { 
        typedef UpperTriMatrix<T,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };

    template <class T, int cs, int rs, bool rm, bool fort>
    struct MCopyHelper<T,UnitLowerTri,cs,rs,rm,fort>
    { 
        typedef SmallLowerTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,UnitLowerTri,UNKNOWN,rs,rm,fort>
    { 
        typedef SmallLowerTriMatrix<T,rs,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,UnitLowerTri,cs,UNKNOWN,rm,fort>
    { 
        typedef SmallLowerTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,UnitLowerTri,UNKNOWN,UNKNOWN,rm,fort>
    { 
        typedef LowerTriMatrix<T,UnitDiag,rm?RowMajor:ColMajor,
                fort?FortranStyle:CStyle> type; 
    };


    // Defined in TMV_CopyU.h
    template <class M1, class M2> 
    static void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);
    template <class M1, class M2> 
    static void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

    // Defiend below
    template <class M1, class M2> 
    static void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2> 
    static void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);

    // Defined in TMV_SwapU.h
    template <class M1, class M2>
    static void Swap(
        BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

    // Defined in TMV_Norm.h
    template <class M>
    static typename M::float_type DoNormF(const BaseMatrix_Tri<M>& m);

    // Defined in TMV_NormU.h
    template <class M>
    static typename M::real_type DoNormSq(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::float_type DoNormSq(
        const BaseMatrix_Tri<M>& m, const typename M::float_type scale);
    template <class M>
    static typename M::float_type DoNorm1(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::float_type DoNormInf(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::value_type DoSumElements(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::float_type DoSumAbsElements(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::real_type DoSumAbs2Elements(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::float_type DoMaxAbsElement(const BaseMatrix_Tri<M>& m);
    template <class M>
    static typename M::real_type DoMaxAbs2Element(const BaseMatrix_Tri<M>& m);

    // Defined in TMV_TriMatrixIO.h
    template <class M>
    static void WriteCompact(std::ostream& os, const BaseMatrix_Tri<M>& m);
    template <class M>
    static void WriteCompact(
        std::ostream& os,
        const BaseMatrix_Tri<M>& m, typename M::float_type thresh) ;
    template <class M>
    static void Read(std::istream& is, BaseMatrix_Tri_Mutable<M>& m);

    // Defined below:
    template <class M>
    static void SetZero(BaseMatrix_Tri_Mutable<M>& m);
    template <class M, class T>
    static void SetAllTo(BaseMatrix_Tri_Mutable<M>& m, const T& val);
    template <class M, class T>
    static void AddToAll(BaseMatrix_Tri_Mutable<M>& m, const T& val);
    template <class M, class RT>
    static void Clip(BaseMatrix_Tri_Mutable<M>& m, const RT& thresh);
    template <class M, class F>
    static void ApplyToAll(BaseMatrix_Tri_Mutable<M>& m, const F& f);
    template <class M>
    static void ConjugateSelf(BaseMatrix_Tri_Mutable<M>& m);
    template <class M>
    static void InvertSelf(BaseMatrix_Tri_Mutable<M>& m);

    template <class M>
    class BaseMatrix_Tri : 
        public BaseMatrix_Calc<M>
    {
    public:
        enum { _colsize = Traits<M>::_size };
        enum { _rowsize = Traits<M>::_size };
        enum { _size = Traits<M>::_size };
        enum { _shape = Traits<M>::_shape };
        enum { _unit = Traits<M>::_unit };
        enum { _unknowndiag = Traits<M>::_unknowndiag };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _stor = Traits<M>::_stor };
        enum { _conj = Traits<M>::_conj };
        enum { _stepi = Traits<M>::_stepi };
        enum { _stepj = Traits<M>::_stepj };
        enum { _diagstep = Traits<M>::_diagstep };

        typedef M type;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;
        typedef typename Traits<M>::inverse_type inverse_type;

        typedef typename Traits<M>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<M>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<M>::const_diag_type const_diag_type;
        typedef typename Traits<M>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<M>::const_subtrimatrix_type 
            const_subtrimatrix_type;
        typedef typename Traits<M>::const_subtrimatrix_step_type 
            const_subtrimatrix_step_type;
        typedef typename Traits<M>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<M>::const_submatrix_step_type 
            const_submatrix_step_type;
        typedef typename Traits<M>::const_subvector_type const_subvector_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_xdview_type const_xdview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;

        typedef typename Traits<M>::const_offdiag_type const_offdiag_type;
        typedef typename Traits<M>::const_unitdiag_type const_unitdiag_type;
        typedef typename Traits<M>::const_nonunitdiag_type 
            const_nonunitdiag_type;
        typedef typename Traits<M>::const_unknowndiag_type 
            const_unknowndiag_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;
        typedef typename Traits<M>::nonconst_type nonconst_type;

        // Derived values:
        enum { _upper = ShapeTraits<Shape(Traits<M>::_shape)>::upper };
        enum { _lower = ShapeTraits<Shape(Traits<M>::_shape)>::lower };
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;
        typedef typename Traits<value_type>::complex_type complex_type;


        //
        // Constructor
        //

        BaseMatrix_Tri() {}
        BaseMatrix_Tri(const BaseMatrix_Tri<M>&) {}
        ~BaseMatrix_Tri() {}

    private:
        void operator=(const BaseMatrix_Tri<M>&);
    public:

        // MJ Temporary -- remove this.
        bool divIsSet() const { return false; }

        //
        // Access 
        //

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
                ( i<0 ?  size()+i : size()-i ),
                diagstep());
        }

        const_diag_sub_type get_diag(int i, int j1, int j2) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        const_row_sub_type row(int i, int j1, int j2) const
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return get_row(i,j1,j2);
        }

        const_col_sub_type col(int j, int i1, int i2) const
        {
            CheckColIndex<_fort>(j,rowsize());
            CheckRowRange<_fort>(i1,i2,colsize());
            return get_col(j,i1,i2);
        }

        const_diag_type diag() const
        {
            CheckDiagTri<_shape>(isunit(),0);
            return const_diag_type(cptr(),size(),diagstep()); 
        }

        const_diag_sub_type diag(int i) const
        {
            CheckDiagIndex<_fort>(i,colsize(),rowsize());
            CheckDiagTri<_shape>(isunit(),i);
            return get_diag(i);
        }

        const_diag_sub_type diag(int i, int j1, int j2) const
        {
            CheckDiagIndex<_fort>(i,j1,j2,colsize(),rowsize());
            CheckDiagTri<_shape>(isunit(),i);
            return get_diag(i,j1,j2);
        }


        //
        // Functions
        //

        value_type sumElements() const
        { return tmv::DoSumElements(cView()); }

        float_type sumAbsElements() const
        { return tmv::DoSumAbsElements(cView()); }

        real_type sumAbs2Elements() const
        { return tmv::DoSumAbs2Elements(cView()); }

        float_type maxAbsElement() const
        { return tmv::DoMaxAbsElement(cView()); }

        real_type maxAbs2Element() const
        { return tmv::DoMaxAbs2Element(cView()); }

        real_type normSq() const
        { return tmv::DoNormSq(cView()); }

        float_type normSq(const float_type scale) const
        { return tmv::DoNormSq(cView(),scale); }

        float_type normF() const 
        { return tmv::DoNormF(cView()); }

        float_type norm() const
        { return normF(); }

        float_type norm1() const
        { return tmv::DoNorm1(cView()); }

        float_type normInf() const
        { return tmv::DoNormInf(cView()); }



        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        const_subtrimatrix_type cSubTriMatrix(int i1, int i2) const
        {
            return const_subtrimatrix_type(
                cptr()+i1*diagstep(), i2-i1, isunit(), stepi(), stepj());
        }

        const_subtrimatrix_step_type cSubTriMatrix(
            int i1, int i2, int istep) const
        {
            return const_subtrimatrix_step_type(
                cptr()+i1*diagstep(), (i2-i1)/istep, isunit(), 
                istep*stepi(), istep*stepj());
        }

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


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        const_submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            // Check each corner of submatrix:
            CheckTri<_shape>(isunit(),i1,j1);
            if (j2 > j1) CheckTri<_shape>(isunit(),i1,j2-1);
            if (i2 > i1) CheckTri<_shape>(isunit(),i2-1,j1);
            if (i2 > i1 && j2 > j1) CheckTri<_shape>(isunit(),i2-1,j2-1);
            return cSubMatrix(i1,i2,j1,j2);
        }

        const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            CheckRowRange<_fort>(i1,i2,istep,colsize());
            CheckColRange<_fort>(j1,j2,jstep,rowsize());
            // Check each corner of submatrix:
            CheckTri<_shape>(isunit(),i1,j1);
            if (j2 > j1) CheckTri<_shape>(isunit(),i1,j2-jstep);
            if (i2 > i1) CheckTri<_shape>(isunit(),i2-istep,j1);
            if (i2 > i1 && j2 > j1) 
                CheckTri<_shape>(isunit(),i2-istep,j2-jstep);
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        const_subvector_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            CheckMatSubVector<_fort>(i,j,istep,jstep,s,colsize(),rowsize());
            // Check each end of subvector:
            CheckTri<_shape>(isunit(),i,j);
            if (s) CheckTri<_shape>(isunit(),i+(s-1)*istep,j+(s-1)*jstep);
            return cSubVector(i,j,istep,jstep,s);
        }

        const_subtrimatrix_type subTriMatrix(int i1, int i2) const
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubTriMatrix(i1,i2);
        }

        const_subtrimatrix_step_type subTriMatrix(
            int i1, int i2, int istep) const
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubTriMatrix(i1,i2,istep);
        }


        //
        // Views
        //

        const_view_type view() const
        { return const_view_type(cptr(),size(),isunit(),stepi(),stepj()); }

        const_cview_type cView() const
        { return view(); }

        const_fview_type fView() const
        { return view(); }

        const_xview_type xView() const
        { return view(); }

        // For Tri Matrices, we add this xdView option that also turns
        // the known D into UnknownDiag
        const_xdview_type xdView() const
        { return view(); }

        const_cmview_type cmView() const
        {
            TMVAssert(iscm() && "Called cmView on non-ColMajor matrix");
            return view(); 
        }

        const_rmview_type rmView() const
        {
            TMVAssert(isrm() && "Called rmView on non-RowMajor matrix");
            return view(); 
        }

        const_view_type constView() const
        { return view(); }

        const_transpose_type transpose() const
        { return const_transpose_type(cptr(),size(),isunit(),stepj(),stepi()); }

        const_conjugate_type conjugate() const
        { return const_conjugate_type(cptr(),size(),isunit(),stepi(),stepj()); }

        const_adjoint_type adjoint() const
        { return const_adjoint_type(cptr(),size(),isunit(),stepj(),stepi()); }

        const_offdiag_type offDiag(int noff) const
        {
            CheckOffDiag(noff,size());
            return const_offdiag_type(
                cptr()+noff*(this->isupper()?stepj():stepi()),
                size()-noff,false,stepi(),stepj()); 
        }

        const_offdiag_type offDiag() const
        {
            CheckOffDiag(size());
            return const_offdiag_type(
                cptr()+(this->isupper()?stepj():stepi()),
                size()-1,false,stepi(),stepj()); 
        }

        const_unitdiag_type viewAsUnitDiag() const
        { return const_unitdiag_type(cptr(),size(),true,stepi(),stepj()); }

        const_nonunitdiag_type viewAsNonUnitDiag() const
        {
            TMVAssert(!isunit());
            return const_nonunitdiag_type(cptr(),size(),false,stepi(),stepj()); 
        }

        const_unknowndiag_type viewAsUnknownDiag() const
        { 
            return const_unknowndiag_type(
                cptr(),size(),isunit(),stepi(),stepj()); 
        }

        const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()), size(), isunit(),
                M::isreal ? stepi() : 2*stepi(),
                M::isreal ? stepj() : 2*stepj());
        }

        const_imagpart_type imagPart() const
        {
            TMVStaticAssert(M::iscomplex);
            TMVAssert(!isunit());
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1, size(), false,
                2*stepi(), 2*stepj());
        }

        const_nonconj_type nonConj() const
        { return const_nonconj_type(cptr(),size(),isunit(),stepi(),stepj()); }

        nonconst_type nonConst() const
        { 
            return nonconst_type(
                const_cast<value_type*>(cptr()),
                size(),isunit(),stepi(),stepj()); 
        }


        //
        // I/O
        //

        void writeCompact(std::ostream& os) const
        { tmv::WriteCompact(os,cView()); }
        void writeCompact(std::ostream& os, float_type thresh) const
        { tmv::WriteCompact(os,cView(),thresh); }



        //
        // Auxilliary routines
        //

        template <class M2>
        void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::Copy(mat(),m2.mat()); 
        }

        template <class M2>
        void newAssignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::NoAliasCopy(mat(),m2.mat()); 
        }

        // Need to do BaseMatrix_Tri_Mutable separately to make sure we
        // check for UnknownDiag correctly.
        // UnknownDiag -> UnitDiag is naively disallowed by the StaticAssert
        // in the above methods.
        template <class M2>
        void assignTo(BaseMatrix_Tri_Mutable<M2>& m2) const
        {
            TMVStaticAssert(!M2::_unit || _unit || _unknowndiag);
            TMVAssert(!m2.isunit() || isunit());
            tmv::Copy(mat(),m2.mat()); 
        }

        template <class M2>
        void newAssignTo(BaseMatrix_Tri_Mutable<M2>& m2) const
        {
            TMVStaticAssert(!M2::_unit || _unit || _unknowndiag);
            TMVAssert(!m2.isunit() || isunit());
            tmv::NoAliasCopy(mat(),m2.mat()); 
        }

        const type& mat() const
        { return *static_cast<const type*>(this); }

        int diagstep() const 
        { return _diagstep == UNKNOWN ? stepi() + stepj() : _diagstep; }
        bool isconj() const { return _conj; }
        bool isunit() const { return mat().isunit(); }
        bool isupper() const { return _upper; }
        bool iscm() const { return mat().iscm(); }
        bool isrm() const { return mat().isrm(); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref, get_row and get_col from BaseMatrix.

        size_t colsize() const { return mat().size(); }
        size_t rowsize() const { return mat().size(); }
        size_t size() const { return mat().size(); }
        int stepi() const { return mat().stepi(); }
        int stepj() const { return mat().stepj(); }

        const value_type* cptr() const { return mat().cptr(); }

    }; // BaseMatrix_Tri

    template <class M>
    class BaseMatrix_Tri_Mutable : 
        public BaseMatrix_Tri<M>,
        public BaseMatrix_Mutable<M>
    {
    public:
        enum { _colsize = Traits<M>::_size };
        enum { _rowsize = Traits<M>::_size };
        enum { _size = Traits<M>::_size };
        enum { _shape = Traits<M>::_shape };
        enum { _unit = Traits<M>::_unit };
        enum { _unknowndiag = Traits<M>::_unknowndiag };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _stor = Traits<M>::_stor };
        enum { _conj = Traits<M>::_conj };
        enum { _stepi = Traits<M>::_stepi };
        enum { _stepj = Traits<M>::_stepj };
        enum { _diagstep = Traits<M>::_diagstep };

        typedef M type;
        typedef BaseMatrix_Tri<M> base_tri;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;

        typedef typename Traits<M>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<M>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<M>::const_diag_type const_diag_type;
        typedef typename Traits<M>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<M>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<M>::const_submatrix_step_type 
            const_submatrix_step_type;
        typedef typename Traits<M>::const_subvector_type const_subvector_type;
        typedef typename Traits<M>::const_subtrimatrix_type 
            const_subtrimatrix_type;
        typedef typename Traits<M>::const_subtrimatrix_step_type 
            const_subtrimatrix_step_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_xdview_type const_xdview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;

        typedef typename Traits<M>::const_offdiag_type const_offdiag_type;
        typedef typename Traits<M>::const_unitdiag_type const_unitdiag_type;
        typedef typename Traits<M>::const_nonunitdiag_type 
            const_nonunitdiag_type;
        typedef typename Traits<M>::const_unknowndiag_type 
            const_unknowndiag_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;

        typedef typename Traits<M>::row_sub_type row_sub_type;
        typedef typename Traits<M>::col_sub_type col_sub_type;
        typedef typename Traits<M>::diag_type diag_type;
        typedef typename Traits<M>::diag_sub_type diag_sub_type;
        typedef typename Traits<M>::submatrix_type submatrix_type;
        typedef typename Traits<M>::submatrix_step_type submatrix_step_type;
        typedef typename Traits<M>::subvector_type subvector_type;
        typedef typename Traits<M>::subtrimatrix_type subtrimatrix_type;
        typedef typename Traits<M>::subtrimatrix_step_type 
            subtrimatrix_step_type;

        typedef typename Traits<M>::view_type view_type;
        typedef typename Traits<M>::cview_type cview_type;
        typedef typename Traits<M>::fview_type fview_type;
        typedef typename Traits<M>::xview_type xview_type;
        typedef typename Traits<M>::xdview_type xdview_type;
        typedef typename Traits<M>::cmview_type cmview_type;
        typedef typename Traits<M>::rmview_type rmview_type;
        typedef typename Traits<M>::transpose_type transpose_type;
        typedef typename Traits<M>::conjugate_type conjugate_type;
        typedef typename Traits<M>::adjoint_type adjoint_type;

        typedef typename Traits<M>::offdiag_type offdiag_type;
        typedef typename Traits<M>::unitdiag_type unitdiag_type;
        typedef typename Traits<M>::nonunitdiag_type nonunitdiag_type;
        typedef typename Traits<M>::unknowndiag_type unknowndiag_type;
        typedef typename Traits<M>::realpart_type realpart_type;
        typedef typename Traits<M>::imagpart_type imagpart_type;
        typedef typename Traits<M>::nonconj_type nonconj_type;

        typedef typename Traits<M>::reference reference;

        // Derived values:
        enum { _upper = ShapeTraits<Shape(Traits<M>::_shape)>::upper };
        enum { _lower = ShapeTraits<Shape(Traits<M>::_shape)>::lower };
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;
        typedef typename Traits<value_type>::complex_type complex_type;


        //
        // Constructor
        //

        BaseMatrix_Tri_Mutable() {}
        BaseMatrix_Tri_Mutable(const BaseMatrix_Tri_Mutable<M>&) {}
        ~BaseMatrix_Tri_Mutable() {}


        //
        // Access 
        //

        reference operator()(int i, int j)
        {
            CheckRowIndex<_fort>(i,size());
            CheckColIndex<_fort>(j,size());
            // MJ: We don't actually want this next check.
            // Let TriRef do the check if the reference is assigned to.
            //CheckTri<_shape>(isunit(),i,j);
            return ref(i,j);
        }

        row_sub_type get_row(int i, int j1, int j2) 
        { return row_sub_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

        col_sub_type get_col(int j, int i1, int i2) 
        { return col_sub_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

        diag_sub_type get_diag(int i) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ?  size()+i : size()-i ),
                diagstep());
        }

        diag_sub_type get_diag(int i, int j1, int j2) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        row_sub_type row(int i, int j1, int j2) 
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return get_row(i,j1,j2);
        }

        col_sub_type col(int j, int i1, int i2) 
        {
            CheckColIndex<_fort>(j,rowsize());
            CheckRowRange<_fort>(i1,i2,colsize());
            return get_col(j,i1,i2);
        }

        diag_type diag() 
        {
            CheckDiagTri<_shape>(isunit(),0);
            return diag_type(ptr(),size(),diagstep()); 
        }

        diag_sub_type diag(int i) 
        {
            CheckDiagIndex<_fort>(i,colsize(),rowsize());
            CheckDiagTri<_shape>(isunit(),i);
            return get_diag(i);
        }

        diag_sub_type diag(int i, int j1, int j2) 
        {
            CheckDiagIndex<_fort>(i,j1,j2,colsize(),rowsize());
            CheckDiagTri<_shape>(isunit(),i);
            return get_diag(i,j1,j2);
        }

        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        value_type operator()(int i, int j) const
        { return base_tri::operator()(i,j); }

        const_row_sub_type get_row(int i, int j1, int j2) const
        { return base_tri::get_row(i,j1,j2); }
        const_col_sub_type get_col(int j, int i1, int i2) const
        { return base_tri::get_col(j,i1,i2); }
        const_diag_sub_type get_diag(int i) const
        { return base_tri::get_diag(i); }
        const_diag_sub_type get_diag(int i, int j1, int j2) const
        { return base_tri::get_diag(i,j1,j2); }

        const_row_sub_type row(int i, int j1, int j2) const
        { return base_tri::row(i,j1,j2); }
        const_col_sub_type col(int j, int i1, int i2) const
        { return base_tri::col(j,i1,i2); }
        const_diag_type diag() const
        { return base_tri::diag(); }
        const_diag_sub_type diag(int i) const
        { return base_tri::diag(i); }
        const_diag_sub_type diag(int i, int j1, int j2) const
        { return base_tri::diag(i,j1,j2); }


        //
        // Op =
        //

        type& operator=(BaseMatrix_Tri_Mutable<M>& m2) 
        {
            TMVAssert(size() == m2.size());
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_size,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_size,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        type& operator=(const BaseMatrix_Tri<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_size,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_size,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVStaticAssert(!M2::_unit || _unit || _unknowndiag);
            TMVAssert(!m2.isunit() || isunit());
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        type& operator=(const BaseMatrix_Diag<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_size,M2::_size>::same));
            TMVAssert(size() == m2.size());
            m2.diag().assignTo(diag());
            offDiag().setZero();
            return mat(); 
        }

        type& operator=(const value_type x)
        {
            TMVAssert(colsize() == rowsize());
            setToIdentity(x);
            return mat();
        }


        //
        // Modifying Functions
        //

        // These are required for all BaseMatrix_Mutable
        type& setZero()
        { tmv::SetZero(mat()); return mat(); }

        type& setAllTo(value_type val)
        { tmv::SetAllTo(mat(),val); return mat(); }

        type& addToAll(value_type val)
        { tmv::AddToAll(mat(),val); return mat(); }

        type& clip(float_type thresh)
        { tmv::Clip(mat(),thresh); return mat(); }

        template <class F>
        type& applyToAll(const F& f)
        { tmv::ApplyToAll(mat(),f); return mat(); }

        type& conjugateSelf()
        { tmv::ConjugateSelf(mat()); return mat(); }

        // Some more that are added for Tri shape:
        type& invertSelf() 
        { tmv::InvertSelf(mat()); return mat(); }

        type& setToIdentity(const value_type x=value_type(1))
        {
            offDiag().setZero();
            if(!isunit()) diag().setAllTo(x);
            return mat(); 
        }


        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        subtrimatrix_type cSubTriMatrix(int i1, int i2) 
        {
            return subtrimatrix_type(
                ptr()+i1*diagstep(), i2-i1, isunit(), stepi(), stepj()); 
        }

        subtrimatrix_step_type cSubTriMatrix(int i1, int i2, int istep) 
        {
            return subtrimatrix_step_type(
                ptr()+i1*diagstep(), (i2-i1)/istep, isunit(),
                istep*stepi(), istep*stepj());
        }

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

        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        submatrix_type subMatrix(int i1, int i2, int j1, int j2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            // Check each corner of submatrix:
            CheckTri<_shape>(isunit(),i1,j1);
            if (j2 > j1) CheckTri<_shape>(isunit(),i1,j2-1);
            if (i2 > i1) CheckTri<_shape>(isunit(),i2-1,j1);
            if (i2 > i1 && j2 > j1) CheckTri<_shape>(isunit(),i2-1,j2-1);
            return cSubMatrix(i1,i2,j1,j2);
        }

        submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            CheckRowRange<_fort>(i1,i2,istep,colsize());
            CheckColRange<_fort>(j1,j2,jstep,rowsize());
            // Check each corner of submatrix:
            CheckTri<_shape>(isunit(),i1,j1);
            if (j2 > j1) CheckTri<_shape>(isunit(),i1,j2-jstep);
            if (i2 > i1) CheckTri<_shape>(isunit(),i2-istep,j1);
            if (i2 > i1 && j2 > j1) 
                CheckTri<_shape>(isunit(),i2-istep,j2-jstep);
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        subvector_type subVector(
            int i, int j, int istep, int jstep, int s) 
        {
            CheckMatSubVector<_fort>(i,j,istep,jstep,s,colsize(),rowsize());
            // Check each end of subvector:
            CheckTri<_shape>(isunit(),i,j);
            if (s) CheckTri<_shape>(isunit(),i+(s-1)*istep,j+(s-1)*jstep);
            return cSubVector(i,j,istep,jstep,s);
        }

        subtrimatrix_type subTriMatrix(int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubTriMatrix(i1,i2);
        }

        subtrimatrix_step_type subTriMatrix(int i1, int i2, int istep) 
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubTriMatrix(i1,i2,istep);
        }


        // Repeat the const versions:
        const_submatrix_type cSubMatrix(
            int i1, int i2, int j1, int j2) const
        { return base_tri::cSubMatrix(i1,i2,j1,j2); }
        const_submatrix_step_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return base_tri::cSubMatrix(i1,i2,j1,j2,istep,jstep); }
        const_subvector_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return base_tri::cSubVector(i,j,istep,jstep,s); }
        const_subtrimatrix_type cSubTriMatrix(int i1, int i2) const
        { return base_tri::cSubTriMatrix(i1,i2); }
        const_subtrimatrix_step_type cSubTriMatrix(
            int i1, int i2, int istep) const
        { return base_tri::cSubTriMatrix(i1,i2,istep); }

        const_submatrix_type subMatrix(
            int i1, int i2, int j1, int j2) const
        { return base_tri::subMatrix(i1,i2,j1,j2); }
        const_submatrix_step_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        { return base_tri::subMatrix(i1,i2,j1,j2,istep,jstep); }
        const_subvector_type subVector(
            int i, int j, int istep, int jstep, int s) const
        { return base_tri::subVector(i,j,istep,jstep,s); }
        const_subtrimatrix_type subTriMatrix(int i1, int i2) const
        { return base_tri::subTriMatrix(i1,i2); }
        const_subtrimatrix_step_type subTriMatrix(
            int i1, int i2, int istep) const
        { return base_tri::subTriMatrix(i1,i2,istep); }


        //
        // Views
        //

        view_type view() 
        { return view_type(ptr(),size(),isunit(),stepi(),stepj()); }

        cview_type cView() 
        { return view(); }

        fview_type fView() 
        { return view(); }

        xview_type xView() 
        { return view(); }

        xdview_type xdView() 
        { return view(); }

        cmview_type cmView() 
        {
            TMVAssert(iscm() && "Called cmView on non-ColMajor matrix");
            return view(); 
        }

        rmview_type rmView() 
        {
            TMVAssert(isrm() && "Called rmView on non-RowMajor matrix");
            return view(); 
        }

        transpose_type transpose() 
        { return transpose_type(ptr(),size(),isunit(),stepj(),stepi()); }

        conjugate_type conjugate() 
        { return conjugate_type(ptr(),size(),isunit(),stepi(),stepj()); }

        adjoint_type adjoint() 
        { return adjoint_type(ptr(),size(),isunit(),stepj(),stepi()); }

        offdiag_type offDiag(int noff) 
        {
            CheckOffDiag(noff,size());
            return offdiag_type(
                ptr()+noff*(this->isupper()?stepj():stepi()),
                size()-noff,false,stepi(),stepj()); 
        }

        offdiag_type offDiag()
        {
            CheckOffDiag(size());
            return offdiag_type(
                ptr()+(this->isupper()?stepj():stepi()),
                size()-1,false,stepi(),stepj()); 
        }

        unitdiag_type viewAsUnitDiag()
        { return unitdiag_type(ptr(),size(),true,stepi(),stepj()); }

        nonunitdiag_type viewAsNonUnitDiag()
        {
            TMVAssert(!isunit());
            return nonunitdiag_type(ptr(),size(),false,stepi(),stepj()); 
        }

        unknowndiag_type viewAsUnknownDiag()
        { return unknowndiag_type(ptr(),size(),isunit(),stepi(),stepj()); }

        realpart_type realPart() 
        {
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()), size(), isunit(),
                M::isreal ? stepi() : 2*stepi(),
                M::isreal ? stepj() : 2*stepj());
        }

        imagpart_type imagPart() 
        {
            TMVStaticAssert(M::iscomplex);
            TMVAssert(!isunit());
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1, size(), false,
                2*stepi(), 2*stepj());
        }

        nonconj_type nonConj()
        { return nonconj_type(ptr(),size(),isunit(),stepi(),stepj()); }


        // Repeat the const versions:
        const_view_type view() const
        { return base_tri::view(); }
        const_cview_type cView() const
        { return base_tri::cView(); }
        const_fview_type fView() const
        { return base_tri::fView(); }
        const_xview_type xView() const
        { return base_tri::xView(); }
        const_xdview_type xdView() const
        { return base_tri::xdView(); }
        const_cmview_type cmView() const
        { return base_tri::cmView(); }
        const_rmview_type rmView() const
        { return base_tri::rmView(); }
        const_transpose_type transpose() const
        { return base_tri::transpose(); }
        const_conjugate_type conjugate() const
        { return base_tri::conjugate(); }
        const_adjoint_type adjoint() const
        { return base_tri::adjoint(); }
        const_offdiag_type offDiag(int noff) const
        { return base_tri::offDiag(noff); }
        const_offdiag_type offDiag() const
        { return base_tri::offDiag(); }
        const_unitdiag_type viewAsUnitDiag() const
        { return base_tri::viewAsUnitDiag(); }
        const_nonunitdiag_type viewAsNonUnitDiag() const
        { return base_tri::viewAsNonUnitDiag(); }
        const_unknowndiag_type viewAsUnknownDiag() const
        { return base_tri::viewAsUnknownDiag(); }
        const_realpart_type realPart() const
        { return base_tri::realPart(); }
        const_imagpart_type imagPart() const
        { return base_tri::imagPart(); }
        const_nonconj_type nonConj() const
        { return base_tri::nonConj(); }


        //
        // I/O
        //

        void read(std::istream& is)
        {
            cview_type mcv = cView();
            tmv::Read(is,mcv); 
        }


        //
        // Auxilliary routines
        //

        const type& mat() const
        { return *static_cast<const type*>(this); }
        type& mat()
        { return *static_cast<type*>(this); }

        int diagstep() const
        { return _diagstep == UNKNOWN ? stepi() + stepj() : _diagstep; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref and cptr from above.

        size_t colsize() const { return mat().size(); }
        size_t rowsize() const { return mat().size(); }
        size_t size() const { return mat().size(); }
        bool isunit() const { return mat().isunit(); }
        int stepi() const { return mat().stepi(); }
        int stepj() const { return mat().stepj(); }
        bool iscm() const { return mat().iscm(); }
        bool isrm() const { return mat().isrm(); }

        value_type* ptr() { return mat().ptr(); }
        reference ref(int i, int j) { return mat().ref(i,j); }
        value_type cref(int i, int j) { return mat().cref(i,j); }

    }; // BaseMatrix_Tri_Mutable

    //
    // setZero
    //

    template <int algo, class M> 
    struct SetZeroU_Helper;

    // algo 1: LowerTri -- Transpose
    template <class M> 
    struct SetZeroU_Helper<1,M> 
    {
        static void call(M& m) 
        { m.transpose().setZero(); } 
    };

    // algo 2: RowMajor NonUnit
    template <class M> 
    struct SetZeroU_Helper<2,M> 
    {
        static void call(M& m) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i,n).setZero();
        } 
    };

    // algo 3: ColMajor NonUnit
    template <class M> 
    struct SetZeroU_Helper<3,M> 
    {
        static void call(M& m) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j+1).setZero();
        } 
    };

    // algo 4: RowMajor Unit
    template <class M> 
    struct SetZeroU_Helper<4,M> 
    {
        static void call(M& m) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i+1,n).setZero();
        } 
    };

    // algo 5: ColMajor Unit
    template <class M> 
    struct SetZeroU_Helper<5,M> 
    {
        static void call(M& m) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j).setZero();
        } 
    };

    // algo 6: UnknownDiag
    template <class M> 
    struct SetZeroU_Helper<6,M> 
    {
        static void call(M& m) 
        {
            if (m.isunit()) {
                typename M::unitdiag_type mu = m.viewAsUnitDiag();
                SetZero(mu);
            } else {
                typename M::nonunitdiag_type mnu = m.viewAsNonUnitDiag();
                SetZero(mnu);
            }
        }
    };

    template <class M>
    static void SetZero(BaseMatrix_Tri_Mutable<M>& m)
    {
        const int algo = 
            M::_lower ? 1 :
            M::_unknowndiag ? 6 :
            M::_unit ? ( M::_rowmajor ? 4 : 5 ) :
            ( M::_rowmajor ? 2 : 3 );
        SetZeroU_Helper<algo,M>::call(m.mat());
    }

    //
    // setAllTo
    //

    template <int algo, class M, class T> 
    struct SetAllToU_Helper;

    // algo 1: LowerTri -- Transpose
    template <class M, class T> 
    struct SetAllToU_Helper<1,M,T>
    { 
        static void call(M& m, const T& val) 
        { m.transpose().setAllTo(val); } 
    };

    // algo 2: RowMajor NonUnit
    template <class M, class T> 
    struct SetAllToU_Helper<2,M,T>
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i,n).setAllTo(val);
        } 
    };

    // algo 3: ColMajor
    template <class M, class T> 
    struct SetAllToU_Helper<3,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j+1).setAllTo(val);
        } 
    };

    // algo 4: RowMajor Unit
    template <class M, class T> 
    struct SetAllToU_Helper<4,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i+1,n).setAllTo(val);
        } 
    };

    // algo 5: ColMajor Unit
    template <class M, class T> 
    struct SetAllToU_Helper<5,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j).setAllTo(val);
        } 
    };

    // algo 6: UnknownDiag
    template <class M, class T> 
    struct SetAllToU_Helper<6,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            if (m.isunit()) {
                typename M::unitdiag_type mu = m.viewAsUnitDiag();
                SetAllTo(mu,val);
            } else {
                typename M::nonunitdiag_type mnu = m.viewAsNonUnitDiag();
                SetAllTo(mnu,val);
            }
        }
    };

    template <class M, class T>
    static void SetAllTo(BaseMatrix_Tri_Mutable<M>& m, const T& val)
    {
        TMVAssert(!m.isunit() || val == T(1));
        const int algo = 
            M::_lower ? 1 :
            M::_unknowndiag ? 6 :
            M::_unit ? ( M::_rowmajor ? 4 : 5 ) :
            ( M::_rowmajor ? 2 : 3 );
        SetAllToU_Helper<algo,M,T>::call(m.mat(),val);
    }

    //
    // addToAll
    //

    template <int algo, class M, class T> 
    struct AddToAllU_Helper;

    // algo 1: LowerTri -- Transpose
    template <class M, class T> 
    struct AddToAllU_Helper<1,M,T> 
    { 
        static void call(M& m, const T& val) 
        { m.transpose().addToAll(val); } 
    };

    // algo 2: RowMajor NonUnit
    template <class M, class T> 
    struct AddToAllU_Helper<2,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i,n).addToAll(val);
        } 
    };

    // algo 3: ColMajor
    template <class M, class T> 
    struct AddToAllU_Helper<3,M,T> 
    {
        static void call(M& m, const T& val) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j+1).addToAll(val);
        } 
    };

    template <class M, class T>
    static void AddToAll(BaseMatrix_Tri_Mutable<M>& m, const T& val)
    {
        TMVStaticAssert(!M::_unit);
        TMVAssert(!m.isunit());
        const int algo = 
            M::_lower ? 1 :
            M::_rowmajor ? 2 : 3;
        AddToAllU_Helper<algo,M,T>::call(m.mat(),val);
    }

    //
    // clip
    //

    template <int algo, class M, class RT> 
    struct ClipU_Helper;

    // algo 1: LowerTri -- Transpose
    template <class M, class RT> 
    struct ClipU_Helper<1,M,RT>
    {
        static void call(M& m, const RT& thresh) 
        { m.transpose().clip(thresh); } 
    };

    // algo 2: RowMajor NonUnit
    template <class M, class RT> 
    struct ClipU_Helper<2,M,RT> 
    {
        static void call(M& m, const RT& thresh) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i,n).clip(thresh);
        } 
    };

    // algo 3: ColMajor
    template <class M, class RT> 
    struct ClipU_Helper<3,M,RT> 
    {
        static void call(M& m, const RT& thresh) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j+1).clip(thresh);
        } 
    };

    // algo 4: RowMajor Unit
    template <class M, class RT> 
    struct ClipU_Helper<4,M,RT> 
    {
        static void call(M& m, const RT& thresh) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i+1,n).clip(thresh);
        } 
    };

    // algo 5: ColMajor Unit
    template <class M, class RT> 
    struct ClipU_Helper<5,M,RT> 
    {
        static void call(M& m, const RT& thresh) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j).clip(thresh);
        } 
    };

    // algo 6: UnknownDiag
    template <class M, class RT> 
    struct ClipU_Helper<6,M,RT> 
    {
        static void call(M& m, const RT& thresh) 
        {
            if (m.isunit()) {
                typename M::unitdiag_type mu = m.viewAsUnitDiag();
                Clip(mu,thresh);
            } else {
                typename M::nonunitdiag_type mnu = m.viewAsNonUnitDiag();
                Clip(mnu,thresh);
            }
        }
    };

    template <class M, class RT>
    static void Clip(BaseMatrix_Tri_Mutable<M>& m, const RT& thresh)
    {
        const int algo = 
            M::_lower ? 1 :
            M::_unknowndiag ? 6 :
            M::_unit ? ( M::_rowmajor ? 4 : 5 ) :
            ( M::_rowmajor ? 2 : 3 );
        ClipU_Helper<algo,M,RT>::call(m.mat(),thresh);
    }

    //
    // applyToAll
    //

    template <int algo, class M, class F> 
    struct ApplyToAllU_Helper;

    // algo 1: LowerTri -- Transpose
    template <class M, class F> 
    struct ApplyToAllU_Helper<1,M,F>
    { 
        static void call(M& m, const F& f) 
        { m.transpose().applyToAll(f); } 
    };

    // algo 2: RowMajor NonUnit
    template <class M, class F> 
    struct ApplyToAllU_Helper<2,M,F> 
    {
        static void call(M& m, const F& f) 
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i,n).applyToAll(f);
        } 
    };

    // algo 3: ColMajor
    template <class M, class F> 
    struct ApplyToAllU_Helper<3,M,F> 
    {
        static void call(M& m, const F& f) 
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j+1).applyToAll(f);
        } 
    };

    template <class M, class F>
    static void ApplyToAll(BaseMatrix_Tri_Mutable<M>& m, const F& f)
    {
        TMVStaticAssert(!M::_unit);
        TMVAssert(!m.isunit());
        const int algo = 
            M::_lower ? 1 :
            M::_rowmajor ? 2 : 3;
        ApplyToAllU_Helper<algo,M,F>::call(m.mat(),f);
    }

    //
    // conjugateSelf
    //

    template <int algo, class M> 
    struct ConjugateU_Helper;

    // algo 0: Not complex, nothing to do
    template <class M>
    struct ConjugateU_Helper<0,M>
    { static void call(M& ) {} };

    // algo 1: LowerTri -- Transpose
    template <class M>
    struct ConjugateU_Helper<1,M>
    {
        static void call(M& m)
        { m.transpose().conjugateSelf(); }
    };

    // algo 2: RowMajor NonUnit
    template <class M>
    struct ConjugateU_Helper<2,M> 
    {
        static void call(M& m)
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i,n).conjugateSelf();
        } 
    };

    // algo 3: ColMajor
    template <class M>
    struct ConjugateU_Helper<3,M> 
    {
        static void call(M& m)
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j+1).conjugateSelf();
        } 
    };

    // algo 4: RowMajor Unit
    template <class M> 
    struct ConjugateU_Helper<4,M> 
    {
        static void call(M& m)
        {
            const int n = m.size();
            for(int i=0;i<n;++i) m.get_row(i,i+1,n).conjugateSelf();
        } 
    };

    // algo 5: ColMajor Unit
    template <class M> 
    struct ConjugateU_Helper<5,M> 
    {
        static void call(M& m)
        {
            const int n = m.size();
            for(int j=0;j<n;++j) m.get_col(j,0,j).conjugateSelf();
        } 
    };

    // algo 6: UnknownDiag
    template <class M> 
    struct ConjugateU_Helper<6,M> 
    {
        static void call(M& m)
        {
            if (m.isunit()) {
                typename M::unitdiag_type mu = m.viewAsUnitDiag();
                ConjugateSelf(mu);
            } else {
                typename M::nonunitdiag_type mnu = m.viewAsNonUnitDiag();
                ConjugateSelf(mnu);
            }
        } 
    };

    template <class M>
    static void ConjugateSelf(BaseMatrix_Tri_Mutable<M>& m)
    {
        const int algo = 
            M::_lower ? 1 :
            M::_unknowndiag ? 6 :
            M::_unit ? ( M::_rowmajor ? 4 : 5 ) : 
            ( M::_rowmajor ? 2 : 3 );
        ConjugateU_Helper<algo,M>::call(m.mat());
    }



    //
    // Trace
    //

    // The BaseMatrix Trace call is efficient for composite types, since
    // it avoid calculating all the elements to do the sum.
    // But if we do have the elements calculated, this overloaded 
    // version will be faster:
    // Note that this TriMatrix version needs to check whether the 
    // diagonal elements are all 1 or not.
    
    template <int algo, class M> 
    struct TraceU_Helper;

    // algo 1: Unit Diag
    template <class M>
    struct TraceU_Helper<1,M> 
    {
        static typename M::value_type call(const M& m)
        { return typename M::value_type(m.size()); }
    };

    // algo 2: NonUnit Diag
    template <class M>
    struct TraceU_Helper<2,M> 
    {
        static typename M::value_type call(const M& m)
        { return m.diag().sumElements(); }
    };

    // algo 3: Unknown Diag
    template <class M>
    struct TraceU_Helper<3,M> 
    {
        static typename M::value_type call(const M& m)
        { 
            if (m.isunit()) return TraceU_Helper<1,M>::call(m);
            else return TraceU_Helper<2,M>::call(m);
        }
    };


    template <class M>
    static typename M::value_type DoTrace(const BaseMatrix_Tri<M>& m)
    {
        enum { algo = (
                M::_unit ? 1 :
                M::_unknowndiag ? 3 :
                2 ) };
        return TraceU_Helper<algo,M>(m); 
    }


    //
    // TMV_Text 
    //

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Tri<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Tri< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Tri_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Tri_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
