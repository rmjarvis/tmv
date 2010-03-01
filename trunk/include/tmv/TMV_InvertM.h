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


#ifndef TMV_InvertM_H
#define TMV_InvertM_H

#include "TMV_BaseMatrix.h"

namespace tmv {

    // Defined below:
    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void NoAliasMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void AliasMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2);

    template <class M1, class M2>
    inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void NoAliasMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void AliasMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);

    // Defined in TMV_InvertD.h
    template <class V>
    inline void ElemInvert(BaseVector_Mutable<V>& v);
    template <class V>
    inline void NoAliasElemInvert(BaseVector_Mutable<V>& v);

    // Defined in TMV_InvertU.h
    template <class M>
    inline void InvertSelf(BaseMatrix_Tri_Mutable<M>& m);
    template <class M>
    inline void NoAliasInvertSelf(BaseMatrix_Tri_Mutable<M>& m);

    // Defined in TMV_MultUL.h
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);

    // Defined in TMV_ElemMultVV.h
    template <bool add, int ix, class T, class V1, class V2, class V3>
    inline void ElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <bool add, int ix, class T, class V1, class V2, class V3>
    inline void NoAliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);


    //
    // minv = x * m^-1
    //

    template <class M1>
    static void InvertM_ThrowSingular(const M1& m1)
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Singular Matrix found in Invert\n";
        exit(1);
#else
        throw SingularMatrix<M1>(m1.mat());
#endif
    }

    template <int algo, int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper;

    // algo 0: cs or rs = 0, nothing to do
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<0,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {}
    };

    // algo 1: size = 1x1
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<1,1,1,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            if (m1.cref(0,0) == T1(0)) InvertM_ThrowSingular(m1);
            m2.ref(0,0) = ZProd<false,false>::quot(x , m1.cref(0,0)); 
        }
    };

    // algo 2: size = 2x2
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<2,2,2,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { 
            typedef typename M1::value_type T1;
            T1 det = DetM_Helper<2,2,M1>::call(m1);
            if (det == T1(0)) InvertM_ThrowSingular(m1);
            // Store these in temporaries just in case there are aliases.
            const T1 m1_00 = m1.cref(0,0);
            const T1 m1_01 = m1.cref(0,1);
            const T1 m1_10 = m1.cref(1,0);
            const T1 m1_11 = m1.cref(1,1);
            m2.ref(0,0) = ZProd<false,false>::prod(
                x , ZProd<false,false>::quot(m1_11 , det));
            m2.ref(0,1) = ZProd<false,false>::prod(
                x , ZProd<false,false>::quot(-m1_01 , det));
            m2.ref(1,0) = ZProd<false,false>::prod(
                x , ZProd<false,false>::quot(-m1_10 , det));
            m2.ref(1,1) = ZProd<false,false>::prod(
                x , ZProd<false,false>::quot(m1_00 , det));
        }
    };

    // algo 11: Use Divider // TODO
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<11,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { }
    };

    // algo 12: Calculate LU decomposition on the spot. // TODO
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<12,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { }
    };

    // algo 13: Calculate QR decomposition on the spot. // TODO
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<13,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { }
    };

    // algo 21: m1 is diagonal.  No alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<21,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower || 
                            ShapeTraits<M2::mshape>::upper);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::diag_type m2d = m2.diag();
            m2.setZero();
            NoAliasCopy(m1.diag(),m2d);
            ElemInvert(m2d);
            Scale(x,m2d);
        }
    };

    // algo 22: m1 is diagonal.  Safe for aliases.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<22,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower || 
                            ShapeTraits<M2::mshape>::upper);
            // For alias versions, we could invert the diagonal and then 
            // zero out the non-diagonal portions, but the non-diagonal part
            // is different for each kind of M2, so to make this code
            // completely generic, we use a temporary copy of m1 which we copy
            // into the diagonal after doing m2.setZero().
            // This is also probably more efficient, since zeroing two
            // triangle matrices is slower than zeroing a contiguous 
            // rectangle matrix.  I don't know if the difference is more or
            // less than the copies required.
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M1::const_diag_type::copy_type m1d(m1.size());
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m1d);
            ElemInvert(m1d);
            Scale(x,m1d);
            m2.setZero();
            NoAliasCopy(m1d,m2d);
        }
    };

    // algo 23: m1,m2 are both diagonal.  No alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<23,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m2d);
            ElemInvert(m2d);
            Scale(x,m2d);
        }
    };

    // algo 24: m1,m2 are both diagonal.  With alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<24,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::diag_type m2d = m2.diag();
            AliasCopy(m1.diag(),m2d);
            ElemInvert(m2d);
            Scale(x,m2d);
        }
    };

    // algo 25: m1,m2 are both diagonal, but use temporary anyway.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<25,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M1::const_diag_type::copy_type m1d(m1.size());
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m1d);
            ElemInvert(m1d);
            Scale(x,m1d);
            NoAliasCopy(m1d,m2d);
        }
    };

    // algo 31: m1 is uppertri.  No alias.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<31,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::uppertri_type m2u = m2.upperTri();
            m2.setZero();
            NoAliasCopy(m1,m2u);
            InvertSelf(m2u);
            Scale(x,m2u);
        }
    };

    // algo 32: m1 is uppertri.  With alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<32,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::uppertri_type m2u = m2.upperTri();
            typename M2::lowertri_type::offdiag_type m2l = 
                m2.upperTri().offDiag();
            AliasCopy(m1,m2u);
            InvertSelf(m2u);
            Scale(x,m2u);
            m2l.setZero();
        }
    };

    // algo 33: m1,m2 are both uppertri.  No alias
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<33,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            NoAliasCopy(m1,m2);
            InvertSelf(m2);
            Scale(x,m2);
        }
    };

    // algo 34: m1,m2 are both uppertri.  With alias check
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<34,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            AliasCopy(m1,m2);
            InvertSelf(m2);
            Scale(x,m2);
        }
    };

    // algo 41: m1 is lowertri.  No alias.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<41,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::lowertri_type m2l = m2.lowerTri();
            m2.setZero();
            NoAliasCopy(m1,m2l);
            InvertSelf(m2l);
            Scale(x,m2l);
        }
    };

    // algo 42: m1 is lowertri.  With alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<42,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            typename M2::lowertri_type m2l = m2.lowerTri();
            typename M2::uppertri_type::offdiag_type m2u = 
                m2.upperTri().offDiag();
            AliasCopy(m1,m2l);
            InvertSelf(m2l);
            Scale(x,m2l);
            m2u.setZero();
        }
    };

    // algo 43: m1,m2 are both lowertri.  No alias
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<43,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            NoAliasCopy(m1,m2);
            InvertSelf(m2);
            Scale(x,m2);
        }
    };

    // algo 44: m1,m2 are both lowertri.  With alias check
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<44,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InvertM_ThrowSingular(m1);
            AliasCopy(m1,m2);
            InvertSelf(m2);
            Scale(x,m2);
        }
    };


    // algo -2: No alias check
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-2,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::mshape>::upper;
            const bool lo1 = ShapeTraits<M1::mshape>::lower;
            const bool up2 = ShapeTraits<M2::mshape>::upper;
            const bool lo2 = ShapeTraits<M2::mshape>::lower;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 21 : M2::mdiagstep==1 ? 23 : 25 ) :
                !lo1 ? ( // m1 is uppertri
                    lo2 ? 31 : 33 ) :
                !up1 ? ( // m1 is lowertri
                    up2 ? 41 : 43 ) :
                cs == 2 && rs == 2 ? 2 :
                0;
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<99,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::mshape>::upper;
            const bool lo1 = ShapeTraits<M1::mshape>::lower;
            const bool up2 = ShapeTraits<M2::mshape>::upper;
            const bool lo2 = ShapeTraits<M2::mshape>::lower;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 22 : M2::mdiagstep==1 ? 24 : 25 ) :
                !lo1 ? ( // m1 is uppertri
                    lo2 ? 32 : 34 ) :
                !up1 ? ( // m1 is lowertri
                    up2 ? 42 : 44 ) :
                cs == 2 && rs == 2 ? 2 :
                0;
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-1,cs,rs,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool checkalias = 
                M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
                M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                checkalias ? 99 :
                -2;
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    template <int algo, int ix, class T, class M1, class M2>
    inline void DoMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M2::mrowsize>::same));
        TMVStaticAssert((Sizes<M1::mrowsize,M2::mcolsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::mcolsize,M1::mrowsize>::size;
        const int rs = Sizes<M2::mrowsize,M1::mcolsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        M2v m2v = m2.cView();
        InvertM_Helper<algo,cs,rs,ix,T,M1,M2v>::call(x,m1,m2v);
    }

    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    { DoMakeInverse<-1>(x,m1,m2); }

    template <int ix, class T, class M1, class M2>
    inline void NoAliasMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    { DoMakeInverse<-2>(x,m1,m2); }

    template <int ix, class T, class M1, class M2>
    inline void AliasMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    { DoMakeInverse<99>(x,m1,m2); }



    //
    // mata = (at * a)^-1
    //

    template <class M1>
    static void InverseATA_ThrowSingular(const M1& m1)
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Singular Matrix found in InverseATA\n";
        exit(1);
#else
        throw SingularMatrix<M1>(m1.mat());
#endif
    }

    template <int algo, int cs, int rs, class M1, class M2>
    struct InverseATA_Helper;

    // algo 0: cs or rs = 0, nothing to do
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<0,cs,rs,M1,M2>
    {
        static inline void call(const M1& , M2& )
        {}
    };

    // algo 1: size = 1x1
    template <class M1, class M2>
    struct InverseATA_Helper<1,1,1,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            if (m1.cref(0,0) == T1(0)) InverseATA_ThrowSingular(m1);
            T1 inva00 = ZProd<false,false>::quot(RT(1) , m1.cref(0,0));
            m2.ref(0,0) = ZProd<false,false>::prod(inva00 , inva00);
        }
    };

    // algo 2: size = Nx2
    template <int cs, class M1, class M2>
    struct InverseATA_Helper<2,cs,2,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        { 
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            SmallMatrix<T1,2,2> ata;
            NoAliasMultMM(Scaling<1,RT>(),m1.adjoint(),m1,ata);
            NoAliasMakeInverse(ata,m2);
        }
    };

    // algo 11: Use Divider // TODO
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<11,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        { }
    };

    // algo 12: Calculate LU decomposition on the spot. // TODO
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<12,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        { }
    };

    // algo 13: Calculate QR decomposition on the spot. // TODO
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<13,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        { }
    };

    // algo 21: m1 is diagonal.  No alias check.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<21,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower || 
                            ShapeTraits<M2::mshape>::upper);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::diag_type m2d = m2.diag();
            m2.setZero();
            NoAliasCopy(m1.diag(),m2d);
            ElemInvert(m2d);
            NoAliasElemMultVV<false>(Scaling<1,RT>(),m2d.conjugate(),m2d,m2d);
        }
    };

    // algo 22: m1 is diagonal.  Safe for aliases.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<22,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower || 
                            ShapeTraits<M2::mshape>::upper);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M1::const_diag_type::copy_type m1d(m1.size());
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m1d);
            ElemInvert(m1d);
            NoAliasElemMultVV<false>(Scaling<1,RT>(),m1d.conjugate(),m1d,m1d);
            m2.setZero();
            NoAliasCopy(m1d,m2d);
        }
    };

    // algo 23: m1,m2 are both diagonal.  No alias check.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<23,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m2d);
            ElemInvert(m2d);
            NoAliasElemMultVV<false>(Scaling<1,RT>(),m2d.conjugate(),m2d,m2d);
        }
    };

    // algo 24: m1,m2 are both diagonal.  With alias check.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<24,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::diag_type m2d = m2.diag();
            AliasCopy(m1.diag(),m2d);
            ElemInvert(m2d);
            NoAliasElemMultVV<false>(Scaling<1,RT>(),m2d.conjugate(),m2d,m2d);
        }
    };

    // algo 25: m1,m2 are both diagonal, but use temporary anyway.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<25,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M1::const_diag_type::copy_type m1d(m1.size());
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m1d);
            ElemInvert(m1d);
            NoAliasElemMultVV<false>(Scaling<1,RT>(),m1d.conjugate(),m1d,m1d);
            NoAliasCopy(m1d,m2d);
        }
    };

    // algo 31: m1 is uppertri.  No alias.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<31,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::real_type RT;
            TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::uppertri_type m2u = m2.upperTri();
            NoAliasCopy(m1,m2u);
            InvertSelf(m2u);
            NoAliasMultMM<false>(Scaling<1,RT>(),m2u,m2u.adjoint(),m2);
        }
    };

    // algo 32: m1 is uppertri.  With alias check.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<32,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::real_type RT;
            TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::uppertri_type m2u = m2.upperTri();
            AliasCopy(m1,m2u);
            InvertSelf(m2u);
            NoAliasMultMM<false>(Scaling<1,RT>(),m2u,m2u.adjoint(),m2);
        }
    };

    // algo 41: m1 is lowertri.  No alias.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<41,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::lowertri_type m2l = m2.lowerTri();
            NoAliasCopy(m1,m2l);
            InvertSelf(m2l);
            NoAliasMultMM<false>(Scaling<1,RT>(),m2l,m2l.adjoint(),m2);
        }
    };

    // algo 42: m1 is lowertri.  With alias check.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<42,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M1::mshape>::lower);
            TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
            TMVStaticAssert(ShapeTraits<M2::mshape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular(m1);
            typename M2::lowertri_type m2l = m2.lowerTri();
            AliasCopy(m1,m2l);
            InvertSelf(m2l);
            NoAliasMultMM<false>(Scaling<1,RT>(),m2l,m2l.adjoint(),m2);
        }
    };

    // algo -2: No alias check
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<-2,cs,rs,M1,M2> 
    {
        static inline void call(const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::mshape>::upper;
            const bool lo1 = ShapeTraits<M1::mshape>::lower;
            const bool up2 = ShapeTraits<M2::mshape>::upper;
            const bool lo2 = ShapeTraits<M2::mshape>::lower;
            TMVStaticAssert((up2 && lo2) || (!up1 && !lo1));
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 21 : M2::mdiagstep==1 ? 23 : 25 ) :
                !lo1 ? 31 : // m1 is uppertri
                !up1 ? 41 : // m1 is lowertri
                rs == 2 ? 2 :
                0;
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<99,cs,rs,M1,M2> 
    {
        static inline void call(const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::mshape>::upper;
            const bool lo1 = ShapeTraits<M1::mshape>::lower;
            const bool up2 = ShapeTraits<M2::mshape>::upper;
            const bool lo2 = ShapeTraits<M2::mshape>::lower;
            TMVStaticAssert((up2 && lo2) || (!up1 && !lo1));
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 22 : M2::mdiagstep==1 ? 24 : 25 ) :
                !lo1 ? 32 : // m1 is uppertri
                !up1 ? 42 : // m1 is lowertri
                rs == 2 ? 2 :
                0;
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<-1,cs,rs,M1,M2> 
    {
        static inline void call(const M1& m1, M2& m2)
        {
            const bool checkalias = 
                M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
                M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                checkalias ? 99 :
                -2;
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <int algo, class M1, class M2>
    inline void DoMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M2::mrowsize>::same));
        TMVStaticAssert((Sizes<M1::mrowsize,M2::mcolsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::mcolsize,M1::mrowsize>::size;
        const int rs = Sizes<M2::mrowsize,M1::mcolsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        M2v m2v = m2.cView();
        InverseATA_Helper<algo,cs,rs,M1,M2v>::call(m1,m2v);
    }

    template <class M1, class M2>
    inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { DoMakeInverseATA<-1>(m1,m2); }

    template <class M1, class M2>
    inline void NoAliasMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { DoMakeInverseATA<-2>(m1,m2); }

    template <class M1, class M2>
    inline void AliasMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { DoMakeInverseATA<99>(m1,m2); }

} // namespace tmv

#endif 
