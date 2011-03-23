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

#ifdef PRINTALGO_InvM
#include <iostream>
#endif

namespace tmv {

    //
    // minv = x * m^-1
    //

    static inline void InvertM_ThrowSingular()
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Singular Matrix found\n";
        exit(1);
#else
        throw Singular("Matrix found\n");
#endif
    }

    template <int algo, int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper;

    // algo 0: cs or rs = 0, nothing to do
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<0,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {}
    };

    // algo 1: size = 1x1
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<1,1,1,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 1: M,N,cs,rs = "<<M<<','<<N<<','<<
                1<<','<<1<<std::endl;
#endif
            typedef typename M1::value_type T1;
            if (m1.cref(0,0) == T1(0)) InvertM_ThrowSingular();
            m2.ref(0,0) = ZProd<false,false>::quot(x , m1.cref(0,0)); 
        }
    };

    // algo 2: size = 2x2
    template <int ix, class T, class M1, class M2>
    struct InvertM_Helper<2,2,2,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 2: M,N,cs,rs = "<<M<<','<<N<<','<<
                2<<','<<2<<std::endl;
#endif
            typedef typename M1::value_type T1;
            T1 det = m1.det();
            if (det == T1(0)) InvertM_ThrowSingular();
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

    // algo 11: Use Divider 
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<11,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 11: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1.divIsSet = "<<m1.divIsSet()<<std::endl;
            std::cout<<"m1.getDivType = "<<TMV_Text(m1.getDivType())<<std::endl;
#endif
            m1.setDiv();
            m1.getDiv()->makeInverse(m2);
            m1.doneDiv();
            Scale(x,m2);
        }
    };

    // algo 12: Calculate LU decomposition on the spot.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<12,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 12: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            // This way is slightly faster than going through m1.lud()
            // since it skips the temporary LU matrix.
            Copy(m1,m2);
            Permutation P(m1.rowsize());
            LU_Decompose(m2,P);
            LU_Inverse(m2,P);
            Scale(x,m2);
        }
    };

    // algo 13: Calculate QR decomposition on the spot.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<13,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 13: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
#if 0
            QRD<M1> qrd(m1);
            qrd.makeInverse(m2);
            Scale(x,m2);
#endif
        }
    };

    // algo 14: Figure out whether to use LU or QR
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<14,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 14: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            if (m1.isSquare())
                InvertM_Helper<12,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
            else
                InvertM_Helper<13,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo 21: m1 is diagonal.  No alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<21,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 21: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower || 
                            ShapeTraits<M2::_shape>::upper);
            typename M2::diag_type m2d = m2.diag();
            m2.setZero();
            NoAliasCopy(m1.diag(),m2d);
            DiagMatrixViewOf(m2d).invertSelf();
            Scale(x,m2d);
        }
    };

    // algo 22: m1 is diagonal.  Safe for aliases.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<22,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 22: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower || 
                            ShapeTraits<M2::_shape>::upper);
            // For alias versions, we could invert the diagonal and then 
            // zero out the non-diagonal portions, but the non-diagonal part
            // is different for each kind of M2, so to make this code
            // completely generic, we use a temporary copy of m1 which we copy
            // into the diagonal after doing m2.setZero().
            // This is also probably more efficient, since zeroing two
            // triangle matrices is slower than zeroing a contiguous 
            // rectangle matrix.  I don't know if the difference is more or
            // less than the copies required.
            typename M1::copy_type m1c = m1;
            m1c.invertSelf();
            typename M1::copy_type::diag_type m1cd = m1c.diag();
            typename M2::diag_type m2d = m2.diag();
            Scale(x,m1cd);
            m2.setZero();
            NoAliasCopy(m1cd,m2d);
        }
    };

    // algo 23: m1,m2 are both diagonal.  No alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<23,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 23: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            typename M2::diag_type m2d = m2.diag();
            NoAliasCopy(m1.diag(),m2d);
            m2.invertSelf();
            Scale(x,m2d);
        }
    };

    // algo 24: m1,m2 are both diagonal.  With alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<24,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 24: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            typename M2::diag_type m2d = m2.diag();
            AliasCopy(m1.diag(),m2d);
            m2.invertSelf();
            Scale(x,m2d);
        }
    };

    // algo 25: m1,m2 are both diagonal, but use temporary anyway.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<25,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 25: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            typename M1::copy_type m1c = m1;
            m1c.invertSelf();
            typename M1::copy_type::diag_type m1cd = m1c.diag();
            typename M2::diag_type m2d = m2.diag();
            Scale(x,m1cd);
            NoAliasCopy(m1cd,m2d);
        }
    };

    // algo 31: m1 is uppertri.  No alias.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<31,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 31: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            typename M2::uppertri_type m2u = m2.upperTri();
            m2.setZero();
            NoAliasCopy(m1,m2u);
            m2u.invertSelf();
            Scale(x,m2u);
        }
    };

    // algo 32: m1 is uppertri.  With alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<32,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int N = m1.rowsize();
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            std::cout<<"InvM algo 32: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            typename M2::uppertri_type m2u = m2.upperTri();
            AliasCopy(m1,m2u);
            m2u.invertSelf();
            Scale(x,m2u);
            if (N > 1) m2.lowerTri().offDiag().setZero();
        }
    };

    // algo 33: m1,m2 are both uppertri.  No alias
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<33,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 33: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            NoAliasCopy(m1,m2);
            m2.invertSelf();
            Scale(x,m2);
        }
    };

    // algo 34: m1,m2 are both uppertri.  With alias check
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<34,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 34: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            AliasCopy(m1,m2);
            m2.invertSelf();
            Scale(x,m2);
        }
    };

    // algo 41: m1 is lowertri.  No alias.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<41,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 41: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            typename M2::lowertri_type m2l = m2.lowerTri();
            m2.setZero();
            NoAliasCopy(m1,m2l);
            m2l.invertSelf();
            Scale(x,m2l);
        }
    };

    // algo 42: m1 is lowertri.  With alias check.
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<42,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int N = m1.rowsize();
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            std::cout<<"InvM algo 42: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            typename M2::lowertri_type m2l = m2.lowerTri();
            AliasCopy(m1,m2l);
            m2l.invertSelf();
            Scale(x,m2l);
            if (N > 1) m2.upperTri().offDiag().setZero();
        }
    };

    // algo 43: m1,m2 are both lowertri.  No alias
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<43,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 43: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            NoAliasCopy(m1,m2);
            m2.invertSelf();
            Scale(x,m2);
        }
    };

    // algo 44: m1,m2 are both lowertri.  With alias check
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<44,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvM algo 44: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            AliasCopy(m1,m2);
            m2.invertSelf();
            Scale(x,m2);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<99,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 22 : M2::_diagstep==1 ? 24 : 25 ) :
                !lo1 ? ( // m1 is uppertri
                    lo2 ? 32 : 34 ) :
                !up1 ? ( // m1 is lowertri
                    up2 ? 42 : 44 ) :
                cs == 2 && rs == 2 ? 2 :
                M1::_hasdivider ? 11 :
                cs == UNKNOWN || rs == UNKNOWN ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"AliasCheck InvertM: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-3,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 21 : M2::_diagstep==1 ? 23 : 25 ) :
                !lo1 ? ( // m1 is uppertri
                    lo2 ? 31 : 33 ) :
                !up1 ? ( // m1 is lowertri
                    up2 ? 41 : 43 ) :
                cs == 2 && rs == 2 ? 2 :
                M1::_hasdivider ? 11 :
                cs == UNKNOWN || rs == UNKNOWN ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"Inline InvertM: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -2: No alias
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-2,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { InvertM_Helper<-3,cs,rs,ix,T,M1,M2>::call(x,m1,m2); }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct InvertM_Helper<-1,cs,rs,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool checkalias = 
                M1::_colsize == UNKNOWN && M1::_rowsize == UNKNOWN &&
                M2::_colsize == UNKNOWN && M2::_rowsize == UNKNOWN;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                checkalias ? 99 :
                -3;
            InvertM_Helper<algo,cs,rs,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    template <int ix, class T, class M1, class M2>
    static inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InvertM_Helper<-1,cs,rs,ix,T,M1,M2v>::call(x,m1.mat(),m2v);
    }
    template <int ix, class T, class M1, class M2>
    static inline void NoAliasMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InvertM_Helper<-3,cs,rs,ix,T,M1,M2v>::call(x,m1.mat(),m2v);
    }
    template <int ix, class T, class M1, class M2>
    static inline void AliasMakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InvertM_Helper<99,cs,rs,ix,T,M1,M2v>::call(x,m1.mat(),m2v);
    }


    //
    // mata = (at * a)^-1
    //

    static inline void InverseATA_ThrowSingular()
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Singular Matrix found\n";
        exit(1);
#else
        throw Singular("Matrix found\n");
#endif
    }

    template <int algo, int cs, int rs, class M1, class M2>
    struct InverseATA_Helper;

    // algo 0: cs or rs = 0, nothing to do
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<0,cs,rs,M1,M2>
    {
        static void call(const M1& , M2& )
        {}
    };

    // algo 1: size = 1x1
    template <class M1, class M2>
    struct InverseATA_Helper<1,1,1,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 1: M,N,cs,rs = "<<M<<','<<N<<','<<
                1<<','<<1<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            if (m1.cref(0,0) == T1(0)) InverseATA_ThrowSingular();
            T1 inva00 = ZProd<false,false>::quot(RT(1) , m1.cref(0,0));
            m2.ref(0,0) = ZProd<false,false>::prod(inva00 , inva00);
        }
    };

    // algo 2: size = Nx2
    template <int cs, class M1, class M2>
    struct InverseATA_Helper<2,cs,2,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 2: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<2<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M1::real_type RT;
            SmallMatrix<T1,2,2> ata;
            NoAliasMultMM<false>(Scaling<1,RT>(),m1.adjoint(),m1,ata);
            NoAliasMakeInverse(Scaling<1,RT>(),ata,m2);
        }
    };

    // algo 11: Use Divider
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<11,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 11: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            m1.setDiv();
            m1.getDiv()->makeInverseATA(m2);
            m1.doneDiv();
        }
    };

    // algo 12: Calculate LU decomposition on the spot. 
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<12,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 12: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            m1.lud().makeInverseATA(m2);
        }
    };

    // algo 13: Calculate QR decomposition on the spot. 
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<13,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 13: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
#if 0
            QRD<M1> qrd(m1);
            qrd.makeInverseATA(m2);
#endif
        }
    };

    // algo 14: Figure out whether to use LU or QR
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<14,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 14: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            if (m1.isSqurae())
                InverseATA_Helper<12,cs,rs,M1,M2>::call(m1,m2);
            else
                InverseATA_Helper<13,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo 21: m1 is diagonal.  No alias check.
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<21,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 21: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower || 
                            ShapeTraits<M2::_shape>::upper);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 22: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower || 
                            ShapeTraits<M2::_shape>::upper);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 23: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 24: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 25: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 31: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 32: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(!ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 41: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
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
        static void call(const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"InvATA algo 42: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::real_type RT;
            TMVStaticAssert(!ShapeTraits<M1::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M1::_shape>::lower);
            TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
            TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
            if (m1.isSingular()) InverseATA_ThrowSingular();
            typename M2::lowertri_type m2l = m2.lowerTri();
            AliasCopy(m1,m2l);
            InvertSelf(m2l);
            NoAliasMultMM<false>(Scaling<1,RT>(),m2l,m2l.adjoint(),m2);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<99,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            TMVStaticAssert((up2 && lo2) || (!up1 && !lo1));
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 22 : M2::_diagstep==1 ? 24 : 25 ) :
                !lo1 ? 32 : // m1 is uppertri
                !up1 ? 42 : // m1 is lowertri
                rs == 2 ? 2 :
                M1::_hasdivider ? 11 :
                cs == UNKNOWN || rs == UNKNOWN ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"AliasCheck InverseATA: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<-3,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            TMVStaticAssert((up2 && lo2) || (!up1 && !lo1));
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up1 && !lo1 ? ( // m1 is diag
                    (up2 || lo2) ? 21 : M2::_diagstep==1 ? 23 : 25 ) :
                !lo1 ? 31 : // m1 is uppertri
                !up1 ? 41 : // m1 is lowertri
                rs == 2 ? 2 :
                M1::_hasdivider ? 11 :
                cs == UNKNOWN || rs == UNKNOWN ? 14 :
                cs == rs ? 12 : 
                13;
#ifdef PRINTALGO_InvM
            const int M = m1.colsize();
            const int N = m1.rowsize();
            std::cout<<"Inline InverseATA: M,N,cs,rs = "<<M<<','<<N<<','<<
                cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: No alias
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<-2,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        { InverseATA_Helper<-3,cs,rs,M1,M2>::call(m1,m2); }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, class M1, class M2>
    struct InverseATA_Helper<-1,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const bool checkalias = 
                M1::_colsize == UNKNOWN && M1::_rowsize == UNKNOWN &&
                M2::_colsize == UNKNOWN && M2::_rowsize == UNKNOWN;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                checkalias ? 99 :
                -3;
            InverseATA_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    static inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InverseATA_Helper<-1,cs,rs,M1,M2v>::call(m1.mat(),m2v);
    }
    template <class M1, class M2>
    static inline void NoAliasMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InverseATA_Helper<-3,cs,rs,M1,M2v>::call(m1.mat(),m2v);
    }
    template <class M1, class M2>
    static inline void AliasMakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVAssert(m1.colsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        const int cs = Sizes<M2::_colsize,M1::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_colsize>::size;
        // Don't make a view for m1, since we want to make sure we keep 
        // a divider object if one is present.
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        InverseATA_Helper<99,cs,rs,M1,M2v>::call(m1.mat(),m2v);
    }

} // namespace tmv

#endif 
