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


#ifndef TMV_DivM_H
#define TMV_DivM_H

#include "TMV_BaseMatrix.h"

#ifdef PRINTALGO_DivM
#include <iostream>
#endif

namespace tmv {

    //
    // m1 /= m2   or   m1 = m2.inverse() * m1
    // (Most of these routines are the same is m1 is a vector rather than
    //  a matrix, so the same Helper structure is used for both.)
    //

    static inline void DivM_ThrowSingular()
    {
#ifdef TMV_NO_THROW
        std::cerr<<"Singular Matrix found\n";
        exit(1);
#else
        throw Singular("Matrix found\n");
#endif
    }

    template <int algo, int s, int xs, class M1, class M2>
    struct LDivEqM_Helper;

    // algo 0: s=0 or xs=0, so nothing to do
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<0,s,xs,M1,M2>
    { static void call(M1& , const M2& ) {} };

    // algo 1: m2 is 1x1
    template <int xs, class M1, class M2>
    struct LDivEqM_Helper<1,1,xs,M1,M2>
    {
        // This one is different for vector and matrix, so break it out here.
        template <class M1x>
        static void call2(BaseVector_Mutable<M1x>& v1, const M2& m2)
        { v1.ref(0) /= m2.cref(0,0); }
        template <class M1x>
        static void call2(BaseMatrix_Mutable<M1x>& m1, const M2& m2)
        { m1.row(0,0,m1.rowsize()) /= m2.cref(0,0); }
        template <class M1x>
        static void call2(BaseMatrix_Rec_Mutable<M1x>& m1, const M2& m2)
        { m1.row(0) /= m2.cref(0,0); }
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 1: s,xs = "<<1<<','<<xs<<std::endl;
#endif
            typedef typename M2::value_type T2;
            if (m2.cref(0,0) == T2(0)) DivM_ThrowSingular();
            call2(m1,m2);
        }
    };

    // algo 2: m2 is 2x2
    template <int xs, class M1, class M2>
    struct LDivEqM_Helper<2,2,xs,M1,M2>
    {
        // Several choices within this selection
        template <int algo2, int dummy>
        struct Helper2;

        template <int dummy>
        struct Helper2<1,dummy>
        {
            // M1 is a vector, so do the whole calculation directly
            template <class V1x>
            static void call(BaseVector_Mutable<V1x>& v1, const M2& m2)
            {
                typedef typename V1x::value_type T1;
                typedef typename M2::value_type T2;
                T2 det = DetM_Helper<2,2,M2>::call(m2);
                if (det == T2(0)) DivM_ThrowSingular();
                T1 a = (v1.cref(0) * m2.cref(1,1) - 
                        v1.cref(1) * m2.cref(0,1))/det;
                T1 b = (v1.cref(1) * m2.cref(0,0) -
                        v1.cref(0) * m2.cref(1,0))/det;
                v1.ref(0) = a;
                v1.ref(1) = b;
            }
            // M1 is a matrix, but only one column. 
            template <class M1x>
            static void call(BaseMatrix_Mutable<M1x>& m1, const M2& m2)
            {
                typedef typename M1x::col_type M1c;
                M1c m1c = m1.col(0);
                LDivEqM_Helper<2,2,1,M1c,M2>::call(m1c,m2);
            }
        };
        // M1 is a matrix, and more than one column.
        template <int dummy>
        struct Helper2<2,dummy>
        {
            static void call(M1& m1, const M2& m2)
            {
                typedef typename M1::value_type T1;
                typedef typename M2::value_type T2;
                SmallMatrix<T2,2,2> m2inv = m2.inverse();
                // TODO: I should write a special in-place mutiply 
                // in the MultMM helper, but right now a statement like
                // m2.transpose() *= m1inv.transpose() results in a temporary.
                // Likewise with a 3x3 matrix.
                // For now, just write a simple loop.
                const int K = xs == UNKNOWN ? m1.rowsize() : xs;
                for(int k=0;k<K;++k) {
                    T1 a = (m1.cref(0,k) * m2inv.cref(0,0) + 
                            m1.cref(1,k) * m2inv.cref(0,1));
                    T1 b = (m1.cref(0,k) * m2inv.cref(1,0) +
                            m1.cref(1,k) * m2inv.cref(1,1));
                    m1.ref(0,k) = a;
                    m1.ref(1,k) = b;
                }
            }
        };
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 2: s,xs = "<<2<<','<<xs<<std::endl;
#endif
            const int algo2 = xs==1 ? 1 : 2;
            Helper2<algo2,1>::call(m1,m2);
        }
    };

    // algo 11: Use Divider 
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<11,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 11: s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
#endif
            m2.setDiv();
            m2.getDiv()->solveInPlace(m1);
            m2.doneDiv();
        }
    };

    // algo 12: Calculate LU decomposition on the spot.
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<12,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 12: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.lud().solveInPlace(m1);
        } 
    };

    // algo 20: Direct transpose
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<20,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 20: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            typedef typename M2::const_transpose_type M2t;
            M2t m2t = m2.transpose();
            LDivEqM_Helper<2,s,xs,M1,M2t>::call(m1,m2t);
        }
    };

    // algo 21: Use Divider (transpose)
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<21,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 21: s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
#endif
            m2.setDiv();
            m2.getDiv()->solveTransposeInPlace(m1);
            m2.doneDiv();
        }
    };

    // algo 22: Calculate LU decomposition on the spot (transpose)
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<22,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDivEq algo 22: s,xs = "<<
                s<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.lud().solveTransposeInPlace(m1);
        } 
    };

    // algo -2: Figure out which algorithm to use for transpose solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<-2,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            // Possible algorithms to choose from:
            //
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar division
            //  2 = s == 2: do division directly
            //
            // 20 = Direct transpose m2
            // 21 = m2 has a divider object.
            // 22 = Do LU decomposition on the spot.

            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                s == 2 ? 20 :
                M2::_hasdivider ? 21 :
                22;
#ifdef PRINTALGO_DivM
            std::cout<<"Inline LDivEq transpose\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Figure out which algorithm to use
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<-1,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            // Possible algorithms to choose from:
            //
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar division
            //  2 = s == 2: do division directly
            //
            // 11 = m2 has a divider object.  Use that.
            // 12 = Do LU decomposition on the spot.

            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                s == 2 ? 2 :
                M2::_hasdivider ? 11 :
                12;
#ifdef PRINTALGO_DivM
            std::cout<<"Inline LDivEq\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    //
    // v1 /= m2
    //
    template <class V1, class M2>
    static inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        // m2 shouldn't be a triangle or diagonal matrix, since they 
        // have their own overrides of this function.  Check to make sure...
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);

        TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<V1::_size,M2::_rowsize>::same));
        TMVAssert(v1.size() == m2.colsize());
        TMVAssert(v1.size() == m2.rowsize());
        const int s1 = Sizes<M2::_colsize,M2::_rowsize>::size;
        const int s = Sizes<V1::_size,s1>::size;
        typedef typename V1::cview_type V1v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        LDivEqM_Helper<-1,s,1,V1v,M2>::call(v1v,m2.mat());
    }
    template <class V1, class M2>
    static inline void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    { LDivEq(v1,m2); }
    template <class V1, class M2>
    static inline void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    { LDivEq(v1,m2); }

    //
    // m1 /= m2
    //
    
    template <class M1, class M2>
    static inline void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);

        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m2.rowsize());
        const int s1 = Sizes<M2::_colsize,M2::_rowsize>::size;
        const int s = Sizes<M1::_colsize,s1>::size;
        const int xs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        LDivEqM_Helper<-1,s,xs,M1v,M2>::call(m1v,m2.mat());
    }
    template <class M1, class M2>
    static inline void NoAliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    { LDivEq(m1,m2); }
    template <class M1, class M2>
    static inline void AliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    { LDivEq(m1,m2); }

    //
    // v1 %= m2
    //
    
    template <class V1, class M2>
    static inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);

        TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<V1::_size,M2::_rowsize>::same));
        TMVAssert(v1.size() == m2.colsize());
        TMVAssert(v1.size() == m2.rowsize());
        const int s1 = Sizes<M2::_colsize,M2::_rowsize>::size;
        const int s = Sizes<V1::_size,s1>::size;
        typedef typename V1::cview_type V1v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        LDivEqM_Helper<-2,s,1,V1v,M2>::call(v1v,m2.mat());
    }
    template <class V1, class M2>
    static inline void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    { RDivEq(v1,m2); }
    template <class V1, class M2>
    static inline void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    { RDivEq(v1,m2); }

    //
    // m1 %= m2
    // -> mt1 /= m2t
    //
    
    template <class M1, class M2>
    static inline void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);

        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int s1 = Sizes<M2::_colsize,M2::_rowsize>::size;
        const int s = Sizes<M1::_rowsize,s1>::size;
        const int xs = M1::_colsize;
        typedef typename M1::transpose_type::cview_type M1t;
        M1t m1t = m1.transpose().cView();
        LDivEqM_Helper<-2,s,xs,M1t,M2>::call(m1t,m2.mat());
    }
    template <class M1, class M2>
    static inline void NoAliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    { RDivEq(m1,m2); }
    template <class M1, class M2>
    static inline void AliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    { RDivEq(m1,m2); }



    //
    // m3 = m1 / m2
    // 
 
    template <int algo, int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper;

    // algo 0: Nothing to do
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<0,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& )
        {}
    };

    // algo 1: m2 is 1x1
    template <int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<1,1,1,xs,ix,T,M1,M2,M3>
    {
        // This one is different for vector and matrix, so break it out here.
        template <class M1x, class M3x>
        static void call2(
            const Scaling<ix,T>& x, const BaseVector<M1x>& v1, 
            const M2& m2, BaseVector_Mutable<M3x>& v3)
        { v3.ref(0) = x * v1.cref(0) / m2.cref(0,0); }
        template <class M1x, class M3x>
        static void call2(
            const Scaling<ix,T>& x, const BaseMatrix_Calc<M1x>& m1,
            const M2& m2, BaseMatrix_Mutable<M3x>& m3)
        {
            m3.row(0,0,m1.rowsize()) =
                x * m1.row(0,0,m1.rowsize()) / m2.cref(0,0); 
        }
        template <class M1x, class M3x>
        static void call2(
            const Scaling<ix,T>& x, const BaseMatrix_Calc<M1x>& m1,
            const M2& m2, BaseMatrix_Rec_Mutable<M3x>& m3)
        { m3.row(0) = x * m1.row(0,0,m1.rowsize()) / m2.cref(0,0); }
        template <class M1x, class M3x>
        static void call2(
            const Scaling<ix,T>& x, const BaseMatrix_Rec<M1x>& m1,
            const M2& m2, BaseMatrix_Rec_Mutable<M3x>& m3)
        { m3.row(0) = x * m1.row(0) / m2.cref(0,0); }
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 1: cs,rs,xs = "<<
                1<<','<<1<<','<<xs<<std::endl;
#endif
            typedef typename M2::value_type T2;
            if (m2.cref(0,0) == T2(0)) DivM_ThrowSingular();
            call2(x,m1,m2,m3);
        }
    };

    // algo 2: m2 is 2x2
    template <int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<2,2,2,xs,ix,T,M1,M2,M3>
    {
        template <int algo2, int dummy>
        struct Helper2;

        template <int dummy>
        struct Helper2<1,dummy>
        {
            // M1,M3 are vectors
            template <class V1x, class V3x>
            static void call(
                const Scaling<ix,T>& x, const BaseVector<V1x>& v1, 
                const M2& m2, BaseVector_Mutable<V3x>& v3)
            {
                typedef typename V3x::value_type T3;
                typedef typename M2::value_type T2;
                T2 det = DetM_Helper<2,2,M2>::call(m2);
                if (det == T2(0)) DivM_ThrowSingular();
                v3.ref(0) = x*(v1.cref(0) * m2.cref(1,1) -
                             v1.cref(1) * m2.cref(0,1))/det;
                v3.ref(1) = x*(v1.cref(1) * m2.cref(0,0) -
                             v1.cref(0) * m2.cref(1,0))/det;
            }
            // M1,m3 are matrices, but only one column. 
            template <class M1x, class M3x>
            static void call(
                const Scaling<ix,T>& x, const BaseMatrix_Calc<M1x>& m1, 
                const M2& m2, BaseMatrix_Mutable<M3x>& m3)
            {
                typedef typename M1x::const_col_type M1c;
                typedef typename M3x::col_type M3c;
                M1c m1c = m1.col(0);
                M3c m3c = m3.col(0);
                LDivM_Helper<2,2,2,1,ix,T,M1c,M2,M3c>::call(x,m1c,m2,m3c);
            }
        };
        // M1 is a matrix, and more than one column.
        template <int dummy>
        struct Helper2<2,dummy>
        {
            static void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                typedef typename M3::value_type T3;
                typedef typename M2::value_type T2;
                SmallMatrix<T2,2,2> m2inv = m2.inverse();
                m3 = x * m2inv * m1;
            }
        };
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 2: cs,rs,xs = "<<
                2<<','<<2<<','<<xs<<std::endl;
#endif
            const int algo2 = xs==1 ? 1 : 2;
            Helper2<algo2,1>::call(x,m1,m2,m3);
        }
    };

    // algo 11: Use Divider
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<11,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 11: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
#endif
            m2.setDiv();
            m2.getDiv()->solve(m1,m3);
            m2.doneDiv();
            Scale(x,m3);
        }
    };

    // algo 12: Calculate LU decomposition on the spot.
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<12,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 12: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.lud().solve(m1,m3);
            Scale(x,m3);
        }
    };

    // algo 13: Calculate QR decomposition on the spot.
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<13,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 13: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
#if 0
            QRD<M2> qrd(m2);
            qrd.solve(m1,m3);
            Scale(x,m3);
#endif
        }
    };

    // algo 14: Figure out whether to use LU or QR
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<14,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 14: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSquare())
                LDivM_Helper<12,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else
                LDivM_Helper<13,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 20: Direct transpose
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<20,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 20: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            typedef typename M2::const_transpose_type M2t;
            M2t m2t = m2.transpose();
            LDivM_Helper<-1,cs,rs,xs,ix,T,M1,M2t,M3>::call(x,m1,m2t,m3);
        }
    };

    // algo 21: Use Divider
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<21,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 21: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
#endif
            m2.setDiv();
            m2.getDiv()->solveTranspose(m1,m3);
            m2.doneDiv();
            Scale(x,m3);
        }
    };

    // algo 22: Calculate LU decomposition on the spot.
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<22,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 22: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.lud().solveTranspose(m1,m3);
            Scale(x,m3);
        }
    };

    // algo 23: Calculate QR decomposition on the spot.
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<23,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 23: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
#if 0
            QRD<M2> qrd(m2);
            qrd.solveTranspose(m1,m3);
            Scale(x,m3);
#endif
        }
    };

    // algo 24: Figure out whether to use LU or QR
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<24,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DivM
            std::cout<<"LDiv algo 24: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSquare())
                LDivM_Helper<22,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else
                LDivM_Helper<23,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -2: Figure out which algorithm to use for transpose solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<-2,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms to choose from:
            //
            //  0 = cs==0 or rs==0, so nothing to do
            //  1 = cs==1, rs==1: reduces to trivial scalar division
            //
            // 20 = Direct transpose
            // 21 = m2 has a divider object.  Use that.
            // 22 = Do LU decomposition on the spot.
            // 23 = Do QR decomposition on the spot.
            // 24 = Unknown whether square or not, find out.

            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                cs == 2 && rs == 2 ? 20 :
                M2::_hasdivider ? 21 :
                cs == UNKNOWN || rs == UNKNOWN ? 24 :
                cs == rs ? 22 :
                23;
#ifdef PRINTALGO_DivM
            std::cout<<"Inline LDiv\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };
   
    // algo -1: Figure out which algorithm to use
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<-1,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms to choose from:
            //
            //  0 = cs==0 or rs==0, so nothing to do
            //  1 = cs==1, rs==1: reduces to trivial scalar division
            //  2 = cs==2, rs==2: do division directly
            //
            // 11 = m2 has a divider object.  Use that.
            // 12 = Do LU decomposition on the spot.
            // 13 = Do QR decomposition on the spot.
            // 14 = Unknown whether square or not, find out.
            //      (This shouldn't usually happen.  In normal usage,
            //       matrices without divider objects know their sizes
            //       at compile time.)

            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                cs == 2 && rs == 2 ? 2 :
                M2::_hasdivider ? 11 :
                cs == UNKNOWN || rs == UNKNOWN ? 14 :
                cs == rs ? 12 :
                13;
#ifdef PRINTALGO_DivM
            std::cout<<"Inline LDiv\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };
   
    //
    // v3 = x * v1 / m2   or   v3 = x * m2.inverse() * v1
    //
    template <int ix, class T, class V1, class M2, class V3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename V3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);

        TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M2::_rowsize>::same));
        TMVAssert(v1.size() == m2.colsize());
        TMVAssert(v3.size() == m2.rowsize());
        const int cs = Sizes<M2::_colsize,V1::_size>::size;
        const int rs = Sizes<M2::_rowsize,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        LDivM_Helper<-1,cs,rs,1,ix,T,V1v,M2,V3v>::call(x,v1v,m2.mat(),v3v);
    }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }

    //
    // m3 = x * m1 / m2   or   m3 = x * m2.inverse() * m1
    //
    template <int ix, class T, class M1, class M2, class M3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename M3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
        TMVStaticAssert(ShapeTraits<M3::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M3::_shape>::lower);

        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m3.colsize() == m2.rowsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        const int cs = Sizes<M2::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M2::_rowsize,M3::_colsize>::size;
        const int xs = Sizes<M3::_rowsize,M1::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        LDivM_Helper<-1,cs,rs,xs,ix,T,M1v,M2,M3v>::call(x,m1v,m2.mat(),m3v);
    }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { LDiv(x,m1,m2,m3); }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { LDiv(x,m1,m2,m3); }

    //
    // v3 = x * v1 % m2   or   v3 = x * v1 * m2.inverse()
    //
    template <int ix, class T, class V1, class M2, class V3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename V3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);

        TMVStaticAssert((Sizes<V1::_size,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<V3::_size,M2::_colsize>::same));
        TMVAssert(v1.size() == m2.rowsize());
        TMVAssert(v3.size() == m2.colsize());
        const int cs = Sizes<M2::_colsize,V3::_size>::size;
        const int rs = Sizes<M2::_rowsize,V1::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        LDivM_Helper<-2,cs,rs,1,ix,T,V1v,M2,V3v>::call(x,v1v,m2.mat(),v3v);
    }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1,m2,v3); }

 
    //
    // m3 = x * m1 % m2   or   m3 = x * m1 * m2.inverse()
    // Switch to m3t = x * m2t.inverse * m1t
    //
    template <int ix, class T, class M1, class M2, class M3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename M3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);
        TMVStaticAssert(ShapeTraits<M2::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M2::_shape>::lower);
        TMVStaticAssert(ShapeTraits<M3::_shape>::upper);
        TMVStaticAssert(ShapeTraits<M3::_shape>::lower);

        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m3.rowsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        const int cs = Sizes<M2::_colsize,M3::_rowsize>::size;
        const int rs = Sizes<M2::_rowsize,M1::_rowsize>::size;
        const int xs = Sizes<M3::_colsize,M1::_colsize>::size;
        typedef typename M1::const_transpose_type::const_cview_type M1t;
        typedef typename M3::transpose_type::cview_type M3t;
        M1t m1t = m1.transpose().cView();
        M3t m3t = m3.transpose().cView();
        LDivM_Helper<-2,cs,rs,xs,ix,T,M1t,M2,M3t>::call(x,m1t,m2.mat(),m3t);
    }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { RDiv(x,m1.mat(),m2,m3); }
    template <int ix, class T, class M1, class M2, class M3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { RDiv(x,m1.mat(),m2,m3); }

 
} // namespace tmv

#endif 
