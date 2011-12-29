

#ifndef TMV_DivM_H
#define TMV_DivM_H

#include "TMV_BaseMatrix.h"
#include "TMV_MultXM_Funcs.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_DivMM_Funcs.h"
#include "TMV_DivVM_Funcs.h" 
#include "TMV_Det.h"

#ifdef PRINTALGO_DIVM
#include <iostream>
#endif

namespace tmv {

    //
    // m1 /= m2   or   m1 = m2.inverse() * m1
    // (Most of these routines are the same is m1 is a vector rather than
    //  a matrix, so the same Helper structure is used for both.)
    //

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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 1: s,xs = "<<1<<','<<xs<<std::endl;
#endif
            typedef typename M2::value_type T2;
            if (m2.cref(0,0) == T2(0)) ThrowSingular("1x1 Matrix");
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
                if (det == T2(0)) ThrowSingular("2x2 Matrix");
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
                const int K = xs == TMV_UNKNOWN ? m1.rowsize() : xs;
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 2: s,xs = "<<2<<','<<xs<<std::endl;
#endif
            const int algo2 = xs==1 ? 1 : 2;
            Helper2<algo2,1>::call(m1,m2);
        }
    };

    // algo 3: m2 is small enough to do inverse directly and multiply.
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<3,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 3: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            m1 = m2.inverse().calc() * m1.copy();
        }
    };

    // algo 5: Direct transpose
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<5,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 5: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            typedef typename M2::const_transpose_type M2t;
            M2t m2t = m2.transpose();
            LDivEqM_Helper<-3,s,xs,M1,M2t>::call(m1,m2t);
        }
    };

    // algo 6: Direct transpose, with alias 
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<6,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 6: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            typedef typename M2::const_transpose_type M2t;
            M2t m2t = m2.transpose();
            LDivEqM_Helper<99,s,xs,M1,M2t>::call(m1,m2t);
        }
    };

    // algo 11: Use Divider 
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<11,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 11: s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
            std::cout<<"divIsInPlace = "<<m2.divIsInPlace()<<std::endl;
            std::cout<<"divIsSaved = "<<m2.divIsSaved()<<std::endl;
            std::cout<<"divType = "<<TMV_Text(m2.getDivType())<<std::endl;
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 12: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.lud().solveInPlace(m1);
        } 
    };

    // algo 21: Use Divider (transpose)
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<21,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 21: s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
            std::cout<<"divIsInPlace = "<<m2.divIsInPlace()<<std::endl;
            std::cout<<"divIsSaved = "<<m2.divIsSaved()<<std::endl;
            std::cout<<"divType = "<<TMV_Text(m2.getDivType())<<std::endl;
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 22: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.lud().solveTransposeInPlace(m1);
        } 
    };

    // algo 31: m2 is diag, invert it directly.
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<31,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 31: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            Scaling<1,typename M1::real_type> one;
            typename M1::noalias_type m1na = m1.noAlias();
            MultMM<false>(one,m2.inverse().calc(),m1,m1na);
        } 
    };

    // algo 32: m2 is diag, invert it directly, with alias check.
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<32,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 32: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            Scaling<1,typename M1::real_type> one;
            MultMM<false>(one,m2.inverse().calc(),m1,m1);
        } 
    };

    // algo 33: m1 is a vector
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<33,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 33: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            Scaling<1,typename M1::real_type> one;
            typename M1::noalias_type m1na = m1.noAlias();
            ElemDivVV(one,m1,m2.diag(),m1na);
        } 
    };

    // algo 34: m1 is a vector, with alias check
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<34,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 34: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            Scaling<1,typename M1::real_type> one;
            ElemDivVV(one,m1,m2.diag(),m1);
        } 
    };

    // algo 35: m1 is diagonal
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<35,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 35: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            Scaling<1,typename M1::real_type> one;
            typename M1::diag_type::noalias_type m1d = m1.diag().noAlias();
            ElemDivVV(one,m1.diag(),m2.diag(),m1d);
        } 
    };

    // algo 36: m1 is diagonal, with alias check
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<36,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 36: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            Scaling<1,typename M1::real_type> one;
            typename M1::diag_type m1d = m1.diag();
            ElemDivVV(one,m1.diag(),m2.diag(),m1d);
        } 
    };

    // algo 41: m2 is triangular
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<41,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 41: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            typename M1::noalias_type m1na = m1.noAlias();
            TriLDivEq(m1na,m2);
        } 
    };

    // algo 42: m2 is triangular, with alias check
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<42,s,xs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq algo 42: s,xs = "<<s<<','<<xs<<std::endl;
#endif
            TriLDivEq(m1,m2);
        } 
    };

    // algo 98: Check aliases for transpose solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<98,s,xs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const bool vec1 = ShapeTraits<M1::_shape>::vector;
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec1 ? 34 :
                    (!up1 && !lo1) ? 36 :
                    32 ) : 
                (!lo2 || !up2) ? 6 : // m2 is triangular
                (s != TMV_UNKNOWN && s <= 4) ? 6 :
                M2::_hasdivider ? 21 :
                22;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq transpose, alias:\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    // algo 99: Check aliases for regular solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<99,s,xs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const bool vec1 = ShapeTraits<M1::_shape>::vector;
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec1 ? 34 :
                    (!up1 && !lo1) ? 36 :
                    32 ) : 
                (!lo2 || !up2) ? 42 : // m2 is triangular
                s == 2 ? 2 :
                (s != TMV_UNKNOWN && s <= 4) ? 3 :
                M2::_hasdivider ? 11 :
                12;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq alias:\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    // algo -4: Figure out which algorithm to use for transpose solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<-4,s,xs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const bool vec1 = ShapeTraits<M1::_shape>::vector;
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec1 ? 33 :
                    (!up1 && !lo1) ? 35 :
                    31 ) : 
                (!lo2 || !up2) ? 5 : // m2 is triangular
                (s != TMV_UNKNOWN && s <= 4) ? 5 :
                M2::_hasdivider ? 21 :
                22;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq transpose:\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Figure out which algorithm to use for regular solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<-3,s,xs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const bool vec1 = ShapeTraits<M1::_shape>::vector;
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;

            // Possible algorithms to choose from:
            //
            //  0 = s==0, so nothing to do
            //  1 = s==1: reduces to trivial scalar division
            //  2 = s==2: do division directly
            //  3 = m2 = m1.inverse() * m2
            //  5 = Direct transpose
            //  6 = Direct transpose with alias check
            //
            // 11-12: regular solution
            // 11 = m2 has a divider object.  Use that.
            // 12 = Do LU decomposition on the spot.
            //
            // 21-22: transpose solution
            // 21 = m2 has a divider object.  Use that.
            // 22 = Do LU decomposition on the spot.
            //
            // 31 = m2 is diagonal
            // 32 = m2 is diagonal, with alias check
            // 33 = m2 is diagonal, m1 is a vector
            // 34 = m2 is diagonal, m1 is a vector, with alias check
            // 35 = m1,m2 are both diagonal
            // 36 = m1,m2 are both diagonal, with alias check
            //
            // 41 = m2 is triangular 
            // 42 = m2 is triangular, with alias check

            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec1 ? 33 :
                    (!up1 && !lo1) ? 35 :
                    31 ) : 
                (!lo2 || !up2) ? 41 : // m2 is triangular
                s == 2 ? 2 :
                (s != TMV_UNKNOWN && s <= 4) ? 3 :
                M2::_hasdivider ? 11 :
                12;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDivEq: regular solution\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,xs = "<<s<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: Check for aliases? for regular solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<-2,s,xs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                M1::_checkalias ? 98 :
                -4;
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases? for regular solution
    template <int s, int xs, class M1, class M2>
    struct LDivEqM_Helper<-1,s,xs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                M1::_checkalias ? 99 :
                -3;
            LDivEqM_Helper<algo,s,xs,M1,M2>::call(m1,m2);
        }
    };

    //
    // v1 /= m2
    //
    template <class V1, class M2>
    inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        // m2 shouldn't be a triangle or diagonal matrix, since they 
        // have their own overrides of this function.  Check to make sure...

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

    //
    // m1 /= m2
    //

    template <class M1, class M2>
    inline void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);

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

    //
    // v1 %= m2
    //

    template <class V1, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);

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

    //
    // m1 %= m2
    // -> mt1 /= m2t
    //

    template <class M1, class M2>
    inline void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);

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



    //
    // m3 = m1 / m2, or v3 = v1 / m2
    // 

    template <int algo, int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper;

    // algo 0: Nothing to do
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<0,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const M1& , const M2& , M3& )
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 1: cs,rs,xs = "<<
                1<<','<<1<<','<<xs<<std::endl;
#endif
            typedef typename M2::value_type T2;
            if (m2.cref(0,0) == T2(0)) ThrowSingular("1x1 Matrix");
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
                if (det == T2(0)) ThrowSingular("2x2 Matrix");
                // v1, v3 might be aliased, so store temporaries:
                T3 a = x*(v1.cref(0) * m2.cref(1,1) - 
                          v1.cref(1) * m2.cref(0,1))/det;
                T3 b = x*(v1.cref(1) * m2.cref(0,0) -
                          v1.cref(0) * m2.cref(1,0))/det;
                v3.ref(0) = a;
                v3.ref(1) = b;
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 2: cs,rs,xs = "<<
                2<<','<<2<<','<<xs<<std::endl;
#endif
            const int algo2 = xs==1 ? 1 : 2;
            Helper2<algo2,1>::call(x,m1,m2,m3);
        }
    };

    // algo 3: m2 is small enough to do inverse directly and multiply.
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<3,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 3: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            m3 = x * m2.inverse().calc() * m1;
        }
    };

    // algo 5: Direct transpose
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<5,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 5: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            typedef typename M2::const_transpose_type M2t;
            M2t m2t = m2.transpose();
            LDivM_Helper<-3,cs,rs,xs,ix,T,M1,M2t,M3>::call(x,m1,m2t,m3);
        }
    };

    // algo 6: Direct transpose, with alias
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<6,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 6: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            typedef typename M2::const_transpose_type M2t;
            M2t m2t = m2.transpose();
            LDivM_Helper<99,cs,rs,xs,ix,T,M1,M2t,M3>::call(x,m1,m2t,m3);
        }
    };

    // algo 11: Use Divider
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<11,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 11: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
            std::cout<<"divIsInPlace = "<<m2.divIsInPlace()<<std::endl;
            std::cout<<"divIsSaved = "<<m2.divIsSaved()<<std::endl;
            std::cout<<"divType = "<<TMV_Text(m2.getDivType())<<std::endl;
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
#ifdef PRINTALGO_DIVM
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 13: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.qrd().solve(m1,m3);
            Scale(x,m3);
        }
    };

    // algo 14: Figure out whether to use LU or QR
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<14,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 14: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSquare())
                LDivM_Helper<12,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else
                LDivM_Helper<13,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 21: Use Divider
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<21,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 21: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"divIsSet = "<<m2.divIsSet()<<std::endl;
            std::cout<<"divIsInPlace = "<<m2.divIsInPlace()<<std::endl;
            std::cout<<"divIsSaved = "<<m2.divIsSaved()<<std::endl;
            std::cout<<"divType = "<<TMV_Text(m2.getDivType())<<std::endl;
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
#ifdef PRINTALGO_DIVM
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
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 23: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M2::real_type>::isinteger);
            m2.qrd().solveTranspose(m1,m3);
            Scale(x,m3);
        }
    };

    // algo 24: Figure out whether to use LU or QR
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<24,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 24: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSquare())
                LDivM_Helper<22,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else
                LDivM_Helper<23,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // 31-38: m2 is diagonal
    // algo 31: Convert to m3 = x * m2.inverse() * m1
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<31,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 31: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<false>(x,m2.inverse().calc(),m1,m3na);
        }
    };

    // algo 32: Convert to m3 = x * m2.inverse() * m1, with alias
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<32,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 32: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            MultMM<false>(x,m2.inverse().calc(),m1,m3);
        }
    };

    // algo 33: m1,m3 are vectors, with alias check
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<33,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 33: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            typename M3::noalias_type m3na = m3.noAlias();
            ElemDivVV(x,m1,m2.diag(),m3na);
        }
    };

    // algo 34: m1,m3 are vectors, with alias check
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<34,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 34: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            ElemDivVV(x,m1,m2.diag(),m3);
        }
    };

    // algo 35: All diagonal
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<35,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 35: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            typename M3::diag_type::noalias_type m3d = m3.diag().noAlias();
            ElemDivVV(x,m1.diag(),m2.diag(),m3d);
        }
    };

    // algo 36: All diagonal, with alias check
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<36,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 36: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            typename M3::diag_type m3d = m3.diag();
            ElemDivVV(x,m1.diag(),m2.diag(),m3d);
        }
    };

    // algo 37: m1 is diagonal
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<37,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 37: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("DiagMatrix");
            typedef typename M3::value_type T3;
            const int s = Sizes<cs,rs>::size;
            typedef typename VCopyHelper<T3,s>::type V;
            V v(m2.size());
            ElemDivVV(x,m1.diag(),m2.diag(),v);
            m3.setZero();
            typename M3::diag_type::noalias_type m3d = m3.diag().noAlias();
            Copy(v,m3d);
        }
    };

    // 41-46: m2 is triangular
    // algo 41: Convert to LDivEq
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<41,cs,rs,xs,ix,T,M1,M2,M3>
    {
        template <class M>
        static void copy(
            const Scaling<ix,T>& x, const BaseMatrix<M>& m1, M3& m3)
        { 
            typename M3::noalias_type m3na = m3.noAlias();
            MultXM<false>(x,m1.mat(),m3na); 
        }

        template <class V>
        static void copy(
            const Scaling<ix,T>& x, const BaseVector<V>& m1, M3& m3)
        {
            typename M3::noalias_type m3na = m3.noAlias();
            MultXV<false>(x,m1.vec(),m3na); 
        }

        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 41: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            copy(x,m1,m3);
            typename M3::noalias_type m3na = m3.noAlias();
            TriLDivEq(m3na,m2);
        }
    };

    // algo 42: Convert to LDivEq with alias check
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<42,cs,rs,xs,ix,T,M1,M2,M3>
    {
        template <class M>
        static void copy(
            const Scaling<ix,T>& x, const BaseMatrix<M>& m1, M3& m3)
        { MultXM<false>(x,m1.mat(),m3); }

        template <class V>
        static void copy(
            const Scaling<ix,T>& x, const BaseVector<V>& m1, M3& m3)
        { MultXV<false>(x,m1.vec(),m3); }

        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv algo 42: cs,rs,xs = "<<
                cs<<','<<rs<<','<<xs<<std::endl;
#endif
            typename M3::noalias_type m3na = m3.noAlias();
            if (SameStorage(m2,m3)) {
                // This is really strange, but we need to check.
                typename M2::copy_type m2c = m2;
                copy(x,m1,m3); 
                TriLDivEq(m3na,m2c);
            } else {
                copy(x,m1,m3); 
                TriLDivEq(m3na,m2);
            }
        }
    };

    // algo 98: Check aliases for transpose solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<98,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const bool up3 = ShapeTraits<M3::_shape>::upper;
            const bool lo3 = ShapeTraits<M3::_shape>::lower;
            const bool vec3 = ShapeTraits<M3::_shape>::vector;

            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec3 ? 34 :
                    (!up3 && !lo3) ? 36 : // all diag
                    (!up1 && !lo1) ? 37 : // m1 diag
                    32 ) :
                !up2 || !lo2 ? 6 : // m2 is triangular
                cs == rs && cs != TMV_UNKNOWN && cs <= 4 ? 6 :
                M2::_hasdivider ? 21 :
                cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ? 24 :
                cs == rs ? 22 :
                23;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv: alias, transpose \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    // algo 99: Check aliases for regular solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<99,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const bool up3 = ShapeTraits<M3::_shape>::upper;
            const bool lo3 = ShapeTraits<M3::_shape>::lower;
            const bool vec3 = ShapeTraits<M3::_shape>::vector;

            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec3 ? 34 :
                    (!up3 && !lo3) ? 36 : // all diag
                    (!up1 && !lo1) ? 37 : // m1 diag
                    32 ) :
                (!lo2 || !up2) ? 42 : // m1 is triangular
                cs == 2 && rs == 2 ? 2 :
                cs == rs && cs != TMV_UNKNOWN && cs <= 4 ? 3 :
                M2::_hasdivider ? 11 :
                cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ? 14 :
                cs == rs ? 12 :
                13;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv: alias\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    // algo -4: Figure out which algorithm to use for transpose solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<-4,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const bool up3 = ShapeTraits<M3::_shape>::upper;
            const bool lo3 = ShapeTraits<M3::_shape>::lower;
            const bool vec3 = ShapeTraits<M3::_shape>::vector;

            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec3 ? 33 :
                    (!up3 && !lo3) ? 35 : // all diag
                    (!up1 && !lo1) ? 37 : // m1 diag
                    31 ) :
                !up2 || !lo2 ? 5 : // m2 is triangular
                cs == rs && cs != TMV_UNKNOWN && cs <= 4 ? 5 :
                M2::_hasdivider ? 21 :
                cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ? 24 :
                cs == rs ? 22 :
                23;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv: transpose\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    // algo -3: Figure out which algorithm to use for regular solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<-3,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool up1 = ShapeTraits<M1::_shape>::upper;
            const bool lo1 = ShapeTraits<M1::_shape>::lower;
            const bool up2 = ShapeTraits<M2::_shape>::upper;
            const bool lo2 = ShapeTraits<M2::_shape>::lower;
            const bool up3 = ShapeTraits<M3::_shape>::upper;
            const bool lo3 = ShapeTraits<M3::_shape>::lower;
            const bool vec3 = ShapeTraits<M3::_shape>::vector;

            // Possible algorithms to choose from:
            //
            //  0 = cs==0 or rs==0, so nothing to do
            //  1 = cs==1, rs==1: reduces to trivial scalar division
            //  2 = cs==2, rs==2: do division directly
            //  3 = invert m2 directly and multiply.
            //  4 = invert m2 directly and multiply, with alias check
            //  5 = Direct transpose
            //  6 = Direct transpose with alias check
            //
            // 11-14: regular solution
            // 11 = m2 has a divider object.  Use that.
            // 12 = Do LU decomposition on the spot.
            // 13 = Do QR decomposition on the spot.
            // 14 = Unknown whether square or not, find out.
            //      (This shouldn't usually happen.  In normal usage,
            //       matrices without divider objects know their sizes
            //       at compile time.)
            //
            // 21-24: transpose solution
            // 21 = m2 has a divider object.  Use that.
            // 22 = Do LU decomposition on the spot.
            // 23 = Do QR decomposition on the spot.
            // 24 = Unknown whether square or not, find out.
            //
            // 31-36: m2 is diagonal
            // 31 = Convert to m3 = x * m2.inverse() * m1
            // 32 = Convert to m3 = x * m2.inverse() * m1, with alias check
            // 33 = v1,v3 are vectors, no alias check
            // 34 = v1,v3 are vectors, with alias check
            // 35 = m1,m3 are diagonal, no alias check
            // 36 = m1,m3 are diagonal, with alias check
            // 37 = m1 is diagonal, no alias check
            //
            // 41-42: m2 is triangular
            // 41 = Turn into LDivEq op
            // 42 = Do alias check

            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                !up2 && !lo2 ? ( // m2 is diag
                    vec3 ? 33 :
                    (!up3 && !lo3) ? 35 : // all diag
                    (!up1 && !lo1) ? 37 : // m1 diag
                    31 ) :
                (!lo2 || !up2) ? 41 : // m2 is triangular
                cs == 2 && rs == 2 ? 2 :
                cs == rs && cs != TMV_UNKNOWN && cs <= 4 ? 3 :
                M2::_hasdivider ? 11 :
                cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ? 14 :
                cs == rs ? 12 :
                13;
#ifdef PRINTALGO_DIVM
            std::cout<<"LDiv: regular\n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,xs = "<<cs<<','<<rs<<','<<xs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
#ifdef PRINTALGO_DIVM
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    // algo -2: Check for aliases? for transpose solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<-2,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                M3::_checkalias ? 98 :
                -4;
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    // algo -1: Check for aliases? for regular solution
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct LDivM_Helper<-1,cs,rs,xs,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                cs == 1 && rs == 1 ? 1 :
                M3::_checkalias ? 99 :
                -3;
            LDivM_Helper<algo,cs,rs,xs,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
        }
    };

    //
    // v3 = x * v1 / m2   or   v3 = x * m2.inverse() * v1
    //
    template <int ix, class T, class V1, class M2, class V3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename V3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);

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

    //
    // m3 = x * m1 / m2   or   m3 = x * m2.inverse() * m1
    //
    template <int ix, class T, class M1, class M2, class M3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename M3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);

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

    //
    // v3 = x * v1 % m2   or   v3 = x * v1 * m2.inverse()
    //
    template <int ix, class T, class V1, class M2, class V3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        typedef typename V1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename V3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);

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


    //
    // m3 = x * m1 % m2   or   m3 = x * m1 * m2.inverse()
    // Switch to m3t = x * m2t.inverse * m1t
    //
    template <int ix, class T, class M1, class M2, class M3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typedef typename M1::real_type RT1;
        typedef typename M2::real_type RT2;
        typedef typename M3::real_type RT3;
        TMVStaticAssert(!Traits<RT1>::isinteger);
        TMVStaticAssert(!Traits<RT2>::isinteger);
        TMVStaticAssert(!Traits<RT3>::isinteger);

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


} // namespace tmv

#endif 
