
#ifndef TMV_MultXM_H
#define TMV_MultXM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Scaling.h"
#include "TMV_MultXV.h"
#include "TMV_CopyM.h"
#include "TMV_SwapM.h"

#ifdef PRINTALGO_XM
#include <iostream>
#include "TMV_MatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_MultXM.cpp
    template <class T1, int C1, class T2>
    void InstMultXM(
        const T2 x, const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstAddMultXM(
        const T2 x, const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2);

    template <class T1, int C1, class T2>
    void InstAliasMultXM(
        const T2 x, const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstAliasAddMultXM(
        const T2 x, const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2);

    //
    // Matrix += x * Matrix
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper;

    // algo 1: trivial: ix == 1, !add, so call Copy (with alias check)
    template <ptrdiff_t cs, ptrdiff_t rs, class T, class M1, class M2>
    struct MultXM_Helper<1,cs,rs,false,1,T,M1,M2>
    {
        static TMV_INLINE void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { 
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 1\n";
#endif
            CopyM_Helper<-3,cs,rs,M1,M2>::call(m1,m2); 
        }
    };

    // algo 101: Same as algo 1, but use algo -3 for Copy
    template <ptrdiff_t cs, ptrdiff_t rs, class T, class M1, class M2>
    struct MultXM_Helper<101,cs,rs,false,1,T,M1,M2>
    {
        static TMV_INLINE void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { 
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 101\n";
#endif
            CopyM_Helper<-1,cs,rs,M1,M2>::call(m1,m2); 
        }
    };

    // algo 201: Same as algo 1, but use algo -2 for Copy
    template <ptrdiff_t cs, ptrdiff_t rs, class T, class M1, class M2>
    struct MultXM_Helper<201,cs,rs,false,1,T,M1,M2>
    {
        static void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 201\n";
#endif
            // Need a quick alias check here, since ExactSameStorage
            // is allowed for MultXM, but not Copy
            if (!SameStorage(m1,m2)) {
                CopyM_Helper<-2,cs,rs,M1,M2>::call(m1,m2); 
            } else {
                TMVAssert(ExactSameStorage(m1,m2));
                Maybe<(M2::_conj != int(M1::_conj))>::conjself(m2);
            }
        }
    };

    // algo 2: Linearize to vector version
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<2,cs,rs,add,ix,T,M1,M2>
    {
        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 2\n";
#endif
            typedef typename M1::const_linearview_type M1l;
            typedef typename M2::linearview_type M2l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            const ptrdiff_t prod = IntTraits2<cs,rs>::prod;
            MultXV_Helper<-4,prod,add,ix,T,M1l,M2l>::call(x,m1l,m2l);
        }
    };

    // algo 11: Loop over columns
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<11,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_type M1c;
            typedef typename M2::col_type M2c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const ptrdiff_t step1 = m1.stepj();
            const ptrdiff_t step2 = m2.stepj();
            IT1 it1 = m1.get_col(0).begin().nonConj();
            IT2 it2 = m2.get_col(0).begin();
            for(;N;--N) {
                MultXV_Helper<-4,cs,add,ix,T,M1c,M2c>::call2(M,x,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 12: Loop over rows
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<12,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 12: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename M2::row_type M2r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const ptrdiff_t step1 = m1.stepi();
            const ptrdiff_t step2 = m2.stepi();
            IT1 it1 = m1.get_row(0).begin().nonConj();
            IT2 it2 = m2.get_row(0).begin();
            for(;M;--M) {
                MultXV_Helper<-4,rs,add,ix,T,M1r,M2r>::call2(N,x,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<15,cs,rs,add,ix,T,M1,M2>
    {
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M,J,N/2>::unroll(x,m1,m2);
                Unroller<I,M,J+N/2,N-N/2>::unroll(x,m1,m2);
            }
        };
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J>
        struct Unroller<I,M,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,1>::unroll(x,m1,m2);
                Unroller<I+M/2,M-M/2,J,1>::unroll(x,m1,m2);
            }
        };
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J>
        struct Unroller<I,M,J,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Maybe<add>::add( 
                    m2.ref(I,J), ZProd<false,false>::prod(x,m1.cref(I,J)) ); 
            }
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,0,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };

        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { 
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 15\n";
#endif
            Unroller<0,cs,0,rs>::unroll(x,m1,m2); 
        }
    };

    // algo 16: Fully unroll by rows
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<16,cs,rs,add,ix,T,M1,M2>
    {
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(x,m1,m2);
                Unroller<I+M/2,M-M/2,J,N>::unroll(x,m1,m2);
            }
        };
        template <ptrdiff_t I, ptrdiff_t J, ptrdiff_t N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(x,m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(x,m1,m2);
            }
        };
        template <ptrdiff_t I, ptrdiff_t J, ptrdiff_t N>
        struct Unroller<I,0,J,N>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Maybe<add>::add( 
                    m2.ref(I,J), ZProd<false,false>::prod(x,m1.cref(I,J)) ); 
            }
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,1,J,0>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };

        static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { 
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 16\n";
#endif
            Unroller<0,cs,0,rs>::unroll(x,m1,m2); 
        }
    };

    // algo 30: Unknown sizes, determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<30,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
#ifdef PRINTALGO_XM
            std::cout<<"XM algo 30\n";
#endif
            if (m1.canLinearize() && m2.canLinearize() &&
                m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
                MultXM_Helper<2,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
            else if ((m1.isrm() && m2.isrm()) || 
                     ( !(m1.iscm() && m2.iscm()) &&
                       (m1.colsize() > m1.rowsize()) ) )
                MultXM_Helper<12,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
            else 
                MultXM_Helper<11,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<90,cs,rs,true,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultXM(xx,m1.xView(),m2.xView()); 
        }
    };
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<90,cs,rs,false,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultXM(xx,m1.xView(),m2.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<91,cs,rs,true,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultXM(xx,m1.xView(),m2.xView()); 
        }
    };
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<91,cs,rs,false,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultXM(xx,m1.xView(),m2.xView()); 
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<97,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXM_Helper<-2,cs,rs,add,ix,T,M1c,M2c>::call(
                TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<197,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXM_Helper<99,cs,rs,add,ix,T,M1c,M2c>::call(
                TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<98,cs,rs,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            // More than just these are non-clobbering.
            // TODO: Add in the steps that work.
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXM_Helper<-2,cs,rs,true,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Need a temporary
                typedef typename M1::copy_type M1c;
                M1c m1c = m1.copy();
                MultXM_Helper<-2,cs,rs,true,ix,T,M1c,M2>::call(x,m1c,m2);
            }
        }
    };
    template <ptrdiff_t cs, ptrdiff_t rs, int ix, class T, class M1, class M2>
    struct MultXM_Helper<98,cs,rs,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXM_Helper<-2,cs,rs,false,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Let Copy handle the aliasing
                Copy(m1,m2);
                Scale(x,m2);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<99,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst =
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -4: No branches or copies 
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-4,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type T2;
            const bool canlin = 
                M1::_canlin && M2::_canlin &&
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_colmajor && M2::_colmajor) );
            const int algo = 
                ( ix == 1 && !add ) ? 1 :
                canlin ? 2 :
                TMV_OPT == 0 ? 11 :
                ( cs != Unknown && rs != Unknown ) ? (
                    ( IntTraits2<cs,rs>::prod <= ptrdiff_t(128/sizeof(T2)) ) ? (
                        ( M1::_rowmajor && M2::_rowmajor ) ? 16 : 15 ) :
                    ( M1::_rowmajor && M2::_rowmajor ) ? 12 : 
                    ( M1::_colmajor && M2::_colmajor ) ? 11 :
                    ( cs > rs ) ? 12 : 11 ) :
                ( M1::_rowmajor && M2::_rowmajor ) ? 12 : 
                11;
#ifdef PRINTALGO_XM
            std::cout<<"XM: algo -4\n";
            std::cout<<"add = "<<add<<"  x = "<<ix<<" "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
#endif
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
#ifdef PRINTALGO_XM
            std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-3,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool canlin = 
                M1::_canlin && M2::_canlin &&
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_colmajor && M2::_colmajor) );
            const int algo = 
                ( ix == 1 && !add ) ? 1 :
                canlin ? 2 :
                TMV_OPT >= 2 && ( cs == Unknown || rs == Unknown ) ? 30 :
                -4;
#ifdef PRINTALGO_XM
            std::cout<<"InlineMultXM: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"cs,rs = "<<cs<<','<<rs<<
                " = "<<m2.colsize()<<','<<m2.rowsize()<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-2,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst =
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                ( ix == 1 && !add ) ? 201 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t cs, ptrdiff_t rs, bool add, int ix, class T, class M1, class M2>
    struct MultXM_Helper<-1,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const bool noclobber = MStepHelper<M1,M2>::same;
            const bool checkalias =
                M2::_checkalias && !noclobber;
            const int algo = 
                ( ix == 1 && !add ) ? 101 :
                checkalias ? 99 : 
                -2;
            MultXM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    template <bool add, int ix, class T, class M1, class M2>
    inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<-1,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    inline void InlineMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<-3,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    inline void InlineAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXM_Helper<98,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

} // namespace tmv

#endif 
