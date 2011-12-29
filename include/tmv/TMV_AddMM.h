

#ifndef TMV_AddMM_H
#define TMV_AddMM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Scaling.h"
#include "TMV_AddVV.h"
#include "TMV_MultXM_Funcs.h"

#ifdef PRINTALGO_AddMM
#include <iostream>
#include "TMV_MatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_AddMM.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMM(
        const T3 x1, const ConstMatrixView<T1,C1>& m1,
        const T3 x2, const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMM(
        const T3 x1, const ConstMatrixView<T1,C1>& m1,
        const T3 x2, const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);


    //
    // Matrix + Matrix
    //

    template <int algo, int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper;

    // algo 0: size == 0, nothing to do
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<0,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& , const M1& , 
            const Scaling<ix2,T2>& , const M2& , M3& )
        {
#ifdef PRINTALGO_AddMM
            std::cout<<"AddMM algo 0\n";
#endif
        }
    };

    // algo 1: Linearize to vector version
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo 1: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::const_linearview_type M1l;
            typedef typename M2::const_linearview_type M2l;
            typedef typename M3::linearview_type M3l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            M3l m3l = m3.linearView();
            InlineAddVV(x1,m1l,x2,m2l,m3l);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<11,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int M = cs == TMV_UNKNOWN ? m2.colsize() : cs;
            int N = rs == TMV_UNKNOWN ? m2.rowsize() : rs;
#ifdef PRINTALGO_AddMM
            std::cout<<"AddMM algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::const_col_type M1c;
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::const_nonconj_type::const_iterator IT2;
            typedef typename M3c::iterator IT3;
            TMVStaticAssert(!M3c::_conj);
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            const int step3 = m3.stepj();
            IT1 it1 = m1.get_col(0).begin().nonConj();
            IT2 it2 = m2.get_col(0).begin().nonConj();
            IT3 it3 = m3.get_col(0).begin();
            for(;N;--N) {
                AddVV_Helper<-4,TMV_UNKNOWN,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
                    M,x1,it1,x2,it2,it3);
                it1.shiftP(step1);
                it2.shiftP(step2);
                it3.shiftP(step3);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<15,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M,J,N/2>::unroll(x1,m1,x2,m2,m3);
                Unroller<I,M,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M/2,J,1>::unroll(x1,m1,x2,m2,m3);
                Unroller<I+M/2,M-M/2,J,1>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo 15: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            Unroller<0,cs,0,rs>::unroll(x1,m1,x2,m2,m3); 
        }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<21,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            int M = cs == TMV_UNKNOWN ? m2.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m2.rowsize() : rs;
#ifdef PRINTALGO_AddMM
            std::cout<<"AddMM algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::const_nonconj_type::const_iterator IT2;
            typedef typename M3r::iterator IT3;
            const int step1 = m1.stepi();
            const int step2 = m2.stepi();
            const int step3 = m3.stepi();
            IT1 it1 = m1.get_row(0).begin().nonConj();
            IT2 it2 = m2.get_row(0).begin().nonConj();
            IT3 it3 = m3.get_row(0).begin();
            for(;M;--M) {
                AddVV_Helper<-4,TMV_UNKNOWN,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
                    N,x1,it1,x2,it2,it3);
                it1.shiftP(step1);
                it2.shiftP(step2);
                it3.shiftP(step3);
            }
        }
    };

    // algo 25: Fully unroll by rows
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<25,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M/2,J,N>::unroll(x1,m1,x2,m2,m3);
                Unroller<I+M/2,M-M/2,J,N>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,1,J,N/2>::unroll(x1,m1,x2,m2,m3);
                Unroller<I,1,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo 25: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            Unroller<0,cs,0,rs>::unroll(x1,m1,x2,m2,m3); 
        }
    };

    // algo 41: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<41,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo 41: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            if (m1.canLinearize() && m2.canLinearize() && m3.canLinearize() &&
                m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj() && 
                m1.stepi() == m3.stepi() && m1.stepj() == m3.stepj())
                AddMM_Helper<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            else if  (
                (m1.isrm() && m2.isrm()) || 
                (m1.isrm() && m3.isrm()) ||
                (m2.isrm() && m3.isrm()) ||
                ( !( (m1.iscm() && m2.iscm()) ||
                     (m1.iscm() && m3.iscm()) ||
                     (m2.iscm() && m3.iscm()) ) && 
                  m1.colsize() > m1.rowsize() ) ) 
                AddMM_Helper<21,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                        x1,m1,x2,m2,m3);
            else 
                AddMM_Helper<11,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<90,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAddMM(xx1,m1.xView(),xx2,m2.xView(),m3.xView());
        }
    };

    // algo 91: Call inst alias
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<91,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAliasAddMM(xx1,m1.xView(),xx2,m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<97,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddMM_Helper<-2,cs,rs,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<197,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddMM_Helper<99,cs,rs,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<98,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_AddMM
            std::cout<<"AddMM Check alias\n";
            std::cout<<"s1,2 = "<<s1<<"  "<<s2<<std::endl;
#endif

            if (!s1 && !s2) {
                // No aliasing 
#ifdef PRINTALGO_AddMM
                std::cout<<"No alias\n";
#endif
                AddMM_Helper<-2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            } else if (!s2) {
#ifdef PRINTALGO_AddMM
                std::cout<<"Do m1 first\n";
#endif
                // Alias with m1 only, do m1 first
                MultXM<false>(x1,m1,m3);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x2,m2,m3na);
            } else if (!s1) {
#ifdef PRINTALGO_AddMM
                std::cout<<"Do m2 first\n";
#endif
                // Alias with m2 only, do m2 first
                MultXM<false>(x2,m2,m3);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x1,m1,m3na);
            } else {
#ifdef PRINTALGO_AddMM
                std::cout<<"Need temporary\n";
#endif
                // Need a temporary
                typename M1::copy_type m1c = m1;
                MultXM<false>(x2,m2,m3);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x1,m1c,m3na);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<99,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                conj ? 197 :
                inst ? 91 :
                98;
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-4,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M3::value_type T3;
            const bool allrm = M1::_rowmajor && M2::_rowmajor && M3::_rowmajor;
            const bool allcm = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool tworm = ( 
                (M1::_rowmajor && M2::_rowmajor) ||
                (M1::_rowmajor && M3::_rowmajor) ||
                (M2::_rowmajor && M3::_rowmajor) );
            const bool twocm = ( 
                (M1::_colmajor && M2::_colmajor) ||
                (M1::_colmajor && M3::_colmajor) ||
                (M2::_colmajor && M3::_colmajor) );
            const bool canlin = 
                M1::_canlin && M2::_canlin && M3::_canlin &&
                ( allrm || allcm );
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                TMV_OPT == 0 ? (allrm ? 21 : 11) :
                canlin ? 1 :
                ( cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T2)) ) ? (
                        tworm ? 25 : 15 ) :
                    tworm ? 21 : 
                    twocm ? 11 :
                    ( cs > rs ) ? 21 : 11 ) :
                tworm ? 21 : 11;
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo -4: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-3,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            const bool allrm = M1::_rowmajor && M2::_rowmajor && M3::_rowmajor;
            const bool allcm = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool canlin = 
                M1::_canlin && M2::_canlin && M3::_canlin &&
                ( allrm || allcm );
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                canlin ? 1 :
                TMV_OPT == 0 ? -4 :
                ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? 41 :
                -4;
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo -3: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                conj ? 97 :
                inst ? 90 :
                -3;
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo -2: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"m1 = "<<m1<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
#ifdef PRINTALGO_AddMM
            std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = 
                M3::_checkalias ? 99 : 
                -2;
#ifdef PRINTALGO_AddMM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"AddMM algo -1: M,N = "<<M<<','<<N<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddMM_Helper<-1,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void InlineAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddMM_Helper<-3,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void InlineAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddMM_Helper<98,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

} // namespace tmv

#endif 
