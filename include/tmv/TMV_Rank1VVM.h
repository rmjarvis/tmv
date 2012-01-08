

#ifndef TMV_Rank1VVM_H
#define TMV_Rank1VVM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseVector.h"
#include "TMV_Prefetch.h"
#include "TMV_Rank1VVM_Funcs.h"
#include "TMV_MultXV.h"

#ifdef PRINTALGO_R1
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#endif

//
// Vector ^ Vector
//

// Should the code use sizes of v1 and v2 to determine whether to 
// apply the scaling to v1 or v2?
#if TMV_OPT >= 2
#define TMV_R1_SCALE
#endif

// UNROLL is the maximum value of cs*rs for which we fully unroll.
// It seems like unrolling is never faster, so I set this to 0, 
// but I'm leaving the code here in case some improvement in the 
// unrolled code changes that.
#define TMV_R1_UNROLL 0

// MIN_COPY_SIZE is the minimum size to copy a vector if its step != 1.
#define TMV_R1_MIN_COPY_SIZE 4

// PREFETCH is the crossover memory size to start using prefetch commands.
// This is undoubtedly a function of the L1 (and L2?) cache size,
// but 2KBytes is probably not too bad for most machines.
// (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_R1_PREFETCH 2048

// The ratio of the time to copy a vector to the time to scale it.  
// ( Or technically it is 1 + this ratio. )
#define TMV_R1_COPY_SCALE_RATIO 4

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
// It doesn't really seem to matter much either way.
#define TMV_R1_ZeroIX (ix == 0)
//#define TMV_R1_ZeroIX (ix != 1)

namespace tmv {

    // Defined in TMV_Rank1VVM.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstRank1Update(
        const T3 x, const ConstVectorView<T1,C1>& v1, 
        const ConstVectorView<T2,C2>& v2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddRank1Update(
        const T3 x, const ConstVectorView<T1,C1>& v1, 
        const ConstVectorView<T2,C2>& v2, MatrixView<T3> m3);
    
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasRank1Update(
        const T3 x, const ConstVectorView<T1,C1>& v1, 
        const ConstVectorView<T2,C2>& v2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddRank1Update(
        const T3 x, const ConstVectorView<T1,C1>& v1, 
        const ConstVectorView<T2,C2>& v2, MatrixView<T3> m3);

    template <int algo, int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper;

    // algo 0: cs or rs = 0, so nothing to do
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<0,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const V1& , const V2& , M3& ) 
        {} 
    };

    // algo 1: cs == 1, so simplifies to a MultXV function
    template <int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<1,1,rs,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 1: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<1<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::row_type M3r;
            M3r m30 = m3.get_row(0);
            typedef typename Traits2<T,typename V1::value_type>::type PT1;
            MultXV_Helper<-4,rs,add,0,PT1,V2,M3r>::call(x*v1.cref(0),v2,m30); 
        }
    };

    // algo 101: same as 1, but use -1 algo
    template <int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<101,1,rs,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 101: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<1<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::row_type M3r;
            M3r m30 = m3.get_row(0);
            typedef typename Traits2<T,typename V1::value_type>::type PT1;
            MultXV_Helper<-1,rs,add,0,PT1,V2,M3r>::call(x*v1.cref(0),v2,m30); 
        }
    };

    // algo 201: same as 1, but use -2 algo
    template <int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<201,1,rs,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 201: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<1<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::row_type M3r;
            M3r m30 = m3.get_row(0);
            typedef typename Traits2<T,typename V1::value_type>::type PT1;
            MultXV_Helper<-2,rs,add,0,PT1,V2,M3r>::call(x*v1.cref(0),v2,m30); 
        }
    };

    // algo 2: rs == 1, so simplifies to a MultXV function
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<2,cs,1,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            std::cout<<"R1 algo 2: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::col_type M3c;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            M3c m30 = m3.get_col(0);
            MultXV_Helper<-4,cs,add,0,PT2,V1,M3c>::call(x*v2.cref(0),v1,m30); 
        }
    };

    // algo 102: same as 2, but use -1 algo
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<102,cs,1,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            std::cout<<"R1 algo 102: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::col_type M3c;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            M3c m30 = m3.get_col(0);
            MultXV_Helper<-1,cs,add,0,PT2,V1,M3c>::call(x*v2.cref(0),v1,m30); 
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<202,cs,1,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            std::cout<<"R1 algo 202: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::col_type M3c;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            M3c m30 = m3.get_col(0);
            MultXV_Helper<-2,cs,add,0,PT2,V1,M3c>::call(x*v2.cref(0),v1,m30); 
        }
    };

    // algo 3: Transpose
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<3,cs,rs,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 3: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::transpose_type M3t;
            M3t m3t = m3.transpose();
            Rank1VVM_Helper<-2,rs,cs,add,ix,T,V2,V1,M3t>::call(
                x,v2,v1,m3t);
        }
    };

    // algo 11: The basic column major loop
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<11,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int N = rs == Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            std::cout<<"R1 algo 12: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename V2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename M3::col_type M3c;

            for(int j=0;j<N;++j) {
                PT2 v2j = ZProd<false,false>::prod(x , v2.cref(j));
                M3c m3j = m3.get_col(j);
                MultXV_Helper<-4,cs,add,0,PT2,V1,M3c>::call(v2j,v1,m3j);
            }
        }
    };

    // algo 12: Column major loop with iterators
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<12,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 12: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

            typedef typename V1::const_nonconj_type::const_iterator IT1;
            const IT1 X = v1.begin().nonConj();

            typedef typename V2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            PT2 Y0;
            IT2 Y = v2.begin().nonConj();
            const bool c2 = V2::_conj;

            typedef typename M3::col_type M3c;
            typedef typename M3c::iterator IT3;

            const int stepj = m3.stepj();
            IT3 A0 = m3.get_col(0).begin();

            for(int j=N;j;--j) {
                Y0 = ZProd<false,c2>::prod(x , *Y++);
                MultXV_Helper<-4,cs,add,0,PT2,V1,M3c>::call2(
                    M,Scaling<0,PT2>(Y0),X,A0);
                A0.shiftP(stepj);
            }
        }
    };

    // algo 13: column major, 4 columns at a time
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<13,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 13: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::isreal);

            if (M) {

                typedef typename V1::value_type T1;
                typedef typename V1::const_iterator IT1;
                T1 X0;

                typedef typename V2::value_type T2;
                typedef typename V2::const_iterator IT2;
                typedef typename Traits2<T,T2>::type PT2;
                PT2 Y0, Y1, Y2, Y3;

                typedef typename M3::value_type T3;
                typedef typename M3::col_type::iterator IT3;

                int N4 = N>>2; // N4 = N/4
                int Nb = N-(N4<<2); // N4 = N%4
                const int stepj = m3.stepj();
                const int stepj_4 = (stepj<<2) - M;
                // step over 4 columns and back to start
                const int stepj_1 = stepj - M;

                IT3 A0 = m3.get_col(0).begin();
                IT3 A1 = A0; A1.shiftP(stepj);
                IT3 A2 = A1; A2.shiftP(stepj);
                IT3 A3 = A2; A3.shiftP(stepj);
                IT2 Y = v2.begin();
                const IT1 X_begin = v1.begin();
                IT1 X = X_begin;

                const bool dopref = M * sizeof(T3) >= TMV_R1_PREFETCH;

                Prefetch_Read(Y.get());
                Prefetch_MultiRead(X.get());
                Prefetch_Write(A0.get());
                Prefetch_Write(A1.get());
                Prefetch_Write(A2.get());
                Prefetch_Write(A3.get());

                int i;

                if (N4) do {
                    Y0 = x * Y[0]; Y1 = x * Y[1]; Y2 = x * Y[2]; Y3 = x * Y[3];
                    Y += 4;
                    X = X_begin;
                    i=M; do {
                        X0 = *X++;
                        Maybe<add>::add(*A0++ , X0 * Y0);
                        Maybe<add>::add(*A1++ , X0 * Y1);
                        Maybe<add>::add(*A2++ , X0 * Y2);
                        Maybe<add>::add(*A3++ , X0 * Y3);
                    } while (--i);
                    A0.shiftP(stepj_4);
                    A1.shiftP(stepj_4);
                    A2.shiftP(stepj_4);
                    A3.shiftP(stepj_4);
                    if (dopref) {
                        Prefetch_Write(A0.get());
                        Prefetch_Write(A1.get());
                        Prefetch_Write(A2.get());
                        Prefetch_Write(A3.get());
                    }
                } while (--N4);
                if (Nb) do {
                    Y0 = x * *Y++;
                    X = X_begin;
                    i=M; do {
                        X0 = *X++;
                        Maybe<add>::add(*A0++ , X0 * Y0);
                    } while (--i);
                    A0.shiftP(stepj_1);
                    if (dopref) {
                        Prefetch_Write(A0.get());
                    }
                } while (--Nb);
            }
        }
    };

    // algo 14: do all columns at once -- rs <= 4, and must be known
    // rs == 1
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<14,cs,1,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        { Rank1VVM_Helper<2,cs,1,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3); }
    };
    // rs == 2
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<14,cs,2,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            int M = cs == Unknown ? m3.colsize() : cs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 14 N==2: M,N,cs,rs,x = "<<M<<','<<2<<
                ','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::isreal);

            typedef typename V1::value_type T1;
            typedef typename V1::const_iterator IT1;
            T1 X0;

            typedef typename V2::value_type T2;
            typedef typename V2::const_iterator IT2;
            typedef typename Traits2<T,T2>::type PT2;
            const PT2 Y0 = x * v2.cref(0);
            const PT2 Y1 = x * v2.cref(1);

            typedef typename M3::value_type T3;
            typedef typename M3::col_type::iterator IT3;
            const int stepj = m3.stepj();
            IT3 A0 = m3.get_col(0).begin();
            IT3 A1 = A0; A1.shiftP(stepj);

            IT1 X = v1.begin();

            if (M) do {
                X0 = *X++;
                Maybe<add>::add(*A0++ , X0 * Y0);
                Maybe<add>::add(*A1++ , X0 * Y1);
            } while (--M);
        }
    };
    // rs == 3
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<14,cs,3,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            int M = cs == Unknown ? m3.colsize() : cs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 14 N==3: M,N,cs,rs,x = "<<M<<','<<3<<
                ','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::isreal);

            typedef typename V1::value_type T1;
            typedef typename V1::const_iterator IT1;
            T1 X0;

            typedef typename V2::value_type T2;
            typedef typename V2::const_iterator IT2;
            typedef typename Traits2<T,T2>::type PT2;
            const PT2 Y0 = x * v2.cref(0);
            const PT2 Y1 = x * v2.cref(1);
            const PT2 Y2 = x * v2.cref(2);

            typedef typename M3::value_type T3;
            typedef typename M3::col_type::iterator IT3;
            const int stepj = m3.stepj();
            IT3 A0 = m3.get_col(0).begin();
            IT3 A1 = A0; A1.shiftP(stepj);
            IT3 A2 = A1; A2.shiftP(stepj);

            IT1 X = v1.begin();

            if (M) do {
                X0 = *X++;
                Maybe<add>::add(*A0++ , X0 * Y0);
                Maybe<add>::add(*A1++ , X0 * Y1);
                Maybe<add>::add(*A2++ , X0 * Y2);
            } while (--M);
        }
    };
    // rs == 4
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<14,cs,4,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            int M = cs == Unknown ? m3.colsize() : cs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 14 N==4: M,N,cs,rs,x = "<<M<<','<<4<<
                ','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::isreal);

            typedef typename V1::value_type T1;
            typedef typename V1::const_iterator IT1;
            T1 X0;

            typedef typename V2::value_type T2;
            typedef typename V2::const_iterator IT2;
            typedef typename Traits2<T,T2>::type PT2;
            const PT2 Y0 = x * v2.cref(0);
            const PT2 Y1 = x * v2.cref(1);
            const PT2 Y2 = x * v2.cref(2);
            const PT2 Y3 = x * v2.cref(3);

            typedef typename M3::value_type T3;
            typedef typename M3::col_type::iterator IT3;
            const int stepj = m3.stepj();
            IT3 A0 = m3.get_col(0).begin();
            IT3 A1 = A0; A1.shiftP(stepj);
            IT3 A2 = A1; A2.shiftP(stepj);
            IT3 A3 = A2; A3.shiftP(stepj);

            IT1 X = v1.begin();

            if (M) do {
                X0 = *X++;
                Maybe<add>::add(*A0++ , X0 * Y0);
                Maybe<add>::add(*A1++ , X0 * Y1);
                Maybe<add>::add(*A2++ , X0 * Y2);
                Maybe<add>::add(*A3++ , X0 * Y3);
            } while (--M);
        }
    };

    // algo 16: fully unroll by columns
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<16,cs,rs,add,ix,T,V1,V2,M3>
    {
        template <int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
            {
                Unroller<J,N/2>::unroll(x,v1,v2,m3);
                Unroller<J+N/2,N-N/2>::unroll(x,v1,v2,m3);
            }
        };
        template <int J>
        struct Unroller<J,1>
        {
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            template <int I, int M>
            struct Unroller2
            {
                static TMV_INLINE void unroll(
                    const V1& v1, const PT2& v2j, M3& m3)
                {
                    Unroller2<I,M/2>::unroll(v1,v2j,m3);
                    Unroller2<I+M/2,M-M/2>::unroll(v1,v2j,m3);
                }
            };
            template <int I>
            struct Unroller2<I,1>
            {
                static TMV_INLINE void unroll(
                    const V1& v1, const PT2& v2j, M3& m3)
                {
                    Maybe<add>::add(
                        m3.ref(I,J) , 
                        ZProd<false,false>::prod(v1.cref(I) , v2j ) ); 
                }
            };
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
            { Unroller2<0,cs>::unroll(v1,x*v2.cref(J),m3); }
        };
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 16: cs,rs,x = "<<cs<<','<<rs<<
                ','<<T(x)<<std::endl;
#endif
            Unroller<0,rs>::unroll(x,v1,v2,m3); 
        }
    };

    // algo 17: column major, 2 columns at a time, complex m3
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<17,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 17: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::iscomplex);

            if (M) {
                typedef typename V1::value_type T1;
                typedef typename V1::const_nonconj_type::const_iterator IT1;
                T1 X0;

                typedef typename V2::value_type T2;
                typedef typename V2::const_nonconj_type::const_iterator IT2;
                typedef typename Traits2<T,T2>::type PT2;
                PT2 Y0, Y1;

                typedef typename M3::value_type T3;
                typedef typename M3::real_type RT3;
                typedef typename M3::col_type::iterator IT3;
                RT3 rA0, iA0, rA1, iA1;

                int N2 = N>>1; // N2 = N/2
                int Nb = N - (N2<<1); // Nb = N%2
                const int stepj = m3.stepj();
                const int stepj_2 = (stepj<<1) - M;
                // step over 2 columns and back to start

                IT3 A0 = m3.get_col(0).begin();
                IT3 A1 = A0; A1.shiftP(stepj);
                IT2 Y = v2.begin().nonConj();
                const IT1 X_begin = v1.begin().nonConj();
                IT1 X = X_begin;

                const bool c1 = V1::_conj;
                const bool c2 = V2::_conj;

                const bool dopref = M * sizeof(T3) >= TMV_R1_PREFETCH;

                Prefetch_Read(Y.get());
                Prefetch_MultiRead(X.get());
                Prefetch_Write(A0.get());
                Prefetch_Write(A1.get());

                int i;

                if (N2) do {
                    Y0 = ZProd<false,c2>::prod(x,Y[0]);
                    Y1 = ZProd<false,c2>::prod(x,Y[1]);
                    Y += 2;
                    X = X_begin;
                    i=M; do {
                        X0 = *X++;
                        rA0 = Maybe<add>::sum( 
                            real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
                        iA0 = Maybe<add>::sum( 
                            imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
                        rA1 = Maybe<add>::sum( 
                            real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
                        iA1 = Maybe<add>::sum( 
                            imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
                        *A0++ = T3(rA0,iA0);
                        *A1++ = T3(rA1,iA1);
                    } while (--i);
                    A0.shiftP(stepj_2);
                    A1.shiftP(stepj_2);
                    if (dopref) {
                        Prefetch_Write(A0.get());
                        Prefetch_Write(A1.get());
                    }
                } while (--N2);
                if (Nb) {
                    Y0 = ZProd<false,c2>::prod(x,*Y);
                    X = X_begin;
                    i=M; do {
                        X0 = *X++;
                        rA0 = Maybe<add>::sum( 
                            real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
                        iA0 = Maybe<add>::sum( 
                            imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
                        *A0++ = T3(rA0,iA0);
                    } while (--i);
                }
            }
        }
    };

    // algo 18: do all columns at once, complex m3: rs <= 4, and must be known
    // rs == 1
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<18,cs,1,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        { Rank1VVM_Helper<2,cs,1,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3); }
    };
    // rs == 2
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<18,cs,2,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            int M = cs == Unknown ? m3.colsize() : cs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 18 N==2: M,N,cs,rs,x = "<<M<<','<<2<<
                ','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::iscomplex);

            if (M) {
                typedef typename V1::value_type T1;
                typedef typename V1::const_nonconj_type::const_iterator IT1;
                T1 X0;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 Y0 = ZProd<false,false>::prod(x,v2.cref(0));
                const PT2 Y1 = ZProd<false,false>::prod(x,v2.cref(1));

                typedef typename M3::value_type T3;
                typedef typename M3::real_type RT3;
                typedef typename M3::col_type::iterator IT3;
                RT3 rA0, iA0, rA1, iA1;

                const int stepj = m3.stepj();
                IT3 A0 = m3.get_col(0).begin();
                IT3 A1 = A0; A1.shiftP(stepj);
                const IT1 X_begin = v1.begin().nonConj();
                IT1 X = X_begin;

                const bool c1 = V1::_conj;

                do {
                    X0 = *X++;
                    rA0 = Maybe<add>::sum( 
                        real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
                    iA0 = Maybe<add>::sum( 
                        imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
                    rA1 = Maybe<add>::sum( 
                        real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
                    iA1 = Maybe<add>::sum( 
                        imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
                    *A0++ = T3(rA0,iA0);
                    *A1++ = T3(rA1,iA1);
                } while (--M);
            }
        }
    };
    // rs == 3
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<18,cs,3,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            int M = cs == Unknown ? m3.colsize() : cs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 18 N==3: M,N,cs,rs,x = "<<M<<','<<3<<
                ','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::iscomplex);

            if (M) {
                typedef typename V1::value_type T1;
                typedef typename V1::const_nonconj_type::const_iterator IT1;
                T1 X0;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 Y0 = ZProd<false,false>::prod(x,v2.cref(0));
                const PT2 Y1 = ZProd<false,false>::prod(x,v2.cref(1));
                const PT2 Y2 = ZProd<false,false>::prod(x,v2.cref(2));

                typedef typename M3::value_type T3;
                typedef typename M3::real_type RT3;
                typedef typename M3::col_type::iterator IT3;
                RT3 rA0, iA0, rA1, iA1, rA2, iA2;

                const int stepj = m3.stepj();
                IT3 A0 = m3.get_col(0).begin();
                IT3 A1 = A0; A1.shiftP(stepj);
                IT3 A2 = A1; A2.shiftP(stepj);
                const IT1 X_begin = v1.begin().nonConj();
                IT1 X = X_begin;

                const bool c1 = V1::_conj;

                do {
                    X0 = *X++;
                    rA0 = Maybe<add>::sum( 
                        real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
                    iA0 = Maybe<add>::sum( 
                        imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
                    rA1 = Maybe<add>::sum( 
                        real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
                    iA1 = Maybe<add>::sum( 
                        imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
                    rA2 = Maybe<add>::sum( 
                        real(*A2) , ZProd<c1,false>::rprod(X0,Y2) );
                    iA2 = Maybe<add>::sum( 
                        imag(*A2) , ZProd<c1,false>::iprod(X0,Y2) );
                    *A0++ = T3(rA0,iA0);
                    *A1++ = T3(rA1,iA1);
                    *A2++ = T3(rA2,iA2);
                } while (--M);
            }
        }
    };
    // rs == 4
    template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<18,cs,4,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            int M = cs == Unknown ? m3.colsize() : cs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 18 N==4: M,N,cs,rs,x = "<<M<<','<<4<<
                ','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::iscomplex);
            if (M) {

                typedef typename V1::value_type T1;
                typedef typename V1::const_nonconj_type::const_iterator IT1;
                T1 X0;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 Y0 = ZProd<false,false>::prod(x,v2.cref(0));
                const PT2 Y1 = ZProd<false,false>::prod(x,v2.cref(1));
                const PT2 Y2 = ZProd<false,false>::prod(x,v2.cref(2));
                const PT2 Y3 = ZProd<false,false>::prod(x,v2.cref(3));

                typedef typename M3::value_type T3;
                typedef typename M3::real_type RT3;
                typedef typename M3::col_type::iterator IT3;
                RT3 rA0, iA0, rA1, iA1;

                const int stepj = m3.stepj();
                IT3 A0 = m3.get_col(0).begin();
                IT3 A1 = A0; A1.shiftP(stepj);
                IT3 A2 = A1; A2.shiftP(stepj);
                IT3 A3 = A2; A3.shiftP(stepj);
                const IT1 X_begin = v1.begin().nonConj();
                IT1 X = X_begin;

                const bool c1 = V1::_conj;

                do {
                    X0 = *X++;
                    rA0 = Maybe<add>::sum( 
                        real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
                    iA0 = Maybe<add>::sum( 
                        imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
                    rA1 = Maybe<add>::sum( 
                        real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
                    iA1 = Maybe<add>::sum( 
                        imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
                    *A0++ = T3(rA0,iA0);
                    *A1++ = T3(rA1,iA1);
                    rA0 = Maybe<add>::sum( 
                        real(*A2) , ZProd<c1,false>::rprod(X0,Y2) );
                    iA0 = Maybe<add>::sum( 
                        imag(*A2) , ZProd<c1,false>::iprod(X0,Y2) );
                    rA1 = Maybe<add>::sum( 
                        real(*A3) , ZProd<c1,false>::rprod(X0,Y3) );
                    iA1 = Maybe<add>::sum( 
                        imag(*A3) , ZProd<c1,false>::iprod(X0,Y3) );
                    *A2++ = T3(rA0,iA0);
                    *A3++ = T3(rA1,iA1);
                } while (--M);
            }
        }
    };

    // algo 26: fully unroll by rows
    // Neither this, nore 16 above seems to be faster for any cs,rs.
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<26,cs,rs,add,ix,T,V1,V2,M3>
    {
        template <int I, int M>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
            {
                Unroller<I,M/2>::unroll(x,v1,v2,m3);
                Unroller<I+M/2,M-M/2>::unroll(x,v1,v2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            typedef typename Traits2<T,typename V1::value_type>::type PT1;
            template <int J, int N>
            struct Unroller2
            {
                static TMV_INLINE void unroll(
                    const PT1& v1i, const V2& v2, M3& m3)
                {
                    Unroller2<J,N/2>::unroll(v1i,v2,m3);
                    Unroller2<J+N/2,N-N/2>::unroll(v1i,v2,m3);
                }
            };
            template <int J>
            struct Unroller2<J,1>
            {
                static TMV_INLINE void unroll(
                    const PT1& v1i, const V2& v2, M3& m3)
                {
                    Maybe<add>::add(
                        m3.ref(I,J),
                        ZProd<false,false>::prod(v1i , v2.cref(J) ) ); 
                }
            };
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
            { Unroller2<0,rs>::unroll(x*v1.cref(I),v2,m3); }
        };
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 26: cs,rs,x = "<<cs<<','<<rs<<
                ','<<T(x)<<std::endl;
#endif
            Unroller<0,cs>::unroll(x,v1,v2,m3); 
        }
    };

    // TODO: Write SSE Rank1Update functions.
#if 0
    // (This isn't actually an SSE version.  Just a starting point.)
#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<31,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 31: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(M3::isreal);
            TMVStaticAssert(V1::isreal);
            TMVStaticAssert(V2::isreal);
            TMVStaticAssert(M3::_colmajor);

            if (M) {

                typedef typename V1::value_type T1;
                typedef typename V1::const_iterator IT1;
                T1 X0;

                typedef typename V2::value_type T2;
                typedef typename V2::const_iterator IT2;
                typedef typename Traits2<T,T2>::type PT2;
                PT2 Y0, Y1, Y2, Y3;

                typedef typename M3::value_type T3;
                typedef typename M3::col_type::iterator IT3;

                int N4 = N>>2; // N4 = N/4
                int Nb = N-(N4<<2); // N4 = N%4
                const int stepj = m3.stepj();
                const int stepj_4 = (stepj<<2) - M;
                // step over 4 columns and back to start
                const int stepj_1 = stepj - M;

                IT3 A0 = m3.get_col(0).begin();
                IT3 A1 = A0; A1.shiftP(stepj);
                IT3 A2 = A1; A2.shiftP(stepj);
                IT3 A3 = A2; A3.shiftP(stepj);
                IT2 Y = v2.begin();
                const IT1 X_begin = v1.begin();
                IT1 X = X_begin;

                const bool dopref = M * sizeof(T3) >= TMV_R1_PREFETCH;

                Prefetch_Read(Y.get());
                Prefetch_MultiRead(X.get());
                Prefetch_Write(A0.get());
                Prefetch_Write(A1.get());
                Prefetch_Write(A2.get());
                Prefetch_Write(A3.get());

                int i;

                if (N4) do {
                    Y0 = x * Y[0]; Y1 = x * Y[1]; Y2 = x * Y[2]; Y3 = x * Y[3];
                    Y += 4;
                    X = X_begin;
                    i=M; do {
                        X0 = *X++;
                        Maybe<add>::add(*A0++ , X0 * Y0);
                        Maybe<add>::add(*A1++ , X0 * Y1);
                        Maybe<add>::add(*A2++ , X0 * Y2);
                        Maybe<add>::add(*A3++ , X0 * Y3);
                    } while (--i);
                    A0.shiftP(stepj_4);
                    A1.shiftP(stepj_4);
                    A2.shiftP(stepj_4);
                    A3.shiftP(stepj_4);
                    if (dopref) {
                        Prefetch_Write(A0.get());
                        Prefetch_Write(A1.get());
                        Prefetch_Write(A2.get());
                        Prefetch_Write(A3.get());
                    }
                } while (--N4);
                if (Nb) do {
                    Y0 = x * *Y++;
                    X = X_begin;
                    i=M; do {
                        X0 = *X++;
                        Maybe<add>::add(*A0++ , X0 * Y0);
                    } while (--i);
                    A0.shiftP(stepj_1);
                    if (dopref) {
                        Prefetch_Write(A0.get());
                    }
                } while (--Nb);
            }
        }
    };
#endif
#endif


    // algo 71: colmajor, ix==0, so might need to copy v1
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<71,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            TMVStaticAssert(TMV_R1_ZeroIX);
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 71: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
#ifdef TMV_R1_SCALE
            const int MM = cs == Unknown ? m3.colsize() : cs;
            const int NN = rs == Unknown ? m3.rowsize() : rs;
            const int algo2 = M3::iscomplex ? 17 : 13;
            if (NN > TMV_R1_COPY_SCALE_RATIO * MM) {
#endif
                Rank1VVM_Helper<81,cs,rs,add,ix,T,V1,V2,M3>::call(
                    x,v1,v2,m3);
#ifdef TMV_R1_SCALE
            } else 
                Rank1VVM_Helper<algo2,cs,rs,add,ix,T,V1,V2,M3>::call(
                    x,v1,v2,m3);
#endif
        }
    };

    template <int ix, class T, class V> class ProdXV;

    // algo 81: copy x*v1
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<81,cs,rs,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 81: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename M3::noalias_type m3na = m3.noAlias();
            Rank1Update<add>(one,ProdXV<ix,T,V1>(x,v1).calc(),v2,m3na);
        }
    };

    // algo 82: copy x*v2
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<82,cs,rs,add,ix,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int N = rs == Unknown ? m3.rowsize() : rs;
            const int M = cs == Unknown ? m3.colsize() : cs;
            std::cout<<"R1 algo 82: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename M3::noalias_type m3na = m3.noAlias();
            Rank1Update<add>(one,v1,ProdXV<ix,T,V2>(x,v2).calc(),m3na);
        }
    };

    // algo 83: copy both
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<83,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
#ifdef PRINTALGO_R1
            std::cout<<"R1 algo 83: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename M3::noalias_type m3na = m3.noAlias();
            if (M > N) {
                Rank1Update<add>(
                    one,v1.copy(),ProdXV<ix,T,V2>(x,v2).calc(),m3na);
            } else {
                Rank1Update<add>(
                    one,ProdXV<ix,T,V1>(x,v1).calc(),v2.copy(),m3na);
            }
        }
    };
    // If ix == 1, don't need the branch
    template <int cs, int rs, bool add, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<83,cs,rs,add,1,T,V1,V2,M3>
    {
        static inline void call(
            const Scaling<1,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
#ifdef PRINTALGO_R1
            const int M = cs == Unknown ? m3.colsize() : cs;
            const int N = rs == Unknown ? m3.rowsize() : rs;
            std::cout<<"R1 algo 83: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typename M3::noalias_type m3na = m3.noAlias();
            Rank1Update<add>(x,v1.copy(),v2.copy(),m3na);
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<90,cs,rs,true,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddRank1Update(xx,v1.xView(),v2.xView(),m3.xView()); 
        }
    };
    template <int cs, int rs, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<90,cs,rs,false,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstRank1Update(xx,v1.xView(),v2.xView(),m3.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <int cs, int rs, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<91,cs,rs,true,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddRank1Update(xx,v1.xView(),v2.xView(),m3.xView()); 
        }
    };
    template <int cs, int rs, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<91,cs,rs,false,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasRank1Update(xx,v1.xView(),v2.xView(),m3.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<97,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename M3::conjugate_type M3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            M3c m3c = m3.conjugate();
            Rank1VVM_Helper<-2,cs,rs,add,ix,T,V1c,V2c,M3c>::call(
                TMV_CONJ(x),v1c,v2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<197,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename M3::conjugate_type M3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            M3c m3c = m3.conjugate();
            Rank1VVM_Helper<99,cs,rs,add,ix,T,V1c,V2c,M3c>::call(
                TMV_CONJ(x),v1c,v2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<98,cs,rs,add,ix,T,V1,V2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const bool s1 = SameStorage(v1,m3);
            const bool s2 = SameStorage(v2,m3);
            if ( !s1 && !s2 ) {
                // No aliasing
                Rank1VVM_Helper<-2,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
            } else if (s1 && !s2) {
                // Use temporary for v1
                Rank1VVM_Helper<81,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
            } else if (s2 && !s1) {
                // Use temporary for v2
                Rank1VVM_Helper<82,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
            } else {
                // Use temporary for both 
                Rank1VVM_Helper<83,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<99,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst =
                (cs == Unknown || cs > 17) &&
                (rs == Unknown || rs > 17) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                conj ? 197 :
                inst ? 91 :
                98;
            Rank1VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<-4,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename M3::value_type T3;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 1 :
                ( rs == 1 ) ? 2 :
                TMV_OPT == 0 ? ( M3::_rowmajor ? 3 : 11 ) :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ?
                ( M3::_rowmajor ? 3 : 11 ) :
                M3::_colmajor ? ( 
                    ( cs == Unknown || rs == Unknown ) ? (
                        ( M3::iscomplex ? 17 : 13 ) ) :
                    ( rs > cs && cs <= 4 && V2::_step == 1 ) ? 3 :
                    ( rs <= 4 ) ? (M3::iscomplex ? 18 : 14) :
                    M3::iscomplex ? 17 : 13 ) :
                M3::_rowmajor ? (
                    ( cs == Unknown || rs == Unknown ) ? 3 :
                    ( cs > rs && rs <= 4 && V1::_step == 1 ) ? (
                        (M3::iscomplex ? 18 : 14) ) :
                    3 ) :
                11;
            Rank1VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(
                x,v1,v2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<-3,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename M3::value_type T3;
            TMVStaticAssert(!M3::_conj);
            // Possible algorithms to choose from:
            //
            //  0 = cs or rs == 0, so nothing to do
            //  1 = cs == 1: reduces to trivial AddVV function
            //  2 = rs == 1: reduces to trivial AddVV function
            //  3 = transpose
            //
            // 11 = simple for loop
            // 12 = use iterators
            // 13 = 4 columns at a time
            // 14 = rs is known and <= 4
            // 16 = fully unroll by cols
            // 17 = complex m3, 2 columns at a time
            // 18 = complex m3, rs is known and <= 4
            //
            // 26 = fully unroll by rows
            //
            // 71 = x != 1, copy v1 if M < N
            //
            // 81 = copy x*v1
            // 82 = copy x*v2
            // 83 = copy both v1 and v2
            //
            //

            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 1 :
                ( rs == 1 ) ? 2 :
                TMV_OPT == 0 ? ( M3::_rowmajor ? 3 : 11 ) :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ?
                ( M3::_rowmajor ? 3 : 11 ) :
                M3::_colmajor ? (
                    ( cs == Unknown || rs == Unknown ) ? (
                        ( V1::_step != 1 ? 81 : 
                          TMV_R1_ZeroIX ? 71 : 
                          M3::iscomplex ? 17 : 13 ) ) :
                    ( rs > cs && cs <= 4 && V2::_step == 1 ) ? 3 :
                    ( cs > TMV_R1_MIN_COPY_SIZE && V1::_step != 1 ) ? 81 :
                    ( TMV_R1_ZeroIX && 
                      rs > IntTraits2<TMV_R1_COPY_SCALE_RATIO,cs>::prod ) ? 71 :
                    ( rs <= 4 ) ? (M3::iscomplex ? 18 : 14) :
                    M3::iscomplex ? 17 : 13 ) :
                M3::_rowmajor ? (
                    ( cs == Unknown || rs == Unknown ) ? 3 :
                    ( cs > rs && rs <= 4 && V1::_step == 1 ) ? (
                        (M3::iscomplex ? 18 : 14) ) :
                    3 ) :
                // nomajor -- don't do anything fancy
                11;
#ifdef PRINTALGO_R1
            std::cout<<"InlineRank1Update_VVM: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
            //std::cout<<"v1 = "<<v1<<std::endl;
            //std::cout<<"v2 = "<<v2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            Rank1VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(
                x,v1,v2,m3);
#ifdef PRINTALGO_R1
            //std::cout<<"m3 => "<<m3<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<-2,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst =
                (cs == Unknown || cs > 17) &&
                (rs == Unknown || rs > 17) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 201 :
                ( rs == 1 ) ? 202 :
                conj ? 97 :
                inst ? 90 :
                -3;
            Rank1VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
    struct Rank1VVM_Helper<-1,cs,rs,add,ix,T,V1,V2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
        {
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 101 :
                ( rs == 1 ) ? 102 :
                M3::_checkalias ? 99 : 
                -2;
            Rank1VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
        }
    };

    template <bool add, int ix, class T, class V1, class V2, class M3>
    inline void Rank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,V2::_size>::same));
        TMVAssert(m3.colsize() == v1.size());
        TMVAssert(m3.rowsize() == v2.size());
        const int cs = Sizes<M3::_colsize,V1::_size>::size;
        const int rs = Sizes<M3::_rowsize,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        Rank1VVM_Helper<-1,cs,rs,add,ix,T,V1v,V2v,M3v>::call(x,v1v,v2v,m3v);
    }

    template <bool add, int ix, class T, class V1, class V2, class M3>
    inline void InlineRank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,V2::_size>::same));
        TMVAssert(m3.colsize() == v1.size());
        TMVAssert(m3.rowsize() == v2.size());
        const int cs = Sizes<M3::_colsize,V1::_size>::size;
        const int rs = Sizes<M3::_rowsize,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        Rank1VVM_Helper<-3,cs,rs,add,ix,T,V1v,V2v,M3v>::call(x,v1v,v2v,m3v);
    }

    template <bool add, int ix, class T, class V1, class V2, class M3>
    inline void InlineAliasRank1Update(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,V2::_size>::same));
        TMVAssert(m3.colsize() == v1.size());
        TMVAssert(m3.rowsize() == v2.size());
        const int cs = Sizes<M3::_colsize,V1::_size>::size;
        const int rs = Sizes<M3::_rowsize,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        Rank1VVM_Helper<98,cs,rs,add,ix,T,V1v,V2v,M3v>::call(x,v1v,v2v,m3v);
    }

} // namespace tmv

#endif 
