

#ifndef TMV_MultMV_H
#define TMV_MultMV_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseVector.h"
#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_MultMV_Funcs.h"
#include "TMV_Prefetch.h"

#ifdef PRINTALGO_MV
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#endif

#ifdef XDEBUG_MV
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_AddVV.h"
#include "tmv/TMV_SumVV.h"
#endif

//
// Matrix * Vector
//

// Check for small (<=4) values of cs or rs 
// This leads to huge speed improvements for such matrices at little
// cost to larger matrices, but it causes a big increase in code size
// since the compiler creates code to handle 4 different known values
// plus the generic code for larger values.
// Hence this requires TMV_OPT = 3
#if TMV_OPT >= 3
#define TMV_MV_SMALL
#endif

// When using an algorithm that does 4 columns at a time or 4 rows at 
// a time, this parameter specifies whether to use specific code for 
// the remaining 1, 2, or 3 columns or rows or generic code.  
// It is mildly expensive in terms of code bloat, but not too bad.  
// And it helps all moderately sized matrices, (e.g. M,N ~ 10-100) 
// unless they happen to have M,N%4 = 0.
// These matrices are pretty common, so we only require TMV_OPT = 2 
// for this.
#if TMV_OPT >= 2
#define TMV_MV_CLEANUP
#endif

// When doing something like w = x * m * v, it is often better to 
// apply x to either v before doing the multiplication, or w after 
// doing it.
// The optimal choice depends on the sizes of v and w, so this parameter
// specifies to determine which is better based on the runtime sizes
// of the vectors.  It is a relatively cheap optimization, but it also
// only really helps when the sizes are very different (eg. 100 x 2).
// There is also a mild increase in code bloat.  So we require TMV_OPT = 2.
#if TMV_OPT >= 2
#define TMV_MV_SCALE
#endif

// UNROLL is the maximum nops to unroll.
// We have other cuts that are based purely on timings in algo -1.
// (Search for "unroll = ".)  This parameter supersedes that value
// to turn off unrolling for larger matrices to reduce code bloat.
#if TMV_OPT >= 3
#define TMV_MV_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_MV_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_MV_UNROLL 9
#else
#define TMV_MV_UNROLL 0
#endif

// The minimum size to copy a vector if its step == Unknown.
#define TMV_MV_COPY_SIZE 4

// PREFETCH is the crossover memory size to start using prefetch commands.
// This is undoubtedly a function of the L1 (and L2?) cache size,
// but 2KBytes is probably not too bad for most machines.
// (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_MV_PREFETCH 2048

// The ratio of the time to copy a vector to the time to scale it.
// ( Or technically it is 1 + this ratio. )
#define TMV_MV_COPY_SCALE_RATIO 4

// The minimum value of cs for which unrolling by columns is
// faster than unrolling by rows (when m1 is colmajor).
// (This seems to be fairly independent of rs.)
// (It is also now almost completely irrelevant, since most of these
//  cases use the non-unrolled algorithm now anyway.)
#define TMV_MV_COL_UNROLL 20

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_MV_ZeroIX (ix==0)
//#define TMV_MV_ZeroIX (ix!=1)

namespace tmv {

    // Defined in TMV_MultMV.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMV(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper;

    // algo 0: cs or rs = 0, so nothing to do
    // Correction: if rs = 0, cs != 0 and !add, then we need to do v3.setZero().
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<0,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const M1& , const V2& , V3& v3) 
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? v3.size() : cs;
            std::cout<<"MV algo 0: M,N,cs,rs,x = "<<M<<','<<0<<
                ','<<cs<<','<<rs<<','<<T(0)<<std::endl;
#endif
            Maybe<(!add && rs == 0 && cs != 0)>::zero(v3); 
        }
    };

    // algo 1: cs == 1, so simplifies to a MultVV function
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<1,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 1: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename V2::value_type T2;
            Maybe<add>::add( 
                v3.ref(0) , ZProd<false,false>::prod(
                    x , MultVV_Helper<-4,rs,M1r,V2>::call(m1.get_row(0),v2)) );
        }
    };

    // algo 101: same as 1, but use -1 algo
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<101,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 101: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename V2::value_type T2;
            Maybe<add>::add( 
                v3.ref(0) , ZProd<false,false>::prod(
                    x , MultVV_Helper<-1,rs,M1r,V2>::call(m1.get_row(0),v2)) );
        }
    };

    // algo 201: same as 1, but use -2 algo
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<201,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 201: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename V2::value_type T2;
            Maybe<add>::add( 
                v3.ref(0) , ZProd<false,false>::prod(
                    x , MultVV_Helper<-2,rs,M1r,V2>::call(m1.get_row(0),v2)) );
        }
    };

    // algo 2: rs == 1, so simplifies to a MultXV function
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<2,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            std::cout<<"MV algo 2: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_type M1c;
            MultXV_Helper<-4,cs,add,0,PT2,M1c,V3>::call( 
                Scaling<0,PT2>(ZProd<false,false>::prod(x,v2.cref(0))) , 
                m1.get_col(0), v3 ); 
        }
    };

    // algo 102: same as 2, but use -1 algo
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<102,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            std::cout<<"MV algo 102: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_type M1c;
            MultXV_Helper<-1,cs,add,0,PT2,M1c,V3>::call( 
                Scaling<0,PT2>(ZProd<false,false>::prod(x,v2.cref(0))) ,
                m1.get_col(0), v3 ); 
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<202,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            std::cout<<"MV algo 202: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_type M1c;
            MultXV_Helper<-2,cs,add,0,PT2,M1c,V3>::call( 
                Scaling<0,PT2>(ZProd<false,false>::prod(x,v2.cref(0))) , 
                m1.get_col(0), v3 ); 
        }
    };

    // algo 11: The basic column major loop
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<11,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

            Maybe<!add>::zero(v3);
            if (N) {
                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                typedef typename M1::const_col_type M1c;
                PT2 Xj;

                const bool c2 = V2::_conj;

                typedef typename M1c::const_nonconj_type::const_iterator IT1;
                typedef typename V2::const_nonconj_type::const_iterator IT2;
                typedef typename V3::iterator IT3;
                IT1 A0j = m1.get_col(0).begin().nonConj();
                IT2 X = v2.begin().nonConj();
                const IT3 Y0 = v3.begin();
                const int Astepj = m1.stepj();

                do {
                    if (*X != T2(0)) {
                        Xj = ZProd<false,c2>::prod(x , *X++);
                        // y += Xj * A.col(j);
                        MultXV_Helper<-4,cs,true,0,PT2,M1c,V3>::call2(
                            M,Scaling<0,PT2>(Xj),A0j,Y0);
                    } else {
                        ++X;
                    }
                    A0j.shiftP(Astepj);
                } while (--N);
            }
        }
    };

    // algo 12: column major, 4 columns at a time
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<12,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void loop_4_cols(
            const int M, const int N,
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVAssert(N%4 == 0);
            TMVAssert(N > 0);

            typedef typename M1::value_type T1;
            typedef typename M1::const_col_type::const_iterator IT1;
            T1 A00, A01, A02, A03;

            typedef typename V2::value_type T2;
            typedef typename V2::const_iterator IT2;
            typedef typename Traits2<T,T2>::type PT2;
            PT2 X0, X1, X2, X3;

            typedef typename V3::value_type T3;
            typedef typename V3::iterator IT3;
            T3 Y0;

            int N_4 = N>>2; // N_4 = N/4
            const int stepj = m1.stepj();
            const int stepj_4 = (stepj<<2) - M;
            // step over 4 columns and back to start

            IT1 A0 = m1.get_col(0).begin();
            IT1 A1 = A0; A1.shiftP(stepj);
            IT1 A2 = A1; A2.shiftP(stepj);
            IT1 A3 = A2; A3.shiftP(stepj);
            IT2 X = v2.begin();
            const IT3 Y_begin = v3.begin();
            IT3 Y = Y_begin;

            const bool dopref = M * sizeof(T1) >= TMV_MV_PREFETCH;

            Prefetch_Read(X.get());
            Prefetch_Read(A0.get());
            Prefetch_Read(A1.get());
            Prefetch_Read(A2.get());
            Prefetch_Read(A3.get());
            Prefetch_MultiWrite(Y.get());

            int i;

            TMVAssert(N_4 > 0);
            do {
                X0 = ZProd<false,false>::prod(x , X[0]); 
                X1 = ZProd<false,false>::prod(x , X[1]); 
                X2 = ZProd<false,false>::prod(x , X[2]); 
                X3 = ZProd<false,false>::prod(x , X[3]);
                X += 4;
                Y = Y_begin;

                i=M; do {
                    Y0 = *Y;
                    A00 = *A0++; A01 = *A1++; A02 = *A2++; A03 = *A3++; 
                    Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
                    *Y++ = Y0;
                } while (--i);
                A0 += stepj_4; A1 += stepj_4; A2 += stepj_4; A3 += stepj_4;
                if (dopref) {
                    Prefetch_Read(A0.get());
                    Prefetch_Read(A1.get());
                    Prefetch_Read(A2.get());
                    Prefetch_Read(A3.get());
                }
            } while (--N_4);
        }
        template <bool addx, class M1x, class V2x>
        static void cleanup(
            const int nb,
            const Scaling<ix,T>& x, const M1x& m1, const V2x& v2, V3& v3)
        {
            TMVAssert(nb == 1 || nb == 2 || nb == 3 || nb == 4);
#ifdef TMV_MV_CLEANUP
            switch (nb) {
              case 1 : 
                   MultMV_Helper<2,cs,1,addx,ix,T,M1x,V2x,V3>::call(
                       x,m1,v2,v3);
                   break;
              case 2 :
                   MultMV_Helper<13,cs,2,addx,ix,T,M1x,V2x,V3>::call(
                       x,m1,v2,v3);
                   break;
              case 3 :
                   MultMV_Helper<13,cs,3,addx,ix,T,M1x,V2x,V3>::call(
                       x,m1,v2,v3);
                   break;
              case 4 :
                   MultMV_Helper<13,cs,4,addx,ix,T,M1x,V2x,V3>::call(
                       x,m1,v2,v3);
            }
#else
            MultMV_Helper<11,cs,Unknown,addx,ix,T,M1x,V2x,V3>::call(
                x,m1,v2,v3);
#endif
        }
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 12: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);
            typedef typename M1::const_colrange_type M1c;
            typedef typename V2::const_subvector_type V2s;
            if (M) {
                if (N > 4) {
                    const int na = ((N>>2)<<2); 
                    const int nb = N-na;
                    Maybe<!add>::zero(v3);
                    loop_4_cols(M,na,x,m1,v2,v3);
                    if (nb) {
                        M1c m1b = m1.cColRange(na,N);
                        V2s v2b = v2.cSubVector(na,N);
                        cleanup<true>(nb,x,m1b,v2b,v3);
                    }
                } else if (N) {
                    cleanup<add>(N,x,m1,v2,v3);
                } else {
                    Maybe<!add>::zero(v3);
                }
            }
        }
    };

    // algo 13: do all columns at once -- rs <= 4, and must be known
    // rs == 2
    template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<13,cs,2,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int M = cs == Unknown ? m1.colsize() : cs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 13 N==2: M,N,cs,rs,x = "<<M<<','<<2<<
                ','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_col_type::const_iterator IT1;
                T1 A00, A01;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 X0 = ZProd<false,false>::prod(x , v2.cref(0));
                const PT2 X1 = ZProd<false,false>::prod(x , v2.cref(1));

                const int stepj = m1.stepj();
                IT1 A0 = m1.get_col(0).begin();
                IT1 A1 = A0; A1.shiftP(stepj);

                typedef typename V3::value_type T3;
                typedef typename V3::iterator IT3;
                T3 Y0;
                IT3 Y = v3.begin();

                do {
                    A00 = *A0++; A01 = *A1++;
                    Maybe<add>::set(Y0,*Y);
                    Maybe<add>::add(
                        Y0 , 
                        ZProd<false,false>::prod(A00,X0) +
                        ZProd<false,false>::prod(A01,X1) );
                    *Y++ = Y0;
                } while (--M);
            }
        }
    };
    // rs == 3
    template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<13,cs,3,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int M = cs == Unknown ? m1.colsize() : cs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 13 N==3: M,N,cs,rs,x = "<<M<<','<<3<<
                ','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_col_type::const_iterator IT1;
                T1 A00, A01, A02;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 X0 = ZProd<false,false>::prod(x , v2.cref(0));
                const PT2 X1 = ZProd<false,false>::prod(x , v2.cref(1));
                const PT2 X2 = ZProd<false,false>::prod(x , v2.cref(2));

                const int stepj = m1.stepj();
                IT1 A0 = m1.get_col(0).begin();
                IT1 A1 = A0; A1.shiftP(stepj);
                IT1 A2 = A1; A2.shiftP(stepj);

                typedef typename V3::value_type T3;
                typedef typename V3::iterator IT3;
                T3 Y0;
                IT3 Y = v3.begin();

                do {
                    A00 = *A0++; A01 = *A1++; A02 = *A2++;
                    Maybe<add>::set(Y0,*Y);
                    Maybe<add>::add(
                        Y0 , 
                        ZProd<false,false>::prod(A00,X0) +
                        ZProd<false,false>::prod(A01,X1) +
                        ZProd<false,false>::prod(A02,X2) );
                    *Y++ = Y0;
                } while (--M);
            }
        }
    };
    // rs == 4
    template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<13,cs,4,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int M = cs == Unknown ? m1.colsize() : cs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 13 N==4: M,N,cs,rs,x = "<<M<<','<<4<<
                ','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_col_type::const_iterator IT1;
                T1 A00, A01, A02, A03;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 X0 = ZProd<false,false>::prod(x , v2.cref(0));
                const PT2 X1 = ZProd<false,false>::prod(x , v2.cref(1));
                const PT2 X2 = ZProd<false,false>::prod(x , v2.cref(2));
                const PT2 X3 = ZProd<false,false>::prod(x , v2.cref(3));

                const int stepj = m1.stepj();
                IT1 A0 = m1.get_col(0).begin();
                IT1 A1 = A0; A1.shiftP(stepj);
                IT1 A2 = A1; A2.shiftP(stepj);
                IT1 A3 = A2; A3.shiftP(stepj);

                typedef typename V3::value_type T3;
                typedef typename V3::iterator IT3;
                T3 Y0;
                IT3 Y = v3.begin();

                do {
                    A00 = *A0++; A01 = *A1++; A02 = *A2++; A03 = *A3++;
                    Maybe<add>::set(Y0,*Y);
                    Maybe<add>::add(
                        Y0 , 
                        ZProd<false,false>::prod(A00,X0) +
                        ZProd<false,false>::prod(A01,X1) +
                        ZProd<false,false>::prod(A02,X2) +
                        ZProd<false,false>::prod(A03,X3) );
                    *Y++ = Y0;
                } while (--M);
            }
        }
    };

    // algo 15: fully unroll by cols, apply x to v3
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<15,cs,rs,add,ix,T,M1,V2,V3>
    {
        typedef typename V2::value_type T2;
        typedef typename V3::value_type T3;
        template <int I, int N, bool first>
        struct Unroller
        {
            static inline void loadx(const V2& v2, T2* X)
            {
                Unroller<I,N/2,first>::loadx(v2,X);
                Unroller<I+N/2,N-N/2,first>::loadx(v2,X);
            }
            static inline void sety(
                const Scaling<ix,T>& x, const T3* Y, V3& v3)
            {
                Unroller<I,N/2,first>::sety(x,Y,v3);
                Unroller<I+N/2,N-N/2,first>::sety(x,Y,v3);
            }
            static inline void calcy(const M1& m1, const T2* X, T3* Y)
            {
                Unroller<I,N/2,first>::calcy(m1,X,Y);
                Unroller<I+N/2,N-N/2,first>::calcy(m1,X,Y);
            }
            static inline void calcyj(
                const int j, const M1& m1, const T2& Xj, T3* Y)
            {
                Unroller<I,N/2,first>::calcyj(j,m1,Xj,Y);
                Unroller<I+N/2,N-N/2,first>::calcyj(j,m1,Xj,Y);
            }
        };
        template <int I, bool first>
        struct Unroller<I,1,first>
        {
            static inline void loadx(const V2& v2, T2* X)
            { X[I] = v2.cref(I); }
            static inline void sety(
                const Scaling<ix,T>& x, const T3* Y, V3& v3)
            { Maybe<add>::add(v3.ref(I) , ZProd<false,false>::prod(x,Y[I])); }
            static inline void calcy(const M1& m1, const T2* X, T3* Y)
                // Note: "I" is really j here...
            { Unroller<0,cs,I==0>::calcyj(I,m1,X[I],Y); }
            static inline void calcyj(
                const int j, const M1& m1, const T2& Xj, T3* Y)
            {
                Maybe<!first>::add(
                    Y[I] , ZProd<false,false>::prod(m1.cref(I,j),Xj)); 
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 15: cs,rs,x = "<<cs<<','<<rs<<
                ','<<T(x)<<std::endl;
#endif
            T2 X[rs];
            T3 Y[cs];
            Unroller<0,rs,false>::loadx(v2,X);
            Unroller<0,rs,false>::calcy(m1,X,Y); 
            Unroller<0,cs,false>::sety(x,Y,v3);
        }
    };

    // algo 16: fully unroll by columnss, apply x to v2
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<16,cs,rs,add,ix,T,M1,V2,V3>
    {
        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename V3::value_type T3;
        template <int I, int N, bool first>
        struct Unroller
        {
            static inline void loadx(
                const Scaling<ix,T>& x, const V2& v2, PT2* X)
            {
                Unroller<I,N/2,first>::loadx(x,v2,X);
                Unroller<I+N/2,N-N/2,first>::loadx(x,v2,X);
            }
            static inline void sety(const T3* Y, V3& v3)
            {
                Unroller<I,N/2,first>::sety(Y,v3);
                Unroller<I+N/2,N-N/2,first>::sety(Y,v3);
            }
            static inline void calcy(const M1& m1, const PT2* X, T3* Y)
            {
                Unroller<I,N/2,first>::calcy(m1,X,Y);
                Unroller<I+N/2,N-N/2,first>::calcy(m1,X,Y);
            }
            static inline void calcyj(
                const int j, const M1& m1, const PT2& Xj, T3* Y)
            {
                Unroller<I,N/2,first>::calcyj(j,m1,Xj,Y);
                Unroller<I+N/2,N-N/2,first>::calcyj(j,m1,Xj,Y);
            }
        };
        template <int I, bool first>
        struct Unroller<I,1,first>
        {
            static inline void loadx(
                const Scaling<ix,T>& x, const V2& v2, PT2* X)
            { X[I] = ZProd<false,false>::prod(x , v2.cref(I)); }
            static inline void sety(const T3* Y, V3& v3)
            { Maybe<add>::add(v3.ref(I) , Y[I]); }
            static inline void calcy(const M1& m1, const PT2* X, T3* Y)
                // Note: "I" is really j here...
            { Unroller<0,cs,I==0>::calcyj(I,m1,X[I],Y); }
            static inline void calcyj(
                const int j, const M1& m1, const PT2& Xj, T3* Y)
            {
                Maybe<!first>::add(
                    Y[I] , ZProd<false,false>::prod(m1.cref(I,j),Xj)); 
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 16: cs,rs,x = "<<cs<<','<<rs<<
                ','<<T(x)<<std::endl;
#endif
            PT2 X[rs];
            T3 Y[cs];
            Unroller<0,rs,false>::loadx(x,v2,X);
            Unroller<0,rs,false>::calcy(m1,X,Y); 
            Unroller<0,cs,false>::sety(Y,v3);
        }
    };

    // algo 17: column major, 2 columns at a time, complex v3
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<17,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void loop_2_cols(
            const int M, const int N,
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVAssert(N%2 == 0);
            TMVAssert(N > 0);

            typedef typename M1::value_type T1;
            typedef typename M1::const_col_type::const_nonconj_type M1c;
            typedef typename M1c::const_iterator IT1;
            T1 A00, A01;

            typedef typename V2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            PT2 X0, X1;

            typedef typename V3::value_type T3;
            typedef typename V3::iterator IT3;
            T3 Y0;

            int N_2 = N>>1; // N_2 = N/2
            const int stepj = m1.stepj();
            const int stepj_2 = (stepj<<1) - M;
            // step over 2 columns and back to start

            IT1 A0 = m1.get_col(0).begin().nonConj();
            IT1 A1 = A0; A1.shiftP(stepj);
            IT2 X = v2.begin().nonConj();
            const IT3 Y_begin = v3.begin();
            IT3 Y = Y_begin;

            const bool c1 = M1::_conj;
            const bool c2 = V2::_conj;

            const bool dopref = M * sizeof(T1) >= TMV_MV_PREFETCH;

            Prefetch_Read(X.get());
            Prefetch_Read(A0.get());
            Prefetch_Read(A1.get());
            Prefetch_MultiWrite(Y.get());

            int i;

            TMVAssert(N_2 > 0);
            do {
                X0 = ZProd<false,c2>::prod(x,X[0]);
                X1 = ZProd<false,c2>::prod(x,X[1]);
                X += 2;
                Y = Y_begin;

                i=M; do {
                    A00 = *A0++; A01 = *A1++;
                    Y0 = *Y;
                    Y0 += ZProd<c1,false>::prod(A00,X0);
                    Y0 += ZProd<c1,false>::prod(A01,X1);
                    *Y++ = Y0;
                } while (--i);
                A0 += stepj_2; A1 += stepj_2;
                if (dopref) {
                    Prefetch_Read(A0.get());
                    Prefetch_Read(A1.get());
                }
            } while (--N_2);
        }
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 17: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);
            if (M) {
                const int na = ((N>>1)<<1);
                const int nb = N-na;
                Maybe<!add>::zero(v3);
                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                typedef typename M1::const_col_type M1c;
                if (na) loop_2_cols(M,na,x,m1,v2,v3);
                if (nb) {
                    MultXV_Helper<-4,cs,true,0,PT2,M1c,V3>::call(
                        Scaling<0,PT2>(ZProd<false,false>::prod(x,v2.cref(na))),
                        m1.get_col(na),v3); 
                }
            }
        }
    };

    // algo 18: do all columns at once, complex v3: rs <= 4, and must be known
    // rs == 2
    template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<18,cs,2,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(V3::iscomplex);
            int M = cs == Unknown ? m1.colsize() : cs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 18 N==2: M,N,cs,rs,x = "<<M<<','<<2<<
                ','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_col_type::const_nonconj_type M1c;
                typedef typename M1c::const_iterator IT1;
                T1 A00, A01;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 X0 = ZProd<false,false>::prod(x,v2.cref(0));
                const PT2 X1 = ZProd<false,false>::prod(x,v2.cref(1));

                typedef typename V3::value_type T3;
                typedef typename V3::iterator IT3;
                T3 Y0;

                const int stepj = m1.stepj();
                IT1 A0 = m1.get_col(0).begin().nonConj();
                IT1 A1 = A0; A1.shiftP(stepj);
                IT3 Y = v3.begin();

                const bool c1 = M1::_conj;

                do {
                    A00 = *A0++; A01 = *A1++;
                    Maybe<add>::set(Y0 , *Y);
                    Maybe<add>::add(Y0 , ZProd<c1,false>::prod(A00,X0));
                    Y0 += ZProd<c1,false>::prod(A01,X1);
                    *Y++ = Y0; 
                } while (--M);
            }
        }
    };
    // rs == 3
    template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<18,cs,3,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int M = cs == Unknown ? m1.colsize() : cs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 18 N==3: M,N,cs,rs,x = "<<M<<','<<3<<
                ','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_col_type::const_nonconj_type M1c;
                typedef typename M1c::const_iterator IT1;
                T1 A00, A01, A02;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 X0 = ZProd<false,false>::prod(x,v2.cref(0));
                const PT2 X1 = ZProd<false,false>::prod(x,v2.cref(1));
                const PT2 X2 = ZProd<false,false>::prod(x,v2.cref(2));

                typedef typename V3::value_type T3;
                typedef typename V3::iterator IT3;
                T3 Y0;

                const int stepj = m1.stepj();
                IT1 A0 = m1.get_col(0).begin().nonConj();
                IT1 A1 = A0; A1.shiftP(stepj);
                IT1 A2 = A1; A2.shiftP(stepj);
                IT3 Y = v3.begin();

                const bool c1 = M1::_conj;

                do {
                    A00 = *A0++; A01 = *A1++; A02 = *A2++;
                    Maybe<add>::set(Y0 , *Y);
                    Maybe<add>::add(Y0 , ZProd<c1,false>::prod(A00,X0));
                    Y0 += ZProd<c1,false>::prod(A01,X1);
                    Y0 += ZProd<c1,false>::prod(A02,X2);
                    *Y++ = Y0; 
                } while (--M);
            }
        }
    };
    // rs == 4
    template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<18,cs,4,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int M = cs == Unknown ? m1.colsize() : cs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 18 N==4: M,N,cs,rs,x = "<<M<<','<<4<<
                ','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_col_type::const_nonconj_type M1c;
                typedef typename M1c::const_iterator IT1;
                T1 A00, A01, A02, A03;

                typedef typename V2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                const PT2 X0 = ZProd<false,false>::prod(x,v2.cref(0));
                const PT2 X1 = ZProd<false,false>::prod(x,v2.cref(1));
                const PT2 X2 = ZProd<false,false>::prod(x,v2.cref(2));
                const PT2 X3 = ZProd<false,false>::prod(x,v2.cref(3));

                typedef typename V3::value_type T3;
                typedef typename V3::iterator IT3;
                T3 Y0;

                const int stepj = m1.stepj();
                IT1 A0 = m1.get_col(0).begin().nonConj();
                IT1 A1 = A0; A1.shiftP(stepj);
                IT1 A2 = A1; A2.shiftP(stepj);
                IT1 A3 = A2; A3.shiftP(stepj);
                IT3 Y = v3.begin();

                const bool c1 = M1::_conj;

                do {
                    A00 = *A0++; A01 = *A1++; A02 = *A2++; A03 = *A3++;
                    Maybe<add>::set(Y0 , *Y);
                    Maybe<add>::add(Y0 , ZProd<c1,false>::prod(A00,X0));
                    Y0 += ZProd<c1,false>::prod(A01,X1);
                    Y0 += ZProd<c1,false>::prod(A02,X2);
                    Y0 += ZProd<c1,false>::prod(A03,X3);
                    *Y++ = Y0; 
                } while (--M);
            }
        }
    };

    // algo 21: The basic row major loop
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<21,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 21: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

            if (M) {
                typedef typename M1::value_type T1;
                typedef typename V2::value_type T2;
                typedef typename Traits2<T1,T2>::type PT;
                typedef typename M1::const_row_type M1r;
                PT Yi;

                typedef typename M1r::const_nonconj_type::const_iterator IT1;
                typedef typename V2::const_nonconj_type::const_iterator IT2;
                typedef typename V3::iterator IT3;

                IT1 Ai0 = m1.get_row(0).begin().nonConj();
                const IT2 X0 = v2.begin().nonConj();
                IT3 Y = v3.begin();
                const int Astepi = m1.stepi();

                do {
                    Yi = MultVV_Helper<-4,rs,M1r,V2>::call2(N,Ai0,X0);
                    Maybe<add>::add(*Y++, ZProd<false,false>::prod(x,Yi));
                    Ai0.shiftP(Astepi);
                } while (--M);
            }
        }
    };

    // algo 22: row major, 4 rows at a time
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<22,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void loop_4_rows(
            const int M, const int N,
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVAssert(M%4 == 0);
            TMVAssert(M>0);

            typedef typename M1::value_type T1;
            typedef typename M1::const_row_type::const_iterator IT1;
            T1 A00, A01, A02, A03;
            T1 A10, A11, A12, A13;
            T1 A20, A21, A22, A23;
            T1 A30, A31, A32, A33;

            typedef typename V2::value_type T2;
            typedef typename V2::const_iterator IT2;
            T2 X0, X1, X2, X3;

            typedef typename V3::value_type T3;
            typedef typename V3::iterator IT3;
            T3 Y0, Y1, Y2, Y3;

            int M_4 = M>>2; // M_4 = M/4
            const int stepi = m1.stepi();
            const int stepi_4 = (stepi<<2) - N;
            // step over 4 rows and back to start

            IT1 A0 = m1.get_row(0).begin();
            IT1 A1 = A0; A1.shiftP(stepi);
            IT1 A2 = A1; A2.shiftP(stepi);
            IT1 A3 = A2; A3.shiftP(stepi);
            const IT2 X_begin = v2.begin();
            IT2 X = X_begin;
            const int N_4 = (N>>2);
            const int Nx = N-(N_4<<2);
            IT3 Y = v3.begin();

            const bool dopref = N * sizeof(T1) >= TMV_MV_PREFETCH;

            Prefetch_MultiRead(X.get());
            Prefetch_Read(A0.get());
            Prefetch_Read(A1.get());
            Prefetch_Read(A2.get());
            Prefetch_Read(A3.get());
            Prefetch_Write(Y.get());

            int j;

            TMVAssert(M_4 > 0);
            do {
                X = X_begin;
                // add ix == 1:   y =    (y+mv)
                // add ix != 1:   y += x*(0+mv)
                // !add ix == 1:  y =    (0+mv)
                // !add ix != 1:  y =  x*(0+mv)
                Y0 = Maybe<add && (ix==1)>::select( Y[0] , T3(0) );
                Y1 = Maybe<add && (ix==1)>::select( Y[1] , T3(0) );
                Y2 = Maybe<add && (ix==1)>::select( Y[2] , T3(0) );
                Y3 = Maybe<add && (ix==1)>::select( Y[3] , T3(0) );
                j=N_4; if (j) do {
                    X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
                    A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; 
                    A0 += 4;
                    A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; 
                    A1 += 4;
                    A20 = A2[0]; A21 = A2[1]; A22 = A2[2]; A23 = A2[3]; 
                    A2 += 4;
                    A30 = A3[0]; A31 = A3[1]; A32 = A3[2]; A33 = A3[3]; 
                    A3 += 4;
                    Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
                    Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
                    Y2 += A20 * X0 + A21 * X1 + A22 * X2 + A23 * X3;
                    Y3 += A30 * X0 + A31 * X1 + A32 * X2 + A33 * X3; 
                } while (--j);
                j=Nx; if (j) do {
                    X0 = *X++; 
                    A00 = *A0++; A10 = *A1++; A20 = *A2++; A30 = *A3++;
                    Y0 += A00 * X0; 
                    Y1 += A10 * X0; 
                    Y2 += A20 * X0; 
                    Y3 += A30 * X0;
                } while (--j);
                A0 += stepi_4; A1 += stepi_4; A2 += stepi_4; A3 += stepi_4;
                if (dopref) {
                    Prefetch_Read(A0.get());
                    Prefetch_Read(A1.get());
                    Prefetch_Read(A2.get());
                    Prefetch_Read(A3.get());
                }
                Maybe<add && (ix!=1)>::add(Y[0],ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(Y[1],ZProd<false,false>::prod(x,Y1));
                Maybe<add && (ix!=1)>::add(Y[2],ZProd<false,false>::prod(x,Y2));
                Maybe<add && (ix!=1)>::add(Y[3],ZProd<false,false>::prod(x,Y3));
                Y += 4;
            } while (--M_4);
        }
        template <class M1x, class V3x>
        static void cleanup(
            const int mb,
            const Scaling<ix,T>& x, const M1x& m1, const V2& v2, V3x& v3)
        {
            TMVAssert(mb == 1 || mb == 2 || mb == 3 || mb == 4);
#ifdef TMV_MV_CLEANUP
            switch (mb) {
              case 1 : 
                   MultMV_Helper<1,1,rs,add,ix,T,M1x,V2,V3x>::call(
                       x,m1,v2,v3);
                   break;
              case 2 :
                   MultMV_Helper<23,2,rs,add,ix,T,M1x,V2,V3x>::call(
                       x,m1,v2,v3);
                   break;
              case 3 :
                   MultMV_Helper<23,3,rs,add,ix,T,M1x,V2,V3x>::call(
                       x,m1,v2,v3);
                   break;
              case 4 :
                   MultMV_Helper<23,4,rs,add,ix,T,M1x,V2,V3x>::call(
                       x,m1,v2,v3);
            }
#else
            MultMV_Helper<21,Unknown,rs,add,ix,T,M1x,V2,V3x>::call(
                x,m1,v2,v3);
#endif
        }
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 22: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);
            if (N) {
                if (M > 4) {
                    const int ma = ((M>>2)<<2);
                    const int mb = M - ma;
                    typedef typename M1::const_rowrange_type M1r;
                    typedef typename V3::subvector_type V3s;
                    loop_4_rows(ma,N,x,m1,v2,v3);
                    if (mb) {
                        M1r m1b = m1.cRowRange(ma,M);
                        V3s v3b = v3.cSubVector(ma,M);
                        cleanup(mb,x,m1b,v2,v3b);
                    }
                } else if (M) {
                    cleanup(M,x,m1,v2,v3);
                }
            } else {
                Maybe<!add>::zero(v3);
            }
        }
    };

    // algo 23: cs <= 4, do all rows at once
    // cs == 2
    template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<23,2,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 23 M==2: M,N,cs,rs,x = "<<2<<','<<N<<
                ','<<2<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);

            if (N) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_row_type::const_iterator IT1;
                T1 A00, A01, A02, A03;
                T1 A10, A11, A12, A13;

                typedef typename V2::value_type T2;
                typedef typename V2::const_iterator IT2;
                T2 X0, X1, X2, X3;

                typedef typename V3::value_type T3;
                T3 Y0, Y1;

                const int stepi = m1.stepi();

                IT1 A0 = m1.get_row(0).begin();
                IT1 A1 = A0; A1.shiftP(stepi);
                IT2 X = v2.begin();
                int N_4 = N>>2;
                int Nx = N-(N_4<<2);

                Y0 = Maybe<add && (ix==1)>::select( v3.cref(0) , T3(0) );
                Y1 = Maybe<add && (ix==1)>::select( v3.cref(1) , T3(0) );
                if (N_4) do {
                    X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
                    A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; 
                    A0 += 4;
                    A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; 
                    A1 += 4;
                    Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
                    Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
                } while (--N_4);
                if (Nx) do {
                    X0 = *X++; A00 = *A0++; A10 = *A1++; 
                    Y0 += A00 * X0; 
                    Y1 += A10 * X0; 
                } while (--Nx);
                Maybe<add && (ix!=1)>::add(
                    v3.ref(0),ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(1),ZProd<false,false>::prod(x,Y1));
            } else 
                Maybe<!add>::zero(v3);
        }
    };
    // cs == 3
    template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<23,3,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 23 M==3: M,N,cs,rs,x = "<<3<<','<<N<<
                ','<<3<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);

            if (N) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_row_type::const_iterator IT1;
                T1 A00, A01, A02, A03;
                T1 A10, A11, A12, A13;
                T1 A20, A21, A22, A23;

                typedef typename V2::value_type T2;
                typedef typename V2::const_iterator IT2;
                T2 X0, X1, X2, X3;

                typedef typename V3::value_type T3;
                T3 Y0, Y1, Y2;

                const int stepi = m1.stepi();

                IT1 A0 = m1.get_row(0).begin();
                IT1 A1 = A0; A1.shiftP(stepi);
                IT1 A2 = A1; A2.shiftP(stepi);
                IT2 X = v2.begin();
                int N_4 = N>>2;
                int Nx = N-(N_4<<2);

                Y0 = Maybe<add && (ix==1)>::select( v3.cref(0) , T3(0) );
                Y1 = Maybe<add && (ix==1)>::select( v3.cref(1) , T3(0) );
                Y2 = Maybe<add && (ix==1)>::select( v3.cref(2) , T3(0) );
                if (N_4) do {
                    X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
                    A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; 
                    A0 += 4;
                    A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; 
                    A1 += 4;
                    A20 = A2[0]; A21 = A2[1]; A22 = A2[2]; A23 = A2[3]; 
                    A2 += 4;
                    Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
                    Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
                    Y2 += A20 * X0 + A21 * X1 + A22 * X2 + A23 * X3;
                } while (--N_4);
                if (Nx) do {
                    X0 = *X++; A00 = *A0++; A10 = *A1++; A20 = *A2++; 
                    Y0 += A00 * X0; 
                    Y1 += A10 * X0; 
                    Y2 += A20 * X0; 
                } while (--Nx);
                Maybe<add && (ix!=1)>::add(
                    v3.ref(0),ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(1),ZProd<false,false>::prod(x,Y1));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(2),ZProd<false,false>::prod(x,Y2));
            } else 
                Maybe<!add>::zero(v3);
        }
    };
    // cs == 4
    template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<23,4,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 23 M==4: M,N,cs,rs,x = "<<4<<','<<N<<
                ','<<4<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::isreal);

            if (N) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_row_type::const_iterator IT1;
                T1 A00, A01, A02, A03;
                T1 A10, A11, A12, A13;
                T1 A20, A21, A22, A23;
                T1 A30, A31, A32, A33;

                typedef typename V2::value_type T2;
                typedef typename V2::const_iterator IT2;
                T2 X0, X1, X2, X3;

                typedef typename V3::value_type T3;
                T3 Y0, Y1, Y2, Y3;

                const int stepi = m1.stepi();

                IT1 A0 = m1.get_row(0).begin();
                IT1 A1 = A0; A1.shiftP(stepi);
                IT1 A2 = A1; A2.shiftP(stepi);
                IT1 A3 = A2; A3.shiftP(stepi);
                IT2 X = v2.begin();
                int N_4 = N>>2;
                int Nx = N-(N_4<<2);

                Y0 = Maybe<add && (ix==1)>::select( v3.cref(0) , T3(0) );
                Y1 = Maybe<add && (ix==1)>::select( v3.cref(1) , T3(0) );
                Y2 = Maybe<add && (ix==1)>::select( v3.cref(2) , T3(0) );
                Y3 = Maybe<add && (ix==1)>::select( v3.cref(3) , T3(0) );
                if (N_4) do {
                    X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
                    A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; 
                    A0 += 4;
                    A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; 
                    A1 += 4;
                    A20 = A2[0]; A21 = A2[1]; A22 = A2[2]; A23 = A2[3]; 
                    A2 += 4;
                    A30 = A3[0]; A31 = A3[1]; A32 = A3[2]; A33 = A3[3]; 
                    A3 += 4;
                    Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
                    Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
                    Y2 += A20 * X0 + A21 * X1 + A22 * X2 + A23 * X3;
                    Y3 += A30 * X0 + A31 * X1 + A32 * X2 + A33 * X3; 
                } while (--N_4);
                if (Nx) do {
                    X0 = *X++; 
                    A00 = *A0++; A10 = *A1++; A20 = *A2++; A30 = *A3++; 
                    Y0 += A00 * X0; 
                    Y1 += A10 * X0; 
                    Y2 += A20 * X0; 
                    Y3 += A30 * X0;
                } while (--Nx);
                Maybe<add && (ix!=1)>::add(
                    v3.ref(0),ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(1),ZProd<false,false>::prod(x,Y1));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(2),ZProd<false,false>::prod(x,Y2));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(3),ZProd<false,false>::prod(x,Y3));
            } else 
                Maybe<!add>::zero(v3);
        }
    };

    // algo 25: fully unroll by rows, apply x to v3
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<25,cs,rs,add,ix,T,M1,V2,V3>
    {
        typedef typename V2::value_type T2;
        typedef typename V3::value_type T3;
        template <int I, int N>
        struct Unroller
        {
            static inline void loadx(const V2& v2, T2* X)
            {
                Unroller<I,N/2>::loadx(v2,X);
                Unroller<I+N/2,N-N/2>::loadx(v2,X);
            }
            static inline void sety(
                const Scaling<ix,T>& x, const T3* Y, V3& v3)
            {
                Unroller<I,N/2>::sety(x,Y,v3);
                Unroller<I+N/2,N-N/2>::sety(x,Y,v3);
            }
            static inline void calcy(const M1& m1, const T2* X, T3* Y)
            {
                Unroller<I,N/2>::calcy(m1,X,Y);
                Unroller<I+N/2,N-N/2>::calcy(m1,X,Y);
            }
            static inline void calcyi(
                const int i, const M1& m1, const T2* X, T3& Yi)
            {
                Unroller<I,N/2>::calcyi(i,m1,X,Yi);
                Unroller<I+N/2,N-N/2>::calcyi(i,m1,X,Yi);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void loadx(const V2& v2, T2* X)
            { X[I] = v2.cref(I); }
            static inline void sety(
                const Scaling<ix,T>& x, const T3* Y, V3& v3)
            { Maybe<add>::add(v3.ref(I) , ZProd<false,false>::prod(x,Y[I])); }
            static inline void calcy(const M1& m1, const T2* X, T3* Y)
            { Unroller<0,rs>::calcyi(I,m1,X,Y[I]); }
            static inline void calcyi(
                const int i, const M1& m1, const T2* X, T3& Yi)
                // Note: "I" is really j here...
            {
                Maybe<(I>0)>::add(
                    Yi , ZProd<false,false>::prod(m1.cref(i,I),X[I])); 
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 25: cs,rs,x = "<<cs<<','<<rs<<
                ','<<T(x)<<std::endl;
#endif
            T2 X[rs];
            T3 Y[cs];
            Unroller<0,rs>::loadx(v2,X);
            Unroller<0,cs>::calcy(m1,X,Y); 
            Unroller<0,cs>::sety(x,Y,v3);
        }
    };

    // algo 26: fully unroll by rows, apply x to v2
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<26,cs,rs,add,ix,T,M1,V2,V3>
    {
        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename V3::value_type T3;
        template <int I, int N>
        struct Unroller
        {
            static inline void loadx(
                const Scaling<ix,T>& x, const V2& v2, PT2* X)
            {
                Unroller<I,N/2>::loadx(x,v2,X);
                Unroller<I+N/2,N-N/2>::loadx(x,v2,X);
            }
            static inline void sety(const T3* Y, V3& v3)
            {
                Unroller<I,N/2>::sety(Y,v3);
                Unroller<I+N/2,N-N/2>::sety(Y,v3);
            }
            static inline void calcy(const M1& m1, const PT2* X, T3* Y)
            {
                Unroller<I,N/2>::calcy(m1,X,Y);
                Unroller<I+N/2,N-N/2>::calcy(m1,X,Y);
            }
            static inline void calcyi(
                const int i, const M1& m1, const PT2* X, T3& Yi)
            {
                Unroller<I,N/2>::calcyi(i,m1,X,Yi);
                Unroller<I+N/2,N-N/2>::calcyi(i,m1,X,Yi);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void loadx(
                const Scaling<ix,T>& x, const V2& v2, PT2* X)
            { X[I] = ZProd<false,false>::prod(x , v2.cref(I)); }
            static inline void sety(const T3* Y, V3& v3)
            { Maybe<add>::add(v3.ref(I) , Y[I]); }
            static inline void calcy(const M1& m1, const PT2* X, T3* Y)
            { Unroller<0,rs>::calcyi(I,m1,X,Y[I]); }
            static inline void calcyi(
                const int i, const M1& m1, const PT2* X, T3& Yi)
                // Note: "I" is really j here...
            {
                Maybe<(I>0)>::add(
                    Yi , ZProd<false,false>::prod(m1.cref(i,I),X[I])); 
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 26: cs,rs,x = "<<cs<<','<<rs<<
                ','<<T(x)<<std::endl;
#endif
            PT2 X[rs];
            T3 Y[cs];
            Unroller<0,rs>::loadx(x,v2,X);
            Unroller<0,cs>::calcy(m1,X,Y); 
            Unroller<0,cs>::sety(Y,v3);
        }
    };

    // algo 27: row major, 2 rows at a time, complex v3
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<27,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void loop_2_rows(
            const int M, const int N,
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVAssert(M%2 == 0);
            TMVAssert(M > 0);

            typedef typename M1::value_type T1;
            typedef typename M1::const_row_type::const_nonconj_type M1r;
            typedef typename M1r::const_iterator IT1;
            T1 A00, A10;

            typedef typename V2::value_type T2;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            T2 X0;

            typedef typename V3::value_type T3;
            typedef typename V3::iterator IT3;
            T3 Y0, Y1;

            int M_2 = M>>1; // M_2 = M/2
            const int stepi = m1.stepi();
            const int stepi_2 = (stepi<<1) - N;
            // step over 2 rows and back to start

            IT1 A0 = m1.get_row(0).begin().nonConj();
            IT1 A1 = A0; A1.shiftP(stepi);
            const IT2 X_begin = v2.begin().nonConj();
            IT2 X = X_begin;
            IT3 Y = v3.begin().nonConj();

            const bool c1 = M1::_conj;
            const bool c2 = V2::_conj;

            const bool dopref = N * sizeof(T1) >= TMV_MV_PREFETCH;

            Prefetch_MultiRead(X.get());
            Prefetch_Read(A0.get());
            Prefetch_Read(A1.get());
            Prefetch_Write(Y.get());

            int j;

            TMVAssert(M_2 > 0);
            do {
                X = X_begin;
                Y0 = Maybe<add && (ix==1)>::select( Y[0] , T3(0) );
                Y1 = Maybe<add && (ix==1)>::select( Y[1] , T3(0) );

                j=N; do {
                    X0 = *X++; A00 = *A0++; A10 = *A1++;
                    Y0 += ZProd<c1,c2>::prod(A00,X0);
                    Y1 += ZProd<c1,c2>::prod(A10,X0);
                } while (--j);
                A0 += stepi_2; A1 += stepi_2; 
                if (dopref) {
                    Prefetch_Read(A0.get());
                    Prefetch_Read(A1.get());
                }
                Maybe<add && (ix!=1)>::add(Y[0],ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(Y[1],ZProd<false,false>::prod(x,Y1));
                Y += 2;
            } while (--M_2);
        }
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 27: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);
            const int ma = ((M>>1)<<1);
            const int mb = M - ma;
            typedef typename M1::const_row_type M1r;
            if (N) {
                if (ma) loop_2_rows(ma,N,x,m1,v2,v3);
                if (mb) {
                    Maybe<add>::add(
                        v3.ref(ma) , ZProd<false,false>::prod(
                            x,
                            MultVV_Helper<-4,rs,M1r,V2>::call(
                                m1.get_row(ma),v2)) 
                    );
                }
            } else 
                Maybe<!add>::zero(v3);
        }
    };

    // algo 28: do all rows at once, complex v3: cs <= 4, and must be known
    // cs == 2
    template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<28,2,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 28: M,N,cs,rs,x = "<<2<<','<<N<<
                ','<<2<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);

            if (N) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_row_type::const_nonconj_type M1r;
                typedef typename M1r::const_iterator IT1;
                T1 A00, A10;

                typedef typename V2::value_type T2;
                typedef typename V2::const_nonconj_type::const_iterator IT2;
                T2 X0;

                typedef typename V3::value_type T3;
                T3 Y0, Y1;

                const int stepi = m1.stepi();

                IT1 A0 = m1.get_row(0).begin().nonConj();
                IT1 A1 = A0; A1.shiftP(stepi);
                IT2 X = v2.begin().nonConj();

                const bool c1 = M1::_conj;
                const bool c2 = V2::_conj;

                Y0 = Maybe<add && (ix==1)>::select(v3.cref(0) , T3(0));
                Y1 = Maybe<add && (ix==1)>::select(v3.cref(1) , T3(0));

                do {
                    X0 = *X++; A00 = *A0++; A10 = *A1++;
                    Y0 += ZProd<c1,c2>::prod(A00,X0);
                    Y1 += ZProd<c1,c2>::prod(A10,X0);
                } while (--N);
                Maybe<add && (ix!=1)>::add(
                    v3.ref(0) , ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(1) , ZProd<false,false>::prod(x,Y1));
            } else 
                Maybe<!add>::zero(v3);
        }
    };
    // cs == 3
    template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<28,3,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 28: M,N,cs,rs,x = "<<3<<','<<N<<
                ','<<3<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);

            if (N) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_row_type::const_nonconj_type M1r;
                typedef typename M1r::const_iterator IT1;
                T1 A00, A10, A20;

                typedef typename V2::value_type T2;
                typedef typename V2::const_nonconj_type::const_iterator IT2;
                T2 X0;

                typedef typename V3::value_type T3;
                T3 Y0, Y1, Y2;

                const int stepi = m1.stepi();

                IT1 A0 = m1.get_row(0).begin().nonConj();
                IT1 A1 = A0; A1.shiftP(stepi);
                IT1 A2 = A1; A2.shiftP(stepi);
                IT2 X = v2.begin().nonConj();

                const bool c1 = M1::_conj;
                const bool c2 = V2::_conj;

                Y0 = Maybe<add && (ix==1)>::select(v3.cref(0) , T3(0));
                Y1 = Maybe<add && (ix==1)>::select(v3.cref(1) , T3(0));
                Y2 = Maybe<add && (ix==1)>::select(v3.cref(2) , T3(0));

                do {
                    X0 = *X++; A00 = *A0++; A10 = *A1++; A20 = *A2++;
                    Y0 += ZProd<c1,c2>::prod(A00,X0);
                    Y1 += ZProd<c1,c2>::prod(A10,X0);
                    Y2 += ZProd<c1,c2>::prod(A20,X0);
                } while (--N);
                Maybe<add && (ix!=1)>::add(
                    v3.ref(0) , ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(1) , ZProd<false,false>::prod(x,Y1));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(2) , ZProd<false,false>::prod(x,Y2));
            } else 
                Maybe<!add>::zero(v3);
        }
    };
    // cs == 4
    template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<28,4,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 28: M,N,cs,rs,x = "<<4<<','<<N<<
                ','<<4<<','<<rs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(V3::iscomplex);

            if (N) {
                typedef typename M1::value_type T1;
                typedef typename M1::const_row_type::const_nonconj_type M1r;
                typedef typename M1r::const_iterator IT1;
                T1 A00, A10, A20, A30;

                typedef typename V2::value_type T2;
                typedef typename V2::const_nonconj_type::const_iterator IT2;
                T2 X0;

                typedef typename V3::value_type T3;
                T3 Y0, Y1, Y2, Y3;

                const int stepi = m1.stepi();

                IT1 A0 = m1.get_row(0).begin().nonConj();
                IT1 A1 = A0; A1.shiftP(stepi);
                IT1 A2 = A1; A2.shiftP(stepi);
                IT1 A3 = A2; A3.shiftP(stepi);
                IT2 X = v2.begin().nonConj();

                const bool c1 = M1::_conj;
                const bool c2 = V2::_conj;

                Y0 = Maybe<add && (ix==1)>::select(v3.cref(0) , T3(0));
                Y1 = Maybe<add && (ix==1)>::select(v3.cref(1) , T3(0));
                Y2 = Maybe<add && (ix==1)>::select(v3.cref(2) , T3(0));
                Y3 = Maybe<add && (ix==1)>::select(v3.cref(3) , T3(0));

                do {
                    X0 = *X++; 
                    A00 = *A0++; A10 = *A1++; A20 = *A2++; A30 = *A3++;
                    Y0 += ZProd<c1,c2>::prod(A00,X0);
                    Y1 += ZProd<c1,c2>::prod(A10,X0);
                    Y2 += ZProd<c1,c2>::prod(A20,X0);
                    Y3 += ZProd<c1,c2>::prod(A30,X0);
                } while (--N);
                Maybe<add && (ix!=1)>::add(
                    v3.ref(0) , ZProd<false,false>::prod(x,Y0));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(1) , ZProd<false,false>::prod(x,Y1));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(2) , ZProd<false,false>::prod(x,Y2));
                Maybe<add && (ix!=1)>::add(
                    v3.ref(3) , ZProd<false,false>::prod(x,Y3));
            } else 
                Maybe<!add>::zero(v3);
        }
    };

    // TODO algo 31, etc. = SSE2, 41, etc. = SSE
    
    // algo 51: colmajor, ix == 0, !add, so might want to use algo 53
    // to do the scaling at the end.
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<51,cs,rs,false,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<0,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 51: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int algo2 = V3::iscomplex ? 17 : 12;
            if (M >= N)
                MultMV_Helper<algo2,cs,rs,false,ix,T,M1,V2,V3>::call(
                    x,m1,v2,v3);
            else 
                MultMV_Helper<53,cs,rs,false,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 52: colmajor, ix==0 && add, so might need to copy v3
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<52,cs,rs,true,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(TMV_MV_ZeroIX);
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 52: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int algo2 = V3::iscomplex ? 17 : 12;
            // If we have a non-trivial scale and the size of 
            // v2 is (significantly) larger than v3 and add = true, 
            // then it is worth doing m1*v2 first and then applying x.

            // On the other hand, if add = false, then we can get around
            // the temporary by multiplying by x after doing m1*v2.
            // This is done in algo 53, but I mention it here to 
            // explain why we have the add=true requirement.

            if (N > TMV_MV_COPY_SCALE_RATIO * M) 
                MultMV_Helper<84,cs,rs,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            else 
                MultMV_Helper<algo2,cs,rs,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 54: colmajor, unknown cs -- check if it is small
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<54,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int algo1 = 
                V3::_step == Unknown ? 86 : 
#ifdef TMV_MV_SCALE
                ( TMV_MV_ZeroIX && add ) ? 52 : 
                ( TMV_MV_ZeroIX && !add ) ? 51 : 
#endif
                V3::iscomplex ? 17 : 12;
            const int algo2 = 
                V2::_step == 1 ? ( V3::iscomplex ? 28 : 23 ) : algo1;

            TMVStaticAssert(cs == Unknown);
            const int M = m1.colsize();
#ifdef PRINTALGO_MV
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 54: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            if (M <= 4) {
                // then it is worth figuring out what M is.
                switch (M) {
                  case 0 :
                       // do nothing
                       break;
                  case 1 :
                       MultMV_Helper<1,1,rs,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                       break;
                  case 2 :
                       MultMV_Helper<algo2,2,rs,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                       break;
                  case 3 :
                       MultMV_Helper<algo2,3,rs,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                       break;
                  case 4 :
                       MultMV_Helper<algo2,4,rs,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                }
            } else 
                MultMV_Helper<algo1,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 53: column major, !add, apply x at the end
    template <int cs, int rs, class T, class M1, class V2, class V3>
    struct MultMV_Helper<53,cs,rs,false,0,T,M1,V2,V3>
    {
        static void call(
            const Scaling<0,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 53: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int algo2 = V3::iscomplex ? 17 : 12;
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            MultMV_Helper<algo2,cs,rs,false,1,RT,M1,V2,V3>::call(one,m1,v2,v3);
            ScaleV_Helper<-3,cs,0,T,V3>::call(x,v3);
        }
    };


    // algo 62: rowmajor, might need to copy v2
    // If we have a non-trivial scale and the size of 
    // v3 is (significantly) larger than v2, then it is worth
    // making x*v2 as a temporary.
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<62,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(TMV_MV_ZeroIX);
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
#ifdef PRINTALGO_MV
            std::cout<<"MV algo 62: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int algo2 = V3::iscomplex ? 27 : 22;

            if (M > TMV_MV_COPY_SCALE_RATIO * N) 
                MultMV_Helper<82,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            else 
                MultMV_Helper<algo2,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 64: rowmajor, unknown rs -- see if it is small
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<64,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int algo1 = 
                V2::_step == Unknown ? 83 : 
#ifdef TMV_MV_SCALE
                TMV_MV_ZeroIX ? 62 :
#endif
                V3::iscomplex ? 27 : 22;
            const int algo2 = V3::iscomplex ? 18 : 13;
            TMVStaticAssert(rs == Unknown);
            const int N = m1.rowsize();
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            std::cout<<"MV algo 64: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            if (N <= 4) {
                // then it is worth figuring out what N is.
                switch (N) {
                  case 0 :
                       // If M > 0, and !add, then we need to do v3.setZero().
                       Maybe<!add>::zero(v3);
                       break;
                  case 1 :
                       MultMV_Helper<2,cs,1,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                       break;
                  case 2 :
                       MultMV_Helper<algo2,cs,2,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                       break;
                  case 3 :
                       MultMV_Helper<algo2,cs,3,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                       break;
                  case 4 :
                       MultMV_Helper<algo2,cs,4,add,ix,T,M1,V2,V3>::call(
                           x,m1,v2,v3);
                }
            } else 
                MultMV_Helper<algo1,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 81: copy v2
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<81,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int N = rs == Unknown ? m1.rowsize() : rs;
            const int M = cs == Unknown ? m1.colsize() : cs;
            std::cout<<"MV algo 81: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typename V3::noalias_type v3na = v3.noAlias();
            MultMV<add>(x,m1,v2.copy(),v3na);
        }
    };

    template <int ix, class T, class V> class ProdXV;

    // algo 82: copy x*v2
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<82,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int N = rs == Unknown ? m1.rowsize() : rs;
            const int M = cs == Unknown ? m1.colsize() : cs;
            std::cout<<"MV algo 82: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename V3::noalias_type v3na = v3.noAlias();
            MultMV<add>(one,m1,ProdXV<ix,T,V2>(x,v2).calc(),v3na);
        }
    };

    // algo 83: Copy v2, figure out where to put x
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<83,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
            if (N >= M) {
                MultMV_Helper<81,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                MultMV_Helper<82,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };
    // If ix == 1, don't need the branch - just go to 81
    template <int cs, int rs, bool add, class T, class M1, class V2, class V3>
    struct MultMV_Helper<83,cs,rs,add,1,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<1,T>& x, const M1& m1, const V2& v2, V3& v3)
        { MultMV_Helper<81,cs,rs,add,1,T,M1,V2,V3>::call(x,m1,v2,v3); }
    };

    // algo 84: v3c = m1*v2, v3 (+)= x*v3c
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<84,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 84: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typename V3::noalias_type v3na = v3.noAlias();
            MultXV<add>(x,(m1*v2).calc(),v3na);
        }
    };

    template <int ix, class T, class M, class V> class ProdMV;

    // algo 85: v3c = x*m1*v2, v3 (+)= v3c
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<85,cs,rs,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_MV
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
            std::cout<<"MV algo 85: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename V3::noalias_type v3na = v3.noAlias();
            MultXV<add>(one,ProdMV<ix,T,M1,V2>(x,m1,v2).calc(),v3na);
        }
    };

    // algo 86: Use temporary for v3, figure out where to put x
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<86,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int M = cs == Unknown ? m1.colsize() : cs;
            const int N = rs == Unknown ? m1.rowsize() : rs;
            if (N >= M) {
                MultMV_Helper<84,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                MultMV_Helper<85,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };
    // If ix == 1, don't need the branch - just go to 84
    template <int cs, int rs, bool add, class T, class M1, class V2, class V3>
    struct MultMV_Helper<86,cs,rs,add,1,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<1,T>& x, const M1& m1, const V2& v2, V3& v3)
        { MultMV_Helper<84,cs,rs,add,1,T,M1,V2,V3>::call(x,m1,v2,v3); }
    };

    // algo 90: call inst
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<90,cs,rs,false,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<90,cs,rs,true,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };

    // algo 91: call inst alias
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<91,cs,rs,false,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<91,cs,rs,true,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<97,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultMV_Helper<-2,cs,rs,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<197,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultMV_Helper<99,cs,rs,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<98,cs,rs,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            if ( !SameStorage(m1,v3) &&
                 !SameStorage(v2,v3) ) {
                // No aliasing
                MultMV_Helper<-2,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else if (SameStorage(m1,v3)) {
                // Use temporary for v3
                MultMV_Helper<86,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                // SameStorage(v2,v3)
                // Use temporary for v2
                MultMV_Helper<83,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<99,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                V3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    template <int cs, int rs, bool cm>
    struct MultMV_Unroll_Helper
    {
        enum { unroll = (
                ( cs == Unknown || rs == Unknown ) ? false :
                IntTraits2<cs,rs>::prod > TMV_MV_UNROLL ? false :
                // These maxunroll values are empirical for my machine,
                // comparing the speed for the fully unrolled algorithm
                // to the regular, non-unrolled algorithm, so they might
                // not be optimal for other machines, or even other
                // compilers, but they should be good enough.
                // And for that matter, I only tested real double, so they
                // might not be best for other data types either.
                cm ? (
                    cs <= 4 ? rs <= 11 :
                    rs == 2 ? cs <= 8 :
                    rs == 3 ? cs <= 6 :
                    rs == 4 ? cs <= 5 :
                    rs == 5 ? cs <= 20 :
                    rs == 6 ? cs <= 10 :
                    rs <= 10 ? cs <= 9 :
                    rs <= 12 ? cs <= 8 :
                    rs <= 19 ? cs <= 7 :
                    rs <= 21 ? cs <= 6 :
                    rs <= 25 ? cs <= 5 :
                    false ) :
                ( // rowmajor
                    rs == 2 ? cs == 2 :
                    rs <= 4 ? cs <= 4 :
                    cs <= 3 ? rs <= 7 :
                    cs <= 4 ? rs <= 11 :
                    cs <= 6 ? rs <= 15 :
                    cs <= 8 ? rs <= 11 :
                    cs <= 9 ? rs <= 8 :
                    cs <= 11 ? rs <= 7 :
                    cs <= 16 ? rs <= 6 : 
                    cs <= 23 ? rs <= 5 :
                    false ) ) };
    };

    // algo -4: No branches or copies
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<-4,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool unroll = 
                MultMV_Unroll_Helper<cs,rs,M1::_colmajor>::unroll;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 1 :
                ( rs == 1 ) ? 2 :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ?
                ( M1::_colmajor ? 11 : 21 ) :
                M1::_colmajor ? (
                    unroll ? (
                        ( cs <= TMV_MV_COL_UNROLL ? 
                          ( (rs<cs) ? 26 : 25 ) :
                          ( (rs<cs) ? 16 : 15 ) ) ) :
                    ( cs != Unknown && cs <= 4 && V2::_step == 1 ) ? (
                        (V3::iscomplex ? 28 : 23) ) :
                    ( !add && TMV_MV_ZeroIX && 
                      rs != Unknown && cs != Unknown && rs > cs ) ? 53 : 
                    ( rs != Unknown && rs <= 4 ) ? (
                        (V3::iscomplex ? 18 : 13) ) :
                    V3::iscomplex ? 17 : 12 ) :
                M1::_rowmajor ? (
                    unroll ? ( (rs<cs) ? 26 : 25 ) :
                    ( rs != Unknown && rs <= 4 && V3::_step == 1 ) ? (
                        (V3::iscomplex ? 18 : 13 ) ) : 
                    ( cs != Unknown && cs <= 4 ) ? (
                        (V3::iscomplex ? 28 : 23) ) :
                    V3::iscomplex ? 27 : 22 ) :
                V2::_step == 1 ? 21 : V3::_step == 1 ? 11 : 21;
            MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<-3,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = cs or rs == 0, so nothing to do
            //  1 = cs == 1: reduces to trivial MultXV function
            //  2 = rs == 1: reduces to trivial MultVV function
            //  3 = call InstMultMV
            //
            // Column Major:
            // 11 = column major, simple for loop
            // 12 = column major, 4 columns at a time
            // 13 = column major, rs is known and <= 4
            // 15 = fully unroll by columns, apply x to v3
            // 16 = fully unroll by columns, apply x to v2
            // 17 = column major, complex v3, 2 columns at a time
            // 18 = column major, complex v3, rs is known and <= 4
            //
            // Row Major:
            // 21 = row major, simple for loop
            // 22 = row major, 4 rows at a time
            // 23 = row major, cs is known and <= 4
            // 25 = fully unroll by rows, apply x to v3
            // 26 = fully unroll by rows, apply x to v2
            // 27 = row major, complex v3, 2 rows at a time
            // 28 = row major, complex v3, cs is known and <= 4
            //
            // Column Major, meta algorithms
            // 51 = column major, ix==0 && !add, so might want algo 53
            // 52 = column major, ix==0 && add, so might need to copy v3
            // 53 = column major, !add, apply x at the end
            // 54 = column major, unknown cs, check if it is small
            //
            // Row Major, meta algorithms
            // 62 = row major, ix==0, so might need to copy v2
            // 64 = row major, unknown rs, check if it is small
            //
            // Copy a vector to new storage:
            // 81 = copy v2
            // 82 = copy x*v2
            // 83 = copy v2, figure out where to put x
            // 84 = temp v3 = m1*v2
            // 85 = temp v3 = x*m1*v2
            // 86 = temp v3, figure out where to put x
            const bool unroll = 
                MultMV_Unroll_Helper<cs,rs,M1::_colmajor>::unroll;
#ifdef TMV_MV_SCALE
            const bool doscale = true;
#else
            const bool doscale = false;
#endif
#ifdef TMV_MV_SMALL
            const bool dosmall = true;
#else
            const bool dosmall = false;
#endif
            const bool V3stepxx = (V3::_step == Unknown);
            const bool V2stepxx = (V2::_step == Unknown);
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : // trivial - nothing to do
                ( cs == 1 ) ? 1 : // trivial - cs = 1
                ( rs == 1 ) ? 2 : // trivial - rs = 1
                TMV_OPT == 0 ? ( M1::_colmajor ? 11 : 21 ) :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ? (
                    ( M1::_colmajor ? 11 : 21 ) ) :

                M1::_colmajor ? 
                ( dosmall && cs == Unknown ) ? 54 :
                cs == Unknown ? (
                    V3stepxx ? 86 : 
                    ( doscale && TMV_MV_ZeroIX && add ) ? 52 :
                    ( TMV_MV_ZeroIX && !add ) ? 51 : 
                    V3::iscomplex ? 17 : 12 ) :
                rs == Unknown ? (
                    ( cs > TMV_MV_COPY_SIZE && V3stepxx ) ? 86 :
                    ( doscale && TMV_MV_ZeroIX && add ) ? 52 : 
                    ( TMV_MV_ZeroIX && !add ) ? 51 : 
                    V3::iscomplex ? 17 : 12 ) :
                unroll ? (
                    ( cs <= TMV_MV_COL_UNROLL ? 
                      ( (rs<cs) ? 26 : 25 ) :
                      ( (rs<cs) ? 16 : 15 ) ) ) :
                ( cs <= 4 && V2::_step == 1 ) ? V3::iscomplex ? 28 : 23 :
                ( cs > TMV_MV_COPY_SIZE && V3stepxx ) ? (
                    rs > cs ? 84 : 85 ) :
                ( add && TMV_MV_ZeroIX &&
                  rs > IntTraits2<TMV_MV_COPY_SCALE_RATIO,cs>::prod ) ? 84 :
                ( !add && TMV_MV_ZeroIX && rs > cs ) ? 53 : 
                rs <= 4 ? (V3::iscomplex ? 18 : 13) :
                V3::iscomplex ? 17 : 12 :

                M1::_rowmajor ? 
                ( dosmall && rs == Unknown ) ? 64 :
                rs == Unknown ? (
                    V2stepxx ? 83 :
                    ( doscale && TMV_MV_ZeroIX ) ? 62 :
                    V3::iscomplex ? 27 : 22 ) :
                cs == Unknown ? (
                    ( rs > TMV_MV_COPY_SIZE && V2stepxx ) ? 83 :
                    ( doscale && TMV_MV_ZeroIX ) ? 62 :
                    V3::iscomplex ? 27 : 22 ) :
                unroll ? ( (rs<cs) ? 26 : 25 ) :
                ( rs <= 4 && V3::_step == 1 ) ? V3::iscomplex ? 18 : 13 : 
                ( rs > TMV_MV_COPY_SIZE && V2stepxx ) ? (
                    rs > cs ? 81 : 82 ) :
                ( doscale && TMV_MV_ZeroIX &&
                  cs > IntTraits2<TMV_MV_COPY_SCALE_RATIO,rs>::prod ) ? 82 :
                cs <= 4 ? V3::iscomplex ? 28 : 23 :
                V3::iscomplex ? 27 : 22 :

                // nomajor -- don't do anything fancy
                V2::_step == 1 ? 21 : V3::_step == 1 ? 11 : 21;
#ifdef PRINTALGO_MV
            std::cout<<"InlineMultMV: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<std::endl;
            std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
#endif
#ifdef XDEBUG_MV
            typedef typename V3::real_type RT;
            typedef typename V3::value_type T3;
            Matrix<T3> m1c = m1;
            Vector<T3> v2c = v2;
            Vector<T3> v3i = v3;
            Vector<T3> v3c = v3;
            for(int i=0; i<v3.size(); ++i) {
                Maybe<add>::add(
                    v3c(i), ZProd<false,false>::prod(
                        x,MultVV(m1c.row(i),v2c)));
            }
#endif
            MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#ifdef XDEBUG_MV
            if (Norm(v3-v3c) > 1.e-3*(Norm(m1c)*Norm(v2c)+(add?Norm(v3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"v2 = "<<v2c<<std::endl;
                std::cout<<"v3 = "<<v3i<<std::endl;
                std::cout<<"v3 => "<<v3<<std::endl;
                std::cout<<"Correct v3 = "<<v3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<-2,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 201 :
                ( rs == 1 ) ? 202 :
                V3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultMV_Helper<-1,cs,rs,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 : 
                ( cs == 1 ) ? 101 :
                ( rs == 1 ) ? 102 :
                V3::_checkalias ? 99 : 
                -2;
            MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(m1.colsize() == v3.size());
        TMVAssert(m1.rowsize() == v2.size());
        typedef typename M1::value_type T1;
        const int cs = Sizes<M1::_colsize,V3::_size>::size;
        const int rs = Sizes<M1::_rowsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultMV_Helper<-1,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(m1.colsize() == v3.size());
        TMVAssert(m1.rowsize() == v2.size());
        typedef typename M1::value_type T1;
        const int cs = Sizes<M1::_colsize,V3::_size>::size;
        const int rs = Sizes<M1::_rowsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultMV_Helper<-3,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineAliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same));
        TMVAssert(m1.colsize() == v3.size());
        TMVAssert(m1.rowsize() == v2.size());
        typedef typename M1::value_type T1;
        const int cs = Sizes<M1::_colsize,V3::_size>::size;
        const int rs = Sizes<M1::_rowsize,V2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultMV_Helper<98,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

} // namespace tmv

#undef TMV_MV_SMALL
#undef TMV_MV_CLEANUP
#undef TMV_MV_SCALE

#undef TMV_MV_UNROLL
#undef TMV_MV_COPY_SIZE
#undef TMV_MV_PREFETCH
#undef TMV_MV_COPY_SCALE_RATIO
#undef TMV_MV_COL_UNROLL
#undef TMV_MV_ZeroIX

#endif 
