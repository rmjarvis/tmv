

#ifndef TMV_MultMM_H
#define TMV_MultMM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_MultMV.h"
#include "TMV_Rank1VVM.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_MultXM_Funcs.h"
#include "TMV_Prefetch.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef PRINTALGO_MM
#include <iostream>
#endif

#ifdef XDEBUG_MM
#include <iostream>
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMM.h"
#endif

// Check for small (<=3) values of cs, rs, or xs
// This leads to significant speed improvements for such matrices at little
// cost to larger matrices, but the large increase in code size and the fact
// that it only benefits a few particular sizes of matrices mean that
// we require TMV_OPT = 3 for it.
#if TMV_OPT >= 3
#define TMV_MM_OPT_SMALL
#endif

// Select the cleanup code for the edge that doesn't fit into KB-sized
// blocks according to the exact value of K.  This can have a significant 
// increase in the speed for moderately sized matrices where the edges 
// are not a negligible fraction of the calculation. So we only 
// require TMV_OPT = 2
#if TMV_OPT >= 2
#define TMV_MM_OPT_CLEANUP
#endif

// For very large matrices, a recursive block algorithm is faster than
// the simpler looping block algorithm.  The dividing line is governed
// by the parameter TMV_MM_MIN_RECURSIVE (see below).  It is not all that 
// expensive with respect to code bloat, and it is significantly faster for 
// large matrices, so we only require TMV_OPT = 2.
#if TMV_OPT >= 2
#define TMV_MM_USE_RECURSIVE_BLOCK
#endif

// For extremely large matrices, it is worth using a recursive Winograd
// algorithm where only 7 matrix multiplies are needed to do the 8
// submatrix calculations in:
// [ A B ] * [ E F ] = [ I J ]
// [ C D ]   [ G H ]   [ K L ]
// This trick is repeated recursively until the minimum size is less than
// MIN_WINOGRAD (see below), and then we call a more standard algorithm.
// It is only used for extremely large matrices.  Since most users 
// won't normally use such large matrices, we only do this for OPT = 3.
#if TMV_OPT >= 3
#define TMV_MM_USE_WINOGRAD
#endif

// The algorithms for large matrices usually involve making temporary
// matrices with the data stored in blocks.  Therefore, to be good
// C++ programmers, we put these into try/catch blocks to catch any
// possible bad_alloc throws.  We have an algorithm specially designed
// for the post-catch calculation, which divides M,N in half and does
// the four sections separately.  Since bad_alloc's are pretty rare and
// this adds a non-trivial amount of code bloat (although not huge by 
// any stretch), I require TMV_OPT=3.
#if TMV_OPT >= 3
#define TMV_MM_OPT_BAD_ALLOC
#endif

// There are a number of values used in the algorithm selection
// that are either arbitrary or empirical.
// So I put them all here to make them easier to change and to  
// track down in the code.

// UNROLL is the maximum nops to unroll.
// This doesn't seem to be necessary.  The non-unrolling versions are
// very fast even for very small matrices.
#define TMV_MM_UNROLL 0

// The block size to use for large matrices.  The value should
// be chosen such that three NxN sized matrices all fit into the 
// L1 cache.  
// For algo 63, I have special block sizes chosen rather than use this,
// since I find that for this algorithm, it is better to use smaller
// blocks.  Also there are special choices when using SSE.
// So this is just for the regular recursion when not copying (algo 61).
#define TMV_MM_BLOCK_SIZE 64

// PREFETCH is the crossover memory size to start using prefetch commands.
// This is undoubtedly a function of the L1 (and L2?) cache size,
// but 2KBytes is probably not too bad for most machines.
// (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_MM_PREFETCH 2048

// The minimum size to use a recursive Winograd algorithm
// (This is repeated in TMV_MultMM_Winograd.h where it is also used.)
#ifdef TMV_MM_USE_RECURSIVE_BLOCK
#define TMV_MM_MIN_WINOGRAD 2048
#else
#define TMV_MM_MIN_WINOGRAD 1024
#endif

// The minimum value of Mb*Nb*Kb*Kb to use recursive algorithm.
// This formula for the crossover between algorithms 63 and 64 is
// purely empirical on my MacBook (with an Intel Core 2 Duo processor),
// so it might not even be the right thing to parametrize on other
// machines.  On the other hand, the difference between the two algorithms
// isn't _that_ extreme.  At most I've seen about a factor of 2, so
// if I use the less optimal algorithm for some sizes, it's not that
// terrible a mistake.
#define TMV_MM_MIN_RECURSIVE 16*1024

// The minimum value of (M*N*K / 16^3) to use multiple threads.
// (We also require that M or N > 64 so we have something to split.)
// The computer that I did most of the testing on only has 2 processors,
// so this optimization might be particular to having 2 threads.
// TODO: Investigate the timings when more than 2 threads are available.
#define TMV_MM_OPENMP_THRESH 64
// There are some sizes where this selection formula doesn't quite work,
// but it seems to be not too bad.  But maybe this deserves a bit more
// work to find a better formula, since the ratio of time for omp to
// not using omp is not at all monotonic in M*N*K.  So perhaps a 
// different formula would be better.

//                    no omp    omp    omp/not    M*N*K/16^3
//  16 x  16 x  16    1.275    9.503    7.453       1
//  16 x  16 x  64    1.224    4.716    3.852       4
//  64 x  16 x  16    1.168    2.891    2.475       4
//  16 x  64 x  16    1.122    2.752    2.452       4
//  16 x  16 x 256    1.259    2.158    1.714      16
//  64 x  64 x  16    1.070    1.303    1.217      16
//  16 x 256 x  16    1.126    1.250    1.110      16
// 256 x  16 x  16    1.186    1.186    1.000      16
//  64 x  16 x  64    1.145    1.096    0.957      16
//  16 x  64 x  64    1.146    1.051    0.917      16

// ^^^^ My formula says don't use omp
// vvvv My formula says use omp

//  64 x  64 x  64    0.655    0.679    1.036      64
//  64 x  16 x 256    1.201    1.213    1.009      64
//  16 x  64 x 256    1.191    1.195    1.003      64
// 256 x  16 x 256    0.566    0.551    0.973     256
//  64 x 256 x  64    0.465    0.412    0.886     256
//  64 x  64 x 256    0.619    0.546    0.882     256
// 256 x  64 x  64    0.463    0.406    0.876     256
// 256 x  64 x 256    0.431    0.334    0.774    1024
//  64 x 256 x  16    1.055    0.836    0.792      64
// 256 x  64 x  16    1.098    0.840    0.765      64
//  64 x 256 x 256    0.444    0.317    0.713    1024
// 256 x 256 x  16    1.080    0.722    0.668     256
// 256 x 256 x  64    0.429    0.271    0.631    1024
// 256 x 256 x 256    0.401    0.239    0.596    4096
//  16 x 256 x  64    1.124    0.539    0.479      64
// 256 x  16 x  64    1.143    0.537    0.469      64
//  16 x 256 x 256    1.143    0.528    0.461     256

// The minimum value of (M*N*K / 16^3) to use copying algorithm
// where each block is copied to temporary storage in a block structure.
#define TMV_MM_MIN_COPY_ALGO 4

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_MM_ZeroIX (ix==0)
//#define TMV_MM_ZeroIX (ix!=1)

namespace tmv {

    // Defined in TMV_MultMM.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstMatrixView<T1,C1>& m1, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);

    //
    // Matrix * Matrix
    //

    template <int algo, int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper;

    // algo 0: Trivial, nothing to do.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<0,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: Trivial, just m3.setZero();
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<1,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { m3.setZero(); }
    };

    // algo 2: cs == 1, so reduces to MultMV
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 2: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::row_type M3r;
            M1r m1r = m1.get_row(0);
            M2t m2t = m2.transpose();
            M3r m3r = m3.get_row(0);
            MultMV_Helper<-4,rs,xs,add,ix,T,M2t,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 102: same as 2, but use -1 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<102,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 102: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::row_type M3r;
            M1r m1r = m1.get_row(0);
            M2t m2t = m2.transpose();
            M3r m3r = m3.get_row(0);
            MultMV_Helper<-1,rs,xs,add,ix,T,M2t,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<202,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 202: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<
                ','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::row_type M3r;
            M1r m1r = m1.get_row(0);
            M2t m2t = m2.transpose();
            M3r m3r = m3.get_row(0);
            MultMV_Helper<-2,rs,xs,add,ix,T,M2t,M1r,M3r>::call(x,m2t,m1r,m3r);
        }
    };

    // algo 3: rs == 1, so reduces to MultMV
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<3,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 3: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultMV_Helper<-4,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 103: same as 3, but use -1 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<103,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 103: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultMV_Helper<-1,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 203: same as 3, but use -2 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<203,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 203: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<
                ','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultMV_Helper<-2,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 4: xs == 1, so reduces to Rank1Update
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<4,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"MM algo 4: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_type M1c;
            typedef typename M2::const_row_type M2r;
            M1c m1c = m1.get_col(0);
            M2r m2r = m2.get_row(0);
            Rank1VVM_Helper<-4,cs,rs,add,ix,T,M1c,M2r,M3>::call(x,m1c,m2r,m3);
        }
    };

    // algo 104: same as 4, but use -1 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<104,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"MM algo 104: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_type M1c;
            typedef typename M2::const_row_type M2r;
            M1c m1c = m1.get_col(0);
            M2r m2r = m2.get_row(0);
            Rank1VVM_Helper<-1,cs,rs,add,ix,T,M1c,M2r,M3>::call(x,m1c,m2r,m3);
        }
    };

    // algo 204: same as 4, but use -2 algo
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<204,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"MM algo 204: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<
                ','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_type M1c;
            typedef typename M2::const_row_type M2r;
            M1c m1c = m1.get_col(0);
            M2r m2r = m2.get_row(0);
            Rank1VVM_Helper<-2,cs,rs,add,ix,T,M1c,M2r,M3>::call(x,m1c,m2r,m3);
        }
    };

    // algo 5: **R, Transpose to get the colmajor version
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<5,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 5: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultMM_Helper<-3,rs,cs,xs,add,ix,T,M2t,M1t,M3t>::call(
                x,m2t,m1t,m3t); 
        }
    };

    // algo 405: same as 5 but go back to -4, rather than -3.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<405,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 405: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultMM_Helper<-4,rs,cs,xs,add,ix,T,M2t,M1t,M3t>::call(
                x,m2t,m1t,m3t); 
        }
    };

    // algo 11: CCC -- Loop over N
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<11,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // This algorithm is the same as for RCC, so just call that one,
            // rather than duplicate code.
            MultMM_Helper<21,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 12: CCC -- Do product in 4x2 blocks of m2
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<12,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_iterator IT1;
        typedef typename M2::value_type T2;
        typedef typename M2::const_col_type::const_iterator IT2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename M3::value_type T3;
        typedef typename M3::col_type::iterator IT3;
        typedef typename Traits2<T1,PT2>::type PT3;

        static void loop_42(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVAssert(N%2 == 0);
            TMVAssert(N >= 2);

            const int N_2 = N>>1; // N_2 = N/2
            const int K_4 = K>>2; // K_4 = K/4 
            const int Ka = K_4<<2; 
            const int Kc = K-Ka;

            T1 a0, a1, a2, a3;
            const int Astepi = m1.stepi();
            const int Astepj = m1.stepj();
            const int Astepj_4 = (Astepj<<2) - M*Astepi;
            // step over 4 cols and back to start
            const int Astepj_1 = Astepj - M*Astepi;
            // step over 1 col and back to start
            const int Astart = -K*Astepj;
            const int Astartx = -Ka*Astepj;
            // step from the end all the way back to the start

            PT2 b00, b10, b20, b30;
            PT2 b01, b11, b21, b31;
            const int Bstepi = m2.stepi();
            const int Bstepj = m2.stepj();
            const int Bstepj_2 = (Bstepj<<1) - K*Bstepi;
            // step over 2 cols and back to start

            PT3 c00, c01;
            const int Cstepj = m3.stepj();
            const int Cstepj_2 = (Cstepj<<1);
            // step over 2 cols

            IT1 A0 = m1.get_col(0).begin();
            IT1 A1 = A0; A1.shiftP(Astepj);
            IT1 A2 = A1; A2.shiftP(Astepj);
            IT1 A3 = A2; A3.shiftP(Astepj);

            IT2 B0 = m2.get_col(0).begin();
            IT2 B1 = B0; B1.shiftP(Bstepj);

            IT3 C0 = m3.get_col(0).begin();
            IT3 C1 = C0; C1.shiftP(Cstepj);

            const bool dopref = K * sizeof(T1) >= TMV_MM_PREFETCH;

            Prefetch_MultiRead(A0.get());
            Prefetch_MultiRead(A1.get());
            Prefetch_Read(B0.get());
            Prefetch_Read(B1.get());
            Prefetch_MultiWrite(C0.get());
            Prefetch_MultiWrite(C1.get());

            int j,i,k;

            j = N_2; do {
                k = K_4; if (k) do {
                    b00 = x * *B0++; b01 = x * *B1++;
                    b10 = x * *B0++; b11 = x * *B1++;
                    b20 = x * *B0++; b21 = x * *B1++;
                    b30 = x * *B0++; b31 = x * *B1++;
                    i = M; if (i) do {
                        c00 = *C0;
                        a0 = *A0++;
                        c00 += a0 * b00;
                        c01 = *C1;
                        c01 += a0 * b01;
                        a1 = *A1++;
                        c00 += a1 * b10;
                        c01 += a1 * b11;
                        a2 = *A2++;
                        c00 += a2 * b20;
                        c01 += a2 * b21;
                        a3 = *A3++;
                        c00 += a3 * b30;
                        *C0++ = c00;
                        c01 += a3 * b31;
                        *C1++ = c01;
                    } while (--i);
                    A0.shiftP(Astepj_4); A1.shiftP(Astepj_4);
                    A2.shiftP(Astepj_4); A3.shiftP(Astepj_4);
                    C0 -= M; C1 -= M; 
                    if (dopref) {
                        Prefetch_Read(A0.get());
                        Prefetch_Read(A1.get());
                        Prefetch_Read(A2.get());
                        Prefetch_Read(A3.get());
                    }
                } while (--k);
                k = Kc; if (k) do {
                    b00 = x * *B0++;   b01 = x * *B1++;
                    i = M; if (i) do {
                        c00 = *C0;
                        a0 = *A0++;
                        c00 += a0 * b00;
                        *C0++ = c00;
                        c01 = *C1;
                        c01 += a0 * b01;
                        *C1++ = c01;
                    } while (--i);
                    A0.shiftP(Astepj_1);
                    C0 -= M; C1 -= M;
                    if (dopref) {
                        Prefetch_Read(A0.get());
                    }
                } while (--k);
                A0.shiftP(Astart); A1.shiftP(Astartx); 
                A2.shiftP(Astartx); A3.shiftP(Astartx);
                B0.shiftP(Bstepj_2); B1.shiftP(Bstepj_2);
                C0.shiftP(Cstepj_2); C1.shiftP(Cstepj_2);
                if (dopref) {
                    Prefetch_Read(B0.get());
                    Prefetch_Read(B1.get());
                    Prefetch_MultiWrite(C0.get());
                    Prefetch_MultiWrite(C1.get());
                }
            } while (--j);
        }
        template <int algo2, int z>
        struct Helper2;

        template <int z>
        struct Helper2<0,z> // algo2 = 0 (cs,rs,xs=0)
        {
            static TMV_INLINE void call(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) 
            {}
        };
        template <int z>
        struct Helper2<2,z> // algo2 = 2 (rs is even)
        {
            static inline void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                TMVStaticAssert(rs != TMV_UNKNOWN);
                TMVStaticAssert(rs != 0);
                TMVStaticAssert(rs%2 == 0);
                const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
                const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
                if (M && K) loop_42(M,rs,K,x,m1,m2,m3); 
            }
        };
        template <int z>
        struct Helper2<3,z> // algo2 = 3 (normal operation)
        {
            static void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
                const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
                const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
                typedef typename M1::const_col_type M1c;
                typedef typename M2::const_col_type M2c;
                typedef typename M2::const_row_sub_type M2r;
                typedef typename M3::col_type M3c;
                typedef typename M3::colrange_type M3a;

                if (M && N && K) {
                    const int na = ((N>>1)<<1);
                    const int nb = N-na;
                    if (na) {
                        loop_42(M,na,K,x,m1,m2,m3);
                    }
                    if (nb) {
                        M2c m2na = m2.get_col(na);
                        M3c m3na = m3.get_col(na);
                        MultMV_Helper<-4,cs,xs,true,ix,T,M1,M2c,M3c>::call(
                            x,m1,m2na,m3na);
                    }
                }
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo2 = 
                ( cs == 0 || rs == 0 || xs == 0 ) ? 0 :
                ( rs == TMV_UNKNOWN ) ? 3 :
                ( rs%2 == 0 ) ? 2 : 3;
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 12: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
            std::cout<<"add = "<<add<<", algo2 = "<<algo2<<std::endl;
#endif
            Maybe<!add>::zero(m3);

            Helper2<algo2,1>::call(x,m1,m2,m3);
        }
    };

    // algo 21: RCC -- Loop over N
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<21,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 21: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            for(int j=0;j<N;++j) {
                // m3.get_col(j) += x * m1 * m2.get_col(j)
                M2c m2j = m2.get_col(j);
                M3c m3j = m3.get_col(j);
                MultMV_Helper<-4,cs,xs,add,ix,T,M1,M2c,M3c>::call(
                    x,m1,m2j,m3j);
            }
        }
    };

    // algo 22: RCC -- Do product in 2x2 blocks of m3
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<22,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        typedef typename M1::value_type T1;
        typedef typename M1::const_row_type::const_iterator IT1;
        typedef typename M2::value_type T2;
        typedef typename M2::const_col_type::const_iterator IT2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename M3::value_type T3;
        typedef typename M3::col_type::iterator IT3;

        static void loop_22(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVAssert(M%2 == 0);
            TMVAssert(N%2 == 0);

            T1 a0, a1;
            const int Astepi = m1.stepi();
            const int Astepj = m1.stepj();
            const int Astepi_2 = (Astepi<<1) - K*Astepj;
            // step over 2 rows and back to start
            const int Astart = -M*Astepi;
            // step from the end all the way back to the start

            PT2 b0, b1;
            const int Bstepj = m2.stepj();
            const int Bstepj_2 = (Bstepj<<1);
            // step over 2 cols (don't need the -K for this one)

            T3 C00, C01, C10, C11;
            T3 c00, c01, c10, c11;
            const int Cstepi = m3.stepi();
            const int Cstepj = m3.stepj();
            const int Cstepj_2 = (Cstepj<<1) - M*Cstepi;
            // step over 2 cols and back to start

            const int M_2 = M>>1; // M_2 = M/2
            const int N_2 = N>>1; // N_2 = N/2
            const int K_4 = (K-1)>>2; // K_4 = (K-1)/4 
            // The K-1 here  ^^^ is to make sure Kc > 0, since the regular
            // loop accesses memory past the end, so we use the Kc loop
            // to make sure that doesn't lead to a seg fault.
            const int Kc = K-(K_4<<2);
            const int Kcm1 = Kc-1;

            IT1 A0 = m1.get_row(0).begin();
            IT1 A1 = A0; A1.shiftP(Astepi);

            IT2 B0 = m2.get_col(0).begin();
            IT2 B1 = B0; B1.shiftP(Bstepj);

            IT3 C0 = m3.get_col(0).begin();
            IT3 C1 = C0; C1.shiftP(Cstepj);

            const bool dopref = K * sizeof(T1) >= TMV_MM_PREFETCH;

            Prefetch_MultiRead(A0.get());
            Prefetch_MultiRead(A1.get());
            Prefetch_MultiRead(B0.get());
            Prefetch_MultiRead(B1.get());
            Prefetch_Write(C0.get());
            Prefetch_Write(C1.get());

            int j,i,k;

            j = N_2; do {
                i = M_2; do {
                    C00 = Maybe<add && (ix==1)>::select(C0[0] , T3(0));
                    C01 = Maybe<add && (ix==1)>::select(C1[0] , T3(0));
                    C10 = Maybe<add && (ix==1)>::select(C0[1] , T3(0));
                    C11 = Maybe<add && (ix==1)>::select(C1[1] , T3(0));
                    b0 = B0[0]; 
                    k = K_4; if (k) do {
                        a0 = A0[0]; c00 = a0 * b0; C00 += c00;
                        a1 = A1[0]; c10 = a1 * b0; C10 += c10;
                        b1 = B1[0]; c01 = a0 * b1; C01 += c01;
                        b0 = B0[1]; c11 = a1 * b1; C11 += c11;
                        a0 = A0[1]; c00 = a0 * b0; C00 += c00;
                        a1 = A1[1]; c10 = a1 * b0; C10 += c10;
                        b1 = B1[1]; c01 = a0 * b1; C01 += c01;
                        b0 = B0[2]; c11 = a1 * b1; C11 += c11;
                        a0 = A0[2]; c00 = a0 * b0; C00 += c00;
                        a1 = A1[2]; c10 = a1 * b0; C10 += c10;
                        b1 = B1[2]; c01 = a0 * b1; C01 += c01;
                        b0 = B0[3]; c11 = a1 * b1; C11 += c11;
                        a0 = A0[3]; c00 = a0 * b0; C00 += c00; A0 += 4;
                        a1 = A1[3]; c10 = a1 * b0; C10 += c10; A1 += 4;
                        b1 = B1[3]; c01 = a0 * b1; C01 += c01; B1 += 4;
                        b0 = B0[4]; c11 = a1 * b1; C11 += c11; B0 += 4;
                    } while (--k);
                    k = Kcm1; if (k) do {
                        a0 = *A0++; c00 = a0 * b0; C00 += c00;
                        a1 = *A1++; c10 = a1 * b0; C10 += c10;
                        b1 = *B1++; c01 = a0 * b1; C01 += c01;
                        b0 = *++B0; c11 = a1 * b1; C11 += c11;
                    } while (--k);
                    a0 = *A0++; c00 = a0 * b0; C00 += c00;
                    a1 = *A1++; c10 = a1 * b0; C10 += c10;
                    b1 = *B1++; c01 = a0 * b1; C01 += c01;
                    ++B0;       c11 = a1 * b1; C11 += c11;
                    Maybe<add && (ix!=1)>::add(C0[0] , x * C00);
                    Maybe<add && (ix!=1)>::add(C0[1] , x * C10); C0 += 2;
                    Maybe<add && (ix!=1)>::add(C1[0] , x * C01);
                    Maybe<add && (ix!=1)>::add(C1[1] , x * C11); C1 += 2;
                    A0.shiftP(Astepi_2); A1.shiftP(Astepi_2);
                    B0 -= K; B1 -= K;
                    if (dopref) {
                        Prefetch_Read(A0.get());
                        Prefetch_Read(A1.get());
                    }
                } while (--i);
                A0.shiftP(Astart); A1.shiftP(Astart);
                B0.shiftP(Bstepj_2); B1.shiftP(Bstepj_2);
                C0.shiftP(Cstepj_2); C1.shiftP(Cstepj_2);
                if (dopref) {
                    Prefetch_MultiRead(B0.get());
                    Prefetch_MultiRead(B1.get());
                    Prefetch_Write(C0.get());
                    Prefetch_Write(C1.get());
                }
            } while (--j);
        }
        template <int algo2, int z>
        struct Helper2;

        template <int z>
        struct Helper2<0,z> // algo2 = 0 (cs,rs=0 or (xs=0 && add))
        {
            static TMV_INLINE void call(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
        };
        template <int z>
        struct Helper2<1,z> // algo2 = 1 (xs=0 and !add)
        {
            static TMV_INLINE void call(
                const Scaling<ix,T>& , const M1& , const M2& , M3& m3) 
            { m3.setZero(); }
        };
        template <int z>
        struct Helper2<2,z> // algo2 = 2 (cs and rs are even)
        {
            static inline void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                TMVStaticAssert(cs != TMV_UNKNOWN);
                TMVStaticAssert(rs != TMV_UNKNOWN);
                TMVStaticAssert(cs != 0);
                TMVStaticAssert(rs != 0);
                TMVStaticAssert(cs%2 == 0);
                TMVStaticAssert(rs%2 == 0);
                const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
                if (K) loop_22(cs,rs,K,x,m1,m2,m3); 
                else if (!add) m3.setZero();
            }
        };
        template <int z>
        struct Helper2<3,z> // algo2 = 3 (normal operation)
        {
            static void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
                const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
                const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
                typedef typename M1::const_row_type M1r;
                typedef typename M2::const_col_type M2c;
                typedef typename M2::const_colrange_type::const_transpose_type M2t;
                typedef typename M3::col_type M3c;
                typedef typename M3::row_sub_type M3r;

                if (M && N && K) {
                    const int na = ((N>>1)<<1);
                    const int nb = N-na;
                    const int ma = ((M>>1)<<1);
                    const int mb = M-ma;
                    if (ma && na) {
                        loop_22(ma,na,K,x,m1,m2,m3);
                    }
                    if (nb) {
                        M2c m2na = m2.get_col(na);
                        M3c m3na = m3.get_col(na);
                        MultMV_Helper<-4,cs,xs,add,ix,T,M1,M2c,M3c>::call(
                            x,m1,m2na,m3na);
                    }
                    if (mb) {
                        M1r m1ma = m1.get_row(ma);
                        M2t m2t = m2.cColRange(0,na).transpose();
                        M3r m3ma = m3.get_row(ma,0,na);
                        const int rsx = rs == TMV_UNKNOWN ? TMV_UNKNOWN : ((rs>>1)<<1);
                        MultMV_Helper<-4,rsx,xs,add,ix,T,M2t,M1r,M3r>::call(
                            x,m2t,m1ma,m3ma);
                    }
                } else if (M && N && !add) { // K == 0 && !add
                    m3.setZero();
                }
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo2 = 
                ( cs == 0 || rs == 0 ) ? 0 :
                ( xs == 0 ) ? ( add ? 0 : 1 ) :
                ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? 3 :
                ( cs%2 == 0 && rs%2 == 0 ) ? 2 : 3;
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 22: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
            std::cout<<"add = "<<add<<", algo2 = "<<algo2<<std::endl;
#endif
            Helper2<algo2,1>::call(x,m1,m2,m3);
        }
    };

    // algo 31: CRC -- Loop over K
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<31,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"MM algo 31: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_col_type M1c;
            typedef typename M2::const_row_type M2r;
            Maybe<!add>::zero(m3);
            for(int k=0;k<K;++k) {
                // m3 += x * m1.get_col(k) * m2.get_row(j)
                M1c m1k = m1.get_col(k);
                M2r m2k = m2.get_row(k);
                Rank1VVM_Helper<-4,cs,rs,true,ix,T,M1c,M2r,M3>::call(
                    x,m1k,m2k,m3);
            }
        }
    };

    // algo 32: CRC -- Do product in 4x2 blocks of m2
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<32,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_iterator IT1;
        typedef typename M2::value_type T2;
        typedef typename M2::const_row_type::const_iterator IT2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename M3::value_type T3;
        typedef typename M3::col_type::iterator IT3;
        typedef typename Traits2<T1,PT2>::type PT3;

        static void loop_42(
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVAssert(K%4 == 0);
            TMVAssert(K >= 4);

            const int N_2 = N>>1; // N_2 = N/2
            const int K_4 = K>>2; // K_4 = K/4 
            const int Na = (N_2<<1);
            const int Nc = N-Na;

            T1 a0, a1, a2, a3;
            const int Astepj = m1.stepj();
            const int Astepj_4 = (Astepj<<2);
            // step over 4 cols

            PT2 b00, b10, b20, b30;
            PT2 b01, b11, b21, b31;
            const int Bstepi = m2.stepi();
            const int Bstepj = m2.stepj();
            const int Bstepi_4 = (Bstepi<<2) - N*Bstepj;
            // step over 4 rows and back to start

            PT3 c00, c01;
            const int Cstepi = m3.stepi();
            const int Cstepj = m3.stepj();
            const int Cstepj_2 = (Cstepj<<1) - M*Cstepi;
            // step over 2 cols and back to start
            const int Cstepj_1 = Cstepj - M*Cstepi;
            // step over 1 col and back to start
            const int Cstart = -N*Cstepj;
            const int Cstartx = -Na*Cstepj;
            // step from the end all the way back to the start

            IT1 A0 = m1.get_col(0).begin();
            IT1 A1 = A0; A1.shiftP(Astepj);
            IT1 A2 = A1; A2.shiftP(Astepj);
            IT1 A3 = A2; A3.shiftP(Astepj);

            IT2 B0 = m2.get_row(0).begin();
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT2 B3 = B2; B3.shiftP(Bstepi);

            IT3 C0 = m3.get_col(0).begin();
            IT3 C1 = C0; C1.shiftP(Cstepj);

            const bool dopref = K * sizeof(T1) >= TMV_MM_PREFETCH;

            Prefetch_Read(A0.get());
            Prefetch_Read(A1.get());
            Prefetch_MultiRead(B0.get());
            Prefetch_MultiRead(B1.get());
            Prefetch_MultiWrite(C0.get());
            Prefetch_MultiWrite(C1.get());

            int j,i,k;

            k = K_4; do {
                j = N_2; if (j) do {
                    b00 = x * *B0++; b01 = x * *B0++;
                    b10 = x * *B1++; b11 = x * *B1++;
                    b20 = x * *B2++; b21 = x * *B2++;
                    b30 = x * *B3++; b31 = x * *B3++;
                    i = M; if (i) do {
                        c00 = *C0;
                        a0 = *A0++;
                        c00 += a0 * b00;
                        c01 = *C1;
                        c01 += a0 * b01;
                        a1 = *A1++;
                        c00 += a1 * b10;
                        c01 += a1 * b11;
                        a2 = *A2++;
                        c00 += a2 * b20;
                        c01 += a2 * b21;
                        a3 = *A3++;
                        c00 += a3 * b30;
                        *C0++ = c00;
                        c01 += a3 * b31;
                        *C1++ = c01;
                    } while (--i);
                    C0.shiftP(Cstepj_2); C1.shiftP(Cstepj_2);
                    A0 -= M; A1 -= M; A2 -= M; A3 -= M; 
                } while (--j);
                if (Nc) {
                    b00 = x * *B0++;
                    b10 = x * *B1++;
                    b20 = x * *B2++;
                    b30 = x * *B3++;
                    i = M; if (i) do {
                        c00 = *C0;
                        a0 = *A0++;
                        c00 += a0 * b00;
                        a1 = *A1++;
                        c00 += a1 * b10;
                        a2 = *A2++;
                        c00 += a2 * b20;
                        a3 = *A3++;
                        c00 += a3 * b30;
                        *C0++ = c00;
                    } while (--i);
                    C0.shiftP(Cstepj_1);
                    A0 -= M; A1 -= M; A2 -= M; A3 -= M; 
                }
                C0.shiftP(Cstart); C1.shiftP(Cstartx);
                A0.shiftP(Astepj_4); A1.shiftP(Astepj_4); 
                A2.shiftP(Astepj_4); A3.shiftP(Astepj_4); 
                B0.shiftP(Bstepi_4); B1.shiftP(Bstepi_4);
                B2.shiftP(Bstepi_4); B3.shiftP(Bstepi_4);
                if (dopref) {
                    Prefetch_Read(A0.get());
                    Prefetch_Read(A1.get());
                    Prefetch_Read(B0.get());
                    Prefetch_Read(B1.get());
                }
            } while (--k);
        }
        template <int algo2, int z>
        struct Helper2;

        template <int z>
        struct Helper2<0,z> // algo2 = 0 (cs,rs,xs=0)
        {
            static TMV_INLINE void call(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
        };
        template <int z>
        struct Helper2<2,z> // algo2 = 2 (rs is even)
        {
            static TMV_INLINE void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                TMVStaticAssert(xs != TMV_UNKNOWN);
                TMVStaticAssert(xs != 0);
                TMVStaticAssert(xs%4 == 0);
                const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
                const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
                if (M && N) loop_42(M,N,xs,x,m1,m2,m3); 
            }
        };
        template <int z>
        struct Helper2<3,z> // algo2 = 3 (normal operation)
        {
            static void call(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
                const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
                const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
                typedef typename M1::const_colrange_type M1c;
                typedef typename M2::const_rowrange_type M2r;

                if (M && N && K) {
                    const int ka = ((K>>2)<<2);
                    const int kb = K-ka;
                    if (ka) {
                        loop_42(M,N,ka,x,m1,m2,m3);
                    }
                    if (kb) {
                        M1c m1c = m1.cColRange(ka,K);
                        M2r m2r = m2.cRowRange(ka,K);
                        const int xsx = 
                            xs == TMV_UNKNOWN ? TMV_UNKNOWN : (xs - ((xs>>2)<<2));
                        MultMM_Helper<31,cs,rs,xsx,true,ix,T,M1c,M2r,M3>::call(
                            x,m1c,m2r,m3);
                    }
                }
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo2 = 
                ( cs == 0 || rs == 0 || xs == 0 ) ? 0 :
                ( xs == TMV_UNKNOWN ) ? 3 :
                ( xs%4 == 0 ) ? 2 : 3;
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 32: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
            std::cout<<"add = "<<add<<", algo2 = "<<algo2<<std::endl;
#endif
            Maybe<!add>::zero(m3);

            Helper2<algo2,1>::call(x,m1,m2,m3);
        }
    };

    // The following algorithms are used for larger matrices that don't 
    // fully fit into the L1 cache, and so memory movement is the 
    // bottleneck.  There are a number of approaches to tackle this issue,
    // but they all usually involve splitting the matrices into blocks,
    // which do fit into the L1 cache.  Then the trick is to multiply
    // these blocks in the correct order to minimize the number of 
    // cache misses.
    // 
    // The optimized BLAS libraries that you can get for most systems 
    // know the exact size of the L1, L2 and L3 caches for the system
    // and move the memory around to optimally take advantage of these
    // sizes.  Since TMV tries to be as portable as possible, we can't
    // be quite as specific in terms of what size to expect, since these
    // numbers change a lot from one system to another.  
    //
    // So we have a single parameter, TMV_MM_BLOCK_SIZE, to define what is a 
    // good choice for the block size.  64 is a good choice for many systems,
    // but if you want to mess around with this number, you might be able to 
    // slightly improve the speed of matrix multiplication with a different 
    // value.
    //
    // For algo 63 (see description below), we actually find that a smaller
    // block is better.  We hard code its blocks to be 16x16 (for MxN) and
    // either 16, 32 or 64 for K depending on the value type in the matrix.
    // This isn't a parameter, because we have special optimized routines for
    // multiplying these blocks that are specialized to these particular 
    // values.
    //
    // Anyway, since we don't know the L1 and L2/L3 cache sizes, we chose
    // algorithms that try to be somewhat "cache oblivious".  This is
    // a term that was apparently coined by Charles Leiserson, and 
    // was applied to matrix multiplication by Harald Prokop in his
    // masters thesis, titled "Cache-Oblivious Algorithms".
    //
    // Basically, the algorithms can be made "oblivious" to the cache size by
    // recursively breaking up the matrix into smaller chunks until you get
    // down to a single block, and then call a kernel function for that
    // block.  The recursion procedure means that when you go back to a
    // particular block again, it is as likely as possible to still be in 
    // either the L2 or L3 cache from when you used it the previous time.
    // So you avoid as many loads from the main memory as possible.

    // algo 61: Split problem into smaller blocks recursively
    // This algorithm is just Prokop's divide and conquer algorithm,
    // which is probably the simplest cache-oblivious algorithm there is.
    // Also, we don't actually use this one anymore, but I leave it here,
    // since it is fairly simple and easy to understand.  It is the basis
    // for algo 63, which gets a bit more complicated.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<61,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        enum { ccc = M1::_colmajor && M2::_colmajor && M3::_colmajor };
        enum { rcc = M1::_rowmajor && M2::_colmajor && M3::_colmajor };
        enum { crc = M1::_colmajor && M2::_rowmajor && M3::_colmajor };
        enum { MB = (
                ccc ? TMV_MM_BLOCK_SIZE :
                rcc ? TMV_MM_BLOCK_SIZE :
                TMV_MM_BLOCK_SIZE ) };
        enum { NB = (
                ccc ? TMV_MM_BLOCK_SIZE :
                rcc ? TMV_MM_BLOCK_SIZE :
                TMV_MM_BLOCK_SIZE ) };
        enum { KB = (
                ccc ? TMV_MM_BLOCK_SIZE :
                rcc ? TMV_MM_BLOCK_SIZE :
                TMV_MM_BLOCK_SIZE ) };
        enum { algo2 = (
                ccc ? ( M3::iscomplex ? 11 : 12 ) :
                rcc ? ( M3::iscomplex ? 21 : 22 ) :
                crc ? ( M3::iscomplex ? 31 : 32 ) :
                21 ) };
        enum { lnMB = IntTraits<MB>::log };
        enum { lnNB = IntTraits<NB>::log };
        enum { lnKB = IntTraits<KB>::log };

        template <bool addx, int csx, int rsx, int xsx>
        static void call1(
            const int i1, const int j1, const int k1, 
            const int i2, const int j2, const int k2, 
            const int M, const int N, const int K,
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVAssert(M >= 1);
            TMVAssert(N >= 1);
            TMVAssert(K >= 1);
            if (M == 1 && N == 1 && K == 1)
                call2<addx,csx,rsx,xsx>(i1,j1,k1,i2,j2,k2,x,m1,m2,m3);
            else if (M >= N && M >= K) { 
                // M is largest
                TMVAssert(M >= N && M >= K);
                const int Mx = M>>1; // = M/2
                const int im = i1 + (Mx<<lnMB); // = i1 + Mx * MB
                call1<addx,MB,rsx,xsx>(i1,j1,k1,im,j2,k2,Mx,N,K,x,m1,m2,m3);
                call1<addx,csx,rsx,xsx>(im,j1,k1,i2,j2,k2,M-Mx,N,K,x,m1,m2,m3);
            } else if (N >= M && N >= K) { 
                // N is largest
                TMVAssert(N >= M && N >= K);
                const int Nx = N>>1; // = N/2
                const int jm = j1 + (Nx<<lnNB); // = j1 + Nx * NB
                call1<addx,csx,NB,xsx>(i1,j1,k1,i2,jm,k2,M,Nx,K,x,m1,m2,m3);
                call1<addx,csx,rsx,xsx>(i1,jm,k1,i2,j2,k2,M,N-Nx,K,x,m1,m2,m3);
            } else { 
                // K is largest
                TMVAssert(K >= M && K >= N);
                const int Kx = K>>1; // = K/2
                const int km = k1 + (Kx<<lnKB); // = k1 + Kx * KB
                call1<addx,csx,rsx,KB>(i1,j1,k1,i2,j2,km,M,N,Kx,x,m1,m2,m3);
                call1<true,csx,rsx,xsx>(i1,j1,km,i2,j2,k2,M,N,K-Kx,x,m1,m2,m3);
            }
        }
        template <bool addx, int csx, int rsx, int xsx>
        static inline void call2(
            const int i1, const int j1, const int k1, 
            const int i2, const int j2, const int k2, 
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_submatrix_type M1s;
            typedef typename M2::const_submatrix_type M2s;
            typedef typename M3::submatrix_type M3s;

            if (i2 == i1 || j2 == j1 || k2 == k1) return;

            TMVAssert(i2 > i1);
            TMVAssert(j2 > j1);
            TMVAssert(k2 > k1);
            TMVAssert(csx == TMV_UNKNOWN || csx == i2-i1);
            TMVAssert(rsx == TMV_UNKNOWN || rsx == j2-j1);
            TMVAssert(xsx == TMV_UNKNOWN || xsx == k2-k1);

            M1s m1s = m1.cSubMatrix(i1,i2,k1,k2);
            M2s m2s = m2.cSubMatrix(k1,k2,j1,j2);
            M3s m3s = m3.cSubMatrix(i1,i2,j1,j2);

            MultMM_Helper<algo2,csx,rsx,xsx,addx,ix,T,M1s,M2s,M3s>::call(
                x,m1s,m2s,m3s);
        }
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
#ifdef PRINTALGO_MM
            std::cout<<"MM algo 61: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(ccc || rcc || crc);
            // These will be the correct value for the final non-block calls.
            const int csx = cs == TMV_UNKNOWN ? TMV_UNKNOWN : (cs-((cs>>lnMB)<<lnMB));
            const int rsx = rs == TMV_UNKNOWN ? TMV_UNKNOWN : (rs-((rs>>lnNB)<<lnNB));
            const int xsx = xs == TMV_UNKNOWN ? TMV_UNKNOWN : (xs-((xs>>lnKB)<<lnKB));
            if (M > MB || N > NB || K > KB) {
                const int Mb = (M>>lnMB)+1; // = M/MB + 1
                const int Nb = (N>>lnNB)+1; // = N/NB + 1
                const int Kb = (K>>lnKB)+1; // = K/KB + 1
                call1<add,csx,rsx,xsx>(0,0,0,M,N,K,Mb,Nb,Kb,x,m1,m2,m3); 
            } else {
                MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            }
        }
    };

    // algo 63: Recursively subdivide m1,m2,m3 just like algo 61, but
    // use temporary memory to copy the blocks to compact storage.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<63,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 63: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            MultMM_RecursiveBlock<add>(x,m1,m2,m3);
        }
    };

    // algo 64: Similar to algo 63 in that we copy the blocks to temporary
    // storage.  However, simply loop over M,N and then do the full K
    // row/column block in one set.  This is faster than algo 63 for
    // smaller matrices.  The dividing line between the two algorithms is 
    // set by TMV_MM_MIN_RECURSIVE.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<64,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 64: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            MultMM_Block<add>(x,m1,m2,m3);
        }
    };

    // algo 66: Split problem into 4 parts then call another algorithm
    // to finish the work.  This is used after a bad_alloc error
    // to split the problem into chunks that will require less temporary
    // memory.  
    //
    // Note: I think I am always catching the bad_alloc's before anything
    // gets written to m3.  But I'm not 100% sure of that.  
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<66,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
#ifdef PRINTALGO_MM
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 66: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            const int Mx = (((M-1)>>5)+1)<<4; // = M/2 rounded up to 16 block
            const int Nx = (((N-1)>>5)+1)<<4; // = N/2 rounded up to 16 block
            typedef typename M1::const_rowrange_type M1r;
            typedef typename M2::const_colrange_type M2c;
            typedef typename M3::submatrix_type M3s;

            TMV_Warning(
                "Caught bad_alloc error in MultMM.\n"
                "Splitting the problem into 4 sub-problems.");
            try {
                M1r m1a = m1.cRowRange(0,Mx);
                M1r m1b = m1.cRowRange(Mx,M);
                M2c m2a = m2.cColRange(0,Nx);
                M2c m2b = m2.cColRange(Nx,N);
                M3s m3aa = m3.cSubMatrix(0,Mx,0,Nx);
                M3s m3ba = m3.cSubMatrix(Mx,M,0,Nx);
                M3s m3ab = m3.cSubMatrix(0,Mx,Nx,N);
                M3s m3bb = m3.cSubMatrix(Mx,M,Nx,N);
                const int xx = TMV_UNKNOWN;

                MultMM_Helper<-2,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
                    x,m1a,m2a,m3aa);
                MultMM_Helper<-2,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
                    x,m1a,m2b,m3ab);
                MultMM_Helper<-2,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
                    x,m1b,m2a,m3ba);
                MultMM_Helper<-2,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
                    x,m1b,m2b,m3bb);
            } catch (std::bad_alloc) {
                // If failed again, use an algorithm that doesn't 
                // allocate memory.
                MultMM_Helper<67,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            }
        }
    };

    // algo 67: Revert to an algorithm that doesn't allocate memory.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<67,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 67: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            const bool ccc = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool rcc = M1::_rowmajor && M2::_colmajor && M3::_colmajor;
            const bool crc = M1::_colmajor && M2::_rowmajor && M3::_colmajor;
            const int algo2 = 
                ccc ? ( M3::iscomplex ? 11 : 12 ) :
                rcc ? ( M3::iscomplex ? 21 : 22 ) :
                crc ? ( M3::iscomplex ? 31 : 32 ) :
                21;

            TMV_Warning(
                "Caught bad_alloc error in MultMM.\n"
                "Using (slower) algorithm that doesn't allocate temporary "
                "memory");
            MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                x,m1,m2,m3);
        }
    };

    // algo 68: Use Alberto and Nicolau's hybrid Winograd algorithm
    // This algorithm is based on the paper by Paolo D’Alberto and
    // Alexandru Nicolau titled "Adaptive Winograd’s Matrix Multiplications"
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<68,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 68: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            MultMM_Winograd<add>(x,m1,m2,m3);
        }
    };

#ifdef _OPENMP
    // algo 69: Split problem into smaller parts with OpenMP for 
    // parallelization.
    // We take a pretty simple approach here, and just split up m2 and m3
    // matrices by columns and let each thread do a single matrix.
    // Then each thread calls algo 63 to calculate its product.
    // Also, we require that all but the last thread has a column width
    // that is a multiple of 16.  This way we get the maximum advantage from
    // our blocking structure while keeping the threads as balanced as 
    // possible.
    //
    // Also, if M is larger than N, then we split the rows of m1 and m3
    // instead.
    //
    // TODO: This way of doing things means that m1 (or m2) gets copied
    // separately for each thread, which is obviously inefficient.
    // So we can improve this by loading up the m1 sub-blocks once
    // here before doing the algo 63 recursion.  But this would require
    // significantly more coding to arrange this correctly, so I haven't
    // done so yet.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<69,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 69: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            MultMM_OpenMP<add>(x,m1,m2,m3);
        }
    };
#endif

    // algo 71: Determine which algorithm to use based on the runtime
    // knowledge of the sizes.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<71,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
#ifdef PRINTALGO_MM
            std::cout<<"MM algo 71: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            const int Mb = (M>>6); // = M/64
            const int Nb = (N>>6); // = N/64
            const int Kb = (K>>6); // = K/64
            const int Mc = M < 16 ? 1 : (M>>4); // = M/16
            const int Nc = N < 16 ? 1 : (N>>4); // = N/16
            const int Kc = K < 16 ? 1 : (K>>4); // = K/16
            const bool twobig = (Mb&&Nb) || (Mb&&Kb) || (Nb&&Kb);

            TMVStaticAssert(!M3::_rowmajor);
            const bool ccc = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool rcc = M1::_rowmajor && M2::_colmajor && M3::_colmajor;
            const bool crc = M1::_colmajor && M2::_rowmajor && M3::_colmajor;
            
            const int algo2 = 
                ccc ? ( M3::iscomplex ? 11 : 
                        cs != TMV_UNKNOWN && xs != TMV_UNKNOWN ? 11 :
                        12 ) :
                rcc ? ( M3::iscomplex ? 21 :
                        cs != TMV_UNKNOWN && xs != TMV_UNKNOWN ? 21 :
                        22 ) :
                crc ? ( M3::iscomplex ? 31 :
                        cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ? 31 :
                        32 ) : 
                11;

            // Put the small matrix option first, so it doesn't have to 
            // go through a bunch of if/else statements.  For large matrices,
            // all these if/else's don't matter for the total time.
            if ( (M < 16 && N < 16 && K < 16) ||
                 (M <= 3 || N <= 3 || K <= 3) ||
                 ( ( M < 16 || N < 16 || K < 16 ) &&
                   ( !twobig || (Mc * Nc * Kc < TMV_MM_MIN_COPY_ALGO) ) ) )
                MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
#ifdef _OPENMP
            else if (!omp_in_parallel() && (Mb || Nb) && 
                     ( Mc * Nc * Kc >= TMV_MM_OPENMP_THRESH ) )
                MultMM_Helper<69,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
            else 
                MultMM_Helper<72,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 72: Same as 71, except 63,64,68 are the only options.
    // (Acutally, now 71 moves on to this one for the last three options,
    //  since it helps reduce the code size a bit when doing openmp.)
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<72,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            std::cout<<"MM algo 72: M,N,K,cs,rs,xs,x = "<<
                m3.colsize()<<','<<m3.rowsize()<<','<<m1.rowsize()<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif

#ifdef TMV_MM_USE_WINOGRAD
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
#endif

#ifdef TMV_MM_USE_WINOGRAD
            if (M >= TMV_MM_MIN_WINOGRAD && N >= TMV_MM_MIN_WINOGRAD && 
                K >= TMV_MM_MIN_WINOGRAD)
                MultMM_Helper<68,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else
#endif
                MultMM_Helper<73,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 73: Same as 72, but no Winograd.
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<73,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            std::cout<<"MM algo 73: M,N,K,cs,rs,xs,x = "<<
                m3.colsize()<<','<<m3.rowsize()<<','<<m1.rowsize()<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif

#ifdef TMV_MM_USE_RECURSIVE_BLOCK
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            const int Mb = (M>>6); // = M/64
            const int Nb = (N>>6); // = N/64
            const int Kb = (K>>6); // = K/64
#endif

            TMVStaticAssert(!M3::_rowmajor);

#ifdef TMV_MM_USE_RECURSIVE_BLOCK
            if (Mb*Nb*Kb*Kb >= TMV_MM_MIN_RECURSIVE)
                MultMM_Helper<63,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else
#endif
                MultMM_Helper<64,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 75: Unknown cs, check if M is small
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<75,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef TMV_MM_OPT_SMALL
            TMVStaticAssert(cs == TMV_UNKNOWN);
            const int M = m3.colsize();
#ifdef PRINTALGO_MM
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 75: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            if (M <= 3) {
                // then it is worth figuring out what M is.
                switch (M) {
                  case 0 :
                       // do nothing
                       break;
                  case 1 :
                       MultMM_Helper<-4,1,rs,xs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 2 :
                       MultMM_Helper<-4,2,rs,xs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 3 :
                       MultMM_Helper<-4,3,rs,xs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                }
            } else 
#endif
                MultMM_Helper<71,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
        }
    };

    // algo 76: Unknown rs, check if N is small
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<76,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef TMV_MM_OPT_SMALL
            TMVStaticAssert(rs == TMV_UNKNOWN);
            const int N = m3.rowsize();
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"MM algo 76: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            if (N <= 3) {
                // then it is worth figuring out what N is.
                switch (N) {
                  case 0 :
                       // do nothing
                       break;
                  case 1 :
                       MultMM_Helper<-4,cs,1,xs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 2 :
                       MultMM_Helper<-4,cs,2,xs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 3 :
                       MultMM_Helper<-4,cs,3,xs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                }
            } else 
#endif
                MultMM_Helper<71,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
        }
    };

    // algo 77: Unknown cs, check if K is small
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<77,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef TMV_MM_OPT_SMALL
            TMVStaticAssert(xs == TMV_UNKNOWN);
            const int K = m1.rowsize();
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"MM algo 77: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            if (K <= 3) {
                // then it is worth figuring out what K is.
                switch (K) {
                  case 0 :
                       Maybe<!add>::zero(m3);
                       break;
                  case 1 :
                       MultMM_Helper<-4,cs,rs,1,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 2 :
                       MultMM_Helper<-4,cs,rs,2,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 3 :
                       MultMM_Helper<-4,cs,rs,3,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                }
            } else 
#endif
                MultMM_Helper<71,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
        }
    };

    template <int ix, class T, class M> class ProdXM;

    // algo 81: copy x*m1
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<81,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs == TMV_UNKNOWN ? m3.colsize() : cs;
            const int K = xs == TMV_UNKNOWN ? m1.rowsize() : xs;
            const int N = rs == TMV_UNKNOWN ? m3.rowsize() : rs;
            std::cout<<"MM algo 81: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            typedef typename Traits2<T,T1>::type PT1;
            const int A = M2::_colmajor || M3::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT1,Rec,cs,xs,A>::type M1c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<add>(one,M1c(ProdXM<ix,T,M1>(x,m1)),m2,m3na);
        }
    };

    // algo 82: copy x*m2
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<82,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs == TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs == TMV_UNKNOWN ? m1.rowsize() : rs;
            std::cout<<"MM algo 82: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::value_type T2;
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            typedef typename Traits2<T,T2>::type PT2;
            const int A = M1::_colmajor && M3::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT2,Rec,xs,rs,A>::type M2c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultMM<add>(one,m1,M2c(ProdXM<ix,T,M2>(x,m2)),m3na);
        }
    };

    // algo 83: Use temporary for m1*m2
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<83,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_MM
            const int M = cs == TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs == TMV_UNKNOWN ? m1.rowsize() : rs;
            std::cout<<"MM algo 83: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<
                ','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const int A = M1::_rowmajor && M2::_rowmajor ? RowMajor : ColMajor;
            typedef typename MCopyHelper<PT3,Rec,cs,rs,A>::type M3c;
            typename M3::noalias_type m3na = m3.noAlias();
            MultXM<add>(x,M3c(m1*m2),m3na);
        }
    };

    // algo 90: call inst
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<90,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<90,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 91: call inst alias
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<91,cs,rs,xs,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<91,cs,rs,xs,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<97,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultMM_Helper<-2,cs,rs,xs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<197,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultMM_Helper<99,cs,rs,xs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<98,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_MM
            std::cout<<"MM algo 98: \n";
            std::cout<<"s1 = "<<s1<<std::endl;
            std::cout<<"s2 = "<<s2<<std::endl;
#endif
            if (!s1 && !s2) {
                // No aliasing
                MultMM_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (s1 && !s2) {
                // Use temporary for m1
                MultMM_Helper<81,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (!s1 && s2) {
                // Use temporary for m2
                MultMM_Helper<82,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // Use temporary for m1*m2
                MultMM_Helper<83,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<99,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                (xs == TMV_UNKNOWN || xs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                M3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<-4,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool ccc = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool rcc = M1::_rowmajor && M2::_colmajor && M3::_colmajor;
            const bool crc = M1::_colmajor && M2::_rowmajor && M3::_colmajor;
            const bool rrc = M1::_rowmajor && M2::_rowmajor && M3::_colmajor;

            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 2 :
                rs == 1 ? 3 :
                xs == 1 ? 4 :
                rs > 1 && rs <= 3 ? 21 :
                cs > 1 && cs <= 3 ? 405 :
                xs > 1 && xs <= 3 ? 31 :
                M3::_rowmajor && !M3::_colmajor ? 405 :
                ccc ? 11 :
                rcc ? 21 :
                crc ? 31 :
                ( M2::_stepi != TMV_UNKNOWN && M3::_stepi != TMV_UNKNOWN ) ? 21 :
                rrc ? 31 : 
                21;
            MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<-3,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;

            // Possible algorithms are:
            //
            // Trivial:
            //  0 = cs or rs == 0, or xs == 0 && add, so nothing to do
            //  1 = xs == 0 && !add, so m3.setZero();
            //  2 = cs == 1: reduces to trivial MultMV function
            //  3 = rs == 1: reduces to trivial MultMV function
            //  4 = xs == 1: reduces to trivial Rank1Update function
            //  5 = transpose
            //
            // CCC: 
            // 11 = CCC, loop over n
            // 12 = CCC, 4x2 blocks in KxN
            //
            // RCC:
            // 21 = RCC, loop over n
            // 22 = RCC, 2x2 blocks in MxN
            //
            // CRC:
            // 31 = CRC, loop over k
            // 32 = CRC, rank-4 updates 
            //
            // Fully Unrolled:
            // I haven't done these, but I don't that there is a 
            // need for them. 
            // The current setup is about as fast or faster (depending on
            // ccc,rcc,etc.) than Eigen's known-size calculation even down to 
            // 2x2 matrices, so I'm not sure that unrolling would speed 
            // things up much if at all.
            //
            // Main drivers for large matrices:
            // 61 = Recurse down to single blocks with divide-and-conquer.
            // 62 = Simpler loop over M,N,then K in blocks.
            // 63 = Same pattern as 61, but copy m1,m2 into a block structure
            // 64 = Same pattern as 62, but copy m1,m2 into a block structure
            // 66 = A simple division into 4 submatrices, 
            //      then call another algo.
            // 68 = Alberto and Nicolau's "hybrid Winograd" algorithm
            // 69 = Parallelize using openmp
            //
            // Runtime algorithm selection
            // 71 = Choose which algorithm based on the (runtime) size
            // 72 = Same as 71, but only large-matrix algorithms
            // 75 = Check if M is small
            // 76 = Check if N is small
            // 77 = Check if K is small
            //
            // Copy matrices to new storage
            // 81 = copy x*m1
            // 82 = copy x*m2
            // 83 = temp m1*m2

            const bool ccc = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool rcc = M1::_rowmajor && M2::_colmajor && M3::_colmajor;
            const bool crc = M1::_colmajor && M2::_rowmajor && M3::_colmajor;
            const bool rrc = M1::_rowmajor && M2::_rowmajor && M3::_colmajor;

            const int Mb = cs == TMV_UNKNOWN ? TMV_UNKNOWN : (cs >> 6);
            const int Nb = rs == TMV_UNKNOWN ? TMV_UNKNOWN : (rs >> 6);
            const int Kb = xs == TMV_UNKNOWN ? TMV_UNKNOWN : (xs >> 6);
            const int Mc = cs == TMV_UNKNOWN ? TMV_UNKNOWN : (cs < 16) ? 1 : (cs>>4);
            const int Nc = rs == TMV_UNKNOWN ? TMV_UNKNOWN : (rs < 16) ? 1 : (rs>>4);
            const int Kc = xs == TMV_UNKNOWN ? TMV_UNKNOWN : (xs < 16) ? 1 : (xs>>4);
            const int McNcKc = IntTraits2<IntTraits2<Mc,Nc>::prod,Kc>::prod;
            const bool twobig = (Mb && Kb) || (Mb && Nb) || (Nb && Kb);


#ifdef _OPENMP
            // This is different from the test we do in algo 71.
            // I'm having trouble with OPENMP for some small known values 
            // of cs,rs,xs.  e.g. 5000x5x5 leads to a compiler error with 
            // g++ 4.3.  So for now, I'm just disabling OPENMP for small 
            // sizes. But really, I should try to figure out what is 
            // going on here.
            //const bool do_openmp = 
            //    (Mb || Nb) && (McNcKc >= TMV_MM_OPENMP_THRESH);
            const bool do_openmp = 
                (Mb && Nb && Kb) && (McNcKc >= TMV_MM_OPENMP_THRESH);
#else
            const bool do_openmp = false;
#endif

#ifdef TMV_MM_USE_WINOGRAD
            const bool do_winograd =
                cs >= TMV_MM_MIN_WINOGRAD && rs >= TMV_MM_MIN_WINOGRAD && 
                xs >= TMV_MM_MIN_WINOGRAD;
#else 
            const bool do_winograd = false;
#endif

#ifdef TMV_MM_USE_RECURSIVE_BLOCK
            const int Kb2 = IntTraits2<Kb,Kb>::prod;
            const int MbNbKb2 = IntTraits2<IntTraits2<Mb,Nb>::prod,Kb2>::prod;
            const bool do_recursive = MbNbKb2 >= TMV_MM_MIN_RECURSIVE;
#else
            const bool do_recursive = false;
#endif
            const bool do_block = 
                (cs >= 16 && rs >= 16 && xs >= 16) ||
                ( twobig && McNcKc >= TMV_MM_MIN_COPY_ALGO );

            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 202 :
                rs == 1 ? 203 :
                xs == 1 ? 204 :
                rs > 1 && rs <= 3 ? 21 :
                cs > 1 && cs <= 3 ? 5 :
                xs > 1 && xs <= 3 ? 31 :
                M3::_rowmajor && !M3::_colmajor ? 5 :
                TMV_OPT == 0 ? ( ccc ? 11 : rcc ? 21 : crc ? 31 : 21 ) :
                // If multiplying Matrix<int> * Matrix<double>
                // (or worse, complex versions of these)
                // then the sophisticated algorithms sometimes have 
                // trouble, so just do the simple naive loops.
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ?
                ( ccc ? 11 : rcc ? 21 : crc ? 31 : 21 ) :
                (cs == TMV_UNKNOWN || rs == TMV_UNKNOWN || xs == TMV_UNKNOWN ) ? (
                    ccc ? ( cs == TMV_UNKNOWN ? 75 : cs <= 3 ? 11 : 71 ) :
                    rcc ? ( xs == TMV_UNKNOWN ? 77 : xs <= 3 ? 21 : 71 ) :
                    crc ? ( cs == TMV_UNKNOWN ? 75 : cs <= 3 ? 31 : 71 ) :
                    71 ) :
                do_openmp ? 69 :
                do_winograd ? 68 :
                do_recursive ? 63 :
                do_block ? 64 :
                // For known sizes, the MultMV calls (algo 11,21) seem 
                // to be faster than the 2x2 algorithms (algo 12,22).  
                // For unknown sizes (cf. selection in algo 71) the 
                // opposite is true.
                ccc ? 11 :
                rcc ? 21 :
                crc ? 31 :
                // When steps are known, direct MultMV is faster than copying.
                // Basically, the non-unit step vectors for MultMV are ok, 
                // since the non-unit step is known and the lengths of the 
                // vectors are not too large.
                ( M2::_stepi != TMV_UNKNOWN && M3::_stepi != TMV_UNKNOWN )  ? 21 :
                // For rrc, 31 seems to be faster than copying for small
                // matrices (which is what we have here).
                rrc ? 31 : 
                !M3::_colmajor ? 83 :
                !(M2::_rowmajor || M2::_colmajor) ? 82 :
                !(M1::_rowmajor || M1::_colmajor) ? 81 :
                // The above should have caught all possibilities.
                // So selection should never get here.
                -999;
            TMVStaticAssert(algo != -999);
#ifdef TMV_MM_OPT_BAD_ALLOC
            // Nothing less than 60 should throw a bad_alloc
            const int algo1 = algo < 60 ? 0 : 66;
#else
            const int algo1 = algo < 60 ? 0 : 67;
#endif
#ifdef PRINTALGO_MM
            const int M = cs==TMV_UNKNOWN ? m3.colsize() : cs;
            const int N = rs==TMV_UNKNOWN ? m3.rowsize() : rs;
            const int K = xs==TMV_UNKNOWN ? m1.rowsize() : xs;
            std::cout<<"InlineMultMM: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"M = "<<M<<"  N = "<<N<<"  K = "<<K<<std::endl;
#ifdef _OPENMP
            std::cout<<"69 ? McNcKc = "<<McNcKc<<
                " >= "<<TMV_MM_OPENMP_THRESH<<std::endl;
#endif
#ifdef TMV_MM_USE_WINOGRAD
            std::cout<<"68 ? all cs,rs,xs >= "<<
                TMV_MM_MIN_WINOGRAD<<std::endl;
#endif
#ifdef TMV_MM_USE_RECURSIVE_BLOCK
            std::cout<<"63 ? MbNbKb2 = "<<MbNbKb2<<
                " >= "<<TMV_MM_MIN_RECURSIVE<<std::endl;
#endif
            std::cout<<"64 ? McNcKc = "<<McNcKc<<
                " >= "<<TMV_MM_MIN_COPY_ALGO<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<"  xs = "<<xs<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_MM
            typedef typename M3::real_type RT;
            typedef typename M3::value_type T3;
            Matrix<T3> m1c = m1;
            Matrix<T3> m2c = m2;
            Matrix<T3> m3i = m3;
            Matrix<T3> m3c = m3;
            for(int j=0; j<m3.colsize(); ++j) {
                typename Matrix<T3>::col_type m3cj = m3c.col(j);
                MultMV<add>(x,m1c,m2c.col(j),m3cj);
            }
#endif
            try {
                MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            } catch (std::bad_alloc) {
                MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
                    x,m1,m2,m3);
            }
#ifdef XDEBUG_MM
            if (Norm(m3-m3c) > 1.e-3*(Norm(m1c)*Norm(m2c)+(add?Norm(m3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"m2 = "<<m2c<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 => "<<m3<<std::endl;
                std::cout<<"Correct m3 = "<<m3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -402: Same as algo -2, but use -4 if no Inst version
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<-402,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                (xs == TMV_UNKNOWN || xs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 2 :
                rs == 1 ? 3 :
                xs == 1 ? 4 :
                inst ? 90 : 
                -4;
            MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                (xs == TMV_UNKNOWN || xs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 202 :
                rs == 1 ? 203 :
                xs == 1 ? 204 :
                M3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMM_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
                ( xs == 0 && !add ) ? 1 :
                cs == 1 ? 102 :
                rs == 1 ? 103 :
                xs == 1 ? 104 :
                M3::_checkalias ? 99 : 
                -2;
            MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        // cs = m3.colsize
        // rs = m3.rowsize
        // xs = "extra size" = m1.rowsize, m2.colsize
        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultMM_Helper<-1,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultMM_Helper<-3,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1,
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs = Sizes<M3::_rowsize,M2::_rowsize>::size;
        const int xs = Sizes<M1::_rowsize,M2::_colsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultMM_Helper<98,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <class M1, int ix, class T, class M2>
    TMV_INLINE void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Rec<M2>& m2)
    { MultMM<false>(x,m1.copy(),m2.mat(),m1.mat()); }

} // namespace tmv

#endif 
