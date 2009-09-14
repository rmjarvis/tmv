///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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


#ifndef TMV_MultMM_H
#define TMV_MultMM_H

#include "TMV_MultMV.h"
#include "TMV_Rank1_VVM.h"
#include "TMV_Matrix.h"
#include "TMV_SmallMatrix.h"
#include "TMV_AddMM.h"
#include "TMV_MultXM.h"

#ifdef _OPENMP
#include "omp.h"
#ifdef PRINTALGO_MM_OMP
#include <fstream>
#endif
#endif

#ifdef PRINTALGO_MM
#include <iostream>
#endif

// Check for small (<=3) values of cs, rs, or xs
// This leads to significant speed improvements for such matrices at little
// cost to larger matrices, but the large increase in code size and the fact
// that it only benefits a few particular sizes of matrices mean that
// we require TMV_OPT = 3 for it.
#if TMV_OPT >= 3
#define TMV_OPT_SMALL
#endif

// Select the cleanup code for the edge that doesn't fit into KB-sized
// blocks according to the exact value of K.  This can have a significant 
// increase in the speed for moderately sized matrices where the edges 
// are not a negligible fraction of the calculation. So we only 
// require TMV_OPT = 2
#if TMV_OPT >= 2
#define TMV_OPT_CLEANUP
#endif

// For very large matrices, a recursive block algorithm is faster than
// the simpler looping block algorithm.  The dividing line is governed
// by the parameter TMV_Q5 (see below).  It is not all that expensive 
// with respect to code bloat, and it is significantly faster for large 
// matrices, so we only require TMV_OPT = 2.
#if TMV_OPT >= 2
#define TMV_USE_RECURSIVE_BLOCK
#endif

// For extremely large matrices, it is worth using a recursive Winograd
// algorithm where only 7 matrix multiplies are needed to do the 8
// submatrix calculations in:
// [ A B ] * [ E F ] = [ I J ]
// [ C D ]   [ G H ]   [ K L ]
// This trick is repeated recursively until the minimum size is less than
// TMV_Q4 (see below), and then we call a more standard algorithm.
// It is only used for extremely large matrices, but it isn't too expensive
// for the code size, so we require TMV_OPT = 2.
#if TMV_OPT >= 2
#define TMV_USE_WINOGRAD
#endif

// The algorithms for large matrices usually involve making temporary
// matrices with the data stored in blocks.  Therefore, to be good
// C++ programmers, we put these into try/catch blocks to catch any
// possible bad_alloc throws.  We have an algorithm specially designed
// for the post-catch calculation, which divides M,N in half and does
// the four sections separately.  For some reason, this is pretty expensive
// with respect to code bloat, although I can't figure out why.  So
// since bad_alloc's are pretty rare, until I figure that out, I require
// TMV_OPT=3.
#if TMV_OPT >= 3
#define TMV_OPT_BAD_ALLOC
#endif

// There are a number of values used in the algorithm selection
// that are either arbitrary or empirical.
// So I put them all here to make them easier to change and to  
// track down in the code.

// Q1 is the maximum value of cs*rs for which a SmallMatrix * SmallVector
// product is fully unrolled. 
// This isn't actually implemented yet.  The non-unrolling versions
// are already pretty fast, probably because the compiler unrolls the
// loops well enough, but at some point I will add manual unrolling
// where this parameter will be relevant.
#if TMV_OPT >= 3
#define TMV_Q1 1024 // It never actually gets this high.
#elif TMV_OPT >= 2
#define TMV_Q1 25
#elif TMV_OPT >= 1
#define TMV_Q1 9
#else
#define TMV_Q1 0
#endif

// Q2 is the block size to use for large matrices.  The value should
// be chosen such that three Q2xQ2 sized matrices all fit into the 
// L1 cache.  
// For algo 63, I have special block sizes chosen rather than Q2,
// since I find that for this algorithm, it is better to use smaller
// blocks.  Also there are special choices when using SSE.
// So Q2 is just for the regular recursion when not copying (algo 61).
#define TMV_Q2 64

// Q3 is the crossover memory size to start using prefetch commands.
// This is undoubtedly a function of the L1 (and L2?) cache size,
// but 2KBytes is probably not too bad for most machines.
// (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_Q3 2048

// Q4 is the minimum size to use a recursive Winograd algorithm
// Note: This needs to be duplicated in TMV_MultMM_Winograd.h
#ifdef TMV_USE_RECURSIVE_BLOCK
#define TMV_Q4 2048
#else
#define TMV_Q4 1024
#endif

// Q5 is the minimum value of Mb*Nb*Kb*Kb to use recursive algorithm.
// This formula for the crossover between algorithms 63 and 64 is
// purely empirical on my MacBook (with an Intel Core 2 Duo processor),
// so it might not even be the right thing to parametrize on other
// machines.  On the other hand, the difference between the two algorithms
// isn't _that_ extreme.  At most I've seen about a factor of 2, so
// if I use the less optimal algorithm for some sizes, it's not that
// terrible a mistake.
#define TMV_Q5 16*1024

// Q6 is the minimum value of (M*N*K / 16^3) to use multiple threads.
// (We also require that M or N > 64 so we have something to split.)
// The computer that I did most of the testing on only has 2 processors,
// so this optimization might be particular to having 2 threads.
// TODO: Investigate the timings when more than 2 threads are available.
#define TMV_Q6 64
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
//  64 x  64 x  64    0.655    0.679    1.036      64
// 256 x  16 x  16    1.186    1.186    1.000      16
//  64 x  16 x  64    1.145    1.096    0.957      16
//  16 x  64 x  64    1.146    1.051    0.917      16

// ^^^^ My formula says don't use omp
// vvvv My formula says use omp

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

// Q7 is the minimum value of (M*N*K / 16^3) to use copying algorithm
// where each block is copied to temporary storage in a block structure.
#define TMV_Q7 4

// ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_ZeroIX (ix==0)
//#define TMV_ZeroIX (ix!=1)

// We split off some of the parts of the calculation to other files.
// First, this helps us compile them into separate object files, so the
// compilation doesn't take forever on a single file.  And second,
// this header file would otherwise be over 6000 lines(!) long.

#include "TMV_MultMM_Kernel.h"
#include "TMV_MultMM_Winograd.h"
#include "TMV_MultMM_Block.h"

namespace tmv {

  //
  // Matrix * Matrix
  //

  template <int algo, int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3> struct MultMM_Helper;

  // algo 0: Trivial, nothing to do.
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<0,cs,rs,xs,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {}
  };

  // algo 1: Trivial, just m3.Zero();
  template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<1,cs,rs,0,false,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { m3.Zero(); }
  };

  // algo 2: cs == 1, so reduces to MultMV
  template <int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<2,1,rs,xs,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 2: M,N,K,cs,rs,xs,x = "<<1<<','<<N<<','<<K<<','<<1<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename M1::const_row_type M1r;
      typedef typename M2::const_transpose_type M2t;
      typedef typename M3::row_type M3r;
      M1r m1r = m1.get_row(0);
      M2t m2t = m2.Transpose();
      M3r m3r = m3.get_row(0);
      MultMV_Helper<-3,rs,xs,add,ix,T,M2t,M1r,M3r>::call(x,m2t,m1r,m3r);
    }
  };

  // algo 2: rs == 1, so reduces to MultMV
  template <int cs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<3,cs,1,xs,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 3: M,N,K,cs,rs,xs,x = "<<M<<','<<1<<','<<K<<','<<cs<<','<<1<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename M2::const_col_type M2c;
      typedef typename M3::col_type M3c;
      M2c m2c = m2.get_col(0);
      M3c m3c = m3.get_col(0);
      MultMV_Helper<-3,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
    }
  };

  // algo 4: xs == 1, so reduces to Rank1Update
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<4,cs,rs,1,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      std::cout<<"algo 4: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<1<<','<<cs<<','<<rs<<','<<1<<','<<T(x)<<std::endl;
#endif
      typedef typename M1::const_col_type M1c;
      typedef typename M2::const_row_type M2r;
      M1c m1c = m1.get_col(0);
      M2r m2r = m2.get_row(0);
      Rank1Update_VVM_Helper<-3,cs,rs,add,ix,T,M1c,M2r,M3>::call(x,m1c,m2r,m3);
    }
  };

  // algo 5: call InstMultMM
  template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<5,cs,rs,xs,false,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 5 !add: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename M3::value_type T3;
      typename M3::xview_type m3x = m3.XView();
      InstMultMM(T3(x),m1.XView(),m2.XView(),m3x); 
    }
  };
  template <int cs, int rs, int xs, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<5,cs,rs,xs,true,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 5 add: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename M3::value_type T3;
      typename M3::xview_type m3x = m3.XView();
      InstAddMultMM(T3(x),m1.XView(),m2.XView(),m3x); 
    }
  };

  // algo 11: CCC -- Loop over N
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<11,cs,rs,xs,add,ix,T,M1,M2,M3>
  {
    static void call(
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

    static void loop_42(const int M, const int N, const int K,
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
      IT1 A1 = A0; A1.ShiftP(Astepj);
      IT1 A2 = A1; A2.ShiftP(Astepj);
      IT1 A3 = A2; A3.ShiftP(Astepj);

      IT2 B0 = m2.get_col(0).begin();
      IT2 B1 = B0; B1.ShiftP(Bstepj);

      IT3 C0 = m3.get_col(0).begin();
      IT3 C1 = C0; C1.ShiftP(Cstepj);

      const bool dopref = K * sizeof(T1) >= TMV_Q3;

      Prefetch_MultiRead(A0.GetP());
      Prefetch_MultiRead(A1.GetP());
      Prefetch_Read(B0.GetP());
      Prefetch_Read(B1.GetP());
      Prefetch_MultiWrite(C0.GetP());
      Prefetch_MultiWrite(C1.GetP());

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
          A0.ShiftP(Astepj_4); A1.ShiftP(Astepj_4);
          A2.ShiftP(Astepj_4); A3.ShiftP(Astepj_4);
          C0 -= M; C1 -= M; 
          if (dopref) {
            Prefetch_Read(A0.GetP());
            Prefetch_Read(A1.GetP());
            Prefetch_Read(A2.GetP());
            Prefetch_Read(A3.GetP());
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
          A0.ShiftP(Astepj_1);
          C0 -= M; C1 -= M;
          if (dopref) {
            Prefetch_Read(A0.GetP());
          }
        } while (--k);
        A0.ShiftP(Astart); A1.ShiftP(Astartx); 
        A2.ShiftP(Astartx); A3.ShiftP(Astartx);
        B0.ShiftP(Bstepj_2); B1.ShiftP(Bstepj_2);
        C0.ShiftP(Cstepj_2); C1.ShiftP(Cstepj_2);
        if (dopref) {
          Prefetch_Read(B0.GetP());
          Prefetch_Read(B1.GetP());
          Prefetch_MultiWrite(C0.GetP());
          Prefetch_MultiWrite(C1.GetP());
        }
      } while (--j);
    }
    template <int algo2, int z> struct Helper2;
    template <int z>
    struct Helper2<0,z> // algo2 = 0 (cs,rs,xs=0)
    {
      static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
    };
    template <int z>
    struct Helper2<2,z> // algo2 = 2 (rs is even)
    {
      static void call(
          const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
      { 
        TMVStaticAssert(rs != UNKNOWN);
        TMVStaticAssert(rs != 0);
        TMVStaticAssert(rs%2 == 0);
        const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
        const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
        if (M && K) loop_42(M,rs,K,x,m1,m2,m3); 
      }
    };
    template <int z>
    struct Helper2<3,z> // algo2 = 3 (normal operation)
    {
      static void call(
          const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
      {
        const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
        const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
        const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
        typedef typename M1::const_col_type M1c;
        typedef typename M2::const_col_type M2c;
        typedef typename M2::const_row_range_type M2r;
        typedef typename M3::col_type M3c;
        typedef typename M3::cols_type M3a;

        if (M && N && K) {
          const int na = ((N>>1)<<1);
          const int nb = N-na;
          if (na) {
            loop_42(M,na,K,x,m1,m2,m3);
          }
          if (nb) {
            M2c m2na = m2.get_col(na);
            M3c m3na = m3.get_col(na);
            MultMV_Helper<-2,cs,xs,true,ix,T,M1,M2c,M3c>::call(x,m1,m2na,m3na);
          }
        }
      }
    };
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      enum { algo2 = (
          ( cs == 0 || rs == 0 || xs == 0 ) ? 0 :
          ( rs == UNKNOWN ) ? 3 :
          ( rs%2 == 0 ) ? 2 : 3 ) };
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 12: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
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
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 21: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      for(int j=0;j<N;++j) {
        // m3.get_col(j) += x * m1 * m2.get_col(j)
        typedef typename M2::const_col_type M2c;
        M2c m2j = m2.get_col(j);
        typedef typename M3::col_type M3c;
        M3c m3j = m3.get_col(j);
        MultMV_Helper<-2,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2j,m3j);
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

    static void loop_22(const int M, const int N, const int K,
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
      IT1 A1 = A0; A1.ShiftP(Astepi);

      IT2 B0 = m2.get_col(0).begin();
      IT2 B1 = B0; B1.ShiftP(Bstepj);

      IT3 C0 = m3.get_col(0).begin();
      IT3 C1 = C0; C1.ShiftP(Cstepj);

      const bool dopref = K * sizeof(T1) >= TMV_Q3;

      Prefetch_MultiRead(A0.GetP());
      Prefetch_MultiRead(A1.GetP());
      Prefetch_MultiRead(B0.GetP());
      Prefetch_MultiRead(B1.GetP());
      Prefetch_Write(C0.GetP());
      Prefetch_Write(C1.GetP());

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
          A0.ShiftP(Astepi_2); A1.ShiftP(Astepi_2);
          B0 -= K; B1 -= K;
          if (dopref) {
            Prefetch_Read(A0.GetP());
            Prefetch_Read(A1.GetP());
          }
        } while (--i);
        A0.ShiftP(Astart); A1.ShiftP(Astart);
        B0.ShiftP(Bstepj_2); B1.ShiftP(Bstepj_2);
        C0.ShiftP(Cstepj_2); C1.ShiftP(Cstepj_2);
        if (dopref) {
          Prefetch_MultiRead(B0.GetP());
          Prefetch_MultiRead(B1.GetP());
          Prefetch_Write(C0.GetP());
          Prefetch_Write(C1.GetP());
        }
      } while (--j);
    }
    template <int algo2, int z> struct Helper2;
    template <int z>
    struct Helper2<0,z> // algo2 = 0 (cs,rs=0 or (xs=0 && add))
    {
      static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
    };
    template <int z>
    struct Helper2<1,z> // algo2 = 1 (xs=0 and !add)
    {
      static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& m3) 
      { m3.Zero(); }
    };
    template <int z>
    struct Helper2<2,z> // algo2 = 2 (cs and rs are even)
    {
      static void call(
          const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
      { 
        TMVStaticAssert(cs != UNKNOWN);
        TMVStaticAssert(rs != UNKNOWN);
        TMVStaticAssert(cs != 0);
        TMVStaticAssert(rs != 0);
        TMVStaticAssert(cs%2 == 0);
        TMVStaticAssert(rs%2 == 0);
        const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
        if (K) loop_22(cs,rs,K,x,m1,m2,m3); 
        else if (!add) m3.Zero();
      }
    };
    template <int z>
    struct Helper2<3,z> // algo2 = 3 (normal operation)
    {
      static void call(
          const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
      {
        const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
        const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
        const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
        typedef typename M1::const_row_type M1r;
        typedef typename M2::const_col_type M2c;
        typedef typename M2::const_cols_type::const_transpose_type M2t;
        typedef typename M3::col_type M3c;
        typedef typename M3::row_range_type M3r;

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
            MultMV_Helper<-2,cs,xs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2na,m3na);
          }
          if (mb) {
            M1r m1ma = m1.get_row(ma);
            M2t m2t = m2.CCols(0,na).Transpose();
            M3r m3ma = m3.get_row(ma,0,na);
            enum { rsx = (rs == UNKNOWN ? UNKNOWN : ((rs>>1)<<1)) };
            MultMV_Helper<-2,rsx,xs,add,ix,T,M2t,M1r,M3r>::call(
                x,m2t,m1ma,m3ma);
          }
        } else if (M && N && !add) { // K == 0 && !add
          m3.Zero();
        }
      }
    };
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      enum { algo2 = (
          ( cs == 0 || rs == 0 ) ? 0 :
          ( xs == 0 ) ? ( add ? 0 : 1 ) :
          ( cs == UNKNOWN || rs == UNKNOWN ) ? 3 :
          ( cs%2 == 0 && rs%2 == 0 ) ? 2 : 3 ) };
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 22: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
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
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      std::cout<<"algo 31: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      Maybe<!add>::zero(m3);
      for(int k=0;k<K;++k) {
        // m3 += x * m1.get_col(k) * m2.get_row(j)
        typedef typename M1::const_col_type M1c;
        M1c m1k = m1.get_col(k);
        typedef typename M2::const_row_type M2r;
        M2r m2k = m2.get_row(k);
        Rank1Update_VVM_Helper<-2,cs,rs,true,ix,T,M1c,M2r,M3>::call(
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

    static void loop_42(const int M, const int N, const int K,
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
      IT1 A1 = A0; A1.ShiftP(Astepj);
      IT1 A2 = A1; A2.ShiftP(Astepj);
      IT1 A3 = A2; A3.ShiftP(Astepj);

      IT2 B0 = m2.get_row(0).begin();
      IT2 B1 = B0; B1.ShiftP(Bstepi);
      IT2 B2 = B1; B2.ShiftP(Bstepi);
      IT2 B3 = B2; B3.ShiftP(Bstepi);

      IT3 C0 = m3.get_col(0).begin();
      IT3 C1 = C0; C1.ShiftP(Cstepj);

      const bool dopref = K * sizeof(T1) >= TMV_Q3;

      Prefetch_Read(A0.GetP());
      Prefetch_Read(A1.GetP());
      Prefetch_MultiRead(B0.GetP());
      Prefetch_MultiRead(B1.GetP());
      Prefetch_MultiWrite(C0.GetP());
      Prefetch_MultiWrite(C1.GetP());

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
          C0.ShiftP(Cstepj_2); C1.ShiftP(Cstepj_2);
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
          C0.ShiftP(Cstepj_1);
          A0 -= M; A1 -= M; A2 -= M; A3 -= M; 
        }
        C0.ShiftP(Cstart); C1.ShiftP(Cstartx);
        A0.ShiftP(Astepj_4); A1.ShiftP(Astepj_4); 
        A2.ShiftP(Astepj_4); A3.ShiftP(Astepj_4); 
        B0.ShiftP(Bstepi_4); B1.ShiftP(Bstepi_4);
        B2.ShiftP(Bstepi_4); B3.ShiftP(Bstepi_4);
        if (dopref) {
          Prefetch_Read(A0.GetP());
          Prefetch_Read(A1.GetP());
          Prefetch_Read(B0.GetP());
          Prefetch_Read(B1.GetP());
        }
      } while (--k);
    }
    template <int algo2, int z> struct Helper2;
    template <int z>
    struct Helper2<0,z> // algo2 = 0 (cs,rs,xs=0)
    {
      static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
    };
    template <int z>
    struct Helper2<2,z> // algo2 = 2 (rs is even)
    {
      static void call(
          const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
      { 
        TMVStaticAssert(xs != UNKNOWN);
        TMVStaticAssert(xs != 0);
        TMVStaticAssert(xs%4 == 0);
        const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
        const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
        if (M && N) loop_42(M,N,xs,x,m1,m2,m3); 
      }
    };
    template <int z>
    struct Helper2<3,z> // algo2 = 3 (normal operation)
    {
      static void call(
          const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
      {
        const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
        const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
        const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
        typedef typename M1::const_cols_type M1c;
        typedef typename M2::const_rows_type M2r;

        if (M && N && K) {
          const int ka = ((K>>2)<<2);
          const int kb = K-ka;
          if (ka) {
            loop_42(M,N,ka,x,m1,m2,m3);
          }
          if (kb) {
            M1c m1c = m1.CCols(ka,K);
            M2r m2r = m2.CRows(ka,K);
            enum { xsx = xs == UNKNOWN ? UNKNOWN : (xs - ((xs>>2)<<2)) };
            MultMM_Helper<31,cs,rs,xsx,true,ix,T,M1c,M2r,M3>::call(
                x,m1c,m2r,m3);
          }
        }
      }
    };
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      enum { algo2 = (
          ( cs == 0 || rs == 0 || xs == 0 ) ? 0 :
          ( xs == UNKNOWN ) ? 3 :
          ( xs%4 == 0 ) ? 2 : 3 ) };
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 32: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
      std::cout<<"add = "<<add<<", algo2 = "<<algo2<<std::endl;
#endif
      Maybe<!add>::zero(m3);

      Helper2<algo2,1>::call(x,m1,m2,m3);
    }
  };

  // algo 51: RRC, No good algorithm here. Copy m2 to column major
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<51,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 51: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      MultMM_Helper<55,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3); 
    }
  };

  // algo 52: **R, Transpose to get the colmajor version
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<52,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 52: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename M1::const_transpose_type M1t;
      typedef typename M2::const_transpose_type M2t;
      typedef typename M3::transpose_type M3t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      M3t m3t = m3.Transpose();
      MultMM_Helper<-1,rs,cs,xs,add,ix,T,M2t,M1t,M3t>(x,m2t,m1t,m3t); 
    }
  };

  // algo 53: m1 is no major -- use temporary (R)
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<53,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 53: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename Traits<T>::real_type RT;
      const Scaling<1,RT> one;
      typedef typename M1::value_type T1;
      typedef typename Traits2<T,T1>::type PT1;
      typedef Matrix<PT1,RowMajor> M1c;
      typedef typename M1c::const_view_type M1cv;
      M1c m1c = x*m1;
      MultMM_Helper<-2,cs,rs,xs,add,1,RT,M1cv,M2,M3>::call(
          one,m1c.View(),m2,m3); 
    }
  };

  // algo 54: m1 is no major -- use temporary (C)
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<54,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 54: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename Traits<T>::real_type RT;
      const Scaling<1,RT> one;
      typedef typename M1::value_type T1;
      typedef typename Traits2<T,T1>::type PT1;
      typedef Matrix<PT1,ColMajor> M1c;
      typedef typename M1c::const_view_type M1cv;
      M1c m1c = x*m1;
      MultMM_Helper<-2,cs,rs,xs,add,1,RT,M1cv,M2,M3>::call(
          one,m1c.View(),m2,m3); 
    }
  };

  // algo 55: m2 is no major -- use temporary (C)
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<55,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 55: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename Traits<T>::real_type RT;
      const Scaling<1,RT> one;
      typedef typename M2::value_type T2;
      typedef typename Traits2<T,T2>::type PT2;
      typedef Matrix<PT2,ColMajor> M2c;
      typedef typename M2c::const_view_type M2cv;
      M2c m2c = x*m2;
      MultMM_Helper<-2,cs,rs,xs,add,1,RT,M1,M2cv,M3>::call(
          one,m1,m2c.View(),m3); 
    }
  };

  // algo 56: m3 is no major -- use temporary (C)
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<56,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_MM
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 56: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typedef typename Traits<T>::real_type RT;
      const Scaling<1,RT> one;
      typedef typename M1::value_type T1;
      typedef typename M2::value_type T2;
      typedef typename Traits2<T1,T2>::type PT3;
      typedef Matrix<PT3,ColMajor> M3c;
      typedef typename M3c::view_type M3cv;
      M3c m3c(M,N);
      M3cv m3cv = m3c.View();
      MultMM_Helper<-2,cs,rs,xs,false,1,RT,M1,M2,M3cv>::call(one,m1,m2,m3cv); 
      Maybe<add>::add(m3 , x * m3cv);
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
  // So we have a single parameter, TMV_Q2, to define what is a good choice
  // for the block size.  64 is a good choice for many systems, but if you
  // want to mess around with this number, you might be able to slightly
  // improve the speed of matrix multiplication with a different value.
  //
  // For algo 63 (see description below), we actually find that a smaller
  // block is better.  We hard code its blocks to be 16x16 (for MxN) and
  // either 16, 32 or 64 for K depending on the value type in the matrix.
  // This isn't a parameter, because we have special optimized routines for
  // multiplying these blocks that are specialized to these particular values.
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
    enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
    enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
    enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
    enum { MB    = ccc ? TMV_Q2 : rcc ? TMV_Q2 : TMV_Q2 };
    enum { NB    = ccc ? TMV_Q2 : rcc ? TMV_Q2 : TMV_Q2 };
    enum { KB    = ccc ? TMV_Q2 : rcc ? TMV_Q2 : TMV_Q2 };
    enum { algo2 = ( 
        ccc ? ( M3::miscomplex ? 11 : 12 ) :
        rcc ? ( M3::miscomplex ? 21 : 22 ) :
        ( M3::miscomplex ? 31 : 32 ) ) };
    enum { lnMB  = IntLog<MB>::log };
    enum { lnNB  = IntLog<NB>::log };
    enum { lnKB  = IntLog<KB>::log };

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
      else if (M >= N && M >= K)  // M is largest
      {
        TMVAssert(M >= N && M >= K);
        const int Mx = M>>1; // = M/2
        const int im = i1 + (Mx<<lnMB); // = i1 + Mx * MB
        call1<addx,MB,rsx,xsx>(i1,j1,k1,im,j2,k2,Mx,N,K,x,m1,m2,m3);
        call1<addx,csx,rsx,xsx>(im,j1,k1,i2,j2,k2,M-Mx,N,K,x,m1,m2,m3);
      }
      else if (N >= M && N >= K) // N is largest
      {
        TMVAssert(N >= M && N >= K);
        const int Nx = N>>1; // = N/2
        const int jm = j1 + (Nx<<lnNB); // = j1 + Nx * NB
        call1<addx,csx,NB,xsx>(i1,j1,k1,i2,jm,k2,M,Nx,K,x,m1,m2,m3);
        call1<addx,csx,rsx,xsx>(i1,jm,k1,i2,j2,k2,M,N-Nx,K,x,m1,m2,m3);
      }
      else // K is largest
      {
        TMVAssert(K >= M && K >= N);
        const int Kx = K>>1; // = K/2
        const int km = k1 + (Kx<<lnKB); // = k1 + Kx * KB
        call1<addx,csx,rsx,KB>(i1,j1,k1,i2,j2,km,M,N,Kx,x,m1,m2,m3);
        call1<true,csx,rsx,xsx>(i1,j1,km,i2,j2,k2,M,N,K-Kx,x,m1,m2,m3);
      }
    }
    template <bool addx, int csx, int rsx, int xsx>
    static void call2(
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
      TMVAssert(csx == UNKNOWN || csx == i2-i1);
      TMVAssert(rsx == UNKNOWN || rsx == j2-j1);
      TMVAssert(xsx == UNKNOWN || xsx == k2-k1);

      M1s m1s = m1.CSubMatrix(i1,i2,k1,k2);
      M2s m2s = m2.CSubMatrix(k1,k2,j1,j2);
      M3s m3s = m3.CSubMatrix(i1,i2,j1,j2);

      MultMM_Helper<algo2,csx,rsx,xsx,addx,ix,T,M1s,M2s,M3s>::call(
          x,m1s,m2s,m3s);
    }
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_MM
      std::cout<<"algo 61: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(ccc || rcc || crc);
      // These will be the correct value for the final non-block calls.
      enum { csx = cs == UNKNOWN ? UNKNOWN : (cs-((cs>>lnMB)<<lnMB)) };
      enum { rsx = rs == UNKNOWN ? UNKNOWN : (rs-((rs>>lnNB)<<lnNB)) };
      enum { xsx = xs == UNKNOWN ? UNKNOWN : (xs-((xs>>lnKB)<<lnKB)) };
      if (M > MB || N > NB || K > KB) {
        const int Mb = (M>>lnMB)+1; // = M/MB + 1
        const int Nb = (N>>lnNB)+1; // = N/NB + 1
        const int Kb = (K>>lnKB)+1; // = K/KB + 1
        call1<add,csx,rsx,xsx>(0,0,0,M,N,K,Mb,Nb,Kb,x,m1,m2,m3); 
      }
      else
        MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  // algo 63: Recursively subdivide m1,m2,m3 just like algo 61, but
  // use temporary memory to copy the blocks to compact storage.
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<63,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 63: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif

      typedef typename M1::value_type T1;
      typedef typename M2::value_type T2;
      typedef typename M3::value_type T3;
      enum { inst = (
#ifdef TMV_INST_MIX
          Traits2<T1,T3>::samebase &&
          Traits2<T2,T3>::samebase &&
#else
          Traits2<T1,T3>::sametype &&
          Traits2<T2,T3>::sametype &&
#endif
          Traits<T1>::isinst &&
          Traits<T2>::isinst &&
          Traits<T3>::isinst) };
      CallRecursiveBlockMultMM<inst,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
          x,m1,m2,m3);
    }
  };

  // algo 64: Similar to algo 63 in that we copy the blocks to temporary
  // storage.  However, simply loop over M,N and then do the full K
  // row/column block in one set.  This is faster than algo 63 for
  // smaller matrices.  The dividing line between the two algorithms is 
  // set by TMV_Q5.
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<64,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 64: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif

      typedef typename M1::value_type T1;
      typedef typename M2::value_type T2;
      typedef typename M3::value_type T3;
      enum { inst = (
#ifdef TMV_INST_MIX
          Traits2<T1,T3>::samebase &&
          Traits2<T2,T3>::samebase &&
#else
          Traits2<T1,T3>::sametype &&
          Traits2<T2,T3>::sametype &&
#endif
          Traits<T1>::isinst &&
          Traits<T2>::isinst &&
          Traits<T3>::isinst) };
      CallBlockMultMM<inst,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
          x,m1,m2,m3);
    }
  };

  // algo 66: Split problem into 4 parts then call another algorithm
  // to finish the work.  This is used after a bad_alloc error
  // to split the problem into chucks that will require less temporary
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
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_MM
      std::cout<<"algo 66: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      const int Mx = (((M-1)>>5)+1)<<4; // = M/2 rounded up to 16 block
      const int Nx = (((N-1)>>5)+1)<<4; // = N/2 rounded up to 16 block
      typedef typename M1::const_rows_type M1r;
      typedef typename M2::const_cols_type M2c;
      typedef typename M3::submatrix_type M3s;

      enum { Mb = cs == UNKNOWN ? UNKNOWN : (cs >> 6) };
      enum { Nb = rs == UNKNOWN ? UNKNOWN : (rs >> 6) };
      enum { Kb = xs == UNKNOWN ? UNKNOWN : (xs >> 6) };
      enum { Mc = cs == UNKNOWN ? UNKNOWN : (cs < 16) ? 1 : (cs>>4) };
      enum { Nc = rs == UNKNOWN ? UNKNOWN : (rs < 16) ? 1 : (rs>>4) };
      enum { Kc = xs == UNKNOWN ? UNKNOWN : (xs < 16) ? 1 : (xs>>4) };
      enum { MbNbKb2 = (
            IntTraits2<IntTraits2<Mb,Nb>::prod,IntTraits2<Kb,Kb>::prod>::prod 
            ) };
      enum { McNcKc = IntTraits2<IntTraits2<Mc,Nc>::prod,Kc>::prod };
      enum { twobig = (Mb && Kb) || (Mb && Nb) || (Nb && Kb) };
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      enum { nosimple = !(ccc || rcc || crc) };
      enum { algo2 = (
          ccc ? 11 : rcc ? 21 : crc ? 31 : rrc ? 31 :
          !M3::mcolmajor ? 56 :
          !(M2::mrowmajor || M2::mcolmajor) ? 55 :
          !(M1::mrowmajor || M1::mcolmajor) ? ( M2::mcolmajor ? 53 : 54 ) :
          -999 ) };
      enum { algo1 = (
          (cs == UNKNOWN || rs == UNKNOWN || xs == UNKNOWN) ? 64 :
          (cs < 16 && rs < 16 && xs < 16) ? algo2 :
          (cs <= 3 || rs <= 3 || xs <= 3) ? algo2 :
#ifdef TMV_USE_RECURSIVE_BLOCK
          (MbNbKb2 >= TMV_Q5) ? 63 :
#endif
          ( (twobig || nosimple) && McNcKc >= TMV_Q7 ) ? 64 :
          algo2 ) };

      try {
        M1r m1a = m1.Rows(0,Mx);
        M1r m1b = m1.Rows(Mx,M);
        M2c m2a = m2.Cols(0,Nx);
        M2c m2b = m2.Cols(Nx,N);
        M3s m3aa = m3.SubMatrix(0,Mx,0,Nx);
        M3s m3ba = m3.SubMatrix(Mx,M,0,Nx);
        M3s m3ab = m3.SubMatrix(0,Mx,Nx,N);
        M3s m3bb = m3.SubMatrix(Mx,M,Nx,N);
        const int xx = UNKNOWN;

        MultMM_Helper<algo1,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
            x,m1a,m2a,m3aa);
        MultMM_Helper<algo1,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
            x,m1a,m2b,m3ab);
        MultMM_Helper<algo1,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
            x,m1b,m2a,m3ba);
        MultMM_Helper<algo1,xx,xx,xs,add,ix,T,M1r,M2c,M3s>::call(
            x,m1b,m2b,m3bb);
      }
      catch (std::bad_alloc)
      {
        TMV_Warning(
            "Caught second bad_alloc error.\n"
            "Using (slower) algorithm that doesn't allocate temporary memory");
        // If failed again, use an algorithm that doesn't allocate memory.
        enum { algo2 = ( 
            ccc ? ( M3::miscomplex ? 11 : 12 ) :
            rcc ? ( M3::miscomplex ? 21 : 22 ) :
            crc ? ( M3::miscomplex ? 31 : 32 ) :
            21 ) };
        MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
      }
    }
  };

  // algo 68: Use Alberto and Nicolau's hybrid Winograd algorithm
  // This algorithm is based on the paper by Paolo D’Alberto and
  // Alexandru Nicolau titled "Adaptive Winograd’s Matrix Multiplications"
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<68,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"algo 68: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      typename M3::xview_type m3x = m3.XView();
      typedef typename M3::value_type T3;
      MultMM_Winograd(add,T3(T(x)),m1.XView(),m2.XView(),m3x);
    }
  };

#ifdef _OPENMP
  // algo 69: Split problem into smaller parts with OpenMP for parallelization
  // We take a pretty simple approach here, and just split up m2 and m3
  // matrices by columns and let each thread do a single matrix.
  // Then each thread calls algo 63 to calculate its product.
  // Also, we require that all but the last thread has a column width
  // that is a multiple of 16.  This way we get the maximum advantage from
  // our blocking structure while keeping the threads as balanced as possible.
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
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_MM
      std::cout<<"algo 69: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif

      enum { Mb = cs == UNKNOWN ? UNKNOWN : (cs >> 6) };
      enum { Nb = rs == UNKNOWN ? UNKNOWN : (rs >> 6) };
      enum { Kb = xs == UNKNOWN ? UNKNOWN : (xs >> 6) };
      enum { Mc = cs == UNKNOWN ? UNKNOWN : (cs < 16) ? 1 : (cs>>4) };
      enum { Nc = rs == UNKNOWN ? UNKNOWN : (rs < 16) ? 1 : (rs>>4) };
      enum { Kc = xs == UNKNOWN ? UNKNOWN : (xs < 16) ? 1 : (xs>>4) };
      enum { MbNbKb2 = (
            IntTraits2<IntTraits2<Mb,Nb>::prod,IntTraits2<Kb,Kb>::prod>::prod 
            ) };
      enum { McNcKc = IntTraits2<IntTraits2<Mc,Nc>::prod,Kc>::prod };
      enum { twobig = (Mb && Kb) || (Mb && Nb) || (Nb && Kb) };
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      enum { nosimple = !(ccc || rcc || crc) };

      // If we are in this function, then we pretty much know that 
      // we want to do one of the large matrix algorithms (63,64,68)
      // for the sub-problems.  So we call algo 72 to determine which
      // one to use.  
      // The algo1 selection here mimics that selection
      // when the sizes are known.
      // It's not perfect, since it uses the unthreaded M,N,K, rather
      // than either M/nthreads or N/nthreads, but it should usually
      // select a pretty good algorithm, and it keeps the 
      // compiler from instantiating the three possible algorithms
      // due to the if statements in algo 72.
      enum { algo1 = (
          (cs == UNKNOWN || rs == UNKNOWN || xs == UNKNOWN) ? 72 :
#ifdef TMV_USE_WINOGRAD
          (cs >= TMV_Q4 && rs >= TMV_Q4 && xs >= TMV_Q4) ? 68 :
#endif
#ifdef TMV_USE_RECURSIVE_BLOCK
          (MbNbKb2 >= TMV_Q5) ? 63 :
#endif
          64 ) };

      bool bad_alloc = false;
#ifdef PRINTALGO_MM_OMP
      std::ofstream fout("omp.out");
#endif
#pragma omp parallel
      {
        try {
          int num_threads = omp_get_num_threads();
          int mythread = omp_get_thread_num();
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
          {
            fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
          }
#endif
          if (num_threads == 1) {
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
            {
              fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
              fout<<"only 1 thread"<<std::endl;
            }
#endif
            MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
          }
          else if (M > N) 
          {
            int Mx = M / num_threads;
            Mx = ((((Mx-1)>>4)+1)<<4); // round up to next multiple of 16
            int i1 = mythread * Mx;
            int i2 = (mythread+1) * Mx;
            if (i2 > M || mythread == num_threads-1) i2 = M;
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
            {
              fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
              fout<<"M > N"<<std::endl;
              fout<<"Mx = "<<Mx<<std::endl;
              fout<<"i1 = "<<i1<<std::endl;
              fout<<"i2 = "<<i2<<std::endl;
            }
#endif
            if (i1 < M)  // Need to make sure, since we rounded up Mx!
            {
              typedef typename M1::const_rows_type M1r;
              typedef typename M3::rows_type M3r;
              enum { csx = UNKNOWN }; 
              M1r m1r = m1.CRows(i1,i2);
              M3r m3r = m3.CRows(i1,i2);

#ifdef PRINTALGO_MM_OMP
#pragma omp critical
              {
                fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
                fout<<"m1r = "<<m1r<<std::endl;
                fout<<"m2 = "<<m2<<std::endl;
                fout<<"m3r = "<<m3r<<std::endl;
              }
#endif
              MultMM_Helper<algo1,csx,rs,xs,add,ix,T,M1r,M2,M3r>::call(
                  x,m1r,m2,m3r);
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
              {
                fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
                fout<<"m3r => "<<m3r<<std::endl;
              }
#endif
            }
          } else {
            int Nx = N / num_threads;
            Nx = ((((Nx-1)>>4)+1)<<4); 
            int j1 = mythread * Nx;
            int j2 = (mythread+1) * Nx;
            if (j2 > N || mythread == num_threads-1) j2 = N;
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
            {
              fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
              fout<<"M <= N"<<std::endl;
              fout<<"Nx = "<<Nx<<std::endl;
              fout<<"j1 = "<<j1<<std::endl;
              fout<<"j2 = "<<j2<<std::endl;
            }
#endif
            if (j1 < N)  
            {
              typedef typename M2::const_cols_type M2c;
              typedef typename M3::cols_type M3c;
              enum { rsx = UNKNOWN }; 
              M2c m2c = m2.CCols(j1,j2);
              M3c m3c = m3.CCols(j1,j2);
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
              {
                fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
                fout<<"m1 = "<<m1<<std::endl;
                fout<<"m2c = "<<m2c<<std::endl;
                fout<<"m3c = "<<m3c<<std::endl;
              }
#endif
              MultMM_Helper<algo1,cs,rsx,xs,add,ix,T,M1,M2c,M3c>::call(
                  x,m1,m2c,m3c);
#ifdef PRINTALGO_MM_OMP
#pragma omp critical
              {
                fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
                fout<<"m3c => "<<m3c<<std::endl;
              }
#endif
            }
          }
        }
        catch (...) 
        // should only be std::bad_alloc, but it's good form to catch
        // everything inside a parallel block
        {
          bad_alloc = true;
        }
      }
#ifdef TMV_OPT_BAD_ALLOC
      enum { algo3 = 66 };
#else
      enum { algo3 = ( 
          ccc ? ( M3::miscomplex ? 11 : 12 ) :
          rcc ? ( M3::miscomplex ? 21 : 22 ) :
          crc ? ( M3::miscomplex ? 31 : 32 ) :
          21 ) };
#endif
      if (bad_alloc)
        MultMM_Helper<algo3,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
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
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_MM
      std::cout<<"algo 71: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      const int Mb = (M>>6); // = M/64
      const int Nb = (N>>6); // = N/64
      const int Kb = (K>>6); // = K/64
      const int Mc = M < 16 ? 1 : (M>>4); // = M/16
      const int Nc = N < 16 ? 1 : (N>>4); // = N/16
      const int Kc = K < 16 ? 1 : (K>>4); // = K/16
      const bool twobig = (Mb&&Nb) || (Mb&&Kb) || (Nb&&Kb);

      TMVStaticAssert(!M3::mrowmajor);
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      enum { nosimple = !(ccc || rcc || crc) };
      enum { algo2 = (
          ccc ? ( M3::miscomplex ? 11 : 12 ) :
          rcc ? ( M3::miscomplex ? 21 : 22 ) :
          crc ? ( M3::miscomplex ? 31 : 32 ) :
          ( M1::mstepj != UNKNOWN && M2::mstepi != UNKNOWN &&
            M3::mstepi != UNKNOWN )  ? 21 :
          rrc ? ( M3::miscomplex ? 31 : 32 ) :
          // Strange -- 31,32 are faster than 51 for the sizes I tested.
          !M3::mcolmajor ? 56 :
          !(M2::mrowmajor || M2::mcolmajor) ? 55 :
          !(M1::mrowmajor || M1::mcolmajor) ? ( M2::mcolmajor ? 53 : 54 ) :
          -999 ) };
      TMVStaticAssert(algo2 != -999);

      // Put the small matrix option first, so it doesn't have to 
      // go through a bunch of if/else statements.  For large matrices,
      // all these if/else's don't matter for the total time.
      if (M < 16 && N < 16 && K < 16)
        MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
      else if (M <= 3 || N <= 3 || K <= 3)
        MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#ifdef _OPENMP
      else if (!omp_in_parallel() && (Mb || Nb) && ( Mc * Nc * Kc >= TMV_Q6 ) )
        MultMM_Helper<69,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
#ifdef TMV_USE_WINOGRAD
      else if (M >= TMV_Q4 && N >= TMV_Q4 && K >= TMV_Q4)
        MultMM_Helper<68,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
#ifdef TMV_USE_RECURSIVE_BLOCK
      else if (Mb*Nb*Kb*Kb >= TMV_Q5)
        MultMM_Helper<63,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
      else if (
          ( M >= 16 && N >= 16 && K >= 16 ) ||
          ( (twobig || nosimple) && (Mc * Nc * Kc >= TMV_Q7) ) )
        MultMM_Helper<64,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
      else
        MultMM_Helper<algo2,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  // algo 72: Same as 71, except 63,64,68 are the only options.
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<72,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
#ifdef PRINTALGO_MM
      std::cout<<"algo 72: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      const int Mb = (M>>6); // = M/64
      const int Nb = (N>>6); // = N/64
      const int Kb = (K>>6); // = K/64
      const int Mc = M < 16 ? 1 : (M>>4); // = M/16
      const int Nc = N < 16 ? 1 : (N>>4); // = N/16
      const int Kc = K < 16 ? 1 : (K>>4); // = K/16
      const bool twobig = (Mb&&Nb) || (Mb&&Kb) || (Nb&&Kb);

      TMVStaticAssert(!M3::mrowmajor);
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      enum { nosimple = !(ccc || rcc || crc) };

      if (false) {}
#ifdef TMV_USE_WINOGRAD
      else if (M >= TMV_Q4 && N >= TMV_Q4 && K >= TMV_Q4)
        MultMM_Helper<68,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
#ifdef TMV_USE_RECURSIVE_BLOCK
      else if (Mb*Nb*Kb*Kb >= TMV_Q5)
        MultMM_Helper<63,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
      else
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
      enum { algo1 = 71 };

#ifdef TMV_OPT_SMALL
      TMVStaticAssert(cs == UNKNOWN);
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      TMVStaticAssert(!M3::mrowmajor);
      TMVStaticAssert(ccc || rcc || crc);
      enum { algo2 = ccc ? 11 : rcc ? 21 : crc ? 31 : 71 };

#ifdef PRINTALGO_MM
      std::cout<<"algo 75: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      if (M <= 3)
      {
        // then it is worth figuring out what M is.
        switch (M)
        {
          case 0 :
            // do nothing
            break;
          case 1 :
            MultMM_Helper<2,1,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            break;
          case 2 :
            MultMM_Helper<algo2,2,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            break;
          case 3 :
            MultMM_Helper<algo2,3,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
      }
      else 
#endif
        MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  // algo 76: Unknown rs, check if N is small
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<76,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      enum { algo1 = 71 };

#ifdef TMV_OPT_SMALL
      TMVStaticAssert(rs == UNKNOWN);
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      TMVStaticAssert(!M3::mrowmajor);
      TMVStaticAssert(ccc || rcc || crc);
      enum { algo2 = ccc ? 11 : rcc ? 21 : crc ? 31 : 71 };

#ifdef PRINTALGO_MM
      std::cout<<"algo 76: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      if (N <= 3)
      {
        // then it is worth figuring out what N is.
        switch (N)
        {
          case 0 :
            // do nothing
            break;
          case 1 :
            MultMM_Helper<3,cs,1,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            break;
          case 2 :
            MultMM_Helper<algo2,cs,2,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            break;
          case 3 :
            MultMM_Helper<algo2,cs,3,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
      }
      else 
#endif
        MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  // algo 77: Unknown cs, check if M is small
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<77,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      enum { algo1 = 71 };

#ifdef TMV_OPT_SMALL
      TMVStaticAssert(xs == UNKNOWN);
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      TMVStaticAssert(!M3::mrowmajor);
      TMVStaticAssert(ccc || rcc || crc);
      enum { algo2 = ccc ? 11 : rcc ? 21 : crc ? 31 : 71 };

#ifdef PRINTALGO_MM
      std::cout<<"algo 77: M,N,K,cs,rs,xs,x = "<<M<<','<<N<<','<<K<<','<<cs<<','<<rs<<','<<xs<<','<<T(x)<<std::endl;
#endif
      if (K <= 3)
      {
        // then it is worth figuring out what K is.
        switch (K)
        {
          case 0 :
            Maybe<!add>::zero(m3);
            break;
          case 1 :
            MultMM_Helper<4,cs,rs,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            break;
          case 2 :
            MultMM_Helper<algo2,cs,rs,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            break;
          case 3 :
            MultMM_Helper<algo2,cs,rs,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
      }
      else 
#endif
        MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  // algo -2: The same as -1, but allow for the possibility of 
  // using InstMultMM.  This is used whenever we recurse back to
  // the start after making a temporary.
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<-2,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M1::value_type T1;
      typedef typename M2::value_type T2;
      typedef typename M3::value_type T3;

      enum { inst = (
          Traits<T1>::isinst &&
          Traits<T2>::isinst &&
          Traits<T3>::isinst &&
#ifdef TMV_INST_MIX
          Traits2<T1,T3>::samebase &&
          Traits2<T2,T3>::samebase &&
#else
          Traits2<T1,T3>::sametype &&
          Traits2<T2,T3>::sametype &&
#endif
          cs == UNKNOWN && rs == UNKNOWN && xs == UNKNOWN ) };

      enum { algo = inst ? 5 : -1 };
      MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, int xs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMM_Helper<-1,cs,rs,xs,add,ix,T,M1,M2,M3> 
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M1::value_type T1;
      typedef typename M2::value_type T2;
      typedef typename M3::value_type T3;

      // Possible algorithms are:
      //
      // Trivial:
      //  0 = cs or rs == 0, or xs == 0 && add, so nothing to do
      //  1 = xs == 0 && !add, so m3.Zero();
      //  2 = cs == 1: reduces to trivial MultMV function
      //  3 = rs == 1: reduces to trivial MultMV function
      //  4 = xs == 1: reduces to trivial Rank1Update function
      //  5 = call InstMultMM (This can happen after copying to a temporary.)
      //
      // CCC: 
      // 11 = CCC, loop over n
      // 12 = CCC, 4x2 blocks in KxN
      //
      // RCC:
      // 21 = RCC, loop over n
      // 22 = RCC, 2x2 blocks in MxN
      // 23 = RCC, 2x2 blocks in MxN, M,N,K known and all compact matrices.
      //      (This is the main kernel for large matrices, so there are
      //       also specialized versions for certain M,N,K, and also
      //       specializations for SSE and (not yet) AltiVec.)
      // 24 = RCC, 2x2 blocks in MxN, all compact matrices, but 
      //      M,N,K don't need to be known at compile time.
      //      Also (and most important really), if using SSE, 
      //      the columns are all aligned on 128 byte boundaries, so 
      //      the SSE commands can be used.
      //
      // CRC:
      // 31 = CRC, loop over k
      // 32 = CRC, rank-4 updates 
      //
      // Fully Unrolled:
      // TODO: I haven't done these, but I'm not sure that there is a need
      // for them.  The current setup is about as fast or faster (depending
      // on ccc,rcc,etc.) than Eigen's known-size calculation even down to 
      // 2x2 matrices, so I'm not sure that these would speed things up 
      // much if at all.
      // 41 = fully unroll CCC, apply x to m2
      // 42 = fully unroll RCC, apply x to m3
      // 43 = fully unroll CRC, apply x to m2
      // 44 = unroll inner K loop, apply x to m2
      // 45 = unroll inner K loop, apply x to m3
      //
      // No good majority:
      // 51 = RRC, no good algorithm, copy m2 to column major.
      // 52 = **R, transpose to **C
      // 53 = m1 is non-major - use temporary (R)
      // 54 = m1 is non-major - use temporary (C)
      // 55 = m2 is non-major - use temporary (C)
      // 56 = m3 is non-major - use temporary (C)
      //
      // Main drivers for large matrices:
      // 61 = Recurse down to single blocks with divide-and-conquer.
      // 62 = Simpler loop over M,N,then K in blocks.
      // 63 = Same pattern as 61, but copy m1,m2 into a block structure
      // 64 = Same pattern as 62, but copy m1,m2 into a block structure
      // 66 = A simple division into 4 submatrices, then call another algo.
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

      enum { ccc = M1::mcolmajor && M2::mcolmajor && M3::mcolmajor };
      enum { rcc = M1::mrowmajor && M2::mcolmajor && M3::mcolmajor };
      enum { crc = M1::mcolmajor && M2::mrowmajor && M3::mcolmajor };
      enum { rrc = M1::mrowmajor && M2::mrowmajor && M3::mcolmajor };
      enum { Mb = cs == UNKNOWN ? UNKNOWN : (cs >> 6) };
      enum { Nb = rs == UNKNOWN ? UNKNOWN : (rs >> 6) };
      enum { Kb = xs == UNKNOWN ? UNKNOWN : (xs >> 6) };
      enum { Mc = cs == UNKNOWN ? UNKNOWN : (cs < 16) ? 1 : (cs>>4) };
      enum { Nc = rs == UNKNOWN ? UNKNOWN : (rs < 16) ? 1 : (rs>>4) };
      enum { Kc = xs == UNKNOWN ? UNKNOWN : (xs < 16) ? 1 : (xs>>4) };
      enum { MbNbKb2 = (
            IntTraits2<IntTraits2<Mb,Nb>::prod,IntTraits2<Kb,Kb>::prod>::prod 
            ) };
      enum { McNcKc = IntTraits2<IntTraits2<Mc,Nc>::prod,Kc>::prod };
      enum { twobig = (Mb && Kb) || (Mb && Nb) || (Nb && Kb) };
      enum { nosimple = !(ccc || rcc || crc) };

#if 0
      enum { algo = 63 };
      enum { algo1 = 21 };
#else
#if TMV_OPT == 0 
      enum { algo = (
          (M3::mrowmajor && !M3::mcolmajor) ? 52 :
          ccc ? 11 :
          crc ? 31 :
          21 ) };
#else
      enum { algo = (
          ( cs == 0 || rs == 0 || (xs == 0 && add) ) ? 0 :
          ( xs == 0 && !add ) ? 1 :
          cs == 1 ? 2 :
          rs == 1 ? 3 :
          xs == 1 ? 4 :
          M3::mrowmajor && !M3::mcolmajor ? 52 :
          (cs == UNKNOWN || rs == UNKNOWN || xs == UNKNOWN ) ? (
            ccc ? ( cs == UNKNOWN ? 75 : cs <= 3 ? 11 : 71 ) :
            rcc ? ( xs == UNKNOWN ? 77 : xs <= 3 ? 21 : 71 ) :
            crc ? ( cs == UNKNOWN ? 75 : cs <= 3 ? 31 : 71 ) :
            71 ) :
#ifdef _OPENMP
          // This is different from the test we do in algo 71.
          // I'm having trouble with OPENMP for some small known values 
          // of cs,rs,xs.  e.g. 5000x5x5 leads to a compiler error with 
          // g++ 4.3.  So for now, I'm just diabling OPENMP for small sizes.
          // But really, I should try to figure out what is going on here.
          //( (Mb || Nb) && (McNcKc >= TMV_Q6) ) ? 69 :
          ( (Mb && Nb && Kb) && (McNcKc >= TMV_Q6) ) ? 69 :
#endif
#ifdef TMV_USE_WINOGRAD
          ( cs >= TMV_Q4 && rs >= TMV_Q4 && xs >= TMV_Q4 ) ? 68 :
#endif
#ifdef TMV_USE_RECURSIVE_BLOCK
          ( MbNbKb2 >= TMV_Q5 ) ? 63 :
#endif
          // TODO: Eigen or BLAS is sometimes faster than TMV for 
          // one small cs,rs or xs, and large values of the other two.
          // I need to investigate what different algorithms might be 
          // better for these cases.  For now calling MultMV or Rank1Update
          // seems to be the best option.
          // In fact, just improving them to include SSE might do the trick.
          ( ccc && (cs <= 3 || rs <= 3 || xs <= 3) ) ? 11 :
          ( rcc && (cs <= 3 || rs <= 3 || xs <= 3) ) ? 21 :
          ( crc && (cs <= 3 || rs <= 3 || xs <= 3) ) ? 31 :
          ( (cs >= 16 && rs >= 16 && xs >= 16) ||
            ( (twobig || nosimple) && McNcKc >= TMV_Q7 ) ) ? 64 :
          // For known sizes, the MultMV calls (algo 11,21) seem to be faster
          // than the 2x2 algorithms (algo 12,22).  For unknown sizes
          // (cf. selection in algo 71) the opposite is true.
          ccc ? 11 :
          rcc ? 21 :
          crc ? 31 :
          // When steps are known, direct MultMV is faster than copying.
          // Basically, the non-unit step vectors for MultMV are ok, 
          // since the non-unit step is known and the lengths of the 
          // vectors are not too large.
          ( M1::mstepj != UNKNOWN && M2::mstepi != UNKNOWN &&
            M3::mstepi != UNKNOWN )  ? 21 :
          rrc ? 31 : 
          !M3::mcolmajor ? 56 :
          !(M2::mrowmajor || M2::mcolmajor) ? 55 :
          !(M1::mrowmajor || M1::mcolmajor) ? ( M2::mcolmajor ? 53 : 54 ) :
          // The above should have caught all possibilities.
          // So selection should never get here.
          -999 ) };
      TMVStaticAssert(algo != -999);
#endif
#ifdef TMV_OPT_BAD_ALLOC
      enum { algo1 = ( algo < 60 ? algo : 66 ) };
#else
      enum { algo1 = (
          algo < 60 ? algo :  // Nothing less than 60 should throw a bad_alloc
          ccc ? ( M3::miscomplex ? 11 : 12 ) :
          rcc ? ( M3::miscomplex ? 21 : 22 ) :
          crc ? ( M3::miscomplex ? 31 : 32 ) :
          21 ) };
#endif
#endif
#ifdef PRINTALGO_MM
      const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
      const int K = xs==UNKNOWN ? int(m1.rowsize()) : xs;
      std::cout<<"InlineMultMM: x = "<<ix<<"  "<<T(x)<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
      std::cout<<"M = "<<M<<"  N = "<<N<<"  K = "<<K<<std::endl;
#ifdef _OPENMP
      std::cout<<"69 ? McNcKc = "<<McNcKc<<" >= Q6 = "<<TMV_Q6<<std::endl;
#endif
#ifdef TMV_USE_WINOGRAD
      std::cout<<"68 ? all cs,rs,xs >= Q4 = "<<TMV_Q4<<std::endl;
#endif
#ifdef TMV_USE_RECURSIVE_BLOCK
      std::cout<<"63 ? MbNbKb2 = "<<MbNbKb2<<" >= Q5 = "<<TMV_Q5<<std::endl;
#endif
      std::cout<<"64 ? McNcKc = "<<McNcKc<<" >= Q7 = "<<TMV_Q7<<std::endl;
      std::cout<<"cs = "<<cs<<"  rs = "<<rs<<"  xs = "<<xs<<std::endl;
      std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
      //std::cout<<"m1 = "<<m1<<std::endl;
      //std::cout<<"m2 = "<<m2<<std::endl;
      //std::cout<<"m3 = "<<m3<<std::endl;
#endif
      try {
        MultMM_Helper<algo,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
            x,m1.mat(),m2.mat(),m3.mat());
      } 
      catch (std::bad_alloc)
      {
        MultMM_Helper<algo1,cs,rs,xs,add,ix,T,M1,M2,M3>::call(
            x,m1.mat(),m2.mat(),m3.mat());
      }
      //std::cout<<"m3 => "<<m3<<std::endl;
    }
  };

  template <bool add, int ix, class T, class M1, class M2, class M3>
  static inline void InlineMultMM(const Scaling<ix,T>& x, 
      const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
      BaseMatrix_Rec_Mutable<M3>& m3)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M2::mrowsize,M3::mrowsize>::same));
    TMVAssert(m1.colsize() == m3.colsize());
    TMVAssert(m1.rowsize() == m2.colsize());
    TMVAssert(m2.rowsize() == m3.rowsize());

    // cs = m3.colsize
    // rs = m3.rowsize
    // xs = "extra size" = m1.rowsize, m2.colsize
    enum { cs = Sizes<M3::mcolsize,M1::mcolsize>::size };
    enum { rs = Sizes<M3::mrowsize,M2::mrowsize>::size };
    enum { xs = Sizes<M1::mrowsize,M2::mcolsize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename M2::const_view_type M2v;
    typedef typename M3::view_type M3v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    M3v m3v = m3.View();
    // This View() bit makes the compiled code a bit smaller, since it then
    // doesn't instantiate separate versions for Matrix and MatrixView.
    // TODO: Make this change to all the other Inline functions.
    MultMM_Helper<-1,cs,rs,xs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
  }

  // Defined in TMV_MultMM.cpp
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3);
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3);

#undef TMV_OPT_SMALL
#undef TMV_OPT_CLEANUP
#undef TMV_USE_RECURSIVE_BLOCK
#undef TMV_USE_WINOGRAD
#undef TMV_OPT_BAD_ALLOC

#undef TMV_Q1
#undef TMV_Q2
#undef TMV_Q3
#undef TMV_Q4
#undef TMV_Q5
#undef TMV_Q6
#undef TMV_Q7
#undef TMV_ZeroIX

} // namespace tmv

#endif 
