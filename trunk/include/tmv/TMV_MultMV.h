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


#ifndef TMV_MultMV_H
#define TMV_MultMV_H

#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_AddVV.h"
#include "TMV_Vector.h"
#include "TMV_SmallVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Prefetch.h"

namespace tmv {

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
#define TMV_OPT_SMALL
#endif

  // When using an algorithm that does 4 columns at a time or 4 rows at a time,
  // this parameter specifies whether to use specific code for the remaining
  // 1, 2, or 3 columns or rows or generic code.  It is mildly expensive
  // in terms of code bloat, but not too bad.  And it helps all moderately
  // sized matrices, (e.g. M,N ~ 10-100) unless they happen to have M,N%4 = 0.
  // These matrices are pretty common, so we only require TMV_OPT = 2 for this.
#if TMV_OPT >= 2
#define TMV_OPT_CLEANUP
#endif

  // When doing something like w = x * m * v, it is often better to 
  // apply x to either v before doing the multiplication, or w after doing it.
  // The optimal choice depends on the sizes of v and w, so this parameter
  // specifies to determine which is better based on the runtime sizes
  // of the vectors.  It is a relatively cheap optimization, but it also
  // only really helps when the sizes are very different (eg. 100 x 2).
  // There is also a mid increase in code bloat.  So we require TMV_OPT = 2.
#if TMV_OPT >= 2
#define TMV_OPT_SCALE
#endif

//#define PRINTALGO_MV

#ifdef PRINTALGO_MV
#include <iostream>
#endif

  // There are a number of values used in the algorithm selection
  // that are either arbitrary or empirical.
  // So I put them all here to make them easier to change and to
  // track down in the code.

  // Q1 is the maximum value of cs*rs for which a SmallMatrix * SmallVector
  // product is fully unrolled. 
  // We have other cuts that are based purely on timings in algo -1.
  // (Search for "unroll = ".)  This parameter supersedes that value
  // to turn off unrolling for larger matrices to reduce code bloat.
#if TMV_OPT >= 3
#define TMV_Q1 1024 // It never actually gets this high.
#elif TMV_OPT >= 2
#define TMV_Q1 25
#elif TMV_OPT >= 1
#define TMV_Q1 9
#else
#define TMV_Q1 0
#endif

  // Q2 is the minimum size to copy a vector if its step != 1.
#define TMV_Q2 4

  // Q3 is the crossover memory size to start using prefetch commands.
  // This is undoubtedly a function of the L1 (and L2?) cache size,
  // but 2KBytes is probably not too bad for most machines.
  // (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_Q3 2048

  // Q4 is the ratio of the time to copy a vector to the time to scale it.
  // ( Or technically Q4 is 1 + this ratio. )
#define TMV_Q4 4

  // Q5 is the minimum value of cs for which unrolling by columns is
  // faster than unrolling by rows (when m1 is colmajor).
  // (This seems to be fairly independent of rs.)
  // (It is also now almost completely irrelevant, since most of these
  //  cases use the non-unrolled algorithm now anyway.)
#define TMV_Q5 20

  // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_ZeroIX (ix==0)
//#define TMV_ZeroIX (ix!=1)

  template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper;

  // algo 0: cs or rs = 0, so nothing to do
  // Correction: if rs = 0, cs != 0 and !add, then we need to do v3.Zero().
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<0,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void call(const Scaling<ix,T>& , const M1& , const V2& , V3& v3) 
    { Maybe<(!add && rs == 0 && cs != 0)>::zero(v3); } 
  };

  // algo 1: cs == 1, so simplifies to a MultVV function
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<1,1,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { 
      typedef typename M1::const_row_type M1r;
      typedef typename V2::value_type T2;
      Maybe<add>::add( v3.ref(0) ,
          x * MultVV_Helper<-1,rs,M1r,V2>::call(m1.get_row(0),v2) );
    }
  };

  // algo 2: rs == 1, so simplifies to an AddVV or MultXV function
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<2,cs,1,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename Traits2<T,typename V2::value_type>::type PT2;
      typedef typename M1::const_col_type M1c;
      AddVV_Helper<-1,cs,add,0,PT2,M1c,V3>::call( 
          Scaling<0,PT2>(x*v2.cref(0)) , m1.get_col(0), v3 ); 
    }
  };

  // algo 3: call InstMultMV
  template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<3,cs,rs,false,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename V3::value_type T3;
      typename V3::xview_type v3x = v3.XView();
      InstMultMV(T3(x),m1.XView(),v2.XView(),v3x);
    }
  };
  template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<3,cs,rs,true,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename V3::value_type T3;
      typename V3::xview_type v3x = v3.XView();
      InstAddMultMV(T3(x),m1.XView(),v2.XView(),v3x);
    }
  };

  // algo 11: The basic column major loop
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<11,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 11: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

      Maybe<!add>::zero(v3);
      if (N) {
        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename M1::const_col_type M1c;
        PT2 Xj;

        enum { c2 = V2::vconj };

        typedef typename M1c::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        IT1 A0j = m1.get_col(0).NonConj().begin();
        IT2 X = v2.NonConj().begin();
        const IT3 Y0 = v3.begin();
        const int Astepj = m1.stepj();

        do {
          if (*X != T2(0)) {
            Xj = ZProd<false,c2>::prod(x , *X++);
            // y += Xj * A.col(j);
            AddVV_Helper<-1,cs,true,0,PT2,M1c,V3>::call2(
                M,Scaling<0,PT2>(Xj),A0j,Y0);
          } else {
            ++X;
          }
          A0j.ShiftP(Astepj);
        } while (--N);
      }
    }
  };

  // algo 12: column major, 4 columns at a time
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<12,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void loop_4_cols(const int M, const int N,
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
      IT1 A1 = A0; A1.ShiftP(stepj);
      IT1 A2 = A1; A2.ShiftP(stepj);
      IT1 A3 = A2; A3.ShiftP(stepj);
      IT2 X = v2.begin();
      const IT3 Y_begin = v3.begin();
      IT3 Y = Y_begin;

      const bool dopref = M * sizeof(T1) >= TMV_Q3;

      Prefetch_Read(X.GetP());
      Prefetch_Read(A0.GetP());
      Prefetch_Read(A1.GetP());
      Prefetch_Read(A2.GetP());
      Prefetch_Read(A3.GetP());
      Prefetch_MultiWrite(Y.GetP());

      int i;

      TMVAssert(N_4 > 0);
      do {
        X0 = x * X[0]; X1 = x * X[1]; X2 = x * X[2]; X3 = x * X[3]; X += 4;
        Y = Y_begin;

        i=M; do {
          Y0 = *Y;
          A00 = *A0++; A01 = *A1++; A02 = *A2++; A03 = *A3++; 
          Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
          *Y++ = Y0;
        } while (--i);
        A0 += stepj_4; A1 += stepj_4; A2 += stepj_4; A3 += stepj_4;
        if (dopref) {
          Prefetch_Read(A0.GetP());
          Prefetch_Read(A1.GetP());
          Prefetch_Read(A2.GetP());
          Prefetch_Read(A3.GetP());
        }
      } while (--N_4);
    }
    template <bool addx, class M1x, class V2x>
    static inline void cleanup(const int nb,
        const Scaling<ix,T>& x, const M1x& m1, const V2x& v2, V3& v3)
    {
      TMVAssert(nb == 1 || nb == 2 || nb == 3 || nb == 4);
#ifdef TMV_OPT_CLEANUP
      switch (nb) {
        case 1 : 
          MultMV_Helper<13,cs,1,addx,ix,T,M1x,V2x,V3>::call(x,m1,v2,v3);
          break;
        case 2 :
          MultMV_Helper<13,cs,2,addx,ix,T,M1x,V2x,V3>::call(x,m1,v2,v3);
          break;
        case 3 :
          MultMV_Helper<13,cs,3,addx,ix,T,M1x,V2x,V3>::call(x,m1,v2,v3);
          break;
        case 4 :
          MultMV_Helper<13,cs,4,addx,ix,T,M1x,V2x,V3>::call(x,m1,v2,v3);
      }
#else
      MultMV_Helper<11,cs,UNKNOWN,addx,ix,T,M1x,V2x,V3>::call(x,m1,v2,v3);
#endif
    }
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 12: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);
      typedef typename M1::const_cols_type M1c;
      typedef typename V2::const_subvector_type V2s;
      if (M) {
        if (N > 4) {
          const int na = ((N>>2)<<2); 
          const int nb = N-na;
          Maybe<!add>::zero(v3);
          loop_4_cols(M,na,x,m1,v2,v3);
          if (nb) {
            M1c m1b = m1.CCols(na,N);
            V2s v2b = v2.CSubVector(na,N);
            cleanup<true>(nb,x,m1b,v2b,v3);
          }
        } else if (N) {
          cleanup<add>(N,x,m1,v2,v3);
        }
      }
    }
  };

  // algo 13: do all columns at once -- rs <= 4, and must be known
  // rs == 1
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<13,cs,1,add,ix,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { MultMV_Helper<2,cs,1,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3); }
  };
  // rs == 2
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<13,cs,2,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 13 N==2: M,N,cs,rs,x = "<<M<<','<<2<<','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);

      if (M) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_iterator IT1;
        T1 A00, A01;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        const PT2 X0 = x * v2.cref(0);
        const PT2 X1 = x * v2.cref(1);

        const int stepj = m1.stepj();
        IT1 A0 = m1.get_col(0).begin();
        IT1 A1 = A0; A1.ShiftP(stepj);

        typedef typename V3::value_type T3;
        typedef typename V3::iterator IT3;
        T3 Y0;
        IT3 Y = v3.begin();

        do {
          A00 = *A0++; A01 = *A1++;
          Maybe<add>::set(Y0,*Y);
          Maybe<add>::add(Y0 , A00 * X0 + A01 * X1);
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
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 13 N==3: M,N,cs,rs,x = "<<M<<','<<3<<','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);

      if (M) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_iterator IT1;
        T1 A00, A01, A02;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        const PT2 X0 = x * v2.cref(0);
        const PT2 X1 = x * v2.cref(1);
        const PT2 X2 = x * v2.cref(2);

        const int stepj = m1.stepj();
        IT1 A0 = m1.get_col(0).begin();
        IT1 A1 = A0; A1.ShiftP(stepj);
        IT1 A2 = A1; A2.ShiftP(stepj);

        typedef typename V3::value_type T3;
        typedef typename V3::iterator IT3;
        T3 Y0;
        IT3 Y = v3.begin();

        do {
          A00 = *A0++; A01 = *A1++; A02 = *A2++;
          Maybe<add>::set(Y0,*Y);
          Maybe<add>::add(Y0 , A00 * X0 + A01 * X1 + A02 * X2);
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
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 13 N==4: M,N,cs,rs,x = "<<M<<','<<4<<','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);

      if (M) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_iterator IT1;
        T1 A00, A01, A02, A03;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        const PT2 X0 = x * v2.cref(0);
        const PT2 X1 = x * v2.cref(1);
        const PT2 X2 = x * v2.cref(2);
        const PT2 X3 = x * v2.cref(3);

        const int stepj = m1.stepj();
        IT1 A0 = m1.get_col(0).begin();
        IT1 A1 = A0; A1.ShiftP(stepj);
        IT1 A2 = A1; A2.ShiftP(stepj);
        IT1 A3 = A2; A3.ShiftP(stepj);

        typedef typename V3::value_type T3;
        typedef typename V3::iterator IT3;
        T3 Y0;
        IT3 Y = v3.begin();

        do {
          A00 = *A0++; A01 = *A1++; A02 = *A2++; A03 = *A3++;
          Maybe<add>::set(Y0,*Y);
          Maybe<add>::add(Y0 , A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3);
          *Y++ = Y0;
        } while (--M);
      }
    }
  };

  // algo 15: column major, 2 columns at a time, complex v3
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<15,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void loop_2_cols(const int M, const int N,
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVAssert(N%2 == 0);
      TMVAssert(N > 0);

      typedef typename M1::value_type T1;
      typedef typename M1::const_col_type::const_nonconj_type::const_iterator IT1;
      T1 A00, A01;

      typedef typename V2::value_type T2;
      typedef typename Traits2<T,T2>::type PT2;
      typedef typename V2::const_nonconj_type::const_iterator IT2;
      PT2 X0, X1;

      typedef typename V3::value_type T3;
      typedef typename V3::real_type RT3;
      typedef typename V3::iterator IT3;
      RT3 rY0, iY0;

      int N_2 = N>>1; // N_2 = N/2
      const int stepj = m1.stepj();
      const int stepj_2 = (stepj<<1) - M;
      // step over 2 columns and back to start

      IT1 A0 = m1.get_col(0).NonConj().begin();
      IT1 A1 = A0; A1.ShiftP(stepj);
      IT2 X = v2.NonConj().begin();
      const IT3 Y_begin = v3.begin();
      IT3 Y = Y_begin;

      enum { c1 = M1::mconj };
      enum { c2 = V2::vconj };

      const bool dopref = M * sizeof(T1) >= TMV_Q3;

      Prefetch_Read(X.GetP());
      Prefetch_Read(A0.GetP());
      Prefetch_Read(A1.GetP());
      Prefetch_MultiWrite(Y.GetP());

      int i;

      TMVAssert(N_2 > 0);
      do {
        X0 = ZProd<false,c2>::prod(x,X[0]);
        X1 = ZProd<false,c2>::prod(x,X[1]);
        X += 2;
        Y = Y_begin;

        i=M; do {
          A00 = *A0++; A01 = *A1++;
          rY0 = real(*Y); iY0 = imag(*Y);
          rY0 += ZProd<c1,false>::rprod(A00,X0);
          rY0 += ZProd<c1,false>::rprod(A01,X1);
          iY0 += ZProd<c1,false>::iprod(A00,X0);
          iY0 += ZProd<c1,false>::iprod(A01,X1);
          *Y++ = T3(rY0,iY0);
        } while (--i);
        A0 += stepj_2; A1 += stepj_2;
        if (dopref) {
          Prefetch_Read(A0.GetP());
          Prefetch_Read(A1.GetP());
        }
      } while (--N_2);
    }
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 15: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);
      if (M) {
        const int na = ((N>>1)<<1);
        const int nb = N-na;
        Maybe<!add>::zero(v3);
        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        typedef typename M1::const_col_type M1c;
        if (na) loop_2_cols(M,na,x,m1,v2,v3);
        if (nb) {
          AddVV_Helper<-1,cs,true,0,PT2,M1c,V3>::call(
              Scaling<0,PT2>(x*v2.cref(na)),m1.get_col(na),v3); 
        }
      }
    }
  };

  // algo 16: do all columns at once, complex v3: rs <= 4, and must be known
  // rs == 1
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<16,cs,1,add,ix,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { MultMV_Helper<2,cs,1,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3); }
  };
  // rs == 2
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<16,cs,2,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVStaticAssert(V3::viscomplex);
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 16 N==2: M,N,cs,rs,x = "<<M<<','<<2<<','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);

      if (M) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_nonconj_type::const_iterator IT1;
        T1 A00, A01;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        enum { c2 = V2::vconj };
        const PT2 X0 = ZProd<false,c2>::prod(x,v2.NonConj().cref(0));
        const PT2 X1 = ZProd<false,c2>::prod(x,v2.NonConj().cref(1));

        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT3;
        typedef typename V3::flatten_type::iterator IT3;
        RT3 rY0, iY0;

        const int stepj = m1.stepj();
        IT1 A0 = m1.get_col(0).NonConj().begin();
        IT1 A1 = A0; A1.ShiftP(stepj);
        IT3 Y = v3.Flatten().begin();

        enum { c1 = M1::mconj };

        do {
          A00 = *A0++; A01 = *A1++;
          Maybe<add>::set(rY0 , *Y);
          Maybe<add>::set(iY0 , Y[1]);
          Maybe<add>::add(rY0 , ZProd<c1,false>::rprod(A00,X0));
          Maybe<add>::add(iY0 , ZProd<c1,false>::iprod(A00,X0));
          rY0 += ZProd<c1,false>::rprod(A01,X1);
          iY0 += ZProd<c1,false>::iprod(A01,X1);
          *Y = rY0; Y[1] = iY0;  Y += 2;
        } while (--M);
      }
    }
  };
  // rs == 3
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<16,cs,3,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 16 N==3: M,N,cs,rs,x = "<<M<<','<<3<<','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);

      if (M) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_nonconj_type::const_iterator IT1;
        T1 A00, A01, A02;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        enum { c2 = V2::vconj };
        const PT2 X0 = ZProd<false,c2>::prod(x,v2.NonConj().cref(0));
        const PT2 X1 = ZProd<false,c2>::prod(x,v2.NonConj().cref(1));
        const PT2 X2 = ZProd<false,c2>::prod(x,v2.NonConj().cref(2));

        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT3;
        typedef typename V3::flatten_type::iterator IT3;
        RT3 rY0, iY0;

        const int stepj = m1.stepj();
        IT1 A0 = m1.get_col(0).NonConj().begin();
        IT1 A1 = A0; A1.ShiftP(stepj);
        IT1 A2 = A1; A2.ShiftP(stepj);
        IT3 Y = v3.Flatten().begin();

        enum { c1 = M1::mconj };

        do {
          A00 = *A0++; A01 = *A1++; A02 = *A2++;
          Maybe<add>::set(rY0 , *Y);
          Maybe<add>::set(iY0 , Y[1]);
          Maybe<add>::add(rY0 , ZProd<c1,false>::rprod(A00,X0));
          Maybe<add>::add(iY0 , ZProd<c1,false>::iprod(A00,X0));
          rY0 += ZProd<c1,false>::rprod(A01,X1);
          iY0 += ZProd<c1,false>::iprod(A01,X1);
          rY0 += ZProd<c1,false>::rprod(A02,X2);
          iY0 += ZProd<c1,false>::iprod(A02,X2);
          *Y = rY0; Y[1] = iY0;  Y += 2;
        } while (--M);
      }
    }
  };
  // rs == 4
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<16,cs,4,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 16 N==4: M,N,cs,rs,x = "<<M<<','<<4<<','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);

      if (M) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_col_type::const_nonconj_type::const_iterator IT1;
        T1 A00, A01, A02, A03;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        enum { c2 = V2::vconj };
        const PT2 X0 = ZProd<false,c2>::prod(x,v2.NonConj().cref(0));
        const PT2 X1 = ZProd<false,c2>::prod(x,v2.NonConj().cref(1));
        const PT2 X2 = ZProd<false,c2>::prod(x,v2.NonConj().cref(2));
        const PT2 X3 = ZProd<false,c2>::prod(x,v2.NonConj().cref(3));

        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT3;
        typedef typename V3::flatten_type::iterator IT3;
        RT3 rY0, iY0;

        const int stepj = m1.stepj();
        IT1 A0 = m1.get_col(0).NonConj().begin();
        IT1 A1 = A0; A1.ShiftP(stepj);
        IT1 A2 = A1; A2.ShiftP(stepj);
        IT1 A3 = A2; A3.ShiftP(stepj);
        IT3 Y = v3.Flatten().begin();

        enum { c1 = M1::mconj };

        do {
          A00 = *A0++; A01 = *A1++; A02 = *A2++; A03 = *A3++;
          Maybe<add>::set(rY0 , *Y);
          Maybe<add>::set(iY0 , Y[1]);
          Maybe<add>::add(rY0 , ZProd<c1,false>::rprod(A00,X0));
          Maybe<add>::add(iY0 , ZProd<c1,false>::iprod(A00,X0));
          rY0 += ZProd<c1,false>::rprod(A01,X1);
          iY0 += ZProd<c1,false>::iprod(A01,X1);
          rY0 += ZProd<c1,false>::rprod(A02,X2);
          iY0 += ZProd<c1,false>::iprod(A02,X2);
          rY0 += ZProd<c1,false>::rprod(A03,X3);
          iY0 += ZProd<c1,false>::iprod(A03,X3);
          *Y = rY0; Y[1] = iY0;  Y += 2;
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
      int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 21: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
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

        IT1 Ai0 = m1.get_row(0).NonConj().begin();
        const IT2 X0 = v2.NonConj().begin();
        IT3 Y = v3.begin();
        const int Astepi = m1.stepi();

        do {
          Yi = MultVV_Helper<-1,rs,M1r,V2>::call2(N,Ai0,X0);
          Maybe<add>::add(*Y++, x * Yi);
          Ai0.ShiftP(Astepi);
        } while (--M);
      }
    }
  };

  // algo 22: row major, 4 rows at a time
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<22,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void loop_4_rows(const int M, const int N,
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
      IT1 A1 = A0; A1.ShiftP(stepi);
      IT1 A2 = A1; A2.ShiftP(stepi);
      IT1 A3 = A2; A3.ShiftP(stepi);
      const IT2 X_begin = v2.begin();
      IT2 X = X_begin;
      const int N_4 = (N>>2);
      const int Nx = N-(N_4<<2);
      IT3 Y = v3.begin();

      const bool dopref = N * sizeof(T1) >= TMV_Q3;

      Prefetch_MultiRead(X.GetP());
      Prefetch_Read(A0.GetP());
      Prefetch_Read(A1.GetP());
      Prefetch_Read(A2.GetP());
      Prefetch_Read(A3.GetP());
      Prefetch_Write(Y.GetP());

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
          A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; A0 += 4;
          A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; A1 += 4;
          A20 = A2[0]; A21 = A2[1]; A22 = A2[2]; A23 = A2[3]; A2 += 4;
          A30 = A3[0]; A31 = A3[1]; A32 = A3[2]; A33 = A3[3]; A3 += 4;
          Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
          Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
          Y2 += A20 * X0 + A21 * X1 + A22 * X2 + A23 * X3;
          Y3 += A30 * X0 + A31 * X1 + A32 * X2 + A33 * X3; 
        } while (--j);
        j=Nx; if (j) do {
          X0 = *X++; A00 = *A0++; A10 = *A1++; A20 = *A2++; A30 = *A3++;
          Y0 += A00 * X0; 
          Y1 += A10 * X0; 
          Y2 += A20 * X0; 
          Y3 += A30 * X0;
        } while (--j);
        A0 += stepi_4; A1 += stepi_4; A2 += stepi_4; A3 += stepi_4;
        if (dopref) {
          Prefetch_Read(A0.GetP());
          Prefetch_Read(A1.GetP());
          Prefetch_Read(A2.GetP());
          Prefetch_Read(A3.GetP());
        }
        Maybe<add && (ix!=1)>::add(Y[0] , x * Y0);
        Maybe<add && (ix!=1)>::add(Y[1] , x * Y1);
        Maybe<add && (ix!=1)>::add(Y[2] , x * Y2);
        Maybe<add && (ix!=1)>::add(Y[3] , x * Y3);
        Y += 4;
      } while (--M_4);
    }
    template <class M1x, class V3x>
    static inline void cleanup(const int mb,
        const Scaling<ix,T>& x, const M1x& m1, const V2& v2, V3x& v3)
    {
      TMVAssert(mb == 1 || mb == 2 || mb == 3 || mb == 4);
#ifdef TMV_OPT_CLEANUP
      switch (mb) {
        case 1 : 
          MultMV_Helper<23,1,rs,add,ix,T,M1x,V2,V3x>::call(x,m1,v2,v3);
          break;
        case 2 :
          MultMV_Helper<23,2,rs,add,ix,T,M1x,V2,V3x>::call(x,m1,v2,v3);
          break;
        case 3 :
          MultMV_Helper<23,3,rs,add,ix,T,M1x,V2,V3x>::call(x,m1,v2,v3);
          break;
        case 4 :
          MultMV_Helper<23,4,rs,add,ix,T,M1x,V2,V3x>::call(x,m1,v2,v3);
      }
#else
      MultMV_Helper<21,UNKNOWN,rs,add,ix,T,M1x,V2,V3x>::call(x,m1,v2,v3);
#endif
    }
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 22: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);
      if (N) {
        if (M > 4) {
          const int ma = ((M>>2)<<2);
          const int mb = M - ma;
          typedef typename M1::const_rows_type M1r;
          typedef typename V3::subvector_type V3s;
          loop_4_rows(ma,N,x,m1,v2,v3);
          if (mb) {
            M1r m1b = m1.CRows(ma,M);
            V3s v3b = v3.CSubVector(ma,M);
            cleanup(mb,x,m1b,v2,v3b);
          }
        } else if (M) {
          cleanup(M,x,m1,v2,v3);
        }
      }
      else Maybe<!add>::zero(v3);
    }
  };

  // algo 23: cs <= 4, do all rows at once
  // cs == 1
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<23,1,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { MultMV_Helper<1,1,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3); }
  };
  // cs == 2
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<23,2,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 23 M==2: M,N,cs,rs,x = "<<2<<','<<N<<','<<2<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);

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
        IT1 A1 = A0; A1.ShiftP(stepi);
        IT2 X = v2.begin();
        int N_4 = N>>2;
        int Nx = N-(N_4<<2);

        Y0 = Maybe<add && (ix==1)>::select( v3.cref(0) , T3(0) );
        Y1 = Maybe<add && (ix==1)>::select( v3.cref(1) , T3(0) );
        if (N_4) do {
          X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
          A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; A0 += 4;
          A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; A1 += 4;
          Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
          Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
        } while (--N_4);
        if (Nx) do {
          X0 = *X++; A00 = *A0++; A10 = *A1++; 
          Y0 += A00 * X0; 
          Y1 += A10 * X0; 
        } while (--Nx);
        Maybe<add && (ix!=1)>::add(v3.ref(0) , x * Y0);
        Maybe<add && (ix!=1)>::add(v3.ref(1) , x * Y1);
      }
      else Maybe<!add>::zero(v3);
    }
  };
  // cs == 3
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<23,3,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 23 M==3: M,N,cs,rs,x = "<<3<<','<<N<<','<<3<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);

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
        IT1 A1 = A0; A1.ShiftP(stepi);
        IT1 A2 = A1; A2.ShiftP(stepi);
        IT2 X = v2.begin();
        int N_4 = N>>2;
        int Nx = N-(N_4<<2);

        Y0 = Maybe<add && (ix==1)>::select( v3.cref(0) , T3(0) );
        Y1 = Maybe<add && (ix==1)>::select( v3.cref(1) , T3(0) );
        Y2 = Maybe<add && (ix==1)>::select( v3.cref(2) , T3(0) );
        if (N_4) do {
          X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
          A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; A0 += 4;
          A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; A1 += 4;
          A20 = A2[0]; A21 = A2[1]; A22 = A2[2]; A23 = A2[3]; A2 += 4;
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
        Maybe<add && (ix!=1)>::add(v3.ref(0) , x * Y0);
        Maybe<add && (ix!=1)>::add(v3.ref(1) , x * Y1);
        Maybe<add && (ix!=1)>::add(v3.ref(2) , x * Y2);
      }
      else Maybe<!add>::zero(v3);
    }
  };
  // cs == 4
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<23,4,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 23 M==4: M,N,cs,rs,x = "<<4<<','<<N<<','<<4<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::visreal);

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
        IT1 A1 = A0; A1.ShiftP(stepi);
        IT1 A2 = A1; A2.ShiftP(stepi);
        IT1 A3 = A2; A3.ShiftP(stepi);
        IT2 X = v2.begin();
        int N_4 = N>>2;
        int Nx = N-(N_4<<2);

        Y0 = Maybe<add && (ix==1)>::select( v3.cref(0) , T3(0) );
        Y1 = Maybe<add && (ix==1)>::select( v3.cref(1) , T3(0) );
        Y2 = Maybe<add && (ix==1)>::select( v3.cref(2) , T3(0) );
        Y3 = Maybe<add && (ix==1)>::select( v3.cref(3) , T3(0) );
        if (N_4) do {
          X0 = X[0]; X1 = X[1]; X2 = X[2]; X3 = X[3]; X += 4;
          A00 = A0[0]; A01 = A0[1]; A02 = A0[2]; A03 = A0[3]; A0 += 4;
          A10 = A1[0]; A11 = A1[1]; A12 = A1[2]; A13 = A1[3]; A1 += 4;
          A20 = A2[0]; A21 = A2[1]; A22 = A2[2]; A23 = A2[3]; A2 += 4;
          A30 = A3[0]; A31 = A3[1]; A32 = A3[2]; A33 = A3[3]; A3 += 4;
          Y0 += A00 * X0 + A01 * X1 + A02 * X2 + A03 * X3;
          Y1 += A10 * X0 + A11 * X1 + A12 * X2 + A13 * X3;
          Y2 += A20 * X0 + A21 * X1 + A22 * X2 + A23 * X3;
          Y3 += A30 * X0 + A31 * X1 + A32 * X2 + A33 * X3; 
        } while (--N_4);
        if (Nx) do {
          X0 = *X++; A00 = *A0++; A10 = *A1++; A20 = *A2++; A30 = *A3++; 
          Y0 += A00 * X0; 
          Y1 += A10 * X0; 
          Y2 += A20 * X0; 
          Y3 += A30 * X0;
        } while (--Nx);
        Maybe<add && (ix!=1)>::add(v3.ref(0) , x * Y0);
        Maybe<add && (ix!=1)>::add(v3.ref(1) , x * Y1);
        Maybe<add && (ix!=1)>::add(v3.ref(2) , x * Y2);
        Maybe<add && (ix!=1)>::add(v3.ref(3) , x * Y3);
      }
      else Maybe<!add>::zero(v3);
    }
  };

  // algo 25: row major, 2 rows at a time, complex v3
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<25,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void loop_2_rows(const int M, const int N,
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVAssert(M%2 == 0);
      TMVAssert(M > 0);

      typedef typename M1::value_type T1;
      typedef typename M1::const_row_type::const_nonconj_type::const_iterator IT1;
      T1 A00, A10;

      typedef typename V2::value_type T2;
      typedef typename V2::const_nonconj_type::const_iterator IT2;
      T2 X0;

      typedef typename V3::value_type T3;
      typedef typename V3::real_type RT3;
      typedef typename V3::iterator IT3;
      RT3 rY0, rY1, iY0, iY1;

      int M_2 = M>>1; // M_2 = M/2
      const int stepi = m1.stepi();
      const int stepi_2 = (stepi<<1) - N;
      // step over 2 rows and back to start

      IT1 A0 = m1.get_row(0).NonConj().begin();
      IT1 A1 = A0; A1.ShiftP(stepi);
      const IT2 X_begin = v2.NonConj().begin();
      IT2 X = X_begin;
      IT3 Y = v3.NonConj().begin();

      enum { c1 = M1::mconj };
      enum { c2 = V2::vconj };

      const bool dopref = N * sizeof(T1) >= TMV_Q3;

      Prefetch_MultiRead(X.GetP());
      Prefetch_Read(A0.GetP());
      Prefetch_Read(A1.GetP());
      Prefetch_Write(Y.GetP());

      int j;

      TMVAssert(M_2 > 0);
      do {
        X = X_begin;
        rY0 = Maybe<add && (ix==1)>::select( real(Y[0]) , RT3(0) );
        iY0 = Maybe<add && (ix==1)>::select( imag(Y[0]) , RT3(0) );
        rY1 = Maybe<add && (ix==1)>::select( real(Y[1]) , RT3(0) );
        iY1 = Maybe<add && (ix==1)>::select( imag(Y[1]) , RT3(0) );

        j=N; do {
          X0 = *X++; A00 = *A0++; A10 = *A1++;
          rY0 += ZProd<c1,c2>::rprod(A00,X0);
          iY0 += ZProd<c1,c2>::iprod(A00,X0);
          rY1 += ZProd<c1,c2>::rprod(A10,X0);
          iY1 += ZProd<c1,c2>::iprod(A10,X0);
        } while (--j);
        A0 += stepi_2; A1 += stepi_2; 
        if (dopref) {
          Prefetch_Read(A0.GetP());
          Prefetch_Read(A1.GetP());
        }
        Maybe<add && (ix!=1)>::add(Y[0] , x * T3(rY0,iY0));
        Maybe<add && (ix!=1)>::add(Y[1] , x * T3(rY1,iY1));
        Y += 2;
      } while (--M_2);
    }
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 25: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);
      const int ma = ((M>>1)<<1);
      const int mb = M - ma;
      typedef typename M1::const_row_type M1r;
      if (N) {
        if (ma) loop_2_rows(ma,N,x,m1,v2,v3);
        if (mb) {
          Maybe<add>::add( v3.ref(ma) ,
              x * MultVV_Helper<-1,rs,M1r,V2>::call(m1.get_row(ma),v2) );
        }
      }
      else Maybe<!add>::zero(v3);
    }
  };

  // algo 26: do all rows at once, complex v3: cs <= 4, and must be known
  // cs == 1
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<26,1,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { MultMV_Helper<1,1,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3); }
  };
  // cs == 2
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<26,2,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 26: M,N,cs,rs,x = "<<2<<','<<N<<','<<2<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);

      if (N) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_row_type::const_nonconj_type::const_iterator IT1;
        T1 A00, A10;

        typedef typename V2::value_type T2;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        T2 X0;

        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT3;
        RT3 rY0, iY0, rY1, iY1;

        const int stepi = m1.stepi();

        IT1 A0 = m1.get_row(0).NonConj().begin();
        IT1 A1 = A0; A1.ShiftP(stepi);
        IT2 X = v2.NonConj().begin();

        enum { c1 = M1::mconj };
        enum { c2 = V2::vconj };

        rY0 = Maybe<add && (ix==1)>::select( real(v3.cref(0)) , RT3(0) );
        iY0 = Maybe<add && (ix==1)>::select( imag(v3.cref(0)) , RT3(0) );
        rY1 = Maybe<add && (ix==1)>::select( real(v3.cref(1)) , RT3(0) );
        iY1 = Maybe<add && (ix==1)>::select( imag(v3.cref(1)) , RT3(0) );

        do {
          X0 = *X++; A00 = *A0++; A10 = *A1++;
          rY0 += ZProd<c1,c2>::rprod(A00,X0);
          iY0 += ZProd<c1,c2>::iprod(A00,X0);
          rY1 += ZProd<c1,c2>::rprod(A10,X0);
          iY1 += ZProd<c1,c2>::iprod(A10,X0);
        } while (--N);
        Maybe<add && (ix!=1)>::add(v3.ref(0) , x * T3(rY0,iY0));
        Maybe<add && (ix!=1)>::add(v3.ref(1) , x * T3(rY1,iY1));
      }
      else Maybe<!add>::zero(v3);
    }
  };
  // cs == 3
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<26,3,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 26: M,N,cs,rs,x = "<<3<<','<<N<<','<<3<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);

      if (N) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_row_type::const_nonconj_type::const_iterator IT1;
        T1 A00, A10, A20;

        typedef typename V2::value_type T2;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        T2 X0;

        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT3;
        RT3 rY0, iY0, rY1, iY1, rY2, iY2;

        const int stepi = m1.stepi();

        IT1 A0 = m1.get_row(0).NonConj().begin();
        IT1 A1 = A0; A1.ShiftP(stepi);
        IT1 A2 = A1; A2.ShiftP(stepi);
        IT2 X = v2.NonConj().begin();

        enum { c1 = M1::mconj };
        enum { c2 = V2::vconj };

        rY0 = Maybe<add && (ix==1)>::select( real(v3.cref(0)) , RT3(0) );
        iY0 = Maybe<add && (ix==1)>::select( imag(v3.cref(0)) , RT3(0) );
        rY1 = Maybe<add && (ix==1)>::select( real(v3.cref(1)) , RT3(0) );
        iY1 = Maybe<add && (ix==1)>::select( imag(v3.cref(1)) , RT3(0) );
        rY2 = Maybe<add && (ix==1)>::select( real(v3.cref(2)) , RT3(0) );
        iY2 = Maybe<add && (ix==1)>::select( imag(v3.cref(2)) , RT3(0) );

        do {
          X0 = *X++; A00 = *A0++; A10 = *A1++; A20 = *A2++;
          rY0 += ZProd<c1,c2>::rprod(A00,X0);
          iY0 += ZProd<c1,c2>::iprod(A00,X0);
          rY1 += ZProd<c1,c2>::rprod(A10,X0);
          iY1 += ZProd<c1,c2>::iprod(A10,X0);
          rY2 += ZProd<c1,c2>::rprod(A20,X0);
          iY2 += ZProd<c1,c2>::iprod(A20,X0);
        } while (--N);
        Maybe<add && (ix!=1)>::add(v3.ref(0) , x * T3(rY0,iY0));
        Maybe<add && (ix!=1)>::add(v3.ref(1) , x * T3(rY1,iY1));
        Maybe<add && (ix!=1)>::add(v3.ref(2) , x * T3(rY2,iY2));
      }
      else Maybe<!add>::zero(v3);
    }
  };
  // cs == 4
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<26,4,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 26: M,N,cs,rs,x = "<<4<<','<<N<<','<<4<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(V3::viscomplex);

      if (N) {
        typedef typename M1::value_type T1;
        typedef typename M1::const_row_type::const_nonconj_type::const_iterator IT1;
        T1 A00, A10, A20, A30;

        typedef typename V2::value_type T2;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        T2 X0;

        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT3;
        RT3 rY0, iY0, rY1, iY1, rY2, iY2, rY3, iY3;

        const int stepi = m1.stepi();

        IT1 A0 = m1.get_row(0).NonConj().begin();
        IT1 A1 = A0; A1.ShiftP(stepi);
        IT1 A2 = A1; A2.ShiftP(stepi);
        IT1 A3 = A2; A3.ShiftP(stepi);
        IT2 X = v2.NonConj().begin();

        enum { c1 = M1::mconj };
        enum { c2 = V2::vconj };

        rY0 = Maybe<add && (ix==1)>::select( real(v3.cref(0)) , RT3(0) );
        iY0 = Maybe<add && (ix==1)>::select( imag(v3.cref(0)) , RT3(0) );
        rY1 = Maybe<add && (ix==1)>::select( real(v3.cref(1)) , RT3(0) );
        iY1 = Maybe<add && (ix==1)>::select( imag(v3.cref(1)) , RT3(0) );
        rY2 = Maybe<add && (ix==1)>::select( real(v3.cref(2)) , RT3(0) );
        iY2 = Maybe<add && (ix==1)>::select( imag(v3.cref(2)) , RT3(0) );
        rY3 = Maybe<add && (ix==1)>::select( real(v3.cref(3)) , RT3(0) );
        iY3 = Maybe<add && (ix==1)>::select( imag(v3.cref(3)) , RT3(0) );

        do {
          X0 = *X++; A00 = *A0++; A10 = *A1++; A20 = *A2++; A30 = *A3++;
          rY0 += ZProd<c1,c2>::rprod(A00,X0);
          iY0 += ZProd<c1,c2>::iprod(A00,X0);
          rY1 += ZProd<c1,c2>::rprod(A10,X0);
          iY1 += ZProd<c1,c2>::iprod(A10,X0);
          rY2 += ZProd<c1,c2>::rprod(A20,X0);
          iY2 += ZProd<c1,c2>::iprod(A20,X0);
          rY3 += ZProd<c1,c2>::rprod(A30,X0);
          iY3 += ZProd<c1,c2>::iprod(A30,X0);
        } while (--N);
        Maybe<add && (ix!=1)>::add(v3.ref(0) , x * T3(rY0,iY0));
        Maybe<add && (ix!=1)>::add(v3.ref(1) , x * T3(rY1,iY1));
        Maybe<add && (ix!=1)>::add(v3.ref(2) , x * T3(rY2,iY2));
        Maybe<add && (ix!=1)>::add(v3.ref(3) , x * T3(rY3,iY3));
      }
      else Maybe<!add>::zero(v3);
    }
  };

  // algo 31: fully unroll by rows, apply x to v3
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<31,cs,rs,add,ix,T,M1,V2,V3>
  {
    typedef typename V2::value_type T2;
    typedef typename V3::value_type T3;
    template <int I, int N>
    struct Unroller
    {
      static inline void loadx(const V2& v2, T2* X)
      {
        Unroller<I,N-1>::loadx(v2,X);
        Unroller<N-1,1>::loadx(v2,X);
      }
      static inline void sety(const Scaling<ix,T>& x, const T3* Y, V3& v3)
      {
        Unroller<I,N-1>::sety(x,Y,v3);
        Unroller<N-1,1>::sety(x,Y,v3);
      }
      static inline void calcy(const M1& m1, const T2* X, T3* Y)
      {
        Unroller<I,N-1>::calcy(m1,X,Y);
        Unroller<N-1,1>::calcy(m1,X,Y);
      }
      static inline void calcyi(const int i, const M1& m1, const T2* X, T3& Yi)
      {
        Unroller<I,N-1>::calcyi(i,m1,X,Yi);
        Unroller<N-1,1>::calcyi(i,m1,X,Yi);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void loadx(const V2& v2, T2* X)
      { X[I] = v2.cref(I); }
      static inline void sety(const Scaling<ix,T>& x, const T3* Y, V3& v3)
      { Maybe<add>::add(v3.ref(I) , x * Y[I]); }
      static inline void calcy(const M1& m1, const T2* X, T3* Y)
      { Unroller<0,rs>::calcyi(I,m1,X,Y[I]); }
      static inline void calcyi(const int i, const M1& m1, const T2* X, T3& Yi)
      // Note: "I" is really j here...
      { Maybe<(I>0)>::add(Yi , ZProd<false,false>::prod(m1.cref(i,I),X[I])); }
    };
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
#ifdef PRINTALGO_MV
      std::cout<<"algo 31: cs,rs,x = "<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      T2 X[rs];
      T3 Y[cs];
      Unroller<0,rs>::loadx(v2,X);
      Unroller<0,cs>::calcy(m1,X,Y); 
      Unroller<0,cs>::sety(x,Y,v3);
    }
  };

  // algo 32: fully unroll by rows, apply x to v2
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<32,cs,rs,add,ix,T,M1,V2,V3>
  {
    typedef typename V2::value_type T2;
    typedef typename Traits2<T,T2>::type PT2;
    typedef typename V3::value_type T3;
    template <int I, int N>
    struct Unroller
    {
      static inline void loadx(const Scaling<ix,T>& x, const V2& v2, PT2* X)
      {
        Unroller<I,N-1>::loadx(x,v2,X);
        Unroller<N-1,1>::loadx(x,v2,X);
      }
      static inline void sety(const T3* Y, V3& v3)
      {
        Unroller<I,N-1>::sety(Y,v3);
        Unroller<N-1,1>::sety(Y,v3);
      }
      static inline void calcy(const M1& m1, const PT2* X, T3* Y)
      {
        Unroller<I,N-1>::calcy(m1,X,Y);
        Unroller<N-1,1>::calcy(m1,X,Y);
      }
      static inline void calcyi(const int i, const M1& m1, const PT2* X, T3& Yi)
      {
        Unroller<I,N-1>::calcyi(i,m1,X,Yi);
        Unroller<N-1,1>::calcyi(i,m1,X,Yi);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void loadx(const Scaling<ix,T>& x, const V2& v2, PT2* X)
      { X[I] = x * v2.cref(I); }
      static inline void sety(const T3* Y, V3& v3)
      { Maybe<add>::add(v3.ref(I) , Y[I]); }
      static inline void calcy(const M1& m1, const PT2* X, T3* Y)
      { Unroller<0,rs>::calcyi(I,m1,X,Y[I]); }
      static inline void calcyi(const int i, const M1& m1, const PT2* X, T3& Yi)
      // Note: "I" is really j here...
      { Maybe<(I>0)>::add(Yi , ZProd<false,false>::prod(m1.cref(i,I),X[I])); }
    };
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
#ifdef PRINTALGO_MV
      std::cout<<"algo 32: cs,rs,x = "<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      PT2 X[rs];
      T3 Y[cs];
      Unroller<0,rs>::loadx(x,v2,X);
      Unroller<0,cs>::calcy(m1,X,Y); 
      Unroller<0,cs>::sety(Y,v3);
    }
  };

  // algo 33: fully unroll by rows, apply x to v3
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<33,cs,rs,add,ix,T,M1,V2,V3>
  {
    typedef typename V2::value_type T2;
    typedef typename V3::value_type T3;
    template <int I, int N, bool first>
    struct Unroller
    {
      static inline void loadx(const V2& v2, T2* X)
      {
        Unroller<I,N-1,first>::loadx(v2,X);
        Unroller<N-1,1,first>::loadx(v2,X);
      }
      static inline void sety(const Scaling<ix,T>& x, const T3* Y, V3& v3)
      {
        Unroller<I,N-1,first>::sety(x,Y,v3);
        Unroller<N-1,1,first>::sety(x,Y,v3);
      }
      static inline void calcy(const M1& m1, const T2* X, T3* Y)
      {
        Unroller<I,N-1,first>::calcy(m1,X,Y);
        Unroller<N-1,1,first>::calcy(m1,X,Y);
      }
      static inline void calcyj(const int j, const M1& m1, const T2& Xj, T3* Y)
      {
        Unroller<I,N-1,first>::calcyj(j,m1,Xj,Y);
        Unroller<N-1,1,first>::calcyj(j,m1,Xj,Y);
      }
    };
    template <int I, bool first>
    struct Unroller<I,1,first>
    {
      static inline void loadx(const V2& v2, T2* X)
      { X[I] = v2.cref(I); }
      static inline void sety(const Scaling<ix,T>& x, const T3* Y, V3& v3)
      { Maybe<add>::add(v3.ref(I) , x * Y[I]); }
      static inline void calcy(const M1& m1, const T2* X, T3* Y)
      // Note: "I" is really j here...
      { Unroller<0,cs,I==0>::calcyj(I,m1,X[I],Y); }
      static inline void calcyj(const int j, const M1& m1, const T2& Xj, T3* Y)
      { Maybe<!first>::add(Y[I] , ZProd<false,false>::prod(m1.cref(I,j),Xj)); }
    };
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
#ifdef PRINTALGO_MV
      std::cout<<"algo 33: cs,rs,x = "<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      T2 X[rs];
      T3 Y[cs];
      Unroller<0,rs,false>::loadx(v2,X);
      Unroller<0,rs,false>::calcy(m1,X,Y); 
      Unroller<0,cs,false>::sety(x,Y,v3);
    }
  };

  // algo 34: fully unroll by columnss, apply x to v2
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<34,cs,rs,add,ix,T,M1,V2,V3>
  {
    typedef typename V2::value_type T2;
    typedef typename Traits2<T,T2>::type PT2;
    typedef typename V3::value_type T3;
    template <int I, int N, bool first>
    struct Unroller
    {
      static inline void loadx(const Scaling<ix,T>& x, const V2& v2, PT2* X)
      {
        Unroller<I,N-1,first>::loadx(x,v2,X);
        Unroller<N-1,1,first>::loadx(x,v2,X);
      }
      static inline void sety(const T3* Y, V3& v3)
      {
        Unroller<I,N-1,first>::sety(Y,v3);
        Unroller<N-1,1,first>::sety(Y,v3);
      }
      static inline void calcy(const M1& m1, const PT2* X, T3* Y)
      {
        Unroller<I,N-1,first>::calcy(m1,X,Y);
        Unroller<N-1,1,first>::calcy(m1,X,Y);
      }
      static inline void calcyj(const int j, const M1& m1, const PT2& Xj, T3* Y)
      {
        Unroller<I,N-1,first>::calcyj(j,m1,Xj,Y);
        Unroller<N-1,1,first>::calcyj(j,m1,Xj,Y);
      }
    };
    template <int I, bool first>
    struct Unroller<I,1,first>
    {
      static inline void loadx(const Scaling<ix,T>& x, const V2& v2, PT2* X)
      { X[I] = x * v2.cref(I); }
      static inline void sety(const T3* Y, V3& v3)
      { Maybe<add>::add(v3.ref(I) , Y[I]); }
      static inline void calcy(const M1& m1, const PT2* X, T3* Y)
      // Note: "I" is really j here...
      { Unroller<0,cs,I==0>::calcyj(I,m1,X[I],Y); }
      static inline void calcyj(const int j, const M1& m1, const PT2& Xj, T3* Y)
      { Maybe<!first>::add(Y[I] , ZProd<false,false>::prod(m1.cref(I,j),Xj)); }
    };
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
#ifdef PRINTALGO_MV
      std::cout<<"algo 34: cs,rs,x = "<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      PT2 X[rs];
      T3 Y[cs];
      Unroller<0,rs,false>::loadx(x,v2,X);
      Unroller<0,rs,false>::calcy(m1,X,Y); 
      Unroller<0,cs,false>::sety(Y,v3);
    }
  };

  // algo 41: colmajor, unknown sizes, so figure out which algo to use
  // I used to select an algorithm based on the runtime knowledge
  // of the sizes, but now I think the current 12 and 15 are always
  // the fastest algorithm.  
  // So now this is only relevaant if ix == 0 and !add to call algo 45.
  // First the generic, trivial one:
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<41,cs,rs,add,ix,T,M1,V2,V3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
#ifdef PRINTALGO_MV
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      std::cout<<"algo 41: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = (V3::viscomplex ? 15 : 12) };
      MultMV_Helper<algo2,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };
  // Specialize ix == 0, !add
  template <int cs, int rs, class T, class M1, class V2, class V3>
  struct MultMV_Helper<41,cs,rs,false,0,T,M1,V2,V3>
  {
    static void call(
        const Scaling<0,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_MV
      std::cout<<"algo 41: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = (V3::viscomplex ? 15 : 12) };
#ifdef TMV_OPT_SCALE
      if (M >= N)
        MultMV_Helper<algo2,cs,rs,false,0,T,M1,V2,V3>::call(x,m1,v2,v3);
      else 
#endif
        MultMV_Helper<45,cs,rs,false,0,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  // algo 42: colmajor, ix==0 && add, so might need to copy v3
  template <int cs, int rs, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<42,cs,rs,true,ix,T,M1,V2,V3> // cs = known
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVStaticAssert(TMV_ZeroIX);
      TMVStaticAssert(cs != UNKNOWN);
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MV
      std::cout<<"algo 42: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      // If we have a non-trivial scale and the size of 
      // v2 is (significantly) larger than v3 and add = true, 
      // then it is worth doing m1*v2 first and then applying x.

      // On the other hand, if add = false, then we can get around
      // the temporary by multiplying by x after doing m1*v2.
      // This is done in algo 45, but I mention it here to 
      // explain why we have the add=true requirement.

#ifdef TMV_OPT_SCALE
      if ( N > TMV_Q4 * M )  
      {
#endif
#ifdef PRINTALGO_MV
        std::cout<<"v3c = m1*v2; v3 (+=) x*v3c\n";
#endif
        typedef typename V3::value_type T3;
        typedef typename Traits2<T1,T2>::type T3c;
        typedef SmallVector<T3c,cs> V3c;
        V3c v3c;
        typedef typename V3c::view_type V3cv;
        typedef typename V3c::const_view_type V3ccv;
        V3cv v3cv = v3c.View();
        V3ccv v3ccv = v3c.View();
        MultMV_Helper<41,cs,rs,false,1,RT,M1,V2,V3cv>::call(one,m1,v2,v3cv);
        AddVV_Helper<-1,cs,true,ix,T,V3ccv,V3>::call(x,v3ccv,v3);
#ifdef TMV_OPT_SCALE
      }
      else MultMV_Helper<41,cs,rs,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
    }
  };
  template <int rs, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<42,UNKNOWN,rs,true,ix,T,M1,V2,V3> // cs = unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVStaticAssert(TMV_ZeroIX);
      const int M = m1.colsize();
      const int cs = UNKNOWN;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MV
      std::cout<<"algo 42: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
#ifdef TMV_OPT_SCALE
      if ( N > TMV_Q4 * M )  
      {
#endif
#ifdef PRINTALGO_MV
        std::cout<<"v3c = m1*v2; v3 (+=) x*v3c\n";
#endif
        typedef typename V3::value_type T3;
        typedef typename Traits2<T1,T2>::type T3c;
        // This is the only difference between the known and unknown versions:
        typedef Vector<T3c> V3c;
        V3c v3c(M);
        typedef typename V3c::view_type V3cv;
        typedef typename V3c::const_view_type V3ccv;
        V3cv v3cv = v3c.View();
        V3ccv v3ccv = v3c.View();
        MultMV_Helper<41,cs,rs,false,1,RT,M1,V2,V3cv>::call(one,m1,v2,v3cv);
        AddVV_Helper<-1,cs,true,ix,T,V3ccv,V3>::call(x,v3ccv,v3);
#ifdef TMV_OPT_SCALE
      }
      else MultMV_Helper<41,cs,rs,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
    }
  };

  // algo 43: colmajor, v3.step() != 1, so copy v3
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<43,cs,rs,add,ix,T,M1,V2,V3> // cs known
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MV
      std::cout<<"algo 43: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      if ( N >= M) 
      {
#ifdef PRINTALGO_MV
        std::cout<<"v3c = m1*v2; v3 (+=) x*v3c\n";
#endif
        typedef typename Traits2<T1,T2>::type T3c;
        typedef SmallVector<T3c,cs> V3c;
        V3c v3c;
        typedef typename V3c::view_type V3cv;
        typedef typename V3c::const_view_type V3ccv;
        V3cv v3cv = v3c.View();
        V3ccv v3ccv = v3c.View();
        MultMV_Helper<41,cs,rs,false,1,RT,M1,V2,V3cv>::call(one,m1,v2,v3cv);
        AddVV_Helper<-1,cs,add,ix,T,V3ccv,V3>::call(x,v3ccv,v3);
      }
      else 
      {
        enum { algo2 = TMV_ZeroIX ? 45 : 41 };
#ifdef PRINTALGO_MV
        std::cout<<"v3c = x*m1*v2; v3 (+=) v3c\n";
#endif
        typedef typename Traits2<T,typename Traits2<T1,T2>::type>::type T3c;
        typedef SmallVector<T3c,cs> V3c;
        V3c v3c;
        typedef typename V3c::view_type V3cv;
        typedef typename V3c::const_view_type V3ccv;
        V3cv v3cv = v3c.View();
        V3ccv v3ccv = v3c.View();
        MultMV_Helper<algo2,cs,rs,false,ix,T,M1,V2,V3cv>::call(x,m1,v2,v3cv);
        AddVV_Helper<-1,cs,add,1,RT,V3ccv,V3>::call(one,v3ccv,v3);
      }
    }
  };
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<43,UNKNOWN,rs,add,ix,T,M1,V2,V3> // cs unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = m1.colsize();
      const int cs = UNKNOWN;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MV
      std::cout<<"algo 43: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      if ( N >= M) 
      {
#ifdef PRINTALGO_MV
        std::cout<<"v3c = m1*v2; v3 (+=) x*v3c\n";
#endif
        typedef typename Traits2<T1,T2>::type T3c;
        typedef Vector<T3c> V3c;
        V3c v3c(M);
        typedef typename V3c::view_type V3cv;
        typedef typename V3c::const_view_type V3ccv;
        V3cv v3cv = v3c.View();
        V3ccv v3ccv = v3c.View();
        MultMV_Helper<41,cs,rs,false,1,RT,M1,V2,V3cv>::call(one,m1,v2,v3cv);
        AddVV_Helper<-1,cs,add,ix,T,V3ccv,V3>::call(x,v3ccv,v3);
      }
      else 
      {
        enum { algo2 = TMV_ZeroIX ? 45 : 41 };
#ifdef PRINTALGO_MV
        std::cout<<"v3c = x*m1*v2; v3 (+=) v3c\n";
#endif
        typedef typename Traits2<T,typename Traits2<T1,T2>::type>::type T3c;
        typedef Vector<T3c> V3c;
        V3c v3c(M);
        typedef typename V3c::view_type V3cv;
        typedef typename V3c::const_view_type V3ccv;
        V3cv v3cv = v3c.View();
        V3ccv v3ccv = v3c.View();
        MultMV_Helper<algo2,cs,rs,false,ix,T,M1,V2,V3cv>::call(x,m1,v2,v3cv);
        AddVV_Helper<-1,cs,add,1,RT,V3ccv,V3>::call(one,v3ccv,v3);
      }
    }
  };

  // algo 44: colmajor, unknown cs -- check if it is small
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<44,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      enum { algo1 = (
          V3::vstep != 1 ? 43 : ( TMV_ZeroIX && add ) ? 42 : 41 ) };
      enum { algo2 = (
          V2::vstep == 1 ? ( V3::viscomplex ? 26 : 23 ) : algo1 ) };

#ifdef TMV_OPT_SMALL
      TMVStaticAssert(cs == UNKNOWN);
      const int M = m1.colsize();
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
#ifdef PRINTALGO_MV
      std::cout<<"algo 44: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      if (M <= 4) 
      {
        // then it is worth figuring out what M is.
        switch (M)
        {
          case 0 :
            // do nothing
            break;
          case 1 :
            MultMV_Helper<1,1,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            break;
          case 2 :
            MultMV_Helper<algo2,2,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            break;
          case 3 :
            MultMV_Helper<algo2,3,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            break;
          case 4 :
            MultMV_Helper<algo2,4,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
      } 
      else MultMV_Helper<algo1,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#else
      MultMV_Helper<algo1,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
    }
  };

  // algo 45: column major, apply x at the end
  template <int cs, int rs, class T, class M1, class V2, class V3>
  struct MultMV_Helper<45,cs,rs,false,0,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<0,T>& x, const M1& m1, const V2& v2, V3& v3)
    { 
      TMVAssert(m1.colsize() < m1.rowsize());
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      typedef typename M1::value_type T1;
      typedef typename V3::value_type T3;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MV
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      std::cout<<"algo 45: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      MultMV_Helper<41,cs,rs,false,1,RT,M1,V2,V3>::call(one,m1,v2,v3);
      MultXV_Helper<-1,cs,0,T,V3>::call(x,v3);
    }
  };


  // algo 51: rowmajor, unknown sizes, so figure out which algo to use
  // This is unnecessary now that I've improved algo 22 to always be
  // as fast or faster than 21, but I'm leaving it hear in case
  // I find another algorithm that is faster for certain sizes.
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<51,cs,rs,add,ix,T,M1,V2,V3> 
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
#ifdef PRINTALGO_MV
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      std::cout<<"algo 51: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = (V3::viscomplex ? 25 : 22) };
      MultMV_Helper<algo2,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  // algo 52: rowmajor, might need to copy vector(s)
  template <int cs, int rs, int ix, bool add, class T, class M1, class V2, class V3>
  struct MultMV_Helper<52,cs,rs,add,ix,T,M1,V2,V3> // rs = known
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVStaticAssert(TMV_ZeroIX);
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
#ifdef PRINTALGO_MV
      std::cout<<"algo 52: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      // If we have a non-trivial scale and the size of 
      // v3 is (significantly) larger than v2, then it is worth
      // making x*v2 as a temporary.

      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef TMV_OPT_SCALE
      if ( M > TMV_Q4 * N ) 
      {
#endif
        typedef typename Traits2<T,T2>::type T2c;
        typedef SmallVector<T2c,rs> V2c;
        V2c v2c;
        typedef typename V2c::view_type V2cv;
        typedef typename V2c::const_view_type V2ccv;
        V2cv v2cv = v2c.View();
        V2ccv v2ccv = v2c.View();
        AddVV_Helper<-1,rs,false,ix,T,V2,V2cv>::call(x,v2,v2cv);
        MultMV_Helper<51,cs,rs,add,1,RT,M1,V2ccv,V3>::call(one,m1,v2ccv,v3);
#ifdef TMV_OPT_SCALE
      }
      else MultMV_Helper<51,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
    }
  };
  template <int cs, int ix, bool add, class T, class M1, class V2, class V3>
  struct MultMV_Helper<52,cs,UNKNOWN,add,ix,T,M1,V2,V3> // rs = unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      TMVStaticAssert(TMV_ZeroIX);
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = m1.rowsize();
      const int rs = UNKNOWN;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
#ifdef PRINTALGO_MV
      std::cout<<"algo 52: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef TMV_OPT_SCALE
      if ( M > TMV_Q4 * N ) 
      {
#endif
        typedef typename Traits2<T,T2>::type T2c;
        typedef Vector<T2c> V2c;
        V2c v2c(N);
        typedef typename V2c::view_type V2cv;
        typedef typename V2c::const_view_type V2ccv;
        V2cv v2cv = v2c.View();
        V2ccv v2ccv = v2c.View();
        AddVV_Helper<-1,rs,false,ix,T,V2,V2cv>::call(x,v2,v2cv);
        MultMV_Helper<51,cs,rs,add,1,RT,M1,V2ccv,V3>::call(one,m1,v2ccv,v3);
#ifdef TMV_OPT_SCALE
      }
      else MultMV_Helper<51,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
    }
  };

  // algo 53: rowmajor, v2.step() != 1, so copy v2
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<53,cs,rs,add,ix,T,M1,V2,V3> // rs = known
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
#ifdef PRINTALGO_MV
      std::cout<<"algo 53: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      // For rowmajor algorithm, if v2 is large, we want to make sure
      // it has unit step, even if this requires a temporary.
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
      if ( M > N ) 
      {
        typedef typename Traits2<T,T2>::type T2c;
        typedef SmallVector<T2c,rs> V2c;
        V2c v2c;
        typedef typename V2c::view_type V2cv;
        typedef typename V2c::const_view_type V2ccv;
        V2cv v2cv = v2c.View();
        V2ccv v2ccv = v2c.View();
        AddVV_Helper<-1,rs,false,ix,T,V2,V2cv>::call(x,v2,v2cv);
        MultMV_Helper<51,cs,rs,add,1,RT,M1,V2ccv,V3>::call(one,m1,v2ccv,v3);
      }
      else 
      {
        typedef SmallVector<T2,rs> V2c;
        V2c v2c;
        typedef typename V2c::view_type V2cv;
        typedef typename V2c::const_view_type V2ccv;
        V2cv v2cv = v2c.View();
        V2ccv v2ccv = v2c.View();
        CopyV_Helper<-1,rs,V2,V2cv>::call(v2,v2cv);
        MultMV_Helper<51,cs,rs,add,ix,T,M1,V2ccv,V3>::call(x,m1,v2ccv,v3);
      }
    }
  };
  template <int cs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<53,cs,UNKNOWN,add,ix,T,M1,V2,V3> // rs = unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = m1.rowsize();
      const int rs = UNKNOWN;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
#ifdef PRINTALGO_MV
      std::cout<<"algo 53: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
      if ( M > N ) 
      {
        typedef typename Traits2<T,T2>::type T2c;
        typedef Vector<T2c> V2c;
        V2c v2c(N);
        typedef typename V2c::view_type V2cv;
        typedef typename V2c::const_view_type V2ccv;
        V2cv v2cv = v2c.View();
        V2ccv v2ccv = v2c.View();
        AddVV_Helper<-1,rs,false,ix,T,V2,V2cv>::call(x,v2,v2cv);
        MultMV_Helper<51,cs,rs,add,1,RT,M1,V2ccv,V3>::call(one,m1,v2ccv,v3);
      }
      else 
      {
        typedef Vector<T2> V2c;
        V2c v2c(N);
        typedef typename V2c::view_type V2cv;
        typedef typename V2c::const_view_type V2ccv;
        V2cv v2cv = v2c.View();
        V2ccv v2ccv = v2c.View();
        CopyV_Helper<-1,rs,V2,V2cv>::call(v2,v2cv);
        MultMV_Helper<51,cs,rs,add,ix,T,M1,V2ccv,V3>::call(x,m1,v2ccv,v3);
      }
    }
  };

  // algo 54: rowmajor, unknown rs -- see if it is small
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<54,cs,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      enum { algo1 = (
          V2::vstep != 1 ? 53 : TMV_ZeroIX ? 52 : 51 ) };
      enum { algo2 = (
          V3::vstep == 1 ? ( V3::viscomplex ? 16 : 13 ) : algo1 ) };
#ifdef TMV_OPT_SMALL
      TMVStaticAssert(rs == UNKNOWN);
      const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
      const int N = m1.rowsize();
      typedef typename M1::value_type T1;
#ifdef PRINTALGO_MV
      std::cout<<"algo 54: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      if (N <= 4) 
      {
        // then it is worth figuring out what N is.
        switch (N)
        {
          case 0 :
            // If M > 0, and !add, then we need to do v3.Zero().
            Maybe<!add>::zero(v3);
            break;
          case 1 :
            MultMV_Helper<2,cs,1,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            break;
          case 2 :
            MultMV_Helper<algo2,cs,2,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            break;
          case 3 :
            MultMV_Helper<algo2,cs,3,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            break;
          case 4 :
            MultMV_Helper<algo2,cs,4,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
      }
      else MultMV_Helper<algo1,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#else
      MultMV_Helper<algo1,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#endif
    }
  };

  // algo -3: The same as -1, but allow for the possibility of using
  // InstMultMV.  This is used when calling MultMV from a degenerate MultMM.
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<-3,cs,rs,add,ix,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef typename V3::value_type T3;
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
          cs == UNKNOWN && rs == UNKNOWN ) };

      enum { algo = inst ? 3 : -1 };
      MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  // algo -2: Determine which algorithm to use, but no if statements, 
  //          and no unrolling.
  // Used by MultMM to avoid if clauses that we already know don't apply.
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<-2,cs,rs,add,ix,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename M1::value_type T1;
      enum { algo = (
          ( rs == 0 || cs == 0 ) ? 0 : // trivial - nothing to do
          ( cs == 1 ) ? 1 : // trivial - cs = 1
          ( rs == 1 ) ? 2 : // trivial - rs = 1
          M1::mcolmajor ? ( // colmajor
            ( cs != UNKNOWN && cs <= 4 && V2::vstep == 1 ) ? (
              (V3::viscomplex ? 26 : 23) ) :
            ( rs != UNKNOWN && rs <= 4 ) ? (
              (V3::viscomplex ? 16 : 13) ) :
            V3::viscomplex ? 15 : 12 ) :
          M1::mrowmajor ? ( // rowmajor
            ( rs != UNKNOWN && rs <= 4 && V3::vstep == 1 ) ? (
              (V3::viscomplex ? 16 : 13 ) ) : 
            ( cs != UNKNOWN && cs <= 4 ) ? (
              (V3::viscomplex ? 26 : 23) ) :
            V3::viscomplex ? 25 : 22 ) :
          // nomajor -- don't do anything fancy
          V2::vstep == 1 ? 21 : V3::vstep == 1 ? 11 : 21 ) };
      MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultMV_Helper<-1,cs,rs,add,ix,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename M1::value_type T1;
      // Possible algorithms to choose from:
      // I spread out the numbers a bit to allow room for more algorithms.
      // Each decade is a set of related algorithms.
      //
      // Trivial:
      //  0 = cs or rs == 0, so nothing to do
      //  1 = cs == 1: reduces to trivial AddVV function
      //  2 = rs == 1: reduces to trivial MultVV function
      //  3 = call InstMultMV
      //
      // Column Major:
      // 11 = column major, simple for loop
      // 12 = column major, 4 columns at a time
      // 13 = column major, rs is known and <= 4
      // 15 = column major, complex v3, 2 columns at a time
      // 16 = column major, complex v3, rs is known and <= 4
      //
      // Row Major:
      // 21 = row major, simple for loop
      // 22 = row major, 4 rows at a time
      // 23 = row major, cs is known and <= 4
      // 25 = row major, complex v3, 2 rows at a time
      // 26 = row major, complex v3, cs is known and <= 4
      //
      // Fully Unrolled:
      // 31 = fully unroll by rows, apply x to v3
      // 32 = fully unroll by rows, apply x to v2
      // 33 = fully unroll by columns, apply x to v3
      // 34 = fully unroll by columns, apply x to v2
      //
      // Column Major, meta algorithms
      // 41 = column major, unknown sizes, so figure out which algo to use
      // 42 = column major, ix==0 && add, so might need to copy v3
      // 43 = column major, v3.step() != 1, so copy v3
      // 44 = column major, unknown cs, check if it is small
      // 45 = column major, !add, apply x at the end
      //
      // Row Major, meta algorithms
      // 51 = row major, unknown sizes, so figure out which algo to use
      // 52 = row major, ix==0, so might need to copy v2
      // 53 = row major, v2.step() != 1, so copy v2
      // 54 = row major, unknown rs, check if it is small
#if TMV_OPT == 0
      enum { algo = M1::mcolmajor ? 11 : 21 };
#else
      enum { unroll = (
          ( cs == UNKNOWN || rs == UNKNOWN ) ? false :
          IntTraits2<cs,rs>::prod > TMV_Q1 ? false :
          // These maxunroll values are empirical for my machine,
          // comparing the speed for the fully unrolled algorithm
          // to the regular, non-unrolled algorithm, so they might
          // not be optimal for other machines, or even other
          // compilers, but they should be good enough.
          // And for that matter, I only tested real double, so they
          // might not be best for other data types either.
          M1::mcolmajor ? (
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
            false ) ) } ;
      enum { algo = (
          ( rs == 0 || cs == 0 ) ? 0 : // trivial - nothing to do
          ( cs == 1 ) ? 1 : // trivial - cs = 1
          ( rs == 1 ) ? 2 : // trivial - rs = 1
          M1::mcolmajor ? ( // colmajor
            cs == UNKNOWN ? 44 :
            rs == UNKNOWN ? (
              V3::vstep != 1 ? 43 : ( TMV_ZeroIX && add ) ? 42 : 41 ) : 
            unroll ? (
              ( cs <= TMV_Q5 ? 
                ( (rs<cs) ? 32 : 31 ) :
                ( (rs<cs) ? 34 : 33 ) ) ) :
            ( cs <= 4 && V2::vstep == 1 ) ? (V3::viscomplex ? 26 : 23) :
            ( rs > TMV_Q2 && V3::vstep != 1 ) ? 43 :
            ( rs > IntTraits2<TMV_Q4,cs>::prod && TMV_ZeroIX && add ) ? 42 :
            ( !add && TMV_ZeroIX && rs > cs ) ? 45 : 
            ( rs <= 4 ) ? (V3::viscomplex ? 16 : 13) :
            V3::viscomplex ? 15 : 12 ) :
          M1::mrowmajor ? ( // rowmajor
            rs == UNKNOWN ? 54 :
            cs == UNKNOWN ? (
              ( V2::vstep != 1 ) ? 53 : ( TMV_ZeroIX ) ? 52 : 51 ) : 
            unroll ? ( (rs<cs) ? 32 : 31 ) :
            ( rs <= 4 && V3::vstep == 1 ) ? (V3::viscomplex ? 16 : 13 ) : 
            ( rs > TMV_Q2 && V2::vstep != 1 ) ? 53 :
            ( cs > IntTraits2<TMV_Q4,rs>::prod && TMV_ZeroIX ) ? 52 :
            ( cs <= 4 ) ? (V3::viscomplex ? 26 : 23) :
            V3::viscomplex ? 25 : 22 ) :
          // nomajor -- don't do anything fancy
          V2::vstep == 1 ? 21 : V3::vstep == 1 ? 11 : 21 ) };
#endif
#ifdef PRINTALGO_MV
      std::cout<<"InlineMultMV: \n";
      std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      std::cout<<"v2 = "<<TypeText(v2)<<std::endl;
      std::cout<<"v3 = "<<TypeText(v3)<<std::endl;
      std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
#endif
      MultMV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  template <bool add, int ix, class T, class M1, class V2, class V3>
  inline void InlineMultMV(const Scaling<ix,T>& x, 
      const BaseMatrix_Rec<M1>& m1, const BaseVector_Calc<V2>& v2, 
      BaseVector_Mutable<V3>& v3)
  {
    //std::cout<<"InlineMultMV "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"add = "<<add<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    //std::cout<<"v3 = "<<TypeText(v3)<<"  "<<v3<<std::endl;
    TMVStaticAssert((Sizes<M1::mcolsize,V3::vsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,V2::vsize>::same));
    TMVAssert(m1.colsize() == v3.size());
    TMVAssert(m1.rowsize() == v2.size());
    typedef typename M1::value_type T1;
    enum { cs = Sizes<M1::mcolsize,V3::vsize>::size };
    enum { rs = Sizes<M1::mrowsize,V2::vsize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename V2::const_view_type V2v;
    typedef typename V3::view_type V3v;
    M1v m1v = m1.View();
    V2v v2v = v2.View();
    V3v v3v = v3.View();
    MultMV_Helper<-1,cs,rs,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    //std::cout<<"End InlineMultMV: v3 = "<<v3<<std::endl;
  }

  // Defined in TMV_MultMV.cpp
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstMultMV(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddMultMV(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

} // namespace tmv

#undef TMV_OPT_SMALL
#undef TMV_OPT_CLEANUP
#undef TMV_OPT_SCALE

#undef TMV_Q1
#undef TMV_Q2
#undef TMV_Q3
#undef TMV_Q4
#undef TMV_Q5
#undef TMV_ZeroIX

#endif 
