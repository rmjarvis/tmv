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
// Boston, MA  02320-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_Rank1Update_VVM_H
#define TMV_Rank1Update_VVM_H

#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_AddVV.h"
#include "TMV_Vector.h"
#include "TMV_SmallVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Prefetch.h"

// TMV_OPT_SCALE determines whether to use the Q4 parameter to 
// determine whether to apply the scaling to v1 or v2.
#if TMV_OPT >= 2
#define TMV_OPT_SCALE
#endif

//#define PRINTALGO

#ifdef PRINTALGO
#include <iostream>
#endif

namespace tmv {

  //
  // Vector ^ Vector
  //

  // There are a number of values used in the algorithm selection
  // that are either arbitrary or empirical.
  // So I put them all here to make them easier to change and to
  // track down in the code.

  // Q1 is the maximum value of cs*rs for which we fully unroll.
  // It seems like unrolling is never faster, so I set this to 0, 
  // but I'm leaving the code here in case some improvement in the 
  // unrolled code changes that.
#define TMV_Q1 0

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

  // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
  // It doesn't really seem to matter much either way.
#define TMV_ZeroIX (ix == 0)
//#define TMV_ZeroIX (ix != 1)

  template <int algo, int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper;

  // algo 0: cs or rs = 0, so nothing to do
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<0,cs,rs,add,ix,T,V1,V2,M3>
  { static void call(const Scaling<ix,T>& , const V1& , const V2& , M3& ) {} };

  // algo 1: cs == 1, so simplifies to an AddVV or MultXV function
  template <int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<1,1,rs,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename M3::row_type M3r;
      M3r m30 = m3.get_row(0);
      typedef typename Traits2<T,typename V1::value_type>::type PT1;
      AddVV_Helper<-1,rs,add,0,PT1,V2,M3r>::call(x*v1.cref(0),v2,m30); 
    }
  };

  // algo 2: rs == 1, so simplifies to an AddVV or MultXV function
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<2,cs,1,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename M3::col_type M3c;
      typedef typename Traits2<T,typename V2::value_type>::type PT2;
      M3c m30 = m3.get_col(0);
      AddVV_Helper<-1,cs,add,0,PT2,V1,M3c>::call(x*v2.cref(0),v1,m30); 
    }
  };

  // algo 3: call InstRank1Update
  template <int cs, int rs, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<3,cs,rs,false,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename M3::value_type T3;
      typename M3::xview_type m3x = m3.XView();
      InstRank1Update(T3(x),v1.XView(),v2.XView(),m3x);
    }
  };
  template <int cs, int rs, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<3,cs,rs,true,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename M3::value_type T3;
      typename M3::xview_type m3x = m3.XView();
      InstAddRank1Update(T3(x),v1.XView(),v2.XView(),m3x);
    }
  };

  // algo 11: The basic column major loop
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<11,cs,rs,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    { 
      const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO
      std::cout<<"algo 11: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

      typedef typename V1::const_nonconj_type::const_iterator IT1;
      const IT1 X = v1.NonConj().begin();

      typedef typename V2::value_type T2;
      typedef typename Traits2<T,T2>::type PT2;
      typedef typename V2::const_nonconj_type::const_iterator IT2;
      PT2 Y0;
      IT2 Y = v2.NonConj().begin();
      enum { c2 = V2::vconj };

      typedef typename M3::col_type M3c;
      typedef typename M3c::iterator IT3;

      const int stepj = m3.stepj();
      IT3 A0 = m3.get_col(0).begin();

      for(int j=N;j;--j) {
        Y0 = ZProd<false,c2>::prod(x , *Y++);
        AddVV_Helper<-1,cs,add,0,PT2,V1,M3c>::call2(M,Scaling<0,PT2>(Y0),X,A0);
        A0.ShiftP(stepj);
      }
    }
  };

  // algo 12: column major, 4 columns at a time
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<12,cs,rs,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO
      std::cout<<"algo 12: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::misreal);

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
        T3 A00, A01, A02, A03;

        int N4 = N>>2; // N4 = N/4
        int Nb = N-(N4<<2); // N4 = N%4
        const int stepj = m3.stepj();
        const int stepj_4 = (stepj<<2) - M;
        // step over 4 columns and back to start
        const int stepj_1 = stepj - M;

        IT3 A0 = m3.get_col(0).begin();
        IT3 A1 = A0; A1.ShiftP(stepj);
        IT3 A2 = A1; A2.ShiftP(stepj);
        IT3 A3 = A2; A3.ShiftP(stepj);
        IT2 Y = v2.begin();
        const IT1 X_begin = v1.begin();
        IT1 X = X_begin;

        const bool dopref = M * sizeof(T3) >= TMV_Q3;

        Prefetch_Read(Y.GetP());
        Prefetch_MultiRead(X.GetP());
        Prefetch_Write(A0.GetP());
        Prefetch_Write(A1.GetP());
        Prefetch_Write(A2.GetP());
        Prefetch_Write(A3.GetP());

        int i;

        if (N4) do {
          Y0 = x * Y[0]; Y1 = x * Y[1]; Y2 = x * Y[2]; Y3 = x * Y[3]; Y += 4;
          X = X_begin;
          i=M; do {
            X0 = *X++;
            Maybe<add>::add(*A0++ , X0 * Y0);
            Maybe<add>::add(*A1++ , X0 * Y1);
            Maybe<add>::add(*A2++ , X0 * Y2);
            Maybe<add>::add(*A3++ , X0 * Y3);
          } while (--i);
          A0.ShiftP(stepj_4);
          A1.ShiftP(stepj_4);
          A2.ShiftP(stepj_4);
          A3.ShiftP(stepj_4);
          if (dopref) {
            Prefetch_Write(A0.GetP());
            Prefetch_Write(A1.GetP());
            Prefetch_Write(A2.GetP());
            Prefetch_Write(A3.GetP());
          }
        } while (--N4);
        if (Nb) do {
          Y0 = x * *Y++;
          X = X_begin;
          i=M; do {
            X0 = *X++;
            Maybe<add>::add(*A0++ , X0 * Y0);
          } while (--i);
          A0.ShiftP(stepj_1);
          if (dopref) {
            Prefetch_Write(A0.GetP());
          }
        } while (--Nb);
      }
    }
  };

  // algo 13: do all columns at once -- rs <= 4, and must be known
  // rs == 1
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<13,cs,1,add,ix,T,V1,V2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    { Rank1Update_VVM_Helper<2,cs,1,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3); }
  };
  // rs == 2
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<13,cs,2,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO
      std::cout<<"algo 13 N==2: M,N,cs,rs,x = "<<M<<','<<2<<','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::misreal);

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
      IT3 A1 = A0; A1.ShiftP(stepj);

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
  struct Rank1Update_VVM_Helper<13,cs,3,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO
      std::cout<<"algo 13 N==3: M,N,cs,rs,x = "<<M<<','<<3<<','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::misreal);

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
      IT3 A1 = A0; A1.ShiftP(stepj);
      IT3 A2 = A1; A2.ShiftP(stepj);

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
  struct Rank1Update_VVM_Helper<13,cs,4,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO
      std::cout<<"algo 13 N==4: M,N,cs,rs,x = "<<M<<','<<4<<','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::misreal);

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
      IT3 A1 = A0; A1.ShiftP(stepj);
      IT3 A2 = A1; A2.ShiftP(stepj);
      IT3 A3 = A2; A3.ShiftP(stepj);

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

  // algo 16: column major, 2 columns at a time, complex m3
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<16,cs,rs,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO
      std::cout<<"algo 16: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::miscomplex);
      TMVAssert(N%2 == 0);

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
        IT3 A1 = A0; A1.ShiftP(stepj);
        IT2 Y = v2.NonConj().begin();
        const IT1 X_begin = v1.NonConj().begin();
        IT1 X = X_begin;

        enum { c1 = V1::vconj };
        enum { c2 = V2::vconj };

        const bool dopref = M * sizeof(T3) >= TMV_Q3;

        Prefetch_Read(Y.GetP());
        Prefetch_MultiRead(X.GetP());
        Prefetch_Write(A0.GetP());
        Prefetch_Write(A1.GetP());

        int i;

        if (N2) do {
          Y0 = ZProd<false,c2>::prod(x,Y[0]);
          Y1 = ZProd<false,c2>::prod(x,Y[1]);
          Y += 2;
          X = X_begin;
          i=M; do {
            X0 = *X++;
            rA0 = Maybe<add>::sum( real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
            iA0 = Maybe<add>::sum( imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
            rA1 = Maybe<add>::sum( real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
            iA1 = Maybe<add>::sum( imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
            *A0++ = T3(rA0,iA0);
            *A1++ = T3(rA1,iA1);
          } while (--i);
          A0.ShiftP(stepj_2);
          A1.ShiftP(stepj_2);
          if (dopref) {
            Prefetch_Write(A0.GetP());
            Prefetch_Write(A1.GetP());
          }
        } while (--N2);
        if (Nb) {
          Y0 = ZProd<false,c2>::prod(x,*Y);
          X = X_begin;
          i=M; do {
            X0 = *X++;
            rA0 = Maybe<add>::sum( real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
            iA0 = Maybe<add>::sum( imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
            *A0++ = T3(rA0,iA0);
          } while (--i);
        }
      }
    }
  };

  // algo 17: do all columns at once, complex m3: rs <= 4, and must be known
  // rs == 1
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<17,cs,1,add,ix,T,V1,V2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    { Rank1Update_VVM_Helper<2,cs,1,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3); }
  };
  // rs == 2
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<17,cs,2,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO
      std::cout<<"algo 17 N==2: M,N,cs,rs,x = "<<M<<','<<2<<','<<cs<<','<<2<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::miscomplex);

      if (M) {
        typedef typename V1::value_type T1;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        T1 X0;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        enum { c2 = V2::vconj };
        const PT2 Y0 = ZProd<false,c2>::prod(x,v2.NonConj().cref(0));
        const PT2 Y1 = ZProd<false,c2>::prod(x,v2.NonConj().cref(1));

        typedef typename M3::value_type T3;
        typedef typename M3::real_type RT3;
        typedef typename M3::col_type::iterator IT3;
        RT3 rA0, iA0, rA1, iA1;

        const int stepj = m3.stepj();
        IT3 A0 = m3.get_col(0).begin();
        IT3 A1 = A0; A1.ShiftP(stepj);
        const IT1 X_begin = v1.NonConj().begin();
        IT1 X = X_begin;

        enum { c1 = V1::vconj };

        do {
          X0 = *X++;
          rA0 = Maybe<add>::sum( real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
          iA0 = Maybe<add>::sum( imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
          rA1 = Maybe<add>::sum( real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
          iA1 = Maybe<add>::sum( imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
          *A0++ = T3(rA0,iA0);
          *A1++ = T3(rA1,iA1);
        } while (--M);
      }
    }
  };
  // rs == 3
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<17,cs,3,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO
      std::cout<<"algo 17 N==3: M,N,cs,rs,x = "<<M<<','<<3<<','<<cs<<','<<3<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::miscomplex);

      if (M) {
        typedef typename V1::value_type T1;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        T1 X0;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        enum { c2 = V2::vconj };
        const PT2 Y0 = ZProd<false,c2>::prod(x,v2.NonConj().cref(0));
        const PT2 Y1 = ZProd<false,c2>::prod(x,v2.NonConj().cref(1));
        const PT2 Y2 = ZProd<false,c2>::prod(x,v2.NonConj().cref(2));

        typedef typename M3::value_type T3;
        typedef typename M3::real_type RT3;
        typedef typename M3::col_type::iterator IT3;
        RT3 rA0, iA0, rA1, iA1, rA2, iA2;

        const int stepj = m3.stepj();
        IT3 A0 = m3.get_col(0).begin();
        IT3 A1 = A0; A1.ShiftP(stepj);
        IT3 A2 = A1; A2.ShiftP(stepj);
        const IT1 X_begin = v1.NonConj().begin();
        IT1 X = X_begin;

        enum { c1 = V1::vconj };

        do {
          X0 = *X++;
          rA0 = Maybe<add>::sum( real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
          iA0 = Maybe<add>::sum( imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
          rA1 = Maybe<add>::sum( real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
          iA1 = Maybe<add>::sum( imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
          rA2 = Maybe<add>::sum( real(*A2) , ZProd<c1,false>::rprod(X0,Y2) );
          iA2 = Maybe<add>::sum( imag(*A2) , ZProd<c1,false>::iprod(X0,Y2) );
          *A0++ = T3(rA0,iA0);
          *A1++ = T3(rA1,iA1);
          *A2++ = T3(rA2,iA2);
        } while (--M);
      }
    }
  };
  // rs == 4
  template <int cs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<17,cs,4,add,ix,T,V1,V2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO
      std::cout<<"algo 17 N==4: M,N,cs,rs,x = "<<M<<','<<4<<','<<cs<<','<<4<<','<<T(x)<<std::endl;
#endif
      TMVStaticAssert(M3::miscomplex);
      if (M) {

        typedef typename V1::value_type T1;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        T1 X0;

        typedef typename V2::value_type T2;
        typedef typename Traits2<T,T2>::type PT2;
        enum { c2 = V2::vconj };
        const PT2 Y0 = ZProd<false,c2>::prod(x,v2.NonConj().cref(0));
        const PT2 Y1 = ZProd<false,c2>::prod(x,v2.NonConj().cref(1));
        const PT2 Y2 = ZProd<false,c2>::prod(x,v2.NonConj().cref(2));
        const PT2 Y3 = ZProd<false,c2>::prod(x,v2.NonConj().cref(3));

        typedef typename M3::value_type T3;
        typedef typename M3::real_type RT3;
        typedef typename M3::col_type::iterator IT3;
        RT3 rA0, iA0, rA1, iA1;

        const int stepj = m3.stepj();
        IT3 A0 = m3.get_col(0).begin();
        IT3 A1 = A0; A1.ShiftP(stepj);
        IT3 A2 = A1; A2.ShiftP(stepj);
        IT3 A3 = A2; A3.ShiftP(stepj);
        const IT1 X_begin = v1.NonConj().begin();
        IT1 X = X_begin;

        enum { c1 = V1::vconj };

        do {
          X0 = *X++;
          rA0 = Maybe<add>::sum( real(*A0) , ZProd<c1,false>::rprod(X0,Y0) );
          iA0 = Maybe<add>::sum( imag(*A0) , ZProd<c1,false>::iprod(X0,Y0) );
          rA1 = Maybe<add>::sum( real(*A1) , ZProd<c1,false>::rprod(X0,Y1) );
          iA1 = Maybe<add>::sum( imag(*A1) , ZProd<c1,false>::iprod(X0,Y1) );
          *A0++ = T3(rA0,iA0);
          *A1++ = T3(rA1,iA1);
          rA0 = Maybe<add>::sum( real(*A2) , ZProd<c1,false>::rprod(X0,Y2) );
          iA0 = Maybe<add>::sum( imag(*A2) , ZProd<c1,false>::iprod(X0,Y2) );
          rA1 = Maybe<add>::sum( real(*A3) , ZProd<c1,false>::rprod(X0,Y3) );
          iA1 = Maybe<add>::sum( imag(*A3) , ZProd<c1,false>::iprod(X0,Y3) );
          *A2++ = T3(rA0,iA0);
          *A3++ = T3(rA1,iA1);
        } while (--M);
      }
    }
  };

  // algo 21: fully unroll by columns
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<21,cs,rs,add,ix,T,V1,V2,M3>
  {
    template <int J, int N>
    struct Unroller
    {
      static inline void unroll(
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
        static inline void unroll(const V1& v1, const PT2& v2j, M3& m3)
        { 
          Unroller2<I,M/2>::unroll(v1,v2j,m3);
          Unroller2<I+M/2,M-M/2>::unroll(v1,v2j,m3);
        }
      };
      template <int I>
      struct Unroller2<I,1>
      {
        static inline void unroll(const V1& v1, const PT2& v2j, M3& m3)
        {
          Maybe<add>::add(m3.ref(I,J) , 
              ZProd<false,false>::prod(v1.cref(I) , v2j ) ); 
        }
      };
      static inline void unroll(
          const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
      { Unroller2<0,cs>::unroll(v1,x*v2.cref(J),m3); }
    };
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
#ifdef PRINTALGO
      std::cout<<"algo 21: cs,rs,x = "<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      Unroller<0,rs>::unroll(x,v1,v2,m3); 
    }
  };

  // algo 22: fully unroll by rows
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<22,cs,rs,add,ix,T,V1,V2,M3>
  {
    template <int I, int M>
    struct Unroller
    {
      static inline void unroll(
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
        static inline void unroll(const PT1& v1i, const V2& v2, M3& m3)
        { 
          Unroller2<J,N/2>::unroll(v1i,v2,m3);
          Unroller2<J+N/2,N-N/2>::unroll(v1i,v2,m3);
        }
      };
      template <int J>
      struct Unroller2<J,1>
      {
        static inline void unroll(const PT1& v1i, const V2& v2, M3& m3)
        {
          Maybe<add>::add(m3.ref(I,J) ,
              ZProd<false,false>::prod(v1i , v2.cref(J) ) ); 
        }
      };
      static inline void unroll(
          const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
      { Unroller2<0,rs>::unroll(x*v1.cref(I),v2,m3); }
    };
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
#ifdef PRINTALGO
      std::cout<<"algo 22: cs,rs,x = "<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      Unroller<0,cs>::unroll(x,v1,v2,m3); 
    }
  };

  // algo 31: colmajor, ix==0, so might need to copy v1
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<31,cs,rs,add,ix,T,V1,V2,M3> // cs known
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX);
      const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
      typedef typename V1::value_type T1;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 31: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = (M3::miscomplex ? 16 : 12 ) };
#ifdef TMV_OPT_SCALE
      if ( N > TMV_Q4 * M) 
      {
#endif
#ifdef PRINTALGO
        std::cout<<"v1c = x * v1; m3 (+=) v1c ^ v2\n";
#endif
        typedef typename Traits2<T,T1>::type T1c;
        typedef SmallVector<T1c,cs> V1c;
        V1c v1c;
        typedef typename V1c::view_type V1cv;
        typedef typename V1c::const_view_type V1ccv;
        V1cv v1cv = v1c.View();
        V1ccv v1ccv = v1c.View();
        AddVV_Helper<-1,cs,false,ix,T,V1,V1cv>::call(x,v1,v1cv);
        Rank1Update_VVM_Helper<algo2,cs,rs,add,1,RT,V1ccv,V2,M3>::call(
            one,v1ccv,v2,m3);
#ifdef TMV_OPT_SCALE
      }
      else Rank1Update_VVM_Helper<algo2,cs,rs,add,ix,T,V1,V2,M3>::call(
            x,v1,v2,m3);
#endif
    }
  };
  template <int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<31,UNKNOWN,rs,add,ix,T,V1,V2,M3> // cs unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      const int M = m3.colsize();
      const int cs = UNKNOWN;
      const int N = (rs == UNKNOWN ? int(m3.rowsize()) : rs);
      typedef typename V1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 31: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = (M3::miscomplex ? 16 : 12 ) };
#ifdef TMV_OPT_SCALE
      if ( N > TMV_Q4 * M) 
      {
#endif
#ifdef PRINTALGO
        std::cout<<"v1c = x * v1; m3 (+=) v1c ^ v2\n";
#endif
        typedef typename Traits2<T,T1>::type T1c;
        typedef Vector<T1c> V1c;
        V1c v1c(M);
        typedef typename V1c::view_type V1cv;
        typedef typename V1c::const_view_type V1ccv;
        V1cv v1cv = v1c.View();
        V1ccv v1ccv = v1c.View();
        AddVV_Helper<-1,cs,false,ix,T,V1,V1cv>::call(x,v1,v1cv);
        Rank1Update_VVM_Helper<algo2,cs,rs,add,1,RT,V1ccv,V2,M3>::call(
            one,v1ccv,v2,m3);
#ifdef TMV_OPT_SCALE
      }
      else Rank1Update_VVM_Helper<algo2,cs,rs,add,ix,T,V1,V2,M3>::call(
            x,v1,v2,m3);
#endif
    }
  };

  // algo 32: colmajor, v1.step() != 1, so copy v1
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<32,cs,rs,add,ix,T,V1,V2,M3> // cs known
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
      typedef typename V1::value_type T1;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 32: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = ( M3::miscomplex ? 16 : 12 ) };
      if ( N >= M) 
      {
#ifdef PRINTALGO
        std::cout<<"v1c = x * v1; m3 (+=) v1c ^ v2\n";
#endif
        typedef typename Traits2<T,T1>::type T1c;
        typedef SmallVector<T1c,cs> V1c;
        V1c v1c;
        typedef typename V1c::view_type V1cv;
        typedef typename V1c::const_view_type V1ccv;
        V1cv v1cv = v1c.View();
        V1ccv v1ccv = v1c.View();
        AddVV_Helper<-1,cs,false,ix,T,V1,V1cv>::call(x,v1,v1cv);
        Rank1Update_VVM_Helper<algo2,cs,rs,add,1,RT,V1ccv,V2,M3>::call(
            one,v1ccv,v2,m3);
      }
      else 
      {
#ifdef PRINTALGO
        std::cout<<"v1c = v1; m3 (+=) x * v1c ^ v2\n";
#endif
        typedef SmallVector<T1,cs> V1c;
        V1c v1c;
        typedef typename V1c::view_type V1cv;
        typedef typename V1c::const_view_type V1ccv;
        V1cv v1cv = v1c.View();
        V1ccv v1ccv = v1c.View();
        CopyV_Helper<-1,cs,V1,V1cv>::call(v1,v1cv);
        Rank1Update_VVM_Helper<algo2,cs,rs,add,ix,T,V1ccv,V2,M3>::call(
            x,v1ccv,v2,m3);
      }
    }
  };
  template <int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<32,UNKNOWN,rs,add,ix,T,V1,V2,M3> // cs unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      const int M = m3.colsize();
      const int cs = UNKNOWN;
      const int N = (rs == UNKNOWN ? int(m3.rowsize()) : rs);
      typedef typename V1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 32: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      enum { algo2 = (M3::miscomplex ? 16 : 12 ) };
      if ( N >= M) 
      {
#ifdef PRINTALGO
        std::cout<<"v1c = x * v1; m3 (+=) v1c ^ v2\n";
#endif
        typedef typename Traits2<T,T1>::type T1c;
        typedef Vector<T1c> V1c;
        V1c v1c(M);
        typedef typename V1c::view_type V1cv;
        typedef typename V1c::const_view_type V1ccv;
        V1cv v1cv = v1c.View();
        V1ccv v1ccv = v1c.View();
        AddVV_Helper<-1,cs,false,ix,T,V1,V1cv>::call(x,v1,v1cv);
        Rank1Update_VVM_Helper<algo2,cs,rs,add,1,RT,V1ccv,V2,M3>::call(
            one,v1ccv,v2,m3);
      }
      else 
      {
#ifdef PRINTALGO
        std::cout<<"v1c = v1; m3 (+=) x * v1c ^ v2\n";
#endif
        typedef Vector<T1> V1c;
        V1c v1c(M);
        typedef typename V1c::view_type V1cv;
        typedef typename V1c::const_view_type V1ccv;
        V1cv v1cv = v1c.View();
        V1ccv v1ccv = v1c.View();
        CopyV_Helper<-1,cs,V1,V1cv>::call(v1,v1cv);
        Rank1Update_VVM_Helper<algo2,cs,rs,add,ix,T,V1ccv,V2,M3>::call(
            x,v1ccv,v2,m3);
      }
    }
  };

  // algo 33: In some cases for small matrices with known cs, rs, 
  // we should loop over rows, rather than over columns, 
  // even though m3 is colmajor.  
  // This algorithm switches to the transpose.
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<33,cs,rs,add,ix,T,V1,V2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
#ifdef PRINTALGO
      const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
      std::cout<<"algo 33: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
      typedef typename M3::transpose_type M3t;
      enum { algo2 = (
            ( cs <= 4 ) ? (M3::miscomplex ? 17 : 13) :
            M3::miscomplex ? 16 : 12 ) };
      M3t m3t = m3.Transpose();
      Rank1Update_VVM_Helper<algo2,rs,cs,add,ix,T,V2,V1,M3t>::call(x,v2,v1,m3t);
    }
  };

  // algo -3: The same as -1, but allow for the possibility of using
  // InstMultMV.  This is used when calling MultMV from a degenerate MultMM.
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<-3,cs,rs,add,ix,T,V1,V2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename V1::value_type T1;
      typedef typename V2::value_type T2;
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
          cs == UNKNOWN && rs == UNKNOWN ) };

      enum { algo = inst ? 3 : -1 };
      Rank1Update_VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
    }
  };

  // algo -2: Determine which algorithm to use, but no if statements, 
  //          and no unrolling.
  // Used by MultMM to avoid if clauses that we already know don't apply.
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<-2,cs,rs,add,ix,T,V1,V2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      enum { algo = (
#if TMV_OPT >= 1
          ( rs == 0 || cs == 0 ) ? 0 : // trivial - nothing to do
          ( cs == 1 ) ? 1 : // trivial - cs = 1
          ( rs == 1 ) ? 2 : // trivial - rs = 1
          M3::mcolmajor ? ( // colmajor
            (cs == UNKNOWN || rs == UNKNOWN ) ? (
                M3::miscomplex ? 16 : 12 ) :
            ( rs > cs && cs <= 4 && V2::vstep == 1 ) ? 33 :
            ( rs <= 4 ) ? (M3::miscomplex ? 17 : 13) :
            M3::miscomplex ? 16 : 12 ) :
#endif
          // nomajor -- don't do anything fancy
          11 ) };
      Rank1Update_VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, bool add, int ix, class T, class V1, class V2, class M3>
  struct Rank1Update_VVM_Helper<-1,cs,rs,add,ix,T,V1,V2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      // Possible algorithms to choose from:
      //
      //  0 = cs or rs == 0, so nothing to do
      //  1 = cs == 1: reduces to trivial AddVV function
      //  2 = rs == 1: reduces to trivial AddVV function
      //  3 = call InstRank1Update
      //
      // 11 = simple for loop
      // 12 = 4 columns at a time
      // 13 = rs is known and <= 4
      // 16 = complex m3, 2 columns at a time
      // 17 = complex m3, rs is known and <= 4
      //
      // These are currently not used.  They aren't faster for and cs,rs.
      // 21 = fully unroll by cols
      // 22 = fully unroll by rows
      //
      // 31 = x != 1, copy v1 if M < N
      // 32 = v1.step() != 1, so copy v1
      // 33 = switch to transposed form
      //

      enum { algo = (
#if TMV_OPT >= 1
          ( rs == 0 || cs == 0 ) ? 0 : // trivial - nothing to do
          ( cs == 1 ) ? 1 : // trivial - cs = 1
          ( rs == 1 ) ? 2 : // trivial - rs = 1
          M3::mcolmajor ? ( // colmajor
            (cs == UNKNOWN || rs == UNKNOWN ) ? (
              ( V1::vstep != 1 ? 32 : 
                TMV_ZeroIX ? 31 : 
                M3::miscomplex ? 16 : 12 ) ) :
            ( rs > cs && cs <= 4 && V2::vstep == 1 ) ? 33 :
            ( cs > TMV_Q2 && V1::vstep != 1 ) ? 32 :
            ( TMV_ZeroIX && rs > IntTraits2<TMV_Q4,cs>::prod ) ? 31 :
            ( rs <= 4 ) ? (M3::miscomplex ? 17 : 13) :
            M3::miscomplex ? 16 : 12 ) :
#endif
          // nomajor -- don't do anything fancy
          11 ) };
#ifdef PRINTALGO
      std::cout<<"InlineRank1Update_VVM: \n";
      std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
      std::cout<<"v1 = "<<TypeText(v1)<<std::endl;
      std::cout<<"v2 = "<<TypeText(v2)<<std::endl;
      std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
      std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
#endif
      Rank1Update_VVM_Helper<algo,cs,rs,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
    }
  };

  template <bool add, int ix, class T, class V1, class V2, class M3>
  inline void InlineRank1Update(const Scaling<ix,T>& x, 
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2, 
      BaseMatrix_Rec_Mutable<M3>& m3)
  {
    TMVStaticAssert(M3::mcolmajor || !M3::mrowmajor); 
    // This should have been handled by CallRank1Update_Helper
    // (And I have to write it this way in case rs == 1, in
    // which case m3 is _both_ colmajor and rowmajor.)

    //std::cout<<"InlineRank1Update "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"v1 = "<<TypeText(v1)<<"  "<<v1<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3<<std::endl;
    TMVStaticAssert((Sizes<M3::mcolsize,V1::vsize>::same));
    TMVStaticAssert((Sizes<M3::mrowsize,V2::vsize>::same));
    TMVAssert(m3.colsize() == v1.size());
    TMVAssert(m3.rowsize() == v2.size());
    enum { cs = Sizes<M3::mcolsize,V1::vsize>::size };
    enum { rs = Sizes<M3::mrowsize,V2::vsize>::size };
    typedef typename V1::const_view_type V1v;
    typedef typename V2::const_view_type V2v;
    typedef typename M3::view_type M3v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    M3v m3v = m3.View();
    Rank1Update_VVM_Helper<-1,cs,rs,add,ix,T,V1v,V2v,M3v>::call(x,v1v,v2v,m3v);
    //std::cout<<"End InlineRank1Update_VVM: m3 = "<<m3<<std::endl;
  }

  // Defined in TMV_Rank1_VVM.cpp
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstRank1Update(const T3 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1, 
      const ConstVectorView<T2,UNKNOWN,C2>& v2, MatrixView<T3> m3);
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddRank1Update(const T3 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1, 
      const ConstVectorView<T2,UNKNOWN,C2>& v2, MatrixView<T3> m3);

} // namespace tmv

#undef TMV_OPT_SCALE

#undef TMV_Q1
#undef TMV_Q2
#undef TMV_Q3
#undef TMV_Q4
#undef TMV_ZeroIX

#endif 
