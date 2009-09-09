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


#ifndef TMV_MultLV_H
#define TMV_MultLV_H

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

//#define PRINTALGO

#ifdef PRINTALGO
#include <iostream>
#endif

  // There are a number of values used in the algorithm selection
  // that are either arbitrary or empirical.
  // So I put them all here to make them easier to change and to
  // track down in the code.

  // Q1 is the maximum value of s*s for which a SmallMatrix * SmallVector
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

  // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
#define TMV_ZeroIX (ix==0)
//#define TMV_ZeroIX (ix!=1)

  template <int algo, int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper;

  // algo 0: s = 0, so nothing to do
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<0,s,add,ix,T,M1,V2,V3>
  { static void call(const Scaling<ix,T>& , const M1& , const V2& , V3& ) {} };

  // algo 1: s == 1, so simplifies to a scalar product
  template <bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<1,1,rs,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      Maybe<add>::add( v3.ref(0) , 
        x * Maybe<!M1::munit>::prod(m1.cref(0,0) , v2.cref(0)) ); 
    }
  };

  // algo 11: The basic column major loop
  // This is designed to work if v2 and v3 are the same storage.
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<11,s,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = (s == UNKNOWN ? m1.size() : s);
      if (N == 0) return;
#ifdef PRINTALGO
      std::cout<<"algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef typename V2::value_type T2;
      typedef typename V3::value_type T3;
      typedef typename Traits2<T,typename V2::value_type>::type PT2;
      typedef typename M1::const_col_range_type M1c;
      typedef typename M1::const_diag_type M1d;
      PT2 Xj;

      enum { unit = M1::munit };
      enum { c1 = M1::mconj };
      enum { c2 = V2::vconj };

      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      typedef typename M1d::const_nonconj_type::const_iterator IT1d;
      typedef typename V2::const_nonconj_type::const_iterator IT2;
      typedef typename V3::iterator IT3;
      IT1 Ajj = A.get_col(N-1,N-1,N).NonConj().begin();
      IT2 X = v2.NonConj().begin() + N-1;
      IT3 Y = v3.begin() + N-1;
      const int Adiagstep = A.diagstep();

      Maybe<add>::add(*Y-- , 
          Maybe<!unit>::template zprod<c1,false>(*Ajj ,
            ZProd<false,c2>::prod(x , *X--)));  
      Ajj.ShiftP(-Adiagstep);

      for(int jj=N-1,len=1; jj; --jj,++len) {
        // (jj = j+1)
        if (*X != T2(0)) {
          Xj = ZProd<false,c2>::prod(x * *X--);
          // y.SubVector(0,j) += x(j) * A.col(j,0,j);
          AddVV_Helper<-1,UNKNOWN,true,0,PT2,M1c,V3>::call2(
              len,Scaling<0,PT2>(Xj),Ajj,Y);
          Maybe<add>::add(*Y-- , Maybe<!unit>::zprod<c1,false>(*Ajj , Xj));  
        } else {
          if (!add) *Y = T3(0);
          --X; --Y;
        }
        Ajj.ShiftP(-Adiagstep);
      }
    }
  };

  // algo 21: The basic row major loop
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<21,s,add,ix,T,M1,V2,V3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { 
      const int N = (s == UNKNOWN ? m1.size() : s);
      if (N == 0) return;
#ifdef PRINTALGO
      std::cout<<"algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef typename Traits2<T1,T2>::type PT;
      typedef typename M1::const_row_range_type M1r;
      PT Yi;

      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      typedef typename V2::const_nonconj_type::const_iterator IT2;
      typedef typename V3::iterator IT3;

      if (N == 0) return;
      enum { unit = M1::munit };

      IT1 Ai0 = m1.get_row(0,0,N).NonConj().begin();
      IT2 X0 = v2.NonConj().begin();
      IT2 X = X0 + N-1;
      IT3 Y = v3.begin() + N-1;
      const int Astepi = m1.stepi();
      int len = Maybe<unit>::sum(-1 , N);

      if (len) do {
        // Yi = A.row(i,0,i+1) * x.SubVector(0,i+1)
        Yi = Maybe<unit>::sum( *X-- , 
            MultVV_Helper<-1,UNKNOWN,M1r,V2>::call2(len,Ai0,X0));
        Maybe<add>::add(*Y--, x * Yi);
        Ai0.ShiftP(-Astepj);
      } while (--len);

      if (unit)
        Maybe<add>::add(*Y , ZProd<false,c2>::prod(x , *X));
    }
  };

  // algo 31: fully unroll by rows, apply x to v3
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<31,s,add,ix,T,M1,V2,V3>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
      {
        Unroller<I+M/2,M-M/2,J,N>::unroll(x,m1,v2,v3);
        Unroller<I,M/2,J,M/2>::unroll(x,m1,v2,v3);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline void unroll(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
      {
        Unroller<I,1,J,N/2>::unroll(x,m1,v2,v3);
        Unroller<I,1,J+N/2,N-N/2>::unroll(x,m1,v2,v3);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    {
      static inline void unroll(
          const Scaling<ix,T>& , const M1& , const V2& , V3& ) {}
    }
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
      { 
        Maybe<add>::add(v3.ref(I) , ZProd<false,false>(x ,
            ZProd<false,false>::prod(m.cref(I,J) , v2.cref(J) )));
      }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    {
      static inline void unroll(
        const Scaling<ix,T>& , const M1& , const V2& , V3& )
    }

    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { Unroller<0,s,0,s>::unroll(x,m1,v2,v3); }
  };

  // algo 43: colmajor, v3.step() != 1, so copy v3
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<43,s,add,ix,T,M1,V2,V3> // s known
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = (rs == UNKNOWN ? m1.rowsize() : rs);
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 43: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

      typedef typename Traits2<T1,T2>::type T3c;
      typedef SmallVector<T3c,s> V3c;
      V3c v3c;
      typedef typename V3c::view_type V3cv;
      typedef typename V3c::const_view_type V3ccv;
      V3cv v3cv = v3c.View();
      V3ccv v3ccv = v3c.View();
      MultLV_Helper<11,s,false,1,RT,M1,V2,V3cv>::call(one,m1,v2,v3cv);
      AddVV_Helper<-1,s,add,ix,T,V3ccv,V3>::call(x,v3ccv,v3);
    }
  };
  template <int rs, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<43,UNKNOWN,rs,add,ix,T,M1,V2,V3> // s unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int s = UNKNOWN;
      const int N = (s == UNKNOWN ? m1.size() : s);
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 43: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

      typedef typename Traits2<T1,T2>::type T3c;
      typedef Vector<T3c> V3c;
      V3c v3c(M);
      typedef typename V3c::view_type V3cv;
      typedef typename V3c::const_view_type V3ccv;
      V3cv v3cv = v3c.View();
      V3ccv v3ccv = v3c.View();
      MultLV_Helper<11,s,false,1,RT,M1,V2,V3cv>::call(one,m1,v2,v3cv);
      AddVV_Helper<-1,s,add,ix,T,V3ccv,V3>::call(x,v3ccv,v3);
    }
  };

  // algo 53: rowmajor, v2.step() != 1, so copy v2
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<53,s,add,ix,T,M1,V2,V3> // rs = known
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = (s == UNKNOWN ? m1.size() : s);
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
#ifdef PRINTALGO
      std::cout<<"algo 53: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
      typedef typename Traits2<T,T2>::type T2c;
      typedef SmallVector<T2c,s> V2c;
      V2c v2c;
      typedef typename V2c::view_type V2cv;
      typedef typename V2c::const_view_type V2ccv;
      V2cv v2cv = v2c.View();
      V2ccv v2ccv = v2c.View();
      AddVV_Helper<-1,s,false,ix,T,V2,V2cv>::call(x,v2,v2cv);
      MultLV_Helper<51,s,add,1,RT,M1,V2ccv,V3>::call(one,m1,v2ccv,v3);
    }
  };
  template <bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<53,UNKNOWN,add,ix,T,M1,V2,V3> // s = unknown
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      const int N = m1.size();
      const int s = UNKNOWN;
      typedef typename M1::value_type T1;
      typedef typename V2::value_type T2;
#ifdef PRINTALGO
      std::cout<<"algo 53: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
      typedef typename Traits2<T,T2>::type T2c;
      typedef Vector<T2c> V2c;
      V2c v2c(N);
      typedef typename V2c::view_type V2cv;
      typedef typename V2c::const_view_type V2ccv;
      V2cv v2cv = v2c.View();
      V2ccv v2ccv = v2c.View();
      AddVV_Helper<-1,s,false,ix,T,V2,V2cv>::call(x,v2,v2cv);
      MultLV_Helper<51,s,add,1,RT,M1,V2ccv,V3>::call(one,m1,v2ccv,v3);
    }
  };

  // algo -2: Determine which algorithm to use, but no if statements, 
  //          and no unrolling.
  // Used by MultMM to avoid if clauses that we already know don't apply.
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<-2,s,add,ix,T,M1,V2,V3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename M1::value_type T1;
      enum { algo = (
          ( s == 0 ) ? 0 : // trivial - nothing to do
          ( s == 1 ) ? 1 : // trivial - s = 1
          M1::mcolmajor ? 11 : 21 ) };
      MultLV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int s, bool add, int ix, class T, class M1, class V2, class V3>
  struct MultLV_Helper<-1,s,add,ix,T,M1,V2,V3>
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
      //  0 = s == 0, so nothing to do
      //  1 = s == 1: reduces to trivial scalar product
      //
      // Column Major:
      // 11 = column major, simple for loop
      //
      // Row Major:
      // 21 = row major, simple for loop
      //
      // Fully Unrolled:
      // 31 = fully unroll by rows
      //
      // Column Major, meta algorithms
      // 43 = column major, v3.step() != 1, so copy v3
      //
      // Row Major, meta algorithms
      // 53 = row major, v2.step() != 1, so copy v2

#if TMV_OPT == 0
      enum { algo = M1::mcolmajor ? 11 : 21 };
#else
      enum { unroll = (
          s == UNKNOWN ? false :
          IntTraits2<s,s>::prod > TMV_Q1 ? false :
          s <= 8 ) };
      enum { algo = (
          ( s == 0 ) ? 0 : // trivial - nothing to do
          ( s == 1 ) ? 1 : // trivial - s = 1
          unroll ? 31 :
          M1::mcolmajor ? ( V3::vstep != 1 ? 43 : 11 ) :
          M1::mrowmajor ? ( V2::vstep != 1 ? 53 : 21 ) :
          V2::vstep == 1 ? 21 : V3::vstep == 1 ? 11 : 21 ) };
#endif
#ifdef PRINTALGO
      std::cout<<"InlineMultLV: \n";
      std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      std::cout<<"v2 = "<<TypeText(v2)<<std::endl;
      std::cout<<"v3 = "<<TypeText(v3)<<std::endl;
      std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
      MultLV_Helper<algo,cs,rs,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
    }
  };

  template <bool add, int ix, class T, class M1, class V2, class V3>
  inline void InlineMultLV(const Scaling<ix,T>& x, 
      const BaseMatrix_Rec<M1>& m1, const BaseVector_Calc<V2>& v2, 
      BaseVector_Mutable<V3>& v3)
  {
    //std::cout<<"InlineMultLV "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"add = "<<add<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    //std::cout<<"v3 = "<<TypeText(v3)<<"  "<<v3<<std::endl;
    TMVStaticAssert((Sizes<M1::msize,V3::vsize>::same));
    TMVStaticAssert((Sizes<M1::msize,V2::vsize>::same));
    TMVStaticAssert((Sizes<V2::vsize,V3::vsize>::same));
    TMVAssert(m1.size() == v3.size());
    TMVAssert(m1.size() == v2.size());
    TMVAssert(v2.size() == v3.size());
    enum { s = Sizes<Sizes<M1::msize,V2::vsize>::size,V3::vsize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename V2::const_view_type V2v;
    typedef typename V3::view_type V3v;
    M1v m1v = m1.View();
    V2v v2v = v2.View();
    V3v v3v = v3.View();
    MultLV_Helper<-1,s,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    //std::cout<<"End InlineMultLV: v3 = "<<v3<<std::endl;
  }

  // Defined in TMV_MultLV.cpp
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstMultLV(const T3 x,
      const ConstLowerTriMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddMultLV(const T3 x,
      const ConstLowerTriMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

} // namespace tmv

#undef TMV_Q1
#undef TMV_Q2
#undef TMV_Q3
#undef TMV_Q4
#undef TMV_Q5
#undef TMV_ZeroIX

#endif 
