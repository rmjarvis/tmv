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


#ifndef TMV_MultVV_H
#define TMV_MultVV_H

#include "TMV_Scaling.h"
#include "TMV_BaseVector.h"

namespace tmv {

  //
  // Vector * Vector
  //

  const int TMV_MultVV_RecurseSize = 16*1024;


  // v1*v2
  template <int algo, int size, class V1, class V2> struct MultVV_Helper;

  // algo 1: simple for loop
  template <int size, class V1, class V2>
  struct MultVV_Helper<1,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      if (n == 0) return PT(0);
      else {
        PT sum(0);
        for(int i=0;i<n;++i) 
          sum += ZProd<false,false>::prod(v1.cref(i) , v2.cref(i));
        return sum;
      }
    }
    static PT call2(int n, IT1 it1, IT2 it2)
    {
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      PT sum(0);
      if (n) do { 
        sum += ZProd<c1,c2>::prod(*it1++ , *it2++); 
      } while (--n);
      return sum;
    }
  };

  // algo 2: 2 at once
  template <int size, class V1, class V2>
  struct MultVV_Helper<2,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      return call2(n,v1.NonConj().begin(),v2.NonConj().begin());
    }
    static PT call2(const int n, IT1 it1, IT2 it2)
    {
      PT sum0(0), sum1(0);
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      if (n_2) {
        do {
          sum0 += ZProd<c1,c2>::prod(it1[0] , it2[0]);
          sum1 += ZProd<c1,c2>::prod(it1[1] , it2[1]);
          it1 += 2;
          it2 += 2;
        } while (--n_2);
        sum0 += sum1;
      }
      if (nb) {
        sum0 += ZProd<c1,c2>::prod((*it1) , (*it2));
      }
      return sum0;
    }
  };

  // algo 3: 4 at once
  template <int size, class V1, class V2>
  struct MultVV_Helper<3,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      return call2(n,v1.NonConj().begin(),v2.NonConj().begin());
    }
    static PT call2(const int n, IT1 it1, IT2 it2)
    {
      PT sum0(0), sum1(0);
      int n_4 = (n>>2);
      int nb = n-(n_4<<2);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      if (n_4) {
        do {
          sum0 += ZProd<c1,c2>::prod(it1[0] , it2[0]);
          sum1 += ZProd<c1,c2>::prod(it1[1] , it2[1]);
          sum0 += ZProd<c1,c2>::prod(it1[2] , it2[2]);
          sum1 += ZProd<c1,c2>::prod(it1[3] , it2[3]);
          it1 += 4;
          it2 += 4;
        } while (--n_4);
        sum0 += sum1;
      }
      if (nb) do {
        sum0 += ZProd<c1,c2>::prod((*it1++) , (*it2++));
      } while (--nb);
      return sum0;
    }
  };

  // algo 4: 8 at once
  template <int size, class V1, class V2>
  struct MultVV_Helper<4,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      return call2(n,v1.NonConj().begin(),v2.NonConj().begin());
    }
    static PT call2(const int n, IT1 it1, IT2 it2)
    {
      PT sum0(0), sum1(0);
      int n_8 = (n>>3);
      int nb = n-(n_8<<3);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      if (n_8) {
        do {
          sum0 += ZProd<c1,c2>::prod(it1[0] , it2[0]);
          sum1 += ZProd<c1,c2>::prod(it1[1] , it2[1]);
          sum0 += ZProd<c1,c2>::prod(it1[2] , it2[2]);
          sum1 += ZProd<c1,c2>::prod(it1[3] , it2[3]);
          sum0 += ZProd<c1,c2>::prod(it1[4] , it2[4]);
          sum1 += ZProd<c1,c2>::prod(it1[5] , it2[5]);
          sum0 += ZProd<c1,c2>::prod(it1[6] , it2[6]);
          sum1 += ZProd<c1,c2>::prod(it1[7] , it2[7]);
          it1 += 8; it2 += 8;
        } while (--n_8);
        sum0 += sum1;
      }
      if (nb) do {
        sum0 += ZProd<c1,c2>::prod((*it1++) , (*it2++));
      } while (--nb);
      return sum0;
    }
  };

  // algo 5: fully unroll
  template <int size, class V1, class V2>
  struct MultVV_Helper<5,size,V1,V2>
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    template <int I, int N>
    struct Unroller
    {
      static inline PT unroll(const V1& v1, const V2& v2)
      {
        return (
            Unroller<I,N/2>::unroll(v1,v2) +
            Unroller<I+N/2,N-N/2>::unroll(v1,v2));
      }
      static inline PT unroll2(const IT1& it1, const IT2& it2)
      {
        return (
            Unroller<I,N/2>::unroll2(it1,it2) +
            Unroller<I+N/2,N-N/2>::unroll2(it1,it2));
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline PT unroll(const V1& v1, const V2& v2)
      { return ZProd<false,false>::prod(v1.cref(I) , v2.cref(I)); }
      static inline PT unroll2(const IT1& it1, const IT2& it2)
      {
        enum { c1 = V1::vconj };
        enum { c2 = V2::vconj };
        return ZProd<c1,c2>::prod(it1[I] , it2[I]); 
      }
    };
    template <int I>
    struct Unroller<I,0>
    {
      static inline PT unroll(const V1& v1, const V2& v2)
      { return PT(0); }
      static inline PT unroll2(const IT1& it1, const IT2& it2)
      { return PT(0); }
    };

    static inline PT call(const V1& v1, const V2& v2)
    { return Unroller<0,size>::unroll(v1,v2); }
    static inline PT call2(const int n, const IT1& it1, const IT2& it2)
    { return Unroller<0,size>::unroll2(it1,it2); }
  };

  // algo 7: complex vectors
  template <int size, class V1, class V2>
  struct MultVV_Helper<7,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename Traits<PT>::real_type RT;
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      return call2(n,v1.NonConj().begin(),v2.NonConj().begin());
    }
    static PT call2(int n, IT1 it1, IT2 it2)
    {
      RT rsum(0), isum(0);
      T1 val1;
      T2 val2;
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n) do {
        val1 = *it1++;
        val2 = *it2++;
        rsum += ZProd<c1,c2>::rprod(val1,val2);
        isum += ZProd<c1,c2>::iprod(val1,val2);
      } while (--n);
      return PT(rsum,isum);
    }
  };

  // algo 8: both complex, both unit, non-conj, 2 at once
  template <int size, class V1, class V2>
  struct MultVV_Helper<8,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename Traits<PT>::real_type RT;
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      return call2(n,v1.NonConj().begin(),v2.NonConj().begin());
    }
    static PT call2(const int n, IT1 it1, IT2 it2)
    {
      RT rsum(0), isum(0);
      T1 val1a, val1b;
      T2 val2a, val2b;
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      if (n_2) do {
        val1a = it1[0]; val2a = it2[0];
        val1b = it1[1]; val2b = it2[1];
        it1 += 2; it2 += 2;
        rsum += (
            ZProd<c1,c2>::rprod(val1a,val2a) +
            ZProd<c1,c2>::rprod(val1b,val2b));
        isum += (
            ZProd<c1,c2>::iprod(val1a,val2a) +
            ZProd<c1,c2>::iprod(val1b,val2b));
      } while (--n_2);
      if (nb) {
        val1a = it1[0]; val2a = it2[0];
        rsum += ZProd<c1,c2>::rprod(val1a,val2a);
        isum += ZProd<c1,c2>::iprod(val1a,val2a);
      }
      return PT(rsum,isum);
    }
  };

  // algo 10: recurse very large vector product
  // This isn't for speed reasons - it's for increased accuracy.
  // For large vectors, the incremental additions can be much smaller
  // than the running sum, so the relative errors can be huge.
  // With the recursive algorithm, the relative error is generally
  // closer to the expected few * epsilon.
  template <int size, class V1, class V2>
  struct MultVV_Helper<10,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V1::const_subvector_type V1s;
    typedef typename V2::const_subvector_type V2s;
    typedef typename Traits<PT>::real_type RT;
    enum { size1 = (size == UNKNOWN ? UNKNOWN : size/2) };
    enum { size2 = (size == UNKNOWN ? UNKNOWN : size-size/2) };
    enum { allunit = (V1::vstep == 1 && V2::vstep == 1) };
    enum { allcomplex = (Traits<T1>::iscomplex && Traits<T2>::iscomplex) };
    enum { algo2 = (
        (sizeof(RT) == 8 && allunit) ? ( 
          ( allcomplex ? 7 : 2 ) ) :
        (sizeof(RT) == 4 && allunit) ? (
          ( allcomplex ? 7 : 3 ) ) :
        (allcomplex && allunit) ? 7 : 
        1 ) };
    static PT call(const V1& v1, const V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      return call2(n,v1.NonConj().begin(),v2.NonConj().begin());
    }
    static PT call2(const int n, const IT1& it1, const IT2& it2)
    {
      if (n == 0) return PT(0);
      else if (n > TMV_MultVV_RecurseSize)
      {
        const int no2 = n/2;
        return (
            MultVV_Helper<10,size1,V1,V2>::call2(no2,it1,it2) + 
            MultVV_Helper<10,size2,V1,V2>::call2(n-no2,it1+no2,it2+no2));
      } 
      else return MultVV_Helper<algo2,size,V1,V2>::call2(n,it1,it2);
    }
  };

  // algo -2: Determine which algorithm to use, but no if statements.
  // (This just foregoes algo 10.)
  template <int size, class V1, class V2>
  struct MultVV_Helper<-2,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename Traits<PT>::real_type RT;
    enum { allunit = (V1::vstep == 1 && V2::vstep == 1) };
    enum { allcomplex = (Traits<T1>::iscomplex && Traits<T2>::iscomplex) };
    enum { algo = (
#if TMV_OPT >= 1
        (size != UNKNOWN && size <= (128/sizeof(PT))) ? 5 :
        (sizeof(RT) == 8 && allunit) ? ( 
          ( allcomplex ? 7 : 2 ) ) :
        (sizeof(RT) == 4 && allunit) ? (
          ( allcomplex ? 7 : 3 ) ) :
        allcomplex ? 7 : 
#endif
        1 ) };
    static PT call(const V1& v1, const V2& v2)
    { return MultVV_Helper<algo,size,V1,V2>::call(v1,v2); }
    static PT call2(const int n, const IT1& it1, const IT2& it2)
    { return MultVV_Helper<algo,size,V1,V2>::call2(n,it1,it2); }
  };

  // algo -1: Determine which algorithm to use.
  // Since MultVV gets called from other routines that might know the
  // size at compile time, but not be in v1 or v2, we move the algorithm
  // selection here, so they can just call algo == -1, rather than
  // duplicating the algo selection everywhere.
  template <int size, class V1, class V2>
  struct MultVV_Helper<-1,size,V1,V2> 
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename Traits2<T1,T2>::type PT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename Traits<PT>::real_type RT;
    enum { allunit = (V1::vstep == 1 && V2::vstep == 1) };
    enum { allcomplex = (Traits<T1>::iscomplex && Traits<T2>::iscomplex) };
    enum { algo = (
#if TMV_OPT >= 1
        (size != UNKNOWN && size <= (128/sizeof(PT))) ? 5 :
#if TMV_OPT >= 2
        (size == UNKNOWN || size > TMV_MultVV_RecurseSize) ? 10 : 
#endif
        (sizeof(RT) == 8 && allunit) ? ( 
          ( allcomplex ? 7 : 2 ) ) :
        (sizeof(RT) == 4 && allunit) ? (
          ( allcomplex ? 7 : 3 ) ) :
        allcomplex ? 7 : 
#endif
        1 ) };
    static PT call(const V1& v1, const V2& v2)
    { return MultVV_Helper<algo,size,V1,V2>::call(v1,v2); }
    static PT call2(const int n, const IT1& it1, const IT2& it2)
    { return MultVV_Helper<algo,size,V1,V2>::call2(n,it1,it2); }
  };

  template <class V1, class V2> 
  static typename Traits2<typename V1::value_type,typename V2::value_type>::type
  InlineMultVV(
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    enum { size = Sizes<V1::vsize,V2::vsize>::size };
    typedef typename V1::const_view_type V1v;
    typedef typename V2::const_view_type V2v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    return MultVV_Helper<-1,size,V1v,V2v>::call(v1v,v2v);
  }

  // Defined in TMV_MultVV.cpp
  template <class T1, bool C1, class T2>
  T2 InstMultVV(const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2>& v2); // BLAS dot, dotu, dotc

  template <class V1, class T2>
  static inline T2 InstMultVV(const V1& v1, const ConstVectorView<T2,UNKNOWN,true>& v2)
  { return TMV_CONJ(InstMultVV(v1.Conjugate(),v2.Conjugate())); }
  template <class T, bool C1>
  static inline std::complex<T> InstMultVV(
      const ConstVectorView<std::complex<T>,UNKNOWN,C1>& v1,
      const ConstVectorView<T> v2)
  { return InstMultVV(v2,v1); }


  //
  // ElementProd functions:
  //

  template <int algo, int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper;

  // algo 1: simple for loop
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<1,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(int n,
        const Scaling<ix1,T1>& x1, IT1 it1, IT2 it2, IT3 it3)
    {
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      if (n) do {
        Maybe<add>::add(*it3++ , 
            ZProd<false,false>::prod(x1 ,
              ZProd<c1,c2>::prod(*it1++ , *it2++)));
      } while (--n);
    }
  };  

  // algo 2: 2 at a time
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<2,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 x, IT2 y, IT3 z)
    {
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
#ifdef PRINTALGO_MD
      std::cout<<"ElemMultVV algo 2:\n";
      std::cout<<"c1,c2 = "<<c1<<','<<c2<<std::endl;
#endif

      if (n_2) do {
        Maybe<add>::add(*z , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(*x , *y)) );
        Maybe<add>::add(z[1] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[1] , y[1]) ));
        x += 2; y += 2; z += 2;
      } while (--n_2);
      if (nb) Maybe<add>::add(*z , 
          ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(*x , *y)) );
    }
  };

  // algo 3: 4 at a time
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<3,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 x, IT2 y, IT3 z)
    {
      int n_4 = (n>>2);
      int nb = n-(n_4<<2);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_4) do {
        Maybe<add>::add(*z , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(*x , *y)));
        Maybe<add>::add(z[1] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[1] , y[1])));
        Maybe<add>::add(z[2] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[2] , y[2])));
        Maybe<add>::add(z[3] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[3] , y[3])));
        x += 4; y += 4; z += 4;
      } while (--n_4);
      if (nb) do {
        Maybe<add>::add(*z++ , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(*x++ , *y++)));
      } while (--nb);
    }
  };

  // algo 4: 8 at a time
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<4,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 x, IT2 y, IT3 z)
    {
      int n_8 = (n>>3);
      int nb = n-(n_8<<3);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_8) do {
        Maybe<add>::add(*z , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(*x , *y)));
        Maybe<add>::add(z[1] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[1] , y[1])));
        Maybe<add>::add(z[2] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[2] , y[2])));
        Maybe<add>::add(z[3] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[3] , y[3])));
        Maybe<add>::add(z[4] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[4] , y[4])));
        Maybe<add>::add(z[5] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[5] , y[5])));
        Maybe<add>::add(z[6] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[6] , y[6])));
        Maybe<add>::add(z[7] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[7] , y[7])));
        x += 8; y += 8; z += 8;
      } while (--n_8);
      if (nb) do {
        Maybe<add>::add(*z++ , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(*x++ , *y++)));
      } while (--nb);
    }
  };

  // algo 5: fully unroll
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<5,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    enum { c1 = V1::vconj };
    enum { c2 = V2::vconj };
    template <int I, int N>
    struct Unroller
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
      {
        Unroller<I,N/2>::unroll(x1,v1,v2,v3);
        Unroller<I+N/2,N-N/2>::unroll(x1,v1,v2,v3);
      }
      static inline void unroll2(
          const Scaling<ix1,T1>& x1, const IT1& x, const IT2& y, const IT3& z)
      {
        Unroller<I,N/2>::unroll2(x1,x,y,z);
        Unroller<I+N/2,N-N/2>::unroll2(x1,x,y,z);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
      { 
        Maybe<add>::add( v3.ref(I) , 
            ZProd<false,false>::prod(x1 ,
              ZProd<false,false>(v1.cref(I) , v2.cref(I)) )); 
      }
      static inline void unroll2(
          const Scaling<ix1,T1>& x1, const IT1& x, const IT2& y, const IT3& z)
      {
        Maybe<add>::add( z[I] , 
            ZProd<false,false>::prod(x1 , ZProd<c1,c2>::prod(x[I],y[I]))); 
      }
    };
    template <int I>
    struct Unroller<I,0>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& , const V1& , const V2& , V3& ) {}
      static inline void unroll2(
          const Scaling<ix1,T1>& , const IT1& , const IT2& , const IT3& ) {}
    };
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    { Unroller<0,size>::unroll(x1,v1,v2,v3); }
    static inline void call2(const int ,
        const Scaling<ix1,T1>& x1, const IT1& x, const IT2& y, const IT3& z)
    { Unroller<0,size>::unroll2(x1,x,y,z); }
  };

  // algo 6: complex vector
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<6,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(int n,
        const Scaling<ix1,T1>& x1, IT1 x, IT2 y, IT3 z)
    {
      typedef typename V1::value_type VT1;
      typedef typename V2::value_type VT2;
      typedef typename V3::real_type RT3;
      typedef typename V3::value_type T3;

      VT1 xval;
      VT2 yval;
      RT3 rz,iz;

      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n) do {
        xval = *x++;
        yval = *y++;
        rz = ZProd<c1,c2>::rprod(xval,yval);
        iz = ZProd<c1,c2>::iprod(xval,yval);
        Maybe<add>::add(*z++ , x1 * T3(rz,iz));
      } while (--n);
    }
  };

  // algo 7: complex vector, 2 at a time
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<7,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 x, IT2 y, IT3 z)
    {
      typedef typename V1::value_type VT1;
      typedef typename V2::value_type VT2;
      typedef typename V3::real_type RT3;
      typedef typename V3::value_type T3;
      typedef typename V3::iterator IT3;

      VT1 xvala, xvalb;
      VT2 yvala, yvalb;
      RT3 rza,iza,rzb,izb;

      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_2) do {
        xvala = x[0]; yvala = y[0];
        rza = ZProd<c1,c2>::rprod(xvala,yvala);
        iza = ZProd<c1,c2>::iprod(xvala,yvala);
        Maybe<add>::add(z[0] , x1 * T3(rza,iza));
        xvalb = x[1]; yvalb = y[1]; x += 2; y += 2;
        rzb = ZProd<c1,c2>::rprod(xvalb,yvalb);
        izb = ZProd<c1,c2>::iprod(xvalb,yvalb);
        Maybe<add>::add(z[1] , x1 * T3(rzb,izb)); z += 2;
      } while (--n_2);
      if (nb) {
        xvala = *x;
        yvala = *y;
        rza = ZProd<c1,c2>::rprod(xvala,yvala);
        iza = ZProd<c1,c2>::iprod(xvala,yvala);
        Maybe<add>::add(*z , x1 * T3(rza,iza));
      }
    }
  };

  // algo 8: complex vectors and complex T
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<8,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x1,v1.NonConj().begin(),v2.NonConj().begin(),v3.begin());
    }
    static void call2(int n,
        const Scaling<ix1,T1>& x1, IT1 x, IT2 y, IT3 z)
    {
      typedef typename V1::value_type VT1;
      typedef typename V2::value_type VT2;
      typedef typename V3::real_type RT3;
      typedef typename V3::value_type T3;

      VT1 xval;
      VT2 yval;
      RT3 rz,iz,rz2,iz2;

      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n) do {
        xval = *x++;
        yval = *y++;
        rz = ZProd<c1,c2>::rprod(xval,yval);
        iz = ZProd<c1,c2>::iprod(xval,yval);
        rz2 = ZProd<false,false>::rprod(x1,T3(rz,iz));
        iz2 = ZProd<false,false>::iprod(x1,T3(rz,iz));
        Maybe<add>::add(*z++ , T3(rz2,iz2));
      } while (--n);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int size, bool add, int ix1, class T1, class V1, class V2, class V3>
  struct ElemMultVV_Helper<-1,size,add,ix1,T1,V1,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    typedef typename V3::real_type RT;
    enum { allunit = (V1::vstep == 1 && V2::vstep == 1 && V3::vstep == 1) };
    enum { algo = (
#if TMV_OPT >= 1
        // Unrolling doesn't ever seem to be faster for this function.
        //(size != UNKNOWN && size <= 32) ? 5 :
        //Traits<T1>::iscomplex ? 8 :
        (allunit && sizeof(RT) == 8) ? V3::viscomplex ? 2 : 3 :
        (allunit && sizeof(RT) == 4) ? V3::viscomplex ? 3 : 4 :
        //V3::viscomplex ? 6 :
#endif
        1 ) };
    static void call(
        const Scaling<ix1,T1>& x1, const V1& v1, const V2& v2, V3& v3)
    { ElemMultVV_Helper<algo,size,add,ix1,T1,V1,V2,V3>::call(x1,v1,v2,v3); }
    static void call2(const int n,
        const Scaling<ix1,T1>& x1, const IT1& x, const IT2& y, IT3& z)
    { ElemMultVV_Helper<algo,size,add,ix1,T1,V1,V2,V3>::call2(n,x1,x,y,z); }
  };

  // v3 (+=) x1 * v1 * v2
  template <bool add, int ix1, class T1, class V1, class V2, class V3>
  inline void InlineElemMultVV(
      const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
      const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
    TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    TMVAssert(v1.size() == v3.size());
    enum { size = Sizes<Sizes<V1::vsize,V2::vsize>::size,V3::vsize>::size };
    typedef typename V1::const_view_type V1v;
    typedef typename V2::const_view_type V2v;
    typedef typename V3::view_type V3v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    V3v v3v = v3.View();
    ElemMultVV_Helper<-1,size,add,ix1,T1,V1v,V2v,V3v>::call(x1,v1v,v2v,v3v);
  }

  // Defined in TMV_MultVV.cpp
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstElemMultVV(const T3 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddElemMultVV(const T3 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1,
      const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

} // namespace tmv

#endif 
