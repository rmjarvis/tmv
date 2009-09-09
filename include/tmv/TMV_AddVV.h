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


#ifndef TMV_AddVV_H
#define TMV_AddVV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

  //
  // Vector += x * Vector
  //

  template <int algo, int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper;

  // algo 0: trivial: ix == 1, !add, so call Copy
  template <int size, class T, class V1, class V2>
  struct AddVV_Helper<0,size,false,1,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<1,T>& , const V1& v1, V2& v2)
    { CopyV_Helper<-1,size,V1,V2>::call(v1,v2); }
    static inline void call2(int n, const Scaling<1,T>& , IT1 it1, IT2 it2)
    { CopyV_Helper<-1,size,V1,V2>::call2(it1,it2); }
  };

  // algo 1: simple for loop
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<1,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      for(int i=0;i<n;++i) 
        Maybe<add>::add(v2.ref(i) , ZProd<false,false>::prod(x , v1.cref(i)));
    }
    static inline void call2(int n, const Scaling<ix,T>& x, IT1 it1, IT2 it2)
    { 
      enum { c1 = V1::vconj };
      for(;n;--n) 
        Maybe<add>::add(*it2++ , ZProd<false,c1>::prod(x , *it1++)); 
    }
  };  

  // algo 2: 2 at a time
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<2,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x,v1.NonConj().begin(),v2.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix,T>& x, IT1 it1, IT2 it2)
    {
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };

      if (n_2) do {
        Maybe<add>::add(it2[0] , ZProd<false,c1>::prod(x , it1[0]));
        Maybe<add>::add(it2[1] , ZProd<false,c1>::prod(x , it1[1]));
        it1 += 2; it2 += 2;
      } while (--n_2);
      if (nb) {
        Maybe<add>::add(*it2 , ZProd<false,c1>::prod(x , *it1)); 
      }
    }
  };

  // algo 3: 4 at a time
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<3,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x,v1.NonConj().begin(),v2.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix,T>& x, IT1 it1, IT2 it2)
    {
      int n_4 = (n>>2);
      int nb = n-(n_4<<2);
      enum { c1 = V1::vconj };

      if (n_4) do {
        Maybe<add>::add(it2[0] , ZProd<false,c1>::prod(x , it1[0]));
        Maybe<add>::add(it2[1] , ZProd<false,c1>::prod(x , it1[1]));
        Maybe<add>::add(it2[2] , ZProd<false,c1>::prod(x , it1[2]));
        Maybe<add>::add(it2[3] , ZProd<false,c1>::prod(x , it1[3]));
        it1 += 4; it2 += 4;
      } while (--n_4);
      if (nb) do {
        Maybe<add>::add(*it2++ , ZProd<false,c1>::prod(x , *it1++)); 
      } while (--nb);
    }
  };

  // algo 4: 8 at a time
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<4,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x,v1.NonConj().begin(),v2.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix,T>& x, IT1 it1, IT2 it2)
    {
      int n_8 = (n>>3);
      int nb = n-(n_8<<3);
      enum { c1 = V1::vconj };

      if (n_8) do {
        Maybe<add>::add(it2[0] , ZProd<false,c1>::prod(x , it1[0]));
        Maybe<add>::add(it2[1] , ZProd<false,c1>::prod(x , it1[1]));
        Maybe<add>::add(it2[2] , ZProd<false,c1>::prod(x , it1[2]));
        Maybe<add>::add(it2[3] , ZProd<false,c1>::prod(x , it1[3]));
        Maybe<add>::add(it2[4] , ZProd<false,c1>::prod(x , it1[4]));
        Maybe<add>::add(it2[5] , ZProd<false,c1>::prod(x , it1[5]));
        Maybe<add>::add(it2[6] , ZProd<false,c1>::prod(x , it1[6]));
        Maybe<add>::add(it2[7] , ZProd<false,c1>::prod(x , it1[7]));
        it1 += 8; it2 += 8;
      } while (--n_8);
      if (nb) do {
        Maybe<add>::add(*it2++ , ZProd<false,c1>::prod(x , *it1++)); 
      } while (--nb);
    }
  };

  // algo 5: fully unroll
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<5,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    template <int I, int N, bool iscomplex>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, const V1& v1, V2& v2)
      {
        Unroller<I,N/2,iscomplex>::unroll(x,v1,v2);
        Unroller<I+N/2,N-N/2,iscomplex>::unroll(x,v1,v2);
      }
      static inline void unroll2(const Scaling<ix,T>& x, 
          const IT1& it1, const IT2& it2)
      {
        Unroller<I,N/2,iscomplex>::unroll2(x,it1,it2);
        Unroller<I+N/2,N-N/2,iscomplex>::unroll2(x,it1,it2);
      }
    };
    template <int I>
    struct Unroller<I,1,false>
    {
      static inline void unroll(const Scaling<ix,T>& x, const V1& v1, V2& v2)
      { Maybe<add>::add(v2.ref(I) , ZProd<false,false>::prod(x , v1.cref(I))); }
      static inline void unroll2(const Scaling<ix,T>& x, 
          const IT1& it1, const IT2& it2)
      {
        enum { c1 = V1::vconj };
        Maybe<add>::add(it2[I] , ZProd<false,c1>::prod(x , it1[I])); 
      }
    };
    template <int I>
    struct Unroller<I,1,true>
    {
      static inline void unroll(const Scaling<ix,T>& x, const V1& v1, V2& v2)
      {
        typedef typename V2::real_type RT;
        typedef typename V2::value_type VT;
        const RT rv = ZProd<false,false>::rprod(x,v1.cref(I));
        const RT iv = ZProd<false,false>::iprod(x,v1.cref(I));
        Maybe<add>::add(v2.ref(I) , VT(rv,iv));
      }
      static inline void unroll2(const Scaling<ix,T>& x, 
          const IT1& it1, const IT2& it2)
      {
        typedef typename V2::real_type RT;
        typedef typename V2::value_type VT;
        const RT rv = ZProd<false,false>::rprod(x,it1[I]);
        const RT iv = ZProd<false,false>::iprod(x,it1[I]);
        Maybe<add>::add(it2[I] , VT(rv,iv));
      }
    };
    template <int I, bool iscomplex>
    struct Unroller<I,0,iscomplex>
    {
      static inline void unroll(const Scaling<ix,T>& , const V1& , V2& ) {}
      static inline void unroll2(const Scaling<ix,T>& , 
          const IT1& , const IT2& ) {}
    };
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    { Unroller<0,size,V2::viscomplex>::unroll(x,v1,v2); }
    static inline void call2(const int , const Scaling<ix,T>& x, 
        const IT1& it1, const IT2& it2)
    { Unroller<0,size,V2::viscomplex>::unroll2(x,it1,it2); }
  };

  // algo 7: complex vector with unit step, convert to real version
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<7,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      typedef typename V1::const_flatten_type V1f;
      typedef typename V2::flatten_type V2f;
      typedef typename V2::real_type RT;
      enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
      enum { algo2 = (
          (size2 != UNKNOWN && size2 <= (128/sizeof(RT))) ? 5 :
          (sizeof(RT) == 8) ? 2 :
          (sizeof(RT) == 4) ? 3 :
          1 ) };
      V1f v1f = v1.Flatten();
      V2f v2f = v2.Flatten();
      AddVV_Helper<algo2,size2,add,ix,T,V1f,V2f>::call(x,v1f,v2f);
    }
    static inline void call2(const int n, const Scaling<ix,T>& x, 
        const IT1& it1, const IT2& it2)
    {
      typedef typename V1::const_flatten_type V1f;
      typedef typename V2::flatten_type V2f;
      typedef typename V1f::const_iterator IT1f;
      typedef typename V2f::iterator IT2f;
      typedef typename V2::real_type RT;
      enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
      const int n2 = size == UNKNOWN ? 2*n : size2;
      enum { algo2 = (
          (size2 != UNKNOWN && size2 <= (128/sizeof(RT))) ? 5 :
          (sizeof(RT) == 8) ? 2 :
          (sizeof(RT) == 4) ? 3 :
          1 ) };
      IT1f it1f = it1.Flatten();
      IT2f it2f = it2.Flatten();
      AddVV_Helper<algo2,size2,add,ix,T,V1f,V2f>::call2(n2,x,it1f,it2f);
    }
  };

  // algo 8: complex vector, but not flatten-able
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<8,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x,v1.NonConj().begin(),v2.begin());
    }
    static inline void call2(int n,
        const Scaling<ix,T>& x, IT1 it1, IT2 it2)
    {
      typedef typename V1::value_type T1;
      typedef typename V2::real_type RT2;
      typedef typename V2::value_type T2;
      T1 val1;
      RT2 rv2,iv2;

      enum { c1 = V1::vconj };

      if (n) do {
        val1 = *it1++;
        rv2 = ZProd<false,c1>::rprod(x,val1);
        iv2 = ZProd<false,c1>::iprod(x,val1);
        Maybe<add>::add(*it2++ , T2(rv2,iv2));
      } while (--n);
    }
  };

  // algo 9: complex vector, 2 at a time
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<9,size,add,ix,T,V1,V2> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v2.size()) : size;
      call2(n,x,v1.NonConj().begin(),v2.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix,T>& x, IT1 it1, IT2 it2)
    {
      typedef typename V1::value_type T1;
      typedef typename V2::real_type RT2;
      typedef typename V2::value_type T2;
      T1 val1a, val1b;
      RT2 rv2a,iv2a,rv2b,iv2b;

      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };

      if (n_2) do {
        val1a = it1[0];
        val1b = it1[1];
        it1 += 2;
        rv2a = ZProd<false,c1>::rprod(x,val1a);
        iv2a = ZProd<false,c1>::iprod(x,val1a);
        rv2b = ZProd<false,c1>::rprod(x,val1b);
        iv2b = ZProd<false,c1>::iprod(x,val1b);
        Maybe<add>::add(it2[0] , T2(rv2a,iv2a));
        Maybe<add>::add(it2[1] , T2(rv2b,iv2b));
        it2 += 2;
      } while (--n_2);
      if (nb) {
        val1a = *it1;
        rv2a = ZProd<false,c1>::rprod(x,val1a);
        iv2a = ZProd<false,c1>::iprod(x,val1a);
        Maybe<add>::add(*it2 , T2(rv2a,iv2a));
      }
    }
  };

  // algo -1: Determine which algorithm to use
  template <int size, bool add, int ix, class T, class V1, class V2>
  struct AddVV_Helper<-1,size,add,ix,T,V1,V2> 
  {
    typedef typename V2::value_type VT;
    typedef typename V2::real_type RT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::iterator IT2;
    enum { allunit = V1::vstep == 1 && V2::vstep == 1 };
    enum { flatten = (
        Traits<T>::isreal && allunit &&
        V1::viscomplex && V2::viscomplex &&
        V1::vconj == int(V2::vconj) ) };
    enum { algo = (
        ( ix == 1 && !add ) ? 0 :
#if TMV_OPT >= 1
        flatten ? 7 :
        ( size != UNKNOWN && size <= (128/sizeof(VT)) ) ? 5 :
        ( sizeof(RT) == 8 && allunit ) ? ( V2::viscomplex ? 8 : 2 ) :
        ( sizeof(RT) == 4 && allunit ) ? ( V2::viscomplex ? 9 : 3 ) :
        V2::viscomplex ? 8 :
#endif
        1 ) };
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    { AddVV_Helper<algo,size,add,ix,T,V1,V2>::call(x,v1,v2); }
    static inline void call2(const int n,
        const Scaling<ix,T>& x, const IT1& it1, const IT2& it2)
    {
      //std::cout<<"AddVV_Helper algo "<<algo<<std::endl;
      AddVV_Helper<algo,size,add,ix,T,V1,V2>::call2(n,x,it1,it2); 
    }
  };

  // v2 += x * v1
  template <int ix, class T, class V1, class V2>
  inline void InlineAddVV(
      const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
      BaseVector_Mutable<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    enum { size = Sizes<V1::vsize,V2::vsize>::size };
    typedef typename V1::const_view_type V1v;
    typedef typename V2::view_type V2v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    AddVV_Helper<-1,size,true,ix,T,V1v,V2v>::call(x,v1v,v2v);
  }

  // Defined in TMV_AddVV.cpp
  template <class T1, bool C1, class T2>  // BLAS axpy
  void InstAddVV(const T2 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1, VectorView<T2> v2);


  //
  // Vector + Vector
  //

  template <int algo, int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2;

  // algo 1: simple for loop
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<1,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      for(int i=0;i<n;++i) 
        v3.ref(i) = (
            ZProd<false,false>::prod(x1 , v1.cref(i)) +
            ZProd<false,false>::prod(x2 , v2.cref(i)));
    }
    static inline void call2(int n,
        const Scaling<ix1,T1>& x1, IT1 it1,
        const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
    {
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };
      for(;n;--n) 
        *it3++ = (
            ZProd<false,c1>::prod(x1 , *it1++) + 
            ZProd<false,c2>::prod(x2 , *it2++)); 
    }
  };  

  // algo 2: 2 at a time
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<2,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),x2,v2.NonConj().begin(),v3.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 it1,
        const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
    {
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_2) do {
        it3[0] = (
            ZProd<false,c1>::prod(x1 , it1[0]) +
            ZProd<false,c2>::prod(x2 , it2[0]));
        it3[1] = (
            ZProd<false,c1>::prod(x1 , it1[1]) +
            ZProd<false,c2>::prod(x2 , it2[1]));
        it1 += 2; it2 += 2; it3 += 2;
      } while (--n_2);
      if (nb) {
        *it3 = (
            ZProd<false,c1>::prod(x1 , *it1) +
            ZProd<false,c2>::prod(x2 , *it2)); 
      }
    }
  };

  // algo 3: 4 at a time
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<3,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),x2,v2.NonConj().begin(),v3.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 it1,
        const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
    {
      int n_4 = (n>>2);
      int nb = n-(n_4<<2);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_4) do {
        it3[0] = (
            ZProd<false,c1>::prod(x1 , it1[0]) +
            ZProd<false,c2>::prod(x2 , it2[0]));
        it3[1] = (
            ZProd<false,c1>::prod(x1 , it1[1]) +
            ZProd<false,c2>::prod(x2 , it2[1]));
        it3[2] = (
            ZProd<false,c1>::prod(x1 , it1[2]) +
            ZProd<false,c2>::prod(x2 , it2[2]));
        it3[3] = (
            ZProd<false,c1>::prod(x1 , it1[3]) +
            ZProd<false,c2>::prod(x2 , it2[3]));
        it1 += 4; it2 += 4; it3 += 4;
      } while (--n_4);
      if (nb) do {
        *it3++ = (
            ZProd<false,c1>::prod(x1 , *it1++) +
            ZProd<false,c2>::prod(x2 , *it2++)); 
      } while (--nb);
    }
  };

  // algo 4: 8 at a time
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<4,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),x2,v2.NonConj().begin(),v3.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 it1,
        const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
    {
      int n_8 = (n>>3);
      int nb = n-(n_8<<3);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_8) do {
        it3[0] = (
            ZProd<false,c1>::prod(x1 , it1[0]) +
            ZProd<false,c2>::prod(x2 , it2[0]));
        it3[1] = (
            ZProd<false,c1>::prod(x1 , it1[1]) +
            ZProd<false,c2>::prod(x2 , it2[1]));
        it3[2] = (
            ZProd<false,c1>::prod(x1 , it1[2]) +
            ZProd<false,c2>::prod(x2 , it2[2]));
        it3[3] = (
            ZProd<false,c1>::prod(x1 , it1[3]) +
            ZProd<false,c2>::prod(x2 , it2[3]));
        it3[4] = (
            ZProd<false,c1>::prod(x1 , it1[4]) +
            ZProd<false,c2>::prod(x2 , it2[4]));
        it3[5] = (
            ZProd<false,c1>::prod(x1 , it1[5]) +
            ZProd<false,c2>::prod(x2 , it2[5]));
        it3[6] = (
            ZProd<false,c1>::prod(x1 , it1[6]) +
            ZProd<false,c2>::prod(x2 , it2[6]));
        it3[7] = (
            ZProd<false,c1>::prod(x1 , it1[7]) +
            ZProd<false,c2>::prod(x2 , it2[7]));
        it1 += 8; it2 += 8; it3 += 8;
      } while (--n_8);
      if (nb) do {
        *it3++ = (
            ZProd<false,c1>::prod(x1 , *it1++) +
            ZProd<false,c2>::prod(x2 , *it2++)); 
      } while (--nb);
    }
  };

  // algo 5: fully unroll
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<5,size,ix1,T1,V1,ix2,T2,V2,V3> // known size, unroll
  {
    template <int I, int N>
    struct Unroller
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const V1& v1,
          const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
      {
        Unroller<I,N/2>::unroll(x1,v1,x2,v2,v3);
        Unroller<I+N/2,N-N/2>::unroll(x1,v1,x2,v2,v3);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const V1& v1,
          const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
      {
        v3.ref(I) = (
            ZProd<false,false>::prod(x1 , v1.cref(I)) +
            ZProd<false,false>::prod(x2 , v2.cref(I))); 
      }
    };
    template <int I>
    struct Unroller<I,0>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const V1& v1,
          const Scaling<ix2,T2>& x2, const V2& v2, V3& v3) {}
    };
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    { Unroller<0,size>::unroll(x1,v1,x2,v2,v3); }
  };

  // algo 7: complex vectors with unit step, convert to real version
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<7,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      typedef typename V1::const_flatten_type V1f;
      typedef typename V2::const_flatten_type V2f;
      typedef typename V3::flatten_type V3f;
      typedef typename V1::real_type RT;
      enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
      enum { algo2 = (
          (size2 != UNKNOWN && size2 <= (128/sizeof(RT))) ? 5 :
          (sizeof(RT) == 8) ? 2 :
          (sizeof(RT) == 4) ? 3 :
          1 ) };
      V1f v1f = v1.Flatten();
      V2f v2f = v2.Flatten();
      V3f v3f = v3.Flatten();
      AddVV_Helper2<algo2,size2,ix1,T1,V1f,ix2,T2,V2f,V3f>::call(
          x1,v1f,x2,v2f,v3f);
    }
  };

  // algo 8: complex vector, but not flatten-able
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<8,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),x2,v2.NonConj().begin(),v3.begin());
    }
    static inline void call2(int n,
        const Scaling<ix1,T1>& x1, IT1 it1,
        const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
    {
      typedef typename V1::value_type VT1;
      typedef typename V2::value_type VT2;
      typedef typename V3::value_type VT3;
      typedef typename V3::real_type RT3;
      VT1 val1;
      VT2 val2;
      RT3 rv3,iv3;

      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n) do {
        val1 = *it1++;
        val2 = *it2++;
        rv3 = (
            ZProd<false,c1>::rprod(x1,val1) +
            ZProd<false,c2>::rprod(x2,val2));
        iv3 = (
            ZProd<false,c1>::iprod(x1,val1) +
            ZProd<false,c2>::iprod(x2,val2));
        *it3++ = VT3(rv3,iv3);
      } while (--n);
    }
  };

  // algo 9: complex vectors, 2 at a time
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<9,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      const int n = size == UNKNOWN ? int(v3.size()) : size;
      call2(n,x1,v1.NonConj().begin(),x2,v2.NonConj().begin(),v3.begin());
    }
    static inline void call2(const int n,
        const Scaling<ix1,T1>& x1, IT1 it1,
        const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
    {
      typedef typename V1::value_type VT1;
      typedef typename V2::value_type VT2;
      typedef typename V3::value_type VT3;
      typedef typename V3::real_type RT3;
      VT1 val1a, val1b;
      VT2 val2a, val2b;
      RT3 rv3a,iv3a,rv3b,iv3b;

      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { c1 = V1::vconj };
      enum { c2 = V2::vconj };

      if (n_2) do {
        val1a = it1[0]; val2a = it2[0];
        val1b = it1[1]; val2b = it2[1];
        it1 += 2; it2 += 2;
        rv3a = (
            ZProd<false,c1>::rprod(x1,val1a) +
            ZProd<false,c2>::rprod(x2,val2a));
        iv3a = (
            ZProd<false,c1>::iprod(x1,val1a) +
            ZProd<false,c2>::iprod(x2,val2a));
        rv3b = (
            ZProd<false,c1>::rprod(x1,val1b) +
            ZProd<false,c2>::rprod(x2,val2b));
        iv3b = (
            ZProd<false,c1>::iprod(x1,val1b) +
            ZProd<false,c2>::iprod(x2,val2b));
        it3[0] = VT3(rv3a,iv3a); it3[1] = VT3(rv3b,iv3b);
        it3 += 2;
      } while (--n_2);
      if (nb) {
        val1a = *it1;
        val2a = *it2;
        rv3a = (
            ZProd<false,c1>::rprod(x1,val1a) +
            ZProd<false,c2>::rprod(x2,val2a));
        iv3a = (
            ZProd<false,c1>::iprod(x1,val1a) +
            ZProd<false,c2>::iprod(x2,val2a));
        *it3 = VT3(rv3a,iv3a);
      }
    }
  };

  // algo -1: Determine which algorithm to use
  template <int size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct AddVV_Helper2<-1,size,ix1,T1,V1,ix2,T2,V2,V3> 
  {
    typedef typename V3::value_type VT;
    typedef typename V2::real_type RT;
    typedef typename V1::const_nonconj_type::const_iterator IT1;
    typedef typename V2::const_nonconj_type::const_iterator IT2;
    typedef typename V3::iterator IT3;
    static inline void call(
        const Scaling<ix1,T1>& x1, const V1& v1,
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      enum { allcomplex = V1::viscomplex && V2::viscomplex && V3::viscomplex };
      enum { allunit = V1::vstep == 1 && V2::vstep == 1 && V3::vstep == 1 };
      enum { flatten = (
          allunit && allcomplex &&
          Traits<T1>::isreal && Traits<T2>::isreal &&
          V1::vconj == int(V3::vconj) && V2::vconj == int(V3::vconj) ) };
      enum { algo = (
#if TMV_OPT >= 1
          flatten ? 7 :
          (size != UNKNOWN && size <= (128/sizeof(VT))) ? 5 :
          (sizeof(RT) == 8 && allunit) ? ( V3::viscomplex ? 8 : 2 ) :
          (sizeof(RT) == 4 && allunit) ? ( 
            ( V3::viscomplex ? ( allcomplex ? 8 : 9 ) : 3 ) ) :
          V3::viscomplex ? 8 :
#endif
          1 ) };
      AddVV_Helper2<algo,size,ix1,T1,V1,ix2,T2,V2,V3>::call(x1,v1,x2,v2,v3);
    }
    static inline void call2(const int n,
        const Scaling<ix1,T1>& x1, const IT1& it1,
        const Scaling<ix2,T2>& x2, const IT2& it2, const IT3& it3)
    {
      enum { allcomplex = V1::viscomplex && V2::viscomplex && V3::viscomplex };
      enum { allunit = V1::vstep == 1 && V2::vstep == 1 && V3::vstep == 1 };
      enum { algo = (
#if TMV_OPT >= 1
          (sizeof(RT) == 8 && allunit) ? ( V3::viscomplex ? 8 : 2 ) :
          (sizeof(RT) == 4 && allunit) ? ( 
            ( V3::viscomplex ? ( allcomplex ? 8 : 9 ) : 3 ) ) :
          V3::viscomplex ? 8 :
#endif
          1 ) };
      AddVV_Helper2<algo,size,ix1,T1,V1,ix2,T2,V2,V3>::call2(
          n,x1,it1,x2,it2,it3);
    } 
  };

  // v3 = x1 * v1 + x2 * v2
  template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  inline void InlineAddVV(
      const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
      const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
      BaseVector_Mutable<V3>& v3)
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
    AddVV_Helper2<-1,size,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(x1,v1v,x2,v2v,v3v);
  }

} // namespace tmv

#endif 
