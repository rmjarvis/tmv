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

#ifndef TMV_SwapV_H
#define TMV_SwapV_H

#include "TMV_BaseVector.h"
#include "TMV_MultXV.h"

namespace tmv {

  //
  // Swap Vectors
  //

  template <int algo, int size, class V1, class V2> struct SwapV_Helper;

  // algo 1: simple for loop
  template <int size, class V1, class V2>
  struct SwapV_Helper<1,size,V1,V2> 
  {
    typedef typename V1::iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      for(int i=0;i<n;++i) TMV_SWAP(v1.ref(i),v2.ref(i));
    }
    static inline void call2(int n, IT1 it1, IT2 it2)
    { for(;n;--n) TMV_SWAP(*it1++,*it2++); }
  };

  // algo 2: 2 at a time
  template <int size, class V1, class V2>
  struct SwapV_Helper<2,size,V1,V2> 
  {
    typedef typename V1::iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      call2(n,v1.begin(),v2.begin());
    }
    static inline void call2(const int n, IT1 it1, IT2 it2)
    {
      typedef typename V1::value_type T1;
      T1 t0, t1;

      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      if (n_2) do {
        t0 = it1[0];
        t1 = it1[1];
        it1[0] = it2[0];
        it1[1] = it2[1]; it1 += 2;
        it2[0] = t0;
        it2[1] = t1; it2 += 2;
      } while (--n_2);
      if (nb) {
        t0 = *it1;
        *it1 = *it2;
        *it2 = t0;
      }
    }
  };

  // algo 3: 4 at a time
  template <int size, class V1, class V2>
  struct SwapV_Helper<3,size,V1,V2> 
  {
    typedef typename V1::iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      call2(n,v1.begin(),v2.begin());
    }
    static inline void call2(const int n, IT1 it1, IT2 it2)
    {
      typedef typename V1::value_type T1;
      T1 t0, t1, t2, t3;

      int n_4 = (n>>2);
      int nb = n-(n_4<<2);

      if (n_4) do {
        t0 = it1[0];
        t1 = it1[1];
        t2 = it1[2];
        t3 = it1[3];
        it1[0] = it2[0];
        it1[1] = it2[1];
        it1[2] = it2[2];
        it1[3] = it2[3]; it1 += 4;
        it2[0] = t0;
        it2[1] = t1;
        it2[2] = t2;
        it2[3] = t3; it2 += 4;
      } while (--n_4);
      if (nb) do {
        t0 = *it1;
        *it1++ = *it2;
        *it2++ = t0;
      } while (--nb);
    }
  };

  // algo 5: fully unroll
  template <int size, class V1, class V2>
  struct SwapV_Helper<5,size,V1,V2>
  {
    template <int I, int N>
    struct Unroller
    {
      static inline void unroll(V1& v1, V2& v2)
      {
        Unroller<I,N/2>::unroll(v1,v2);
        Unroller<I+N/2,N-N/2>::unroll(v1,v2);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void unroll(V1& v1, V2& v2)
      { TMV_SWAP(v1.ref(I),v2.ref(I)); }
    };
    template <int I>
    struct Unroller<I,0>
    { static inline void unroll(V1& v1, V2& v2) {} };
    static inline void call(V1& v1, V2& v2)
    { Unroller<0,size>::unroll(v1,v2); }
  };

  // algo 7: complex vectors with unit step, convert to real version
  template <int size, class V1, class V2>
  struct SwapV_Helper<7,size,V1,V2>
  {
    static inline void call(V1& v1, V2& v2)
    {
      typedef typename V1::flatten_type V1f;
      typedef typename V2::flatten_type V2f;
      typedef typename V1::real_type RT;
      enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
      enum { algo2 = (
          size != UNKNOWN && size <= 32 ? 5 :
          sizeof(RT) == 8 ? 2 :
          sizeof(RT) == 4 ? 3 :
          1 ) };
      V1f v1f = v1.Flatten();
      V2f v2f = v2.Flatten();
      SwapV_Helper<algo2,size2,V1f,V2f>::call(v1f,v2f);
    }
  };

  // algo 8: complex vectors, but not unit step
  template <int size, class V1, class V2>
  struct SwapV_Helper<8,size,V1,V2> 
  {
    typedef typename V1::iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      call2(n,v1.begin(),v2.begin());
    }
    static inline void call2(int n, IT1 it1, IT2 it2)
    {
      typedef typename V1::real_type RT;
      RT t0, t1;
      if (n) do {
        t0 = real(*it1);
        t1 = imag(*it1);
        real(*it1) = real(*it2);
        imag(*it1) = imag(*it2); ++it1;
        real(*it2) = t0;
        imag(*it2) = t1; ++it2;
      } while (--n);
    }
  };

  // algo 9: complex vectors, but v1 is conjugate
  template <int size, class V1, class V2>
  struct SwapV_Helper<9,size,V1,V2> 
  {
    typedef typename V1::iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      call2(n,v1.begin(),v2.begin());
    }
    static inline void call2(int n, IT1 it1, IT2 it2)
    {
      typedef typename IT1::nonconj_type IT1n;
      typedef typename V1::real_type RT;
      RT t0, t1;
      IT1n it1n = it1.NonConj();
      if (n) do {
        t0 = real(*it1n);
        t1 = imag(*it1n);
        real(*it1n) = real(*it2);
        imag(*it1n) = -imag(*it2); ++it1n;
        real(*it2) = t0;
        imag(*it2) = -t1; ++it2;
      } while (--n);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int size, class V1, class V2>
  struct SwapV_Helper<-1,size,V1,V2> 
  {
    typedef typename V1::iterator IT1;
    typedef typename V2::iterator IT2;
    typedef typename V1::real_type RT;
    static inline void call(V1& v1, V2& v2)
    {
      enum { allunit = (V1::vstep == 1 && V2::vstep == 1) };
      enum { algo = (
#if TMV_OPT >= 1
          // Strangely, algo 7 doesn't seem to be faster.
          //(V1::viscomplex && allunit && !V1::vconj) ? 7 :
          size != UNKNOWN && size <= (128/sizeof(RT)) ? 5 :
          (V1::viscomplex) ? (V1::vconj ? 9 : 8) :
          (sizeof(RT) == 8 && allunit) ? 2 :
          (sizeof(RT) == 4 && allunit) ? 3 :
#endif
          1 ) };
      SwapV_Helper<algo,size,V1,V2>::call(v1,v2);
    }
    static inline void call2(int n, IT1 it1, IT2 it2)
    {
      enum { allunit = (V1::vstep == 1 && V2::vstep == 1) };
      enum { algo = (
#if TMV_OPT >= 1
          (V1::viscomplex) ? (V1::vconj ? 9 : 8) :
          (sizeof(RT) == 8 && allunit) ? 2 :
          (sizeof(RT) == 4 && allunit) ? 3 :
#endif
          1 ) };
      SwapV_Helper<algo,size,V1,V2>::call2(n,it1,it2);
    }
  };

  template <class V1, class V2>
  inline void InlineSwap(
      BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    TMVStaticAssert((Traits2<T1,T2>::sametype));
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    enum { size = Sizes<V1::vsize,V2::vsize>::size };
    typedef typename V1::view_type V1v;
    typedef typename V2::view_type V2v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    SwapV_Helper<-1,size,V1v,V2v>::call(v1v,v2v);
  }

  // Defined in TMV_Vector.cpp
  template <class T, bool C>
  void InstSwap(VectorView<T,UNKNOWN,C> v1, VectorView<T> v2); // BLAS swap

  template <bool checkalias, bool conj2, bool inst, class V1, class V2>
  struct CallSwapv;

  template <bool conj2, bool inst, class V1, class V2>
  struct CallSwapv<true,conj2,inst,V1,V2>
  {
    static inline void call(V1& v1, V2& v2)
    {
      if (!SameStorage(v1,v2)) 
      {
        CallSwapv<false,conj2,inst,V1,V2>::call(v1,v2);
      }
      else
      {
        if (v1.step() == v2.step()) 
        {
          if (V1::vconj != int(V2::vconj)) 
          {
            v1.ConjugateSelf();
            v2.ConjugateSelf();
          }
          else {}  // They are equal
        } 
        else if (v1.step()*v2.step() < 0) // then no clobbering
        {
          CallSwapv<false,conj2,inst,V1,V2>::call(v1,v2);
        }
        else // Have to swap with a full temporary
        { 
          typename V1::copy_type v1c = v1.copy();
          v1 = v2;
          v2 = v1c;
        }
      } 
    }
  };
  template <bool inst, class V1, class V2>
  struct CallSwapv<false,true,inst,V1,V2> // conj = true
  {
    static inline void call(V1& v1, V2& v2)
    { 
      typedef typename V1::conjugate_type V1c;
      typedef typename V2::conjugate_type V2c;
      V1c v1c = v1.Conjugate();
      V2c v2c = v2.Conjugate();
      CallSwapv<false,false,inst,V1c,V2c>::call(v1c,v2c);
    }
  };
  template <class V1, class V2>
  struct CallSwapv<false,false,false,V1,V2> // inst = false
  {
    static inline void call(V1& v1, V2& v2)
    { InlineSwap(v1,v2); }
  };
  template <class V1, class V2>
  struct CallSwapv<false,false,true,V1,V2> // inst = true
  {
    static inline void call(V1& v1, V2& v2)
    { InstSwap(v1.XView(),v2.XView()); }
  };

  template <class V1, class V2> 
  inline void DoSwap(BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    enum { checkalias = (
        V1::vsize == UNKNOWN &&
        V2::vsize == UNKNOWN) };
    enum { inst = (
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits2<T1,T2>::sametype &&
        checkalias) };
    CallSwapv<checkalias,V2::vconj,inst,V1,V2>::call(v1.vec(),v2.vec());
  }
  template <class V1, class V2> 
  inline void Swap(BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
  { DoSwap(v1,v2); }


  //
  // ReverseSelf
  //

  template <int algo, int size, class V> struct ReverseV_Helper;

  // algo 1: simple for loop
  template <int size, class V>
  struct ReverseV_Helper<1,size,V>
  {
    static inline void call(V& v)
    { 
      const int n = size == UNKNOWN ? int(v.size()) : size;
      const int no2 = n/2;
      if (no2) 
        for(int i1=0;i1<no2;++i1) v.CSwap(i1,n-i1-1);
    }
  };

  // algo 2: call swap on two halves
  template <int size, class V>
  struct ReverseV_Helper<2,size,V> 
  {
    static inline void call(V& v)
    {
      typedef typename V::value_type T;
      typedef typename V::subvector_type V1;
      typedef typename V::subvector_type::reverse_type V2;
      const int n = size == UNKNOWN ? int(v.size()) : size;
      enum { sizeo2 = size == UNKNOWN ? UNKNOWN : size/2 };
      if (n > 1)  {
        V1 v1 = v.CSubVector(0,n/2);
        V2 v2 = v.CSubVector(n-n/2,n).Reverse();
        enum { algo2 = (
#if TMV_OPT >= 1
            size != UNKNOWN && size <= 64 ? 5 :
            sizeof(T) == 8 ? 2 :
            sizeof(T) == 4 ? 4 :
#endif
            1 ) };
        SwapV_Helper<algo2,sizeo2,V1,V2>::call(v1,v2);
      }
    }
  };

  // algo 5: fully unroll
  template <int size, class V>
  struct ReverseV_Helper<5,size,V> 
  {
    template <int I, int N>
    struct Unroller
    {
      static inline void unroll(V& v)
      {
        Unroller<I,N-1>::unroll(v);
        v.CSwap(N-1,size-N);
      }
    };
    template <int I>
    struct Unroller<I,0>
    { static inline void unroll(V& v) {} };
    static inline void call(V& v)
    { Unroller<0,size/2>::unroll(v); }
  };

  template <class V>
  inline void InlineReverseSelf(BaseVector_Mutable<V>& v)
  {
    const int size = V::vsize;
    enum { algo = (
#if TMV_OPT >= 1
        size != UNKNOWN && size <= 32 ? 5 :
        V::viscomplex ? 2 :
#endif
        1 ) };
    typedef typename V::view_type Vv;
    Vv vv = v.View();
    ReverseV_Helper<algo,size,Vv>::call(vv);
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  void InstReverseSelf(VectorView<T> v);
  template <class T>
  inline void InstReverseSelf(VectorView<T,UNKNOWN,true> v)
  { return InstReverseSelf(v.Conjugate()); }
  
  template <bool inst, class V>
  struct CallReverseSelf // inst = false
  {
    static inline void call(V& v)
    { InlineReverseSelf(v); }
  };
  template <class V>
  struct CallReverseSelf<true,V> // inst = true
  {
    static inline void call(V& v)
    { InstReverseSelf(v.XView()); }
  };

  template <class V>
  inline void ReverseSelf(BaseVector_Mutable<V>& v)
  {
    typedef typename V::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    CallReverseSelf<inst,V>::call(v.vec());
  }


  //
  // ConjugateSelf
  //

  template <int algo, int size, class V> struct ConjugateV_Helper;

  // algo 1: simple for loop
  template <int size, class V>
  struct ConjugateV_Helper<1,size,V> 
  {
    static inline void call(V& v)
    { 
      const int n=v.size();
      for(int i=0;i<n;++i) v.ref(i) = TMV_CONJ(v.cref(i));
    }
  };

  // algo 2: v.Imag() *= -1
  template <int size, class V>
  struct ConjugateV_Helper<2,size,V> 
  {
    static inline void call(V& v)
    { 
      typedef typename V::real_type RT;
      typedef typename V::imagview_type Vi;
      const Scaling<-1,RT> mone;
      Vi vi = v.Imag();
      MultXV_Helper<-1,size,-1,RT,Vi>::call(mone,vi);
    }
  };

  // algo 5: fully unroll
  template <int size, class V>
  struct ConjugateV_Helper<5,size,V>
  {
    template <int I, int N>
    struct Unroller
    {
      static inline void unroll(V& v)
      {
        Unroller<I,N/2>::unroll(v);
        Unroller<I+N/2,N-N/2>::unroll(v);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void unroll(V& v)
      { v.ref(I) = TMV_CONJ(v.cref(I)); }
    };
    template <int I>
    struct Unroller<I,0>
    { static inline void unroll(V& v) {} };
    static inline void call(V& v)
    { Unroller<0,size>::unroll(v); }
  };

  // algo -1: Determine which algorithm to use
  template <int size, class V>
  struct ConjugateV_Helper<-1,size,V>
  {
    static inline void call(V& v)
    {
      enum { algo = (
#if TMV_OPT >= 1
          size != UNKNOWN && size <= 32 ? 5 :
#endif
          2 ) };
      ConjugateV_Helper<algo,size,V>::call(v);
    }
  };

  template <class V>
  inline void InlineConjugateSelf(BaseVector_Mutable<V>& v)
  {
    typedef typename V::view_type Vv;
    Vv vv = v.View();
    ConjugateV_Helper<-1,V::vsize,Vv>::call(vv); 
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  void InstConjugateSelf(VectorView<T> v);
  template <class T>
  inline void InstConjugateSelf(VectorView<T,UNKNOWN,true> v)
  { InstConjugateSelf(v.Conjugate()); }

  template <bool iscomplex, bool inst, class V>
  struct CallConjugateSelfv // iscomplex = false
  { static inline void call(V& v) {} };
  template <class V>
  struct CallConjugateSelfv<true,false,V> // inst = false
  { static inline void call(V& v) { InlineConjugateSelf(v); } };
  template <class V>
  struct CallConjugateSelfv<true,true,V> // inst = true
  { static inline void call(V& v) { InstConjugateSelf(v.XView()); } };

  template <class V>
  inline void ConjugateSelf(BaseVector_Mutable<V>& v)
  {
    typedef typename V::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    CallConjugateSelfv<V::viscomplex,inst,V>::call(v.vec());
  }

} // namespace tmv

#endif
