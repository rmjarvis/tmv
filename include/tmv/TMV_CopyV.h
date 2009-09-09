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

#ifndef TMV_CopyV_H
#define TMV_CopyV_H

#include "TMV_BaseVector.h"
#include <string.h> // for memmove

namespace tmv {

  //
  // Copy Vectors
  //

  template <int algo, int size, class V1, class V2> struct CopyV_Helper;

  // algo 1: simple for loop
  template <int size, class V1, class V2>
  struct CopyV_Helper<1,size,V1,V2>
  {
    typedef typename V1::const_iterator IT1;
    typedef typename V2::iterator IT2;
    static inline void call(const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      for(int i=0;i<n;++i) v2.ref(i) = v1.cref(i); 
    }
    static inline void call2(int n, IT1 it1, IT2 it2)
    { for(;n;--n) *it2++ = *it1++; }
  };

  // algo 2: memmove
  template <int size, class V1, class V2>
  struct CopyV_Helper<2,size,V1,V2>
  {
    static inline void call(const V1& v1, V2& v2)
    {
      const int n = size == UNKNOWN ? int(v1.size()) : size;
      memmove(v2.ptr(),v1.cptr(),n*sizeof(typename V2::value_type));
    }
    static inline void call2(const int n, 
        typename V1::const_iterator it1, typename V2::iterator it2)
    { memmove(it2.GetP(),it1.GetP(),n*sizeof(typename V2::value_type)); }
  };

  // algo 3: std::copy
  template <int size, class V1, class V2>
  struct CopyV_Helper<3,size,V1,V2> 
  {
    static inline void call(const V1& v1, V2& v2)
    { std::copy(v1.begin(),v1.end(),v2.begin()); }
  };

  // algo 4: fully unroll
  template <int size, class V1, class V2>
  struct CopyV_Helper<4,size,V1,V2>
  {
    template <int I, int N>
    struct Unroller
    {
      static inline void dounroll(const V1& v1, V2& v2)
      {
        Unroller<I,N/2>::dounroll(v1,v2);
        Unroller<I+N/2,N-N/2>::dounroll(v1,v2);
        //Unroller<I,N-1>::dounroll(v1,v2);
        //Unroller<I+N-1,1>::dounroll(v1,v2);
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline void dounroll(const V1& v1, V2& v2)
      { v2.ref(I) = v1.cref(I); }
    };
    template <int I>
    struct Unroller<I,0>
    { static inline void dounroll(const V1& v1, V2& v2) {} };
    static inline void call(const V1& v1, V2& v2)
    { Unroller<0,size>::dounroll(v1,v2); }
  };

  // algo -1: Determine which algorithm to use
  template <int size, class V1, class V2>
  struct CopyV_Helper<-1,size,V1,V2>
  {
    static inline void call(const V1& v1, V2& v2)
    {
      typedef typename V1::value_type T1;
      typedef typename V2::value_type T2;
      enum { algo = (
#if TMV_OPT > 0
          // only very small vectors seem to be faster unrolled
          ( size != UNKNOWN && size <= 8 ) ? 4 :
          ( Traits2<T1,T2>::sametype && V1::vconj == int(V2::vconj) &&
            V1::vstep == 1 && V2::vstep == 1 ) ? 2 :
#endif
          1 ) };
      CopyV_Helper<algo,size,V1,V2>::call(v1,v2);
    }
    static inline void call2(const int n,
        typename V1::const_iterator it1, typename V2::iterator it2)
    {
      typedef typename V1::value_type T1;
      typedef typename V2::value_type T2;
      enum { algo = (
#if TMV_OPT > 0
          ( Traits2<T1,T2>::sametype && V1::vconj == int(V2::vconj) &&
            V1::vstep == 1 && V2::vstep == 1 ) ? 2 :
#endif
          1 ) };
      CopyV_Helper<algo,size,V1,V2>::call2(n,it1,it2);
    }
  };

  template <class V1, class V2>
  inline void InlineCopy(
      const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    enum { size = Sizes<V1::vsize,V2::vsize>::size };
    typedef typename V1::const_view_type V1v;
    typedef typename V2::view_type V2v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    CopyV_Helper<-1,size,V1v,V2v>::call(v1v,v2v); 
  }

  // Defined in TMV_Vector.cpp
  template <class T, bool C>
  void InstCopy(const ConstVectorView<T,UNKNOWN,C>& v1,
      VectorView<T> v2); // BLAS copy
  template <class T, bool C>
  inline void InstCopy(const ConstVectorView<T,UNKNOWN,C>& v1,
      VectorView<T,UNKNOWN,true> v2)
  { InstCopy(v1.Conjugate(),v2.Conjugate()); }

  template <bool checkalias, bool inst, class V1, class V2>
  struct CallCopyV  // checkalias = true
  {
    static inline void call(const V1& v1, V2& v2)
    {
      if (!SameStorage(v1,v2)) 
      {
        CallCopyV<false,inst,V1,V2>::call(v1,v2);
      }
      else 
      {
        if (v1.step() == v2.step()) 
        {
          if (V1::vconj != int(V2::vconj)) v2.ConjugateSelf();
          else {}  // They are already equal.
        } 
        else if (v1.step()*v2.step() < 0 || 
            std::abs(v2.step()) < std::abs(v1.step())) 
        { // then no clobering
          CallCopyV<false,inst,V1,V2>::call(v1,v2);
        } 
        else 
        {
          v2 = v1.copy();
        }
      }
    }
  };
  template <class V1, class V2>
  struct CallCopyV<false,false,V1,V2> // inst = false
  {
    static inline void call(const V1& v1, V2& v2)
    { InlineCopy(v1,v2); }
  };
  template <class V1, class V2>
  struct CallCopyV<false,true,V1,V2> // inst = true
  {
    static inline void call(const V1& v1, V2& v2)
    { InstCopy(v1.XView(),v2.XView()); }
  };

  template <class V1, class V2> 
  inline void Copy(const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
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
        checkalias ) };
    CallCopyV<checkalias,inst,V1,V2>::call(v1.vec(),v2.vec());
  }

} // namespace tmv

#endif
