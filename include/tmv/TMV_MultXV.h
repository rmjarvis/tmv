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


#ifndef TMV_MultXV_H
#define TMV_MultXV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"
#include "TMV_CopyV.h"
#include "TMV_AddVV.h"

namespace tmv {

  //
  // Vector *= Scalar
  //

  // v *= x
  template <int algo, int size, int ix, class T, class V>
  struct MultXV_Helper;

  // algo 0: trivial: ix == 1, so nothing to do
  template <int size, class T, class V>
  struct MultXV_Helper<1,size,1,T,V>
  { static inline void call(const Scaling<1,T>& , V& ) {} };

  // algo 1: simple for loop
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<1,size,ix,T,V>
  {
    typedef typename V::iterator IT;
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      call2(n,x,v.begin());
    }
    static inline void call2(int n, const Scaling<ix,T>& x, IT it)
    { if (n) for(;n;--n) *it++ *= x; }
  };

  // algo 2: 2 at a time
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<2,size,ix,T,V>
  {
    typedef typename V::iterator IT;
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      call2(n,x,v.begin());
    }
    static inline void call2(const int n, const Scaling<ix,T>& x, IT it)
    {
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);

      if (n_2) do {
        it[0] *= x;
        it[1] *= x;
        it += 2;
      } while (--n_2);
      if (nb) *it *= x;
    }
  };

  // algo 3: 4 at a time
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<3,size,ix,T,V>
  {
    typedef typename V::iterator IT;
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      call2(n,x,v.begin());
    }
    static inline void call2(const int n, const Scaling<ix,T>& x, IT it)
    {
      int n_4 = (n>>2);
      int nb = n-(n_4<<2);

      if (n_4) do {
        it[0] *= x;
        it[1] *= x;
        it[2] *= x;
        it[3] *= x;
        it += 4;
      } while (--n_4);
      if (nb) do {
        *it++ *= x; 
      } while (--nb);
    }
  };

  // algo 5: fully unroll
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<5,size,ix,T,V> 
  {
    template <int I, int N, bool iscomplex>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, V& v)
      {
        Unroller<I,N/2,iscomplex>::unroll(x,v);
        Unroller<I+N/2,N-N/2,iscomplex>::unroll(x,v);
        //Unroller<I,N-1>::unroll(x,v);
        //Unroller<I+N-1,1>::unroll(x,v);
      }
    };
    template <int I>
    struct Unroller<I,1,false>
    {
      static inline void unroll(const Scaling<ix,T>& x, V& v)
      { v.ref(I) *= x; }
    };
    template <int I>
    struct Unroller<I,1,true>
    {
      static inline void unroll(const Scaling<ix,T>& x, V& v)
      { 
        typedef typename V::real_type RT;
        typedef typename V::value_type VT;
        const RT rv = ZProd<false,false>::rprod(x,v.cref(I));
        const RT iv = ZProd<false,false>::iprod(x,v.cref(I));
        v.ref(I) = VT(rv,iv);
      }
    };
    template <int I, bool iscomplex>
    struct Unroller<I,0,iscomplex>
    { static inline void unroll(const Scaling<ix,T>& , V& ) {} };

    static inline void call(const Scaling<ix,T>& x, V& v)
    { Unroller<0,size,V::viscomplex>::unroll(x,v); }
  };

  // algo 7: complex vector with unit step, convert to real version.
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<7,size,ix,T,V> 
  {
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      typedef typename V::flatten_type VF;
      typedef typename V::real_type RT;
      enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
      enum { algo2 = (
          (size2 != UNKNOWN && size2 <= (128/sizeof(RT))) ? 5 :
          (sizeof(RT) == 8) ? 2 : 
          (sizeof(RT) == 4) ? 3 :
          1 ) };
      VF vf = v.Flatten();
      MultXV_Helper<algo2,size2,ix,T,VF>::call(x,vf);
    }
  };

  // algo 8: complex vector, but not flatten-able
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<8,size,ix,T,V> 
  {
    typedef typename V::iterator IT;
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      call2(n,x,v.begin());
    }
    static inline void call2(int n, const Scaling<ix,T>& x, IT it)
    {
      typedef typename V::value_type VT;
      typedef typename V::real_type RT;
      VT val;
      RT rv,iv;
      enum { cx = false };
      enum { c1 = V::vconj };

      if (n) do {
        val = *it;
        rv = ZProd<cx,c1>::rprod(x,val);
        iv = ZProd<cx,c1>::iprod(x,val);
        *it++ = VT(rv,iv);
      } while (--n);
    }
  };

  // algo 9: complex vector, 2 at a time
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<9,size,ix,T,V> 
  {
    typedef typename V::iterator IT;
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      call2(n,x,v.begin());
    }
    static inline void call2(const int n, const Scaling<ix,T>& x, IT it)
    {
      typedef typename V::value_type VT;
      typedef typename V::real_type RT;
      VT vala, valb;
      RT rva,iva,rvb,ivb;

      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);
      enum { cx = false };
      enum { c1 = V::vconj };

      if (n_2) do {
        vala = it[0]; valb = it[1];
        rva = ZProd<cx,c1>::rprod(x,vala);
        iva = ZProd<cx,c1>::iprod(x,vala);
        rvb = ZProd<cx,c1>::rprod(x,valb);
        ivb = ZProd<cx,c1>::iprod(x,valb);
        it[0] = VT(rva,iva);
        it[1] = VT(rvb,ivb);
        it += 2;
      } while (--n_2);
      if (nb) {
        vala = *it;
        rva = ZProd<cx,c1>::rprod(x,vala);
        iva = ZProd<cx,c1>::iprod(x,vala);
        *it = VT(rva,iva);
      }
    }
  };

  // algo -1: Determine which algorithm to use
  template <int size, int ix, class T, class V>
  struct MultXV_Helper<-1,size,ix,T,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename V::iterator IT;
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      enum { algo = (
          (ix == 1) ? 0 :
#if TMV_OPT >= 1
          (Traits<T>::isreal && V::viscomplex && V::vstep == 1) ? 7 :
          (size != UNKNOWN && size <= (128/sizeof(VT))) ? 5 :
          (sizeof(RT) == 8 && V::vstep == 1) ? ( V::viscomplex ? 8 : 2 ) :
          (sizeof(RT) == 4 && V::vstep == 1) ? ( V::viscomplex ? 9 : 3 ) :
          (V::viscomplex && V::vstep == 1) ? 8 :
#endif
          1 ) };
      MultXV_Helper<algo,size,ix,T,V>::call(x,v);
    }
    static inline void call2(const int n, const Scaling<ix,T>& x, const IT& it)
    {
      enum { algo = (
#if TMV_OPT >= 1
          (sizeof(RT) == 8 && V::vstep == 1) ? ( V::viscomplex ? 8 : 2 ) :
          (sizeof(RT) == 4 && V::vstep == 1) ? ( V::viscomplex ? 9 : 3 ) :
          (V::viscomplex && V::vstep == 1) ? 8 :
#endif
          1 ) };
      MultXV_Helper<algo,size,ix,T,V>::call2(n,x,it);
    }
  };

  template <int ix, class T, class V>
  inline void InlineMultXV(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v)
  {
    enum { size = V::vsize };
    typedef typename V::view_type Vv;
    Vv vv = v.View();
    MultXV_Helper<-1,size,ix,T,Vv>::call(x,vv);
  }

  // Defined in TMV_MultXV.cpp
  template <class T>
  void InstMultXV(const T x, VectorView<T> v);  // BLAS scal

  // TODO: add a conj parameter to CallMultXV
  template <class T>
  inline void InstMultXV(const T x, VectorView<T,UNKNOWN,true> v)
  { InstMultXV(TMV_CONJ(x),v.Conjugate()); }

  //
  // Vector * / Scalar
  // -Vector
  //

  // v2 = x * v1
  // We just use AddVV_Helper with add = false
  template <int ix, class T, class V1, class V2>
  inline void InlineMultXV(const Scaling<ix,T>& x, 
      const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    enum { size = Sizes<V1::vsize,V2::vsize>::size };
    typedef typename V1::const_view_type V1v;
    typedef typename V2::view_type V2v;
    V1v v1v = v1.View();
    V2v v2v = v2.View();
    AddVV_Helper<-1,size,false,ix,T,V1v,V2v>::call(x,v1v,v2v);
  }

} // namespace tmv

#endif 
