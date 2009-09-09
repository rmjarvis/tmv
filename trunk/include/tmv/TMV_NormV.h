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

#ifndef TMV_NormV_H
#define TMV_NormV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

  // 
  // SumElements
  //

  template <int algo, int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper;

  // algo 1: simple for loop
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<1,size,comp,ix,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      ret sum(0);
      for(int i=0;i<n;++i) sum += Component<comp,VT>::f(x * v.cref(i));
      return sum;
    }
  };

  // algo 2: 2 at a time
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<2,size,comp,ix,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      ret sum0(0), sum1(0);
      VT v0, v1;
      typedef typename V::const_iterator IT;
      IT it = v.begin();
      int n_2 = (n>>1);
      const int nb = n-(n_2<<1);

      if (n_2) {
        do {
          v0 = x * *it; v1 = x*it[1]; it += 2;
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);
        } while (--n_2);
        sum0 += sum1;
      }
      if (nb) {
        v0 = x * *it;
        Component<comp,VT>::applyf(v0);
        sum0 += Component<comp,VT>::get(v0);
      }
      return sum0;
    }
  };

  // algo 3: 4 at a time
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<3,size,comp,ix,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      ret sum0(0), sum1(0);
      VT v0, v1;
      typedef typename V::const_iterator IT;
      IT it = v.begin();

      int n_4 = (n>>2);
      int nb = n-(n_4<<2);

      if (n_4) {
        do {
          v0 = x*it[0]; v1 = x*it[1];
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);

          v0 = x*it[2]; v1 = x*it[3]; it += 4;
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);
        } while (--n_4);
        sum0 += sum1;
      }
      if (nb) do {
        v0 = x*(*it++);
        Component<comp,VT>::applyf(v0);
        sum0 += Component<comp,VT>::get(v0);
      } while (--nb);
      return sum0;
    }
  };

  // algo 4: 8 at a time
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<4,size,comp,ix,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    {
      const int n = size == UNKNOWN ? int(v.size()) : size;
      ret sum0(0), sum1(0);
      VT v0, v1;
      typedef typename V::const_iterator IT;
      IT it = v.begin();

      int n_8 = (n>>3);
      int nb = n-(n_8<<3);

      if (n_8) {
        do {
          v0 = x*it[0]; v1 = x*it[1];
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);

          v0 = x*it[2]; v1 = x*it[3];
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);

          v0 = x*it[4]; v1 = x*it[5];
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);

          v0 = x*it[6]; v1 = x*it[7]; it += 8;
          Component<comp,VT>::applyf(v0);
          Component<comp,VT>::applyf(v1);
          sum0 += Component<comp,VT>::get(v0);
          sum1 += Component<comp,VT>::get(v1);
        } while (--n_8);
        sum0 += sum1;
      }
      if (nb) do {
        v0 = x*(*it++);
        Component<comp,VT>::applyf(v0);
        sum0 += Component<comp,VT>::get(v0);
      } while (--nb);
      return sum0;
    }
  };

  // algo 5: fully unroll
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<5,size,comp,ix,V> // known size, unroll
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    template <int I, int N>
    struct Unroller
    {
      static inline ret unroll(const V& v, const Scaling<ix,RT>& x)
      {
        return (
            Unroller<I,N/2>::unroll(v,x) +
            Unroller<I+N/2,N-N/2>::unroll(v,x));
      }
    };
    template <int I>
    struct Unroller<I,1>
    {
      static inline ret unroll(const V& v, const Scaling<ix,RT>& x)
      { return Component<comp,VT>::f(x * v.cref(I)); }
    };
    template <int I>
    struct Unroller<I,0>
    { 
      static inline ret unroll(const V& v, const Scaling<ix,RT>& x)
      { return ret(0); }
    };
    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    { return Unroller<0,size>::unroll(v,x); }
  };

  // algo 7: complex with unit step, convert to real
  // (This only makes sense for comp = NORM and ABS2.)
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<7,size,comp,ix,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    {
      typedef typename V::const_flatten_type Vf;
      Vf vf = v.Flatten();
      enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
      enum { algo2 = (
#if TMV_OPT >= 1
          (size2 != UNKNOWN && size2 <= 64) ? 5 :
          sizeof(RT) == 8 ? 3 :
          sizeof(RT) == 4 ? 4 :
#endif
          1 ) };
      return SumElementsV_Helper<algo2,size2,comp,ix,Vf>::call(vf,x);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int size, COMPType comp, int ix, class V>
  struct SumElementsV_Helper<-1,size,comp,ix,V> 
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<VT>::type ret;

    static inline ret call(const V& v, const Scaling<ix,RT>& x)
    {
      enum { maxunroll = (
          comp == VALUE ? 80 :
          comp == ABS_COMP ? 64 :
          comp == ABS2_COMP ? 64 :
          comp == NORM_COMP ? 64 : 
          /* no other options currently */ 0 ) };
      enum { algo = (
#if TMV_OPT >= 1
          ( ( comp == ABS2_COMP || comp == NORM_COMP ) && 
            ( V::viscomplex && V::vstep == 1) ) ? 7 :
          ( V::vsize != UNKNOWN && V::vsize <= int(maxunroll) ) ? 5 :
          (sizeof(RT) == 8 && V::vstep == 1) ? (V::viscomplex ? 2 : 3) :
          (sizeof(RT) == 4 && V::vstep == 1) ? (V::viscomplex ? 3 : 4) :
#endif
          1 ) };
      return SumElementsV_Helper<algo,size,comp,ix,V>::call(v,x);
    }
  };

  template <class V>
  inline typename V::value_type InlineSumElements(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename V::const_view_type Vv;
    Vv vv = v.View();
    return SumElementsV_Helper<-1,V::vsize,VALUE,1,Vv>::call(
        vv,Scaling<1,RT>());
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  T InstSumElements(const ConstVectorView<T>& v); 

  template <bool conj, bool inst, class V>
  struct CallSumElementsv // inst = false
  {
    static inline typename V::value_type call(const V& v)
    {
      typedef typename V::const_conjugate_type Vc;
      return TMV_CONJ(CallSumElementsv<false,inst,Vc>::call(v.Conjugate()));
    }
  };
  template <class V>
  struct CallSumElementsv<false,false,V>
  {
    static inline typename V::value_type call(const V& v)
    { return InlineSumElements(v); }
  };
  template <class V>
  struct CallSumElementsv<false,true,V>
  {
    static inline typename V::value_type call(const V& v)
    { return InstSumElements(v.XView()); }
  };

  template <class V>
  inline typename V::value_type SumElements(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    return CallSumElementsv<V::vconj,inst,V>::call(v.vec());
  }


  // 
  // SumAbsElements
  //

  // TODO: This (the real version) is one routine where the 
  // BLAS function (dasum) on my computer is significantly 
  // (~30%) faster than TMV.
  // I haven't figured out any way to make the above functions faster though.
  template <class V>
  inline typename V::real_type InlineSumAbsElements(
      const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename V::const_view_type Vv;
    Vv vv = v.View();
    return SumElementsV_Helper<-1,V::vsize,ABS_COMP,1,Vv>::call(
        vv,Scaling<1,RT>());
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  RealType(T) InstSumAbsElements(const ConstVectorView<T>& v);

  template <bool inst, class V>
  struct CallSumAbsElementsv // inst = false
  {
    static inline typename V::real_type call(const V& v)
    { return InlineSumAbsElements(v); }
  };
  template <class V>
  struct CallSumAbsElementsv<true,V>
  {
    static inline typename V::real_type call(const V& v)
    { return InstSumAbsElements(v.XView()); }
  };

  template <class V>
  inline typename V::real_type SumAbsElements(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type T;
    typedef typename V::const_nonconj_type Vn;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    return CallSumAbsElementsv<inst,Vn>::call(v.NonConj());
  }


  // 
  // SumAbs2Elements
  //

  template <class V>
  inline typename V::real_type InlineSumAbs2Elements(
      const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename V::const_view_type Vv;
    Vv vv = v.View();
    return SumElementsV_Helper<-1,V::vsize,ABS2_COMP,1,Vv>::call(
        vv,Scaling<1,RT>());
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  RealType(T) InstSumAbs2Elements(const ConstVectorView<T>& v); // BLAS asum

  template <bool inst, class V>
  struct CallSumAbs2Elementsv // inst = false
  {
    static inline typename V::real_type call(const V& v)
    { return InlineSumAbs2Elements(v); }
  };
  template <class V>
  struct CallSumAbs2Elementsv<true,V>
  {
    static inline typename V::real_type call(const V& v)
    { return InstSumAbs2Elements(v.XView()); }
  };

  template <class V>
  inline typename V::real_type SumAbs2Elements(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type T;
    typedef typename V::const_nonconj_type Vn;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    return CallSumAbs2Elementsv<inst,Vn>::call(v.NonConj());
  }


  // 
  // NormSq
  //

  template <class V>
  inline typename V::real_type InlineNormSq(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename V::const_view_type Vv;
    Vv vv = v.View();
    return SumElementsV_Helper<-1,V::vsize,NORM_COMP,1,Vv>::call(
        vv,Scaling<1,RT>());
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  RealType(T) InstNormSq(const ConstVectorView<T>& v); 

  template <bool inst, class V>
  struct CallNormSqv // inst = false
  {
    static inline typename V::real_type call(const V& v)
    { return InlineNormSq(v); }
  };
  template <class V>
  struct CallNormSqv<true,V>
  {
    static inline typename V::real_type call(const V& v)
    { return InstNormSq(v.XView()); }
  };

  template <class V>
  inline typename V::real_type NormSq(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type T;
    typedef typename V::const_nonconj_type Vn;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    return CallNormSqv<inst,Vn>::call(v.NonConj());
  }


  // 
  // NormSq with scaling
  //

  template <class V>
  inline typename V::real_type InlineNormSq(
      const BaseVector_Calc<V>& v, typename V::real_type scale)
  {
    typedef typename V::value_type VT;
    typedef typename V::real_type RT;
    typedef typename V::const_view_type Vv;
    Vv vv = v.View();
    return SumElementsV_Helper<-1,V::vsize,NORM_COMP,0,Vv>::call(
        vv,Scaling<0,RT>(scale));
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  RealType(T) InstNormSq(const ConstVectorView<T>& v, RealType(T) scale); 

  template <bool inst, class V>
  struct CallNormSq_scalev // inst = false
  {
    static inline typename V::real_type call(const V& v, 
        const typename V::real_type scale)
    { return InlineNormSq(v,scale); }
  };
  template <class V>
  struct CallNormSq_scalev<true,V>
  {
    static inline typename V::real_type call(const V& v,
        const typename V::real_type scale)
    { return InstNormSq(v.XView(),scale); }
  };

  template <class V>
  inline typename V::real_type NormSq(const BaseVector_Calc<V>& v,
      const typename V::real_type scale)
  {
    typedef typename V::value_type T;
    typedef typename V::const_nonconj_type Vn;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    return CallNormSq_scalev<inst,Vn>::call(v.NonConj(),scale);
  }


  //
  // Norm2
  //

  // This helper struct works for either Vector or Matrix "V"
  template <int algo, class V> struct Norm_Helper;

  // algo 1: simple: sqrt(NormSq(v))
  template <class V> 
  struct Norm_Helper<1,V>
  {
    typedef typename V::real_type RT;
    static inline RT call(const V& v)
    { return TMV_SQRT(InlineNormSq(v)); }
  };

  // algo 2: Robust algorithm with checks for overflow and underflow.
  // This one always calls MaxAbsElement and then NormSq.
  // This is inefficient if there are no problems.
  // Since no problems is the usual case, I switched to the below
  // version (algo 3) that calls NormSq first and then redoes it if there
  // are problems.
  template <class V> 
  struct Norm_Helper<2,V>
  {
    typedef typename V::real_type RT;
    static RT call(const V& v)
    {
      const RT eps = Epsilon<RT>();

      // Start with the maximum |v(i)|.  It will tell us how (and if)
      // we need to use a scaling for NormSq().
      RT vmax = v.MaxAbsElement();

      // If vmax = 0, then norm2 = 0:
      if (vmax == RT(0)) return RT(0);

      // If vmax^2 * eps = 0, but vmax != 0, then a naive NormSq()
      // will produce underflow rounding errors.  Find a better scaling.
      // eps is a pure power of 2, so no rounding errors from
      // rescaling by a power of eps.
      else if (vmax * vmax * eps == RT(0)) {
        const RT inveps = RT(1)/eps;
        RT scale = inveps;
        vmax *= scale;
        RT eps2 = eps*eps;
        while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
        return TMV_SQRT(v.NormSq(scale))/scale;
      }

      // If 1/vmax == 0, then vmax is already inf, so no hope of
      // making it more accurate.  (And need to check, since otherwise
      // the next section would cause an infinite loop.)
      else if (RT(1)/vmax == RT(0)) {
        return vmax;
      }

      // If 1/(vmax^2) == 0, then a naive NormSq() will produce 
      // overflow.  Find a better scaling.
      else if (RT(1)/(vmax*vmax) == RT(0)) {
        RT scale = eps;
        vmax *= scale;
        while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
        return TMV_SQRT(v.NormSq(scale))/scale;
      }

      // No problems with overflow or underflow.
      else return TMV_SQRT(v.NormSq());
    }
  };

  // algo 3: Robust algorithm with checks for overflow and underflow.
  // This version is slower if there is a problem, but since
  // there usually isn't a problem, it is generally faster.
  template <class V> 
  struct Norm_Helper<3,V>
  {
    typedef typename V::real_type RT;
    static RT call(const V& v)
    {
      const RT eps = Epsilon<RT>();
      const RT vnormsq = v.NormSq();

      if (vnormsq * eps == RT(0)) {
        // Possible underflow errors:

        // If vmax = 0, then norm2 = 0:
        RT vmax = v.MaxAbsElement();
        if (vmax == RT(0)) return RT(0);

        // If vmax^2 * eps = 0, but vmax != 0, then vnormsq has
        // underflow rounding errors.  Find a better scaling.
        // eps is a pure power of 2, so no rounding errors from
        // rescaling by a power of eps.
        else if (vmax * vmax * eps == RT(0)) {
          const RT inveps = RT(1)/eps;
          RT scale = inveps;
          vmax *= scale;
          RT eps2 = eps*eps;
          while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
          return TMV_SQRT(v.NormSq(scale))/scale;
        }

        else return TMV_SQRT(vnormsq);
      }

      else if (RT(1)/vnormsq == RT(0)) {
        // Possible overflow errors:

        // If 1/vmax == 0, then vmax is already inf, so no hope of
        // making it more accurate.  (And need to check, since otherwise
        // the next section would cause an infinite loop.)
        RT vmax = v.MaxAbsElement();
        if (RT(1)/vmax == RT(0)) {
          return vmax;
        }

        // If 1/(vmax^2) == 0, then vnormsq has overflow errors. 
        // Find a better scaling.
        else if (RT(1)/(vmax*vmax) == RT(0)) {
          RT scale = eps;
          vmax *= scale;
          while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
          return TMV_SQRT(v.NormSq(scale))/scale;
        }

        else return TMV_SQRT(vnormsq);
      }

      // No problems with overflow or underflow.
      else return TMV_SQRT(vnormsq);
    }
  };

  template <class V>
  inline typename V::real_type InlineNorm2(const BaseVector_Calc<V>& v)
  {
    typedef typename V::const_view_type Vv;
    Vv vv = v.View();
#if TMV_OPT == 0
    return Norm_Helper<1,Vv>::call(vv);
#else
    if (v.size() == 0) return typename V::real_type(0);
    else return Norm_Helper<3,Vv>::call(vv);
#endif
  }

  // Defined in TMV_Vector.cpp
  template <class T>
  RealType(T) InstNorm2(const ConstVectorView<T>& v); // BLAS nrm2

  template <bool inst, class V>
  struct CallNorm2v // inst = false
  {
    static inline typename V::real_type call(const V& v)
    { return InlineNorm2(v); }
  };
  template <class V>
  struct CallNorm2v<true,V>
  {
    static inline typename V::real_type call(const V& v)
    { return InstNorm2(v.XView()); }
  };

  template <class V>
  inline typename V::real_type Norm2(const BaseVector_Calc<V>& v)
  {
    typedef typename V::value_type T;
    typedef typename V::const_nonconj_type Vn;
    enum { inst = (
        Traits<T>::isinst &&
        V::vsize == UNKNOWN) };
    return CallNorm2v<inst,Vn>::call(v.NonConj());
  }

} // namespace tmv

#endif
