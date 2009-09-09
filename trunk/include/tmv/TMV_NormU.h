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

#ifndef TMV_NormU_H
#define TMV_NormU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_NormV.h"
#include "TMV_MinMax.h"

namespace tmv {

  // 
  // SumElements
  //

  template <int algo, int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper;

  // algo 1: m1 is unitdiag
  template <int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper<1,s,comp,ix,M1> 
  {
    typedef typename M1::value_type MT;
    typedef typename M1::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<MT>::type ret;

    static inline ret call(const M1& m, const Scaling<ix,RT>& x)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      typedef typename M1::const_offdiag_type Mo;
      Mo mo = m.OffDiag();
      enum { sm1 = IntTraits2<s,-1>::sum };
      return SumElementsU_Helper<-1,sm1,comp,ix,Mo>::call(mo,x) + x*RT(N);
    }
  };

  // algo 2: loop over rows
  template <int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper<2,s,comp,ix,M1> 
  {
    typedef typename M1::value_type MT;
    typedef typename M1::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<MT>::type ret;

    static inline ret call(const M1& m, const Scaling<ix,RT>& x)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      typedef typename M1::const_row_range_type Mr;
      ret sum(0);
      for(int i=0;i<N;++i) {
        Mr mr = m.get_row(i,i,N);
        sum += SumElementsV_Helper<-1,UNKNOWN,comp,ix,Mr>::call(mr,x);
      }
      return sum;
    }
  };

  // algo 3: loop over columns
  template <int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper<3,s,comp,ix,M1> 
  {
    typedef typename M1::value_type MT;
    typedef typename M1::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<MT>::type ret;

    static inline ret call(const M1& m, const Scaling<ix,RT>& x)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      typedef typename M1::const_col_range_type Mc;
      ret sum(0);
      for(int j=0;j<N;++j) {
        Mc mc = m.get_col(j,0,j+1);
        sum += SumElementsV_Helper<-1,UNKNOWN,comp,ix,Mc>::call(mc,x);
      }
      return sum;
    }
  };

  // algo 5: Fully unroll by rows
  template <int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper<5,s,comp,ix,M1>
  {
    typedef typename M1::value_type MT;
    typedef typename M1::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<MT>::type ret;

    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      {
        return (
            Unroller<I,M/2,J,N>::unroll(m,x) +
            Unroller<I+M/2,M-M/2,J,N>::unroll(m,x));
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      {
        return (
            Unroller<I,1,J,N/2>::unroll(m,x) +
            Unroller<I,1,J+N/2,N-N/2>::unroll(m,x));
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    { 
      static inline ret unroll(const M1& , const Scaling<ix,RT>& ) 
      { return ret(0); } 
    };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      { return Component<comp,MT>::f(x * m.cref(I,J)); }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { 
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      { return ret(0); }
    };
    static inline ret call(const M1& m, const Scaling<ix,RT>& x)
    { return Unroller<0,s,0,s>::unroll(m,x); }
  };

  // algo 6: Fully unroll by columns
  template <int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper<6,s,comp,ix,M1>
  {
    typedef typename M1::value_type MT;
    typedef typename M1::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<MT>::type ret;

    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      {
        return (
            Unroller<I,M-(N-N/2),J,N/2>::unroll(m,x) +
            Unroller<I,M,J+N/2,N-N/2>::unroll(m,x));
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,1>
    {
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      {
        return (
            Unroller<I,M/2,J,1>::unroll(m,x) +
            Unroller<I+M/2,M-M/2,J,1>::unroll(m,x));
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,0>
    { 
      static inline ret unroll(const M1& , const Scaling<ix,RT>& ) 
      { return ret(0); } 
    };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      { return Component<comp,MT>::f(x * m.cref(I,J)); }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { 
      static inline ret unroll(const M1& m, const Scaling<ix,RT>& x)
      { return ret(0); }
    };
    static inline ret call(const M1& m, const Scaling<ix,RT>& x)
    { return Unroller<0,s,0,s>::unroll(m,x); }
  };

  // algo -1: Determine which algorithm to use
  template <int s, COMPType comp, int ix, class M1>
  struct SumElementsU_Helper<-1,s,comp,ix,M1>
  {
    typedef typename M1::value_type MT;
    typedef typename M1::real_type RT;
    typedef typename Maybe<comp!=VALUE>::template RealType<MT>::type ret;
    static inline ret call(const M1& m, const Scaling<ix,RT>& x)
    {
      TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
      enum { unroll = (
          s == UNKNOWN ? false :
          // Norm is faster with the regular algorithm except for 
          // very small matrices.
          (s > 3 && comp == NORM_COMP) ? false :
#if TMV_OPT == 1
          s <= 3 ? true :
#elif TMV_OPT == 2
          s <= 5 ? true :
#elif TMV_OPT == 3
          s <= 10 ? true :
#endif
          false ) };
      enum { algo = (
          M1::munit ? 1 :
          unroll ? ( M1::mrowmajor ? 5 : 6 ) :
          M1::mcolmajor ? 3 : 2 ) };
      return SumElementsU_Helper<algo,s,comp,ix,M1>::call(m,x);
    }
  };

  template <class M>
  inline typename M::value_type InlineSumElements(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    typedef typename M::value_type MT;
    typedef typename M::real_type RT;
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return SumElementsU_Helper<-1,M::msize,VALUE,1,Mv>::call(
        mv,Scaling<1,RT>());
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  T InstSumElements(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  T InstSumElements(const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool lower, bool conj, bool inst, class M>
  struct CallSumElementsU;

  template <bool conj, bool inst, class M>
  struct CallSumElementsU<true,conj,inst,M> // lower = true
  {
    static inline typename M::value_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallSumElementsU<false,conj,inst,Mt>::call(mt);
    }
  };
  template <bool inst, class M>
  struct CallSumElementsU<false,true,inst,M> // conj = true
  {
    static inline typename M::value_type call(const M& m)
    {
      typedef typename M::const_conjugate_type Mc;
      return TMV_CONJ(
          CallSumElementsU<false,false,inst,Mc>::call(m.Conjugate()));
    }
  };
  template <class M>
  struct CallSumElementsU<false,false,true,M> // inst = true
  {
    static inline typename M::value_type call(const M& m)
    { return InstSumElements(m.XView()); }
  };
  template <class M>
  struct CallSumElementsU<false,false,false,M> // inst = false
  {
    static inline typename M::value_type call(const M& m)
    { return InlineSumElements(m); }
  };

  template <class M>
  inline typename M::value_type SumElements(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallSumElementsU<lower,M::mconj,inst,M>::call(m.mat());
  }


  // 
  // SumAbsElements
  //

  template <class M>
  inline typename M::real_type InlineSumAbsElements(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    typedef typename M::value_type MT;
    typedef typename M::real_type RT;
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return SumElementsU_Helper<-1,M::msize,ABS_COMP,1,Mv>::call(
        mv,Scaling<1,RT>());
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstSumAbsElements(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstSumAbsElements(
      const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool lower, bool inst, class M>
  struct CallSumAbsElementsU;

  template <bool inst, class M>
  struct CallSumAbsElementsU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallSumAbsElementsU<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallSumAbsElementsU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstSumAbsElements(m.XView()); }
  };
  template <class M>
  struct CallSumAbsElementsU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineSumAbsElements(m); }
  };

  template <class M>
  inline typename M::real_type SumAbsElements(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallSumAbsElementsU<lower,inst,Mn>::call(m.NonConj());
  }


  // 
  // SumAbs2Elements
  //

  template <class M>
  inline typename M::real_type InlineSumAbs2Elements(
      const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    typedef typename M::value_type MT;
    typedef typename M::real_type RT;
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return SumElementsU_Helper<-1,M::msize,ABS2_COMP,1,Mv>::call(
        mv,Scaling<1,RT>());
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstSumAbs2Elements(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstSumAbs2Elements(
      const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool lower, bool inst, class M>
  struct CallSumAbs2ElementsU;

  template <bool inst, class M>
  struct CallSumAbs2ElementsU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallSumAbs2ElementsU<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallSumAbs2ElementsU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstSumAbs2Elements(m.XView()); }
  };
  template <class M>
  struct CallSumAbs2ElementsU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineSumAbs2Elements(m); }
  };

  template <class M>
  inline typename M::real_type SumAbs2Elements(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallSumAbs2ElementsU<lower,inst,Mn>::call(m.NonConj());
  }


  // 
  // NormSq
  //

  template <class M>
  inline typename M::real_type InlineNormSq(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    typedef typename M::value_type MT;
    typedef typename M::real_type RT;
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return SumElementsU_Helper<-1,M::msize,NORM_COMP,1,Mv>::call(
        mv,Scaling<1,RT>());
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstNormSq(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstNormSq(const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool lower, bool inst, class M>
  struct CallNormSqU;

  template <bool inst, class M>
  struct CallNormSqU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallNormSqU<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallNormSqU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstNormSq(m.XView()); }
  };
  template <class M>
  struct CallNormSqU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineNormSq(m); }
  };

  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallNormSqU<lower,inst,Mn>::call(m.NonConj());
  }


  // 
  // NormSq with scaling
  //

  template <class M>
  inline typename M::real_type InlineNormSq(
      const BaseMatrix_Tri<M>& m, typename M::real_type scale)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    typedef typename M::value_type MT;
    typedef typename M::real_type RT;
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return SumElementsU_Helper<-1,M::msize,NORM_COMP,0,Mv>::call(
        mv,Scaling<0,RT>(scale));
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstNormSq(const ConstUpperTriMatrixView<T,D,1>& m, 
      RealType(T) scale); 
  template <class T, DiagType D>
  RealType(T) InstNormSq(const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m,
      RealType(T) scale); 

  template <bool lower, bool inst, class M>
  struct CallNormSq_scaleU;

  template <bool inst, class M>
  struct CallNormSq_scaleU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m,
        const typename M::real_type scale)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallNormSq_scaleU<false,inst,Mt>::call(mt,scale);
    }
  };
  template <class M>
  struct CallNormSq_scaleU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m,
        const typename M::real_type scale)
    { return InstNormSq(m.XView(),scale); }
  };
  template <class M>
  struct CallNormSq_scaleU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m, 
        const typename M::real_type scale)
    { return InlineNormSq(m,scale); }
  };

  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix_Tri<M>& m,
      const typename M::real_type scale)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallNormSq_scaleU<lower,inst,Mn>::call(m.NonConj(),scale);
  }



  //
  // NormF
  //

  // Norm_Helper in TMV_NormV.h works for UpperTriMatrix as well, so no
  // need to repeat that here.
  template <class M>
  inline typename M::real_type InlineNormF(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
#if TMV_OPT == 0
    return Norm_Helper<1,Mv>::call(mv);
#else
    if (m.size() == 0) return typename M::real_type(0);
    else return Norm_Helper<3,Mv>::call(mv);
#endif
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstNormF(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstNormF(const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool lower, bool inst, class M>
  struct CallNormFU;

  template <bool inst, class M>
  struct CallNormFU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallNormFU<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallNormFU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstNormF(m.XView()); }
  };
  template <class M>
  struct CallNormFU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineNormF(m); }
  };

  template <class M>
  inline typename M::real_type NormF(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallNormFU<lower,inst,Mn>::call(m.NonConj());
  }


  // 
  // MaxAbsElement
  //

  template <int algo, int s, class M1>
  struct MaxAbsElementU_Helper;

  // algo 1: m1 is unitdiag
  template <int s, class M1>
  struct MaxAbsElementU_Helper<1,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      typedef typename M1::const_offdiag_type Mo;
      Mo mo = m.OffDiag();
      enum { sm1 = IntTraits2<s,-1>::sum };
      enum { algo2 = ( M1::mrowmajor ? 2 : 3 ) };
      RT temp = MaxAbsElementU_Helper<algo2,sm1,Mo>::call(mo);
      return (temp > RT(1)) ? temp : RT(1);
    }
  };

  // algo 2: loop over rows
  template <int s, class M1>
  struct MaxAbsElementU_Helper<2,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      RT max(0);
      for(int i=0;i<N;++i) {
        RT temp = InlineMaxAbsElement(m.get_row(i,i,N),0);
        if (temp > max) max = temp;
      }
      return max;
    }
  };

  // algo 3: loop over columns
  template <int s, class M1>
  struct MaxAbsElementU_Helper<3,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      RT max(0);
      for(int j=0;j<N;++j) {
        RT temp = InlineMaxAbsElement(m.get_col(j,0,j+1),0);
        if (temp > max) max = temp;
      }
      return max;
    }
  };

  template <class M>
  inline typename M::real_type InlineMaxAbsElement(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    enum { s = M::msize };
    enum { algo = ( M::munit ? 1 : M::mrowmajor ? 2 : 3 ) };
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return MaxAbsElementU_Helper<algo,s,M>::call(mv);
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstMaxAbsElement(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstMaxAbsElement(
      const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m);

  template <bool lower, bool inst, class M>
  struct CallMaxAbsElementU;

  template <bool inst, class M>
  struct CallMaxAbsElementU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallMaxAbsElementU<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallMaxAbsElementU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstMaxAbsElement(m.XView()); }
  };
  template <class M>
  struct CallMaxAbsElementU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineMaxAbsElement(m); }
  };

  template <class M>
  inline typename M::real_type MaxAbsElement(const BaseMatrix_Tri<M>& m)
  {
    TMVAssert(m.size() > 0);
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallMaxAbsElementU<lower,inst,Mn>::call(m.NonConj());
  }


  // 
  // Norm1
  //

  // TODO: Norm1 and NormInf would benefit from an unroller.
  // Unlike with a regular matrix, where small matrices can unroll
  // each column (or row), with a triangle matrix, we lose the knowledge
  // of the length of each column in the for loop, so a full unroller
  // would be able to keep that.
  template <int algo, int s, class M1>
  struct Norm1U_Helper;

  // algo 1: loop over columns
  template <int s, class M1>
  struct Norm1U_Helper<1,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      RT max(0);
      for(int j=0;j<N;++j) {
        // If unit,     temp = 1 + SumAbsElements(m.col(j,0,j)
        // If non-unit, temp =     SumAbsElements(m.col(j,0,j+1)
        RT temp = Maybe<M1::munit>::sum( RT(1) , 
            InlineSumAbsElements(m.get_col(j,0,
                Maybe<M1::munit>::select(j,j+1))));
        if (temp > max) max = temp;
      }
      return max;
    }
  };

  // algo 2: loop over rows
  template <int s, class M1>
  struct Norm1U_Helper<2,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      int N = (s == UNKNOWN ? m.size() : s);
      if (N <= 8) return Norm1U_Helper<1,s,M1>::call(m);

      typedef typename M1::value_type MT;
      typedef typename tmv::Vector<RT>::iterator IT1;
      typedef typename M1::const_row_range_type::const_iterator IT2;

      MT value;

      // If unit,     start with all 1's.
      // If non-unit, start with all 0's.
      tmv::Vector<RT> temp(N,
          Maybe<M1::munit>::select(RT(1),RT(0)) );
      const IT1 begin1 = temp.begin();
      const IT1 end1 = temp.end();

      IT1 it1 = begin1;
      IT2 it2 = m.get_row(0,Maybe<M1::munit>::select(1,0),N).begin();
      // If unit, then we only need to add the offdiag to temp.
      // This is effected by: --N, ++it1, and the above select for it2.
      Maybe<M1::munit>::decrement(N);
      Maybe<M1::munit>::increment(it1);
      int end_step = m.diagstep() - N; // back to the start of next row

      do {
        do {
          value = *it2++;
          Component<ABS_COMP,MT>::applyf(value);
          *it1++ += TMV_REAL(value);
        } while (it1 != end1);
        it2 += end_step++;
        it1 -= (--N);
      } while (N);
      return InlineMaxElement(temp,0);
    }
  };

  template <class M>
  inline typename M::real_type InlineNorm1(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    enum { s = M::msize };
    enum { algo = (
#if TMV_OPT >= 1
        ( M::mrowmajor && (s == UNKNOWN || s > 8) ) ? 2 :
#endif
        1 ) };
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return Norm1U_Helper<algo,s,Mv>::call(mv);
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstNorm1(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstNorm1(const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool lower, bool inst, class M>
  struct CallNorm1U;
  template <bool lower, bool inst, class M>
  struct CallNormInfU;

  template <bool inst, class M>
  struct CallNorm1U<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallNormInfU<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallNorm1U<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstNorm1(m.XView()); }
  };
  template <class M>
  struct CallNorm1U<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineNorm1(m); }
  };

  template <class M>
  inline typename M::real_type Norm1(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallNorm1U<lower,inst,Mn>::call(m.NonConj());
  }

  // 
  // NormInf
  //

  template <int algo, int s, class M1>
  struct NormInfU_Helper;

  // algo 1: loop over rows
  template <int s, class M1>
  struct NormInfU_Helper<1,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      const int N = (s == UNKNOWN ? m.size() : s);
      RT max(0);
      for(int i=0;i<N;++i) {
        // If unit,     temp = 1 + SumAbsElements(m.row(i,0,i)
        // If non-unit, temp =     SumAbsElements(m.row(i,0,i+1)
        RT temp = Maybe<M1::munit>::sum( RT(1) , 
            InlineSumAbsElements(m.get_row(i,
                Maybe<M1::munit>::select(i+1,i),N)));
        if (temp > max) max = temp;
      }
      return max;
    }
  };

  // algo 2: loop over cols
  template <int s, class M1>
  struct NormInfU_Helper<2,s,M1> 
  {
    typedef typename M1::real_type RT;
    static inline RT call(const M1& m)
    {
      int N = (s == UNKNOWN ? m.size() : s);
      if (N <= 8) return NormInfU_Helper<1,s,M1>::call(m);

      typedef typename M1::value_type MT;
      typedef typename tmv::Vector<RT>::iterator IT1;
      typedef typename M1::const_col_range_type::const_iterator IT2;

      MT value;

      // If unit,     start with all 1's.
      // If non-unit, start with all 0's.
      tmv::Vector<RT> temp(N,
          Maybe<M1::munit>::select(RT(1),RT(0)) );
      const IT1 begin1 = temp.begin();

      IT1 it1 = begin1;
      IT2 it2 = m.get_col(Maybe<M1::munit>::select(1,0),0,1).begin();
      // If unit, then we only need to add the offdiag to temp.
      // This is effected by: --N, and the above select for it2.
      Maybe<M1::munit>::decrement(N);
      int end_step = m.stepj()-1; // back to the start of next column

      int M=1, i;
      do {
        i=M; do {
          value = *it2++;
          Component<ABS_COMP,MT>::applyf(value);
          *it1++ += TMV_REAL(value);
        } while (--i);
        it1 -= M++;
        it2 += end_step--;
      } while (--N);
      return InlineMaxElement(temp,0);
    }
  };

  template <class M>
  inline typename M::real_type InlineNormInf(const BaseMatrix_Tri<M>& m)
  {
    TMVStaticAssert(ShapeTraits<M::mshape>::upper);
    enum { s = M::msize };
    enum { algo = (
#if TMV_OPT >= 1
        ( M::mcolmajor && (s == UNKNOWN || s > 8) ) ? 2 :
#endif
        1 ) };
    typedef typename M::const_view_type Mv;
    Mv mv = m.View();
    return NormInfU_Helper<algo,s,Mv>::call(mv);
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D>
  RealType(T) InstNormInf(const ConstUpperTriMatrixView<T,D,1>& m); 
  template <class T, DiagType D>
  RealType(T) InstNormInf(const ConstUpperTriMatrixView<T,D,UNKNOWN,1>& m); 

  template <bool inst, class M>
  struct CallNormInfU<true,inst,M> // lower = true
  {
    static inline typename M::real_type call(const M& m)
    {
      typedef typename M::const_transpose_type Mt;
      Mt mt = m.Transpose();
      return CallNorm1U<false,inst,Mt>::call(mt);
    }
  };
  template <class M>
  struct CallNormInfU<false,true,M> // inst = true
  {
    static inline typename M::real_type call(const M& m)
    { return InstNormInf(m.XView()); }
  };
  template <class M>
  struct CallNormInfU<false,false,M> // inst = false
  {
    static inline typename M::real_type call(const M& m)
    { return InlineNormInf(m); }
  };

  template <class M>
  inline typename M::real_type NormInf(const BaseMatrix_Tri<M>& m)
  {
    typedef typename M::value_type T;
    typedef typename M::const_nonconj_type Mn;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::msize == UNKNOWN) };
    enum { lower = ShapeTraits<M::mshape>::lower };
    return CallNormInfU<lower,inst,Mn>::call(m.NonConj());
  }

} // namespace tmv

#endif
