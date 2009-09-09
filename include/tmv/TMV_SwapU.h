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

//
// This file contains the code for:
//
// DoSwap(m1,m2)
// m.TransposeSelf()
// m.PermuteRows(P)
// m.ReversePermuteRows(P)

#ifndef TMV_SwapU_H
#define TMV_SwapU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultXM.h"

namespace tmv {

  //
  // Swap Matrices
  //

  template <int algo, int s, class M1, class M2>
  struct SwapU_Helper;

  // algo 1: m1 is unitdiag
  template <int s, class M1, class M2>
  struct SwapU_Helper<1,s,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      typedef typename M1::offdiag_type M1o;
      typedef typename M2::offdiag_type M2o;
      M1o m1o = m1.OffDiag();
      M2o m2o = m2.OffDiag();
      enum { sm1 = IntTraits2<s,-1>::sum };
      SwapU_Helper<-1,sm1,M1o,M2o>::call(m1o,m2o);
    }
  };

  // algo 2: Loop over rows
  template <int s, class M1, class M2>
  struct SwapU_Helper<2,s,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      int N = (s == UNKNOWN ? m2.size() : s);
      typedef typename M1::row_range_type M1r;
      typedef typename M2::row_range_type M2r;
      typedef typename M1r::iterator IT1;
      typedef typename M2r::iterator IT2;
      const int step1 = m1.diagstep();
      const int step2 = m2.diagstep();
      IT1 it1 = m1.get_row(0,0,N).begin();
      IT2 it2 = m2.get_row(0,0,N).begin();
      for(;N;--N) {
        SwapV_Helper<-1,UNKNOWN,M1r,M2r>::call2(N,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      } 
    }
  };

  // algo 3: Loop over columns
  template <int s, class M1, class M2>
  struct SwapU_Helper<3,s,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      const int N = (s == UNKNOWN ? m2.size() : s);
      typedef typename M1::col_range_type M1c;
      typedef typename M2::col_range_type M2c;
      typedef typename M1c::iterator IT1;
      typedef typename M2c::iterator IT2;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      IT1 it1 = m1.get_col(0,0,1).begin();
      IT2 it2 = m2.get_col(0,0,1).begin();
      for(int j=0;j<N;++j) {
        SwapV_Helper<-1,UNKNOWN,M1c,M2c>::call2(j+1,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      } 
    }
  };

  // algo 5: Fully unroll by rows
  template <int s, class M1, class M2>
  struct SwapU_Helper<5,s,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,N>::unroll(m1,m2);
        Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(m1,m2);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline void unroll(M1& m1, M2& m2)
      {
        Unroller<I,1,J,N/2>::unroll(m1,m2);
        Unroller<I,1,J+N/2,N-N/2>::unroll(m1,m2);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    { static inline void unroll(M1& , M2& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(M1& m1, M2& m2)
      { TMV_SWAP(m2.ref(I,J) , m1.ref(I,J) ); }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { static inline void unroll(M1& , M2& ) {} };

    static inline void call(M1& m1, M2& m2)
    { Unroller<0,s,0,s>::unroll(m1,m2); }
  };

  // algo 6: Fully unroll by columns
  template <int s, class M1, class M2>
  struct SwapU_Helper<6,s,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(M1& m1, M2& m2)
      {
        Unroller<I,M-(N-N/2),J,N/2>::unroll(m1,m2);
        Unroller<I,M,J+N/2,N-N/2>::unroll(m1,m2);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,1>
    {
      static inline void unroll(M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,1>::unroll(m1,m2);
        Unroller<I+M/2,M-M/2,J,1>::unroll(m1,m2);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,0>
    { static inline void unroll(M1& , M2& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(M1& m1, M2& m2)
      { TMV_SWAP(m2.ref(I,J) , m1.ref(I,J)); }
    };
    template <int I, int J>
    struct Unroller<I,0,J,1>
    { static inline void unroll(M1& , M2& ) {} };

    static inline void call(M1& m1, M2& m2)
    { Unroller<0,s,0,s>::unroll(m1,m2); }
  };

  // algo -1: Determine which algorithm to use
  template <int s, class M1, class M2>
  struct SwapU_Helper<-1,s,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      TMVStaticAssert(M1::munit == int(M2::munit));
      typedef typename M2::value_type T2;
      enum { unit = M1::munit };
      enum { unroll = (
          s == UNKNOWN ? false :
#if TMV_OPT == 1
          s <= 3 ? true :
#elif TMV_OPT == 2
          s <= 5 ? true :
#elif TMV_OPT == 3
          s <= 10 ? true :
#endif
          false ) };
      enum { algo = (
          unit ? 1 :
          unroll ? ( M2::mrowmajor ? 5 : 6 ) :
          M2::mcolmajor ? 3 : 2 ) };
      //std::cout<<"algo = "<<algo<<std::endl;
      SwapU_Helper<algo,s,M1,M2>::call(m1,m2);
    }
  };

  template <class M1, class M2>
  inline void InlineSwap(
      BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
  {
    //std::cout<<"Inline Swap: m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    TMVStaticAssert(ShapeTraits<M1::mshape>::upper);
    TMVStaticAssert(ShapeTraits<M2::mshape>::upper);
    TMVStaticAssert(M1::munit == int(M2::munit));
    TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
    TMVAssert(m1.size() == m2.size());
    enum { s = Sizes<M1::msize,M2::msize>::size };
    typedef typename M1::view_type M1v;
    typedef typename M2::view_type M2v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    SwapU_Helper<-1,s,M1v,M2v>::call(m1v,m2v);
    //std::cout<<"m1 => "<<m1.mat()<<std::endl;
    //std::cout<<"m2 => "<<m2.mat()<<std::endl;
  }

  // Defined in TMV_TriMatrix.cpp
  template <class T, DiagType D, bool C1>
  void InstSwap(UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C1> m1,
      UpperTriMatrixView<T,D> m2); 

  template <bool checkalias, bool lower, bool conj, bool inst, class M1, class M2>
  struct CallSwapU  // checkalias = true
  {
    static inline void call(M1& m1, M2& m2)
    {
      if (!SameStorage(m1,m2)) 
      {
        CallSwapU<false,lower,conj,inst,M1,M2>::call(m1,m2);
      }
      else
      {
        if (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
        {
          if (M1::mconj != int(M2::mconj)) 
          {
            m1.ConjugateSelf();
            m2.ConjugateSelf();
          }
          else {}  // They are equal.
        } 
        else // Have to swap with a full temporary
        { 
          typename M1::copy_type m1c = m1.copy();
          m1 = m2;
          m2 = m1c;
        }
      }
    }
  };
  template <bool conj, bool inst, class M1, class M2>
  struct CallSwapU<false,true,conj,inst,M1,M2> // lower = true
  {
    static inline void call(M1& m1, M2& m2)
    { 
      typedef typename M1::transpose_type M1t;
      typedef typename M2::transpose_type M2t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      CallSwapU<false,false,conj,inst,M1t,M2t>::call(m1t,m2t);
    }
  };
  template <bool inst, class M1, class M2>
  struct CallSwapU<false,false,true,inst,M1,M2> // conj = true
  {
    static inline void call(M1& m1, M2& m2)
    { 
      typedef typename M1::conjugate_type M1c;
      typedef typename M2::conjugate_type M2c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      CallSwapU<false,false,false,inst,M1c,M2c>::call(m1c,m2c);
    }
  };
  template <class M1, class M2>
  struct CallSwapU<false,false,false,true,M1,M2> // inst = true
  {
    static inline void call(M1& m1, M2& m2)
    { InstSwap(m1.XView(),m2.XView()); }
  };
  template <class M1, class M2>
  struct CallSwapU<false,false,false,false,M1,M2> // inst = false
  {
    static inline void call(M1& m1, M2& m2)
    { InlineSwap(m1,m2); }
  };

  template <class M1, class M2> 
  inline void DoSwap(
      BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
  {
    TMVStaticAssert(M1::mshape == int(M2::mshape));
    TMVStaticAssert(M1::munit == int(M2::munit));
    typedef typename M1::value_type T1;
    typedef typename M2::value_type T2;
    enum { checkalias = (
        M1::msize == UNKNOWN && 
        M2::msize == UNKNOWN ) };
    enum { inst = (
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits2<T1,T2>::sametype &&
        (M1::mcolmajor || M1::mrowmajor) &&
        (M2::mcolmajor || M2::mrowmajor) &&
        checkalias ) };
    enum { lower = ShapeTraits<M2::mshape>::lower };
    CallSwapU<checkalias,lower,M2::mconj,inst,M1,M2>::call(m1.mat(),m2.mat());
  }
  template <class M1, class M2> 
  inline void Swap(
      BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
  { DoSwap(m1,m2); }

} // namespace tmv

#endif
