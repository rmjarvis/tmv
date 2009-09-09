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

#ifndef TMV_CopyM_H
#define TMV_CopyM_H

#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

  //
  // Copy Matrices
  //

  template <int algo, int cs, int rs, class M1, class M2>
  struct CopyM_Helper;

  // algo 1: Linearize to vector version
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<1,cs,rs,M1,M2>
  {
    static inline void call(const M1& m1, M2& m2)
    {
      typedef typename M1::const_linearview_type M1l;
      typedef typename M2::linearview_type M2l;
      M1l m1l = m1.LinearView();
      M2l m2l = m2.LinearView();
      enum { cs_rs = IntTraits2<cs,rs>::prod };
      CopyV_Helper<-1,cs_rs,M1l,M2l>::call(m1l,m2l);
    }
  };

  // algo 2: Loop over rows
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<2,cs,rs,M1,M2>
  {
    static inline void call(const M1& m1, M2& m2)
    {
      int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::const_row_type M1r;
      typedef typename M2::row_type M2r;
      typedef typename M1r::const_iterator IT1;
      typedef typename M2r::iterator IT2;
      const int step1 = m1.stepi();
      const int step2 = m2.stepi();
      IT1 it1 = m1.get_row(0).begin();
      IT2 it2 = m2.get_row(0).begin();
      for(;M;--M) {
        CopyV_Helper<-1,rs,M1r,M2r>::call2(N,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 3: Loop over columns
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<3,cs,rs,M1,M2>
  {
    static inline void call(const M1& m1, M2& m2)
    {
      const int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::const_col_type M1c;
      typedef typename M2::col_type M2c;
      typedef typename M1c::const_iterator IT1;
      typedef typename M2c::iterator IT2;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      IT1 it1 = m1.get_col(0).begin();
      IT2 it2 = m2.get_col(0).begin();
      for(;N;--N) {
        CopyV_Helper<-1,cs,M1c,M2c>::call2(M,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<4,cs,rs,M1,M2>
  {
    static inline void call(const M1& m1, M2& m2)
    {
#if TMV_OPT >= 2
      if (m1.CanLinearize() && m2.CanLinearize() &&
          m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
        CopyM_Helper<1,cs,rs,M1,M2>::call(m1,m2);
      else if (m1.isrm() && m2.isrm() || 
          ( !(m1.iscm() && m2.iscm()) &&
            (m1.colsize() > m1.rowsize()) ) )
        CopyM_Helper<2,cs,rs,M1,M2>::call(m1,m2);
      else 
        CopyM_Helper<3,cs,rs,M1,M2>::call(m1,m2);
#else
      enum { algo2 = (
#if TMV_OPT >= 1
          ( M1::mrowmajor && M2::mrowmajor ) ? 2 :
          ( M1::mcolmajor && M2::mcolmajor ) ? 3 :
          ( cs == UNKNOWN || rs == UNKNOWN ) ? 3 :
          ( cs > rs ) ? 2 : 
#endif
          3 ) };
      CopyM_Helper<algo2,cs,rs,M1,M2>::call(m1,m2);
#endif
    }
  };

  // algo 5: Fully unroll by rows
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<5,cs,rs,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(const M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,N>::unroll(m1,m2);
        Unroller<I+M/2,M-M/2,J,N>::unroll(m1,m2);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline void unroll(const M1& m1, M2& m2)
      {
        Unroller<I,1,J,N/2>::unroll(m1,m2);
        Unroller<I,1,J+N/2,N-N/2>::unroll(m1,m2);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    { static inline void unroll(const M1& , M2& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(const M1& m1, M2& m2)
      { m2.ref(I,J) = m1.cref(I,J); }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { static inline void unroll(const M1& , M2& ) {} };

    static inline void call(const M1& m1, M2& m2)
    { Unroller<0,cs,0,rs>::unroll(m1,m2); }
  };

  // algo 6: Fully unroll by columns
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<6,cs,rs,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(const M1& m1, M2& m2)
      {
        Unroller<I,M,J,N/2>::unroll(m1,m2);
        Unroller<I,M,J+N/2,N-N/2>::unroll(m1,m2);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,1>
    {
      static inline void unroll(const M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,1>::unroll(m1,m2);
        Unroller<I+M/2,M-M/2,J,1>::unroll(m1,m2);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,0>
    { static inline void unroll(const M1& , M2& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(const M1& m1, M2& m2)
      { m2.ref(I,J) = m1.cref(I,J); }
    };
    template <int I, int J>
    struct Unroller<I,0,J,1>
    { static inline void unroll(const M1& , M2& ) {} };

    static inline void call(const M1& m1, M2& m2)
    { Unroller<0,cs,0,rs>::unroll(m1,m2); }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, class M1, class M2>
  struct CopyM_Helper<-1,cs,rs,M1,M2>
  {
    static inline void call(const M1& m1, M2& m2)
    {
      typedef typename M2::value_type T2;
      enum { canlin = (
          M1::mcanlin && M2::mcanlin &&
          ( (M1::mrowmajor && M2::mrowmajor) ||
            (M1::mcolmajor && M2::mcolmajor) ) ) };
      enum { algo = (
#if TMV_OPT >= 1
          canlin ? 1 :
          ( cs != UNKNOWN && rs != UNKNOWN ) ? (
            ( IntTraits2<cs,rs>::prod <= (128/sizeof(T2)) ) ? (
              ( M1::mrowmajor && M2::mrowmajor ) ? 5 : 6 ) :
            ( M1::mrowmajor && M2::mrowmajor ) ? 2 : 
            ( M1::mcolmajor && M2::mcolmajor ) ? 3 :
            ( cs > rs ) ? 2 : 3 ) :
#endif
          4 ) };
      //std::cout<<"algo = "<<algo<<std::endl;
      CopyM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
    }
  };

  template <class M1, class M2>
  inline void InlineCopy(
      const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
  {
    //std::cout<<"Inline Copy: m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same));
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    enum { cs = Sizes<M1::mcolsize,M2::mcolsize>::size };
    enum { rs = Sizes<M1::mrowsize,M2::mrowsize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename M2::view_type M2v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    CopyM_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    //std::cout<<"m2 => "<<m2.mat()<<std::endl;
  }

  // Defined in TMV_Matrix.cpp
  template <class T, bool C>
  void InstCopy(const ConstMatrixView<T,UNKNOWN,UNKNOWN,C>& m1,
      MatrixView<T> m2); // LAP lacpy

  template <bool checkalias, bool conj, bool inst, class M1, class M2>
  struct CallCopym  // checkalias = true
  {
    static inline void call(const M1& m1, M2& m2)
    {
      if (!SameStorage(m1,m2)) 
      {
        CallCopym<false,conj,inst,M1,M2>::call(m1,m2);
      }
      else
      {
        if (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) {
          if (M1::mconj != int(M2::mconj)) {
            m2.ConjugateSelf();
          }
          else {}  // They are already equal.
        } else {
          m2 = m1.copy();
        }
      }
    }
  };
  template <bool inst, class M1, class M2>
  struct CallCopym<false,true,inst,M1,M2> // conj = true
  {
    static inline void call(const M1& m1, M2& m2)
    {
      typedef typename M1::const_conjugate_type M1c;
      typedef typename M2::conjugate_type M2c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      CallCopym<false,false,inst,M1c,M2c>::call(m1c,m2c);
    }
  };
  template <class M1, class M2>
  struct CallCopym<false,false,true,M1,M2> // inst = true
  {
    static inline void call(const M1& m1, M2& m2)
    { InstCopy(m1.XView(),m2.XView()); }
  };
  template <class M1, class M2>
  struct CallCopym<false,false,false,M1,M2> // inst = false
  {
    static inline void call(const M1& m1, M2& m2)
    { InlineCopy(m1,m2); }
  };

  template <class M1, class M2> 
  inline void Copy(
      const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
  {
    typedef typename M1::value_type T1;
    typedef typename M2::value_type T2;
    enum { checkalias = (
        M1::mcolsize == UNKNOWN && 
        M2::mcolsize == UNKNOWN &&
        M1::mrowsize == UNKNOWN && 
        M2::mrowsize == UNKNOWN ) };
    enum { inst = (
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits2<T1,T2>::sametype &&
        (M1::mrowmajor || M1::mcolmajor) &&
        (M1::mrowmajor || M2::mcolmajor) &&
        checkalias ) };
    CallCopym<checkalias,M2::mconj,inst,M1,M2>::call(m1.mat(),m2.mat());
  }

} // namespace tmv

#endif
