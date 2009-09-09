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


#ifndef TMV_AddMM_H
#define TMV_AddMM_H

#include "TMV_AddVV.h"
#include "TMV_CopyM.h"

namespace tmv {

  //
  // Matrix += Matrix
  //

  template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper;

  // algo 0: trivial: ix == 1, !add, so call Copy
  template <int cs, int rs, class T, class M1, class M2>
  struct AddMM_Helper<0,cs,rs,false,1,T,M1,M2>
  {
    static inline void call(const Scaling<1,T>& , const M1& m1, M2& m2)
    { CopyM_Helper<-1,cs,rs,M1,M2>::call(m1,m2); }
  };

  // algo 1: Linearize to vector version
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<1,cs,rs,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M1::const_linearview_type M1l;
      typedef typename M2::linearview_type M2l;
      M1l m1l = m1.LinearView();
      M2l m2l = m2.LinearView();
      enum { prod = IntTraits2<cs,rs>::prod };
      AddVV_Helper<-1,prod,add,ix,T,M1l,M2l>::call(x,m1l,m2l);
    }
  };

  // algo 2: Loop over rows
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<2,cs,rs,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::const_row_type M1r;
      typedef typename M2::row_type M2r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      typedef typename M2r::iterator IT2;
      const int step1 = m1.stepi();
      const int step2 = m2.stepi();
      IT1 it1 = m1.get_row(0).NonConj().begin();
      IT2 it2 = m2.get_row(0).begin();
      for(;M;--M) {
        AddVV_Helper<-1,rs,add,ix,T,M1r,M2r>::call2(N,x,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 3: Loop over columns
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<3,cs,rs,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      const int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::const_col_type M1c;
      typedef typename M2::col_type M2c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      typedef typename M2c::iterator IT2;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      IT1 it1 = m1.get_col(0).NonConj().begin();
      IT2 it2 = m2.get_col(0).begin();
      for(;N;--N) {
        AddVV_Helper<-1,cs,add,ix,T,M1c,M2c>::call2(M,x,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<4,cs,rs,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
#if TMV_OPT >= 2
      if (m1.CanLinearize() && m2.CanLinearize() &&
          m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
        AddMM_Helper<1,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
      else if (m1.isrm() && m2.isrm() || 
          ( !(m1.iscm() && m2.iscm()) &&
            (m1.colsize() > m1.rowsize()) ) )
        AddMM_Helper<2,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
      else 
        AddMM_Helper<3,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
#else
      enum { algo2 = (
#if TMV_OPT >= 1
          ( M1::mrowmajor && M2::mrowmajor ) ? 2 :
          ( M1::mcolmajor && M2::mcolmajor ) ? 3 :
          ( cs == UNKNOWN || rs == UNKNOWN ) ? 3 :
          ( cs > rs ) ? 2 : 
#endif
          3 ) };
      AddMM_Helper<algo2,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
#endif
    }
  };

  // algo 5: Fully unroll by rows
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<5,cs,rs,add,ix,T,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,N>::unroll(x,m1,m2);
        Unroller<I+M/2,M-M/2,J,N>::unroll(x,m1,m2);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Unroller<I,1,J,N/2>::unroll(x,m1,m2);
        Unroller<I,1,J+N/2,N-N/2>::unroll(x,m1,m2);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    { static inline void unroll(const Scaling<ix,T>& , const M1& , M2& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Maybe<add>::add( m2.ref(I,J) , 
            ZProd<false,false>::prod(x , m1.cref(I,J)) ); 
      }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { static inline void unroll(const Scaling<ix,T>& , const M1& , M2& ) {} };

    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { Unroller<0,cs,0,rs>::unroll(x,m1,m2); }
  };

  // algo 6: Fully unroll by columns
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<6,cs,rs,add,ix,T,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Unroller<I,M,J,N/2>::unroll(x,m1,m2);
        Unroller<I,M,J+N/2,N-N/2>::unroll(x,m1,m2);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,1>
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,1>::unroll(x,m1,m2);
        Unroller<I+M/2,M-M/2,J,1>::unroll(x,m1,m2);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,0>
    { static inline void unroll(const Scaling<ix,T>& , const M1& , M2& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      { 
        Maybe<add>::add( m2.ref(I,J) , 
            ZProd<false,false>::prod(x , m1.cref(I,J)) ); 
      }
    };
    template <int I, int J>
    struct Unroller<I,0,J,1>
    { static inline void unroll(const Scaling<ix,T>& , const M1& , M2& ) {} };

    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { Unroller<0,cs,0,rs>::unroll(x,m1,m2); }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
  struct AddMM_Helper<-1,cs,rs,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M2::value_type T2;
      enum { canlin = (
          M1::mcanlin && M2::mcanlin &&
          ( (M1::mrowmajor && M2::mrowmajor) ||
            (M1::mcolmajor && M2::mcolmajor) ) ) };
      enum { algo = (
          ( ix == 1 && !add ) ? 0 :
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
      AddMM_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
    }
  };

  // m2 += x * m1
  template <int ix, class T, class M1, class M2>
  inline void InlineAddMM(
      const Scaling<ix,T>& x, const BaseMatrix_Rec<M1>& m1, 
      BaseMatrix_Rec_Mutable<M2>& m2)
  {
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
    AddMM_Helper<-1,cs,rs,true,ix,T,M1v,M2v>::call(x,m1v,m2v);
  }

  // Defined in TMV_AddMM.cpp
  template <class T1, bool C1, class T2>
  void InstAddMM(const T2 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, MatrixView<T2> m2);

  //
  // Matrix + Matrix
  //

  template <int algo, int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2;

  // algo 1: Linearize to vector version
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      typedef typename M1::const_linearview_type M1l;
      typedef typename M2::const_linearview_type M2l;
      typedef typename M3::linearview_type M3l;
      M1l m1l = m1.LinearView();
      M2l m2l = m2.LinearView();
      M3l m3l = m3.LinearView();
      InlineAddVV(x1,m1l,x2,m2l,m3l);
    }
  };

  // algo 2: Loop over rows
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::const_row_type M1r;
      typedef typename M2::const_row_type M2r;
      typedef typename M3::row_type M3r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      typedef typename M2r::const_nonconj_type::const_iterator IT2;
      typedef typename M3r::iterator IT3;
      const int step1 = m1.stepi();
      const int step2 = m2.stepi();
      const int step3 = m3.stepi();
      IT1 it1 = m1.get_row(0).NonConj().begin();
      IT2 it2 = m2.get_row(0).NonConj().begin();
      IT3 it3 = m3.get_row(0).begin();
      for(;M;--M) {
        AddVV_Helper2<-1,UNKNOWN,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
            N,x1,it1,x2,it2,it3);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
        it3.ShiftP(step3);
      }
    }
  };

  // algo 3: Loop over columns
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<3,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      const int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::const_col_type M1c;
      typedef typename M2::const_col_type M2c;
      typedef typename M3::col_type M3c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      typedef typename M2c::const_nonconj_type::const_iterator IT2;
      typedef typename M3c::iterator IT3;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      const int step3 = m3.stepj();
      IT1 it1 = m1.get_col(0).NonConj().begin();
      IT2 it2 = m2.get_col(0).NonConj().begin();
      IT3 it3 = m3.get_col(0).begin();
      for(;N;--N) {
        AddVV_Helper2<-1,UNKNOWN,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
            M,x1,it1,x2,it2,it3);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
        it3.ShiftP(step3);
      }
    }
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<4,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
#if TMV_OPT >= 2
      if (m1.CanLinearize() && m2.CanLinearize() && m3.CanLinearize() &&
          m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj() && 
          m1.stepi() == m3.stepi() && m1.stepj() == m3.stepj())
        AddMM_Helper2<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
      else if ( 
          ( (m1.isrm() && m2.isrm()) || (m1.isrm() && m3.isrm()) ||
            (m2.isrm() && m3.isrm()) ) ||
          ( !( (m1.iscm() && m2.iscm()) || (m1.iscm() && m3.iscm()) ||
               (m2.iscm() && m3.iscm()) ) && m1.colsize() > m1.rowsize()) )
        AddMM_Helper2<2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
      else 
        AddMM_Helper2<3,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
#else
      enum { algo2 = (
#if TMV_OPT >= 1
          ( (M1::mrowmajor && M2::mrowmajor) ||
            (M1::mrowmajor && M3::mrowmajor) ||
            (M2::mrowmajor && M3::mrowmajor) ) ? 2 :
#endif
          3 ) };
      AddMM_Helper2<algo2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
#endif
    }
  };

  // algo 5: Fully unroll by rows
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<5,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      {
        Unroller<I,M/2,J,N>::unroll(x1,m1,x2,m2,m3);
        Unroller<I+M/2,M-M/2,J,N>::unroll(x1,m1,x2,m2,m3);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      {
        Unroller<I,1,J,N/2>::unroll(x1,m1,x2,m2,m3);
        Unroller<I,1,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
    };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
    };
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    { Unroller<0,cs,0,rs>::unroll(x1,m1,x2,m2,m3); }
  };

  // algo 6: Fully unroll by columns
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<6,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      {
        Unroller<I,M,J,N/2>::unroll(x1,m1,x2,m2,m3);
        Unroller<I,M,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,1>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      {
        Unroller<I,M/2,J,1>::unroll(x1,m1,x2,m2,m3);
        Unroller<I+M/2,M-M/2,J,1>::unroll(x1,m1,x2,m2,m3);
      }
    };
    template <int I, int M, int J>
    struct Unroller<I,M,J,0>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
    };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
    };
    template <int I, int J>
    struct Unroller<I,0,J,1>
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
    };
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    { Unroller<0,cs,0,rs>::unroll(x1,m1,x2,m2,m3); }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddMM_Helper2<-1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      typedef typename M3::value_type T3;
      enum { canlin = (
          M1::mcanlin && M2::mcanlin && M3::mcanlin &&
          ( (M1::mrowmajor && M2::mrowmajor && M3::mrowmajor) ||
            (M1::mcolmajor && M2::mcolmajor && M3::mcolmajor) ) ) };
      enum { algo = (
#if TMV_OPT >= 1
          canlin ? 1 :
          ( cs != UNKNOWN && rs != UNKNOWN ) ? (
            ( IntTraits2<cs,rs>::prod <= (128/sizeof(T2)) ) ? (
              ( (M1::mrowmajor && M2::mrowmajor) ||
                (M1::mrowmajor && M3::mrowmajor) ||
                (M2::mrowmajor && M3::mrowmajor) ) ? 5 : 6 ) :
            ( (M1::mrowmajor && M2::mrowmajor) ||
              (M1::mrowmajor && M3::mrowmajor) ||
              (M2::mrowmajor && M3::mrowmajor) ) ? 2 :
            ( (M1::mcolmajor && M2::mcolmajor) ||
              (M1::mcolmajor && M3::mcolmajor) ||
              (M2::mcolmajor && M3::mcolmajor) ) ? 3 :
            ( cs > rs ) ? 2 : 3 ) :
#endif
          4 ) };
      AddMM_Helper2<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
    }
  };

  // m3 = x1 * m1 + x2 * m2
  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void InlineAddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
      const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
      BaseMatrix_Rec_Mutable<M3>& m3)
  {
    //std::cout<<"AddMM: x1 = "<<ix1<<"  "<<T1(x1)<<"  x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3<<std::endl;
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M3::mrowsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M3::mrowsize>::same));
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.colsize() == m3.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.rowsize() == m3.rowsize());
    enum { cs = Sizes<Sizes<M1::mcolsize,M2::mcolsize>::size,M3::mcolsize>::size };
    enum { rs = Sizes<Sizes<M1::mrowsize,M2::mrowsize>::size,M3::mrowsize>::size };
    //std::cout<<"colsize = "<<colsize<<std::endl;
    //std::cout<<"rowsize = "<<rowsize<<std::endl;
    //std::cout<<"rowmajor = "<<rowmajor<<std::endl;
    typedef typename M1::const_view_type M1v;
    typedef typename M2::const_view_type M2v;
    typedef typename M3::view_type M3v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    M3v m3v = m3.View();
    AddMM_Helper2<-1,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(x1,m1v,x2,m2v,m3v);
    //std::cout<<"m3 => "<<m3<<std::endl;
  }

} // namespace tmv

#endif 
