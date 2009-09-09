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


#ifndef TMV_AddUU_H
#define TMV_AddUU_H

#include "TMV_AddVV.h"
#include "TMV_CopyU.h"

namespace tmv {

  //
  // U += U
  //

  template <int algo, int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper;

  // algo 0: trivial: ix == 1, !add, so call Copy
  template <int s, class T, class M1, class M2>
  struct AddUU_Helper<0,s,false,1,T,M1,M2>
  {
    static inline void call(const Scaling<1,T>& , const M1& m1, M2& m2)
    { CopyU_Helper<-1,s,M1,M2>::call(m1,m2); }
  };

  // algo 1: m1 is unitdiag
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<1,s,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M1::const_offdiag_type M1o;
      typedef typename M2::offdiag_type M2o;
      M1o m1o = m1.OffDiag();
      M2o m2o = m2.OffDiag();
      enum { sm1 = IntTraits2<s,-1>::sum };
      AddUU_Helper<-1,sm1,add,ix,T,M1o,M2o>::call(x,m1o,m2o);
      typedef typename M2::diag_type M2d;
      M2d m2d = m2.diag();
      Maybe<add>::addtoall(m2d,T(x));
    }
  };

  // algo 2: Loop over rows
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<2,s,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      int N = (s == UNKNOWN ? m2.size() : s);
      typedef typename M1::const_row_range_type M1r;
      typedef typename M2::row_range_type M2r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      typedef typename M2r::iterator IT2;
      const int step1 = m1.diagstep();
      const int step2 = m2.diagstep();
      IT1 it1 = m1.get_row(0,0,N).NonConj().begin();
      IT2 it2 = m2.get_row(0,0,N).begin();
      for(;N;--N) {
        AddVV_Helper<-1,UNKNOWN,add,ix,T,M1r,M2r>::call2(N,x,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 3: Loop over columns
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<3,s,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      int N = (s == UNKNOWN ? m2.size() : s);
      typedef typename M1::const_col_range_type M1c;
      typedef typename M2::col_range_type M2c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      typedef typename M2c::iterator IT2;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      IT1 it1 = m1.get_col(0,0,1).NonConj().begin();
      IT2 it2 = m2.get_col(0,0,1).begin();
      int M=1;
      for(;N;--N) {
        AddVV_Helper<-1,UNKNOWN,add,ix,T,M1c,M2c>::call2(M++,x,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 5: Fully unroll by rows
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<5,s,add,ix,T,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,N>::unroll(x,m1,m2);
        Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(x,m1,m2);
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
            ZProd<false,false>::prod(x,m1.cref(I,J)) );
      }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { static inline void unroll(const Scaling<ix,T>& , const M1& , M2& ) {} };

    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { Unroller<0,s,0,s>::unroll(x,m1,m2); }
  };

  // algo 6: Fully unroll by columns
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<6,s,add,ix,T,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, const M1& m1, M2& m2)
      {
        Unroller<I,M-(N-N/2),J,N/2>::unroll(x,m1,m2);
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
            ZProd<false,false>::prod(x,m1.cref(I,J)) );
      }
    };
    template <int I, int J>
    struct Unroller<I,0,J,1>
    { static inline void unroll(const Scaling<ix,T>& , const M1& , M2& ) {} };

    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { Unroller<0,s,0,s>::unroll(x,m1,m2); }
  };

  // algo 10: LowerTri, transpose:
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<10,s,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M1::const_transpose_type M1t;
      typedef typename M2::transpose_type M2t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      AddUU_Helper<-1,s,add,ix,T,M1t,M2t>::call(x,m1t,m2t);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int s, bool add, int ix, class T, class M1, class M2>
  struct AddUU_Helper<-1,s,add,ix,T,M1,M2>
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      TMVStaticAssert(!M2::munit || (ix == 1 && M1::munit && !add));
      typedef typename M2::value_type T2;
      enum { unit1 = M1::munit };
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
          !ShapeTraits<M2::mshape>::upper ? 10 :
          (ix == 1 && !add) ? 0 :
          unit1 ? 1 :
          unroll ? ( M2::mrowmajor ? 5 : 6 ) :
          M2::mrowmajor ? 2 : 3 ) };
      AddUU_Helper<algo,s,add,ix,T,M1,M2>::call(x,m1,m2);
    }
  };

  // m2 = x * m1
  template <int ix, class T, class M1, class M2>
  inline void InlineAddMM(
      const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
      BaseMatrix_Tri_Mutable<M2>& m2)
  {
    TMVStaticAssert(ShapeTraits<M1::mshape>::upper == 
        int(ShapeTraits<M2::mshape>::upper));
    TMVStaticAssert(!M2::munit || (ix == 1 && M1::munit));
    TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
    TMVAssert(m1.size() == m2.size());
    enum { s = Sizes<M1::msize,M2::msize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename M2::view_type M2v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    AddUU_Helper<-1,s,true,ix,T,M1v,M2v>::call(x,m1v,m2v);
  }

  // Defined in TMV_AddUU.cpp
  template <class T1, DiagType D1, bool C1, class T2>
  void InstAddMM(const T2 x,
      const ConstUpperTriMatrixView<T1,D1,UNKNOWN,UNKNOWN,C1>& m1, 
      UpperTriMatrixView<T2,NonUnitDiag> m2);
  template <class T1, DiagType D1, bool C1, class T2>
  void InstAddMM(const T2 x,
      const ConstLowerTriMatrixView<T1,D1,UNKNOWN,UNKNOWN,C1>& m1,
      LowerTriMatrixView<T2,NonUnitDiag> m2);


  //
  // U = x * U + x * U
  //

  template <int algo, int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2;

  // algo 1: m1 and/or m2 is unitdiag
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<1,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      typedef typename M1::const_offdiag_type M1o;
      typedef typename M2::const_offdiag_type M2o;
      typedef typename M3::offdiag_type M3o;
      M1o m1o = m1.OffDiag();
      M2o m2o = m2.OffDiag();
      M3o m3o = m3.OffDiag();
      enum { sm1 = IntTraits2<s,-1>::sum };
      AddUU_Helper2<-1,sm1,ix1,T1,M1o,ix2,T2,M2o,M3o>::call(x1,m1o,x2,m2o,m3o);
      // TODO: This could be made slightly more efficient with 
      // a Maybe<> function.
      if (m1.isunit()) m3.diag().SetAllTo(T1(x1));
      else m3.diag() = x1*m1.diag();
      if (m2.isunit()) m3.diag().AddToAll(T2(x2));
      else m3.diag() += x2*m2.diag();
    }
  };

  // algo 2: Loop over rows
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<2,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      int N = (s == UNKNOWN ? m2.size() : s);
      typedef typename M1::const_row_range_type M1r;
      typedef typename M2::const_row_range_type M2r;
      typedef typename M3::row_range_type M3r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      typedef typename M2r::const_nonconj_type::const_iterator IT2;
      typedef typename M3r::iterator IT3;
      const int step1 = m1.diagstep();
      const int step2 = m2.diagstep();
      const int step3 = m3.diagstep();
      IT1 it1 = m1.get_row(0,0,N).NonConj().begin();
      IT2 it2 = m2.get_row(0,0,N).NonConj().begin();
      IT3 it3 = m3.get_row(0,0,N).begin();
      for(;N;--N) {
        AddVV_Helper2<-1,UNKNOWN,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
            N,x1,it1,x2,it2,it3);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
        it3.ShiftP(step3);
      }
    }
  };

  // algo 3: Loop over columns
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<3,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      int N = (s == UNKNOWN ? m2.size() : s);
      typedef typename M1::const_col_range_type M1c;
      typedef typename M2::const_col_range_type M2c;
      typedef typename M3::col_range_type M3c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      typedef typename M2c::const_nonconj_type::const_iterator IT2;
      typedef typename M3c::iterator IT3;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      const int step3 = m3.stepj();
      IT1 it1 = m1.get_col(0,0,1).NonConj().begin();
      IT2 it2 = m2.get_col(0,0,1).NonConj().begin();
      IT3 it3 = m3.get_col(0,0,1).begin();
      int M=1;
      for(;N;--N) {
        AddVV_Helper2<-1,UNKNOWN,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
            M++,x1,it1,x2,it2,it3);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
        it3.ShiftP(step3);
      }
    }
  };

  // algo 5: Fully unroll by rows
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<5,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      {
        Unroller<I,M/2,J,N>::unroll(x1,m1,x2,m2,m3);
        Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(x1,m1,x2,m2,m3);
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
    { Unroller<0,s,0,s>::unroll(x1,m1,x2,m2,m3); }
  };

  // algo 6: Fully unroll by columns
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<6,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(
          const Scaling<ix1,T1>& x1, const M1& m1, 
          const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
      {
        Unroller<I,M-(N-N/2),J,N/2>::unroll(x1,m1,x2,m2,m3);
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
    { Unroller<0,s,0,s>::unroll(x1,m1,x2,m2,m3); }
  };

  // algo 10: LowerTri, transpose:
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<10,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      typedef typename M1::const_transpose_type M1t;
      typedef typename M2::const_transpose_type M2t;
      typedef typename M3::transpose_type M3t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      M3t m3t = m3.Transpose();
      AddUU_Helper2<-1,s,ix1,T1,M1t,ix2,T2,M2t,M3t>::call(x1,m1t,x2,m2t,m3t);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct AddUU_Helper2<-1,s,ix1,T1,M1,ix2,T2,M2,M3>
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1, 
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      typedef typename M3::value_type T3;
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
          !ShapeTraits<M3::mshape>::upper ? 10 :
          ( M1::munit || M2::munit ) ? 1 :
          unroll ? (
            ( (M1::mrowmajor && M2::mrowmajor) ||
              (M1::mrowmajor && M3::mrowmajor) ||
              (M2::mrowmajor && M3::mrowmajor) ) ? 5 : 6 ) :
          ( (M1::mrowmajor && M2::mrowmajor) ||
            (M1::mrowmajor && M3::mrowmajor) ||
            (M2::mrowmajor && M3::mrowmajor) ) ? 2 : 3 ) };
      AddUU_Helper2<algo,s,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
    }
  };

  // m3 = x1 * m1 + x2 * m2
  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void InlineAddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
      const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
      BaseMatrix_Tri_Mutable<M3>& m3)
  {
    TMVStaticAssert(ShapeTraits<M1::mshape>::upper == 
        int(ShapeTraits<M3::mshape>::upper));
    TMVStaticAssert(ShapeTraits<M2::mshape>::upper == 
        int(ShapeTraits<M3::mshape>::upper));
    TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
    TMVStaticAssert((Sizes<M1::msize,M3::msize>::same));
    TMVStaticAssert((Sizes<M1::msize,M3::msize>::same));
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.size() == m3.size());
    enum { s = Sizes<Sizes<M1::msize,M2::msize>::size,M3::msize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename M2::const_view_type M2v;
    typedef typename M3::view_type M3v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    M3v m3v = m3.View();
    AddUU_Helper2<-1,s,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(x1,m1v,x2,m2v,m3v);
  }

} // namespace tmv

#endif 
