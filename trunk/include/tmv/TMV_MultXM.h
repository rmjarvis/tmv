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
// This file defines the basic composite type for a product of a 
// matrix and a scalar.  It also implements the calculation for 
// dense rectangular matrices.

#ifndef TMV_MultXM_H
#define TMV_MultXM_H

#include "TMV_MultXV.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_AddMM.h"

namespace tmv {

  //
  // Matrix *= Scalar 
  //

  template <int algo, int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper;

  // algo 0: trivial: ix == 1, so nothing to do
  template <int cs, int rs, class T, class M1>
  struct MultXM_Helper<1,cs,rs,1,T,M1>
  { static inline void call(const Scaling<1,T>& x, M1& m) {} };

  // algo 1: Linearize to vector version
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<1,cs,rs,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      typedef typename M1::linearview_type Ml;
      Ml ml = m.LinearView();
      MultXV_Helper<-1,Ml::vsize,ix,T,Ml>::call(x,ml);
    }
  };

  // algo 2: Loop over rows
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<2,cs,rs,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      int M = cs == UNKNOWN ? int(m.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      typedef typename M1::row_type Mr;
      typedef typename Mr::iterator IT;
      const int step = m.stepi();
      IT it = m.get_row(0).begin();
      for(;M;--M) {
        MultXV_Helper<-1,rs,ix,T,Mr>::call2(N,x,it);
        it.ShiftP(step);
      }
    }
  };

  // algo 3: Loop over columns
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<3,cs,rs,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      const int M = cs == UNKNOWN ? int(m.colsize()) : cs;
      int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      typedef typename M1::col_type Mc;
      typedef typename Mc::iterator IT;
      const int step = m.stepj();
      IT it = m.get_col(0).begin();
      for(;N;--N) {
        MultXV_Helper<-1,cs,ix,T,Mc>::call2(M,x,it);
        it.ShiftP(step);
      }
    }
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<4,cs,rs,ix,T,M1>
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
#if TMV_OPT >= 2
      if (m.CanLinearize()) 
        MultXM_Helper<1,cs,rs,ix,T,M1>::call(x,m);
      else if ( m.isrm() ||
          (!m.iscm() && (m.colsize() > m.rowsize()) ) )
        MultXM_Helper<2,cs,rs,ix,T,M1>::call(x,m);
      else
        MultXM_Helper<3,cs,rs,ix,T,M1>::call(x,m);
#else
      enum { algo2 = (
#if TMV_OPT >= 1
          M1::mrowmajor ? 2 :
          M1::mcolmajor ? 3 :
          ( cs == UNKNOWN || rs == UNKNOWN ) ? 3 :
          ( cs > rs ) ? 2 :
#endif
          3 ) };
      MultXM_Helper<algo2,cs,rs,ix,T,M1>::call(x,m);
#endif
    }
  };

  // algo 5: Fully unroll by rows
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<5,cs,rs,ix,T,M1>
  {
    template <int I, int M, int J, int N, bool iscomplex>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        Unroller<I,M/2,J,N,iscomplex>::unroll(x,m);
        Unroller<I+M/2,M-M/2,J,N,iscomplex>::unroll(x,m);
      }
    };
    template <int I, int J, int N, bool iscomplex>
    struct Unroller<I,1,J,N,iscomplex>
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        Unroller<I,1,J,N/2,iscomplex>::unroll(x,m);
        Unroller<I,1,J+N/2,N-N/2,iscomplex>::unroll(x,m);
      }
    };
    template <int I, int J, int N, bool iscomplex>
    struct Unroller<I,0,J,N,iscomplex>
    { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1,false>
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      { m.ref(I,J) *= x; }
    };
    template <int I, int J>
    struct Unroller<I,1,J,1,true>
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        typedef typename M1::real_type RT;
        typedef typename M1::value_type VT;
        const RT rm = ZProd<false,M1::mconj>::rprod(x,m.NonConj().cref(I,J));
        const RT im = ZProd<false,M1::mconj>::iprod(x,m.NonConj().cref(I,J));
        m.ref(I,J) = VT(rm,im);
      }
    };
    template <int I, int J, bool iscomplex>
    struct Unroller<I,1,J,0,iscomplex>
    { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };

    static inline void call(const Scaling<ix,T>& x, M1& m)
    { Unroller<0,cs,0,rs,M1::miscomplex>::unroll(x,m); }
  };

  // algo 6: Fully unroll by columns
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<6,cs,rs,ix,T,M1>
  {
    template <int I, int M, int J, int N, bool iscomplex>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        Unroller<I,M,J,N/2,iscomplex>::unroll(x,m);
        Unroller<I,M,J+N/2,N-N/2,iscomplex>::unroll(x,m);
      }
    };
    template <int I, int M, int J, bool iscomplex>
    struct Unroller<I,M,J,1,iscomplex>
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        Unroller<I,M/2,J,1,iscomplex>::unroll(x,m);
        Unroller<I+M/2,M-M/2,J,1,iscomplex>::unroll(x,m);
      }
    };
    template <int I, int M, int J, bool iscomplex>
    struct Unroller<I,M,J,0,iscomplex>
    { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1,false>
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      { m.ref(I,J) *= x; }
    };
    template <int I, int J>
    struct Unroller<I,1,J,1,true>
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        typedef typename M1::real_type RT;
        typedef typename M1::value_type VT;
        const RT rm = ZProd<false,M1::mconj>::rprod(x,m.NonConj().cref(I,J));
        const RT im = ZProd<false,M1::mconj>::iprod(x,m.NonConj().cref(I,J));
        m.ref(I,J) = VT(rm,im);
      }
    };
    template <int I, int J, bool iscomplex>
    struct Unroller<I,0,J,1,iscomplex>
    { static inline void unroll(const Scaling<ix,T>& , M1& ) {} };

    static inline void call(const Scaling<ix,T>& x, M1& m)
    { Unroller<0,cs,0,rs,M1::miscomplex>::unroll(x,m); }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, int ix, class T, class M1>
  struct MultXM_Helper<-1,cs,rs,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      typedef typename M1::value_type T1;
      enum { algo = (
          (ix == 1) ? 0 :
#if TMV_OPT >= 1
          M1::mcanlin ? 1 :
          ( cs != UNKNOWN && rs != UNKNOWN ) ? (
            ( IntTraits2<cs,rs>::prod <= (128/sizeof(T1)) ) ? (
              ( M1::mrowmajor ? 5 : 6 ) ) :
            M1::mrowmajor ? 2 :
            M1::mcolmajor ? 3 :
            ( cs > rs ) ? 2 : 3 ) :
#endif
          4 ) };
      //std::cout<<"InlineMultXM: x = "<<T(x)<<std::endl;
      //std::cout<<"m = "<<TypeText(m)<<"  "<<m<<std::endl;
      //std::cout<<"algo = "<<algo<<std::endl;
      MultXM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
      //std::cout<<"m => "<<m<<std::endl;
    }
  };

  // m *= x
  template <int ix, class T, class M>
  inline void InlineMultXM(const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m)
  { MultXM_Helper<-1,M::mcolsize,M::mrowsize,ix,T,M>::call(x,m.mat()); }

  // Defined in TMV_MultXM.cpp
  template <class T>
  void InstMultXM(const T x, MatrixView<T> m);


  //
  // Matrix = Matrix * Scalar
  //

  // m2 = x * m1
  template <int ix, class T, class M1, class M2>
  inline void InlineMultXM(
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
    AddMM_Helper<-1,cs,rs,false,ix,T,M1v,M2v>::call(x,m1v,m2v);
  }

  // Defined in TMV_MultXM.cpp
  template <class T1, bool C1, class T2>
  void InstMultXM(const T2 x, const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1, 
      MatrixView<T2> m2);

} // namespace tmv

#endif
