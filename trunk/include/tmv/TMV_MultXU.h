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

#ifndef TMV_MultXU_H
#define TMV_MultXU_H

#include "TMV_MultXV.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_AddUU.h"

namespace tmv {

  //
  // Matrix *= Scalar 
  //

  template <int algo, int s, int ix, class T, class M1>
  struct MultXU_Helper;

  // algo 0: trivial: ix == 1, so nothing to do
  template <int s, class T, class M1>
  struct MultXU_Helper<0,s,1,T,M1>
  { static inline void call(const Scaling<1,T>& , M1& ) {} };

  // algo 2: Loop over rows
  template <int s, int ix, class T, class M1>
  struct MultXU_Helper<2,s,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      int N = (s == UNKNOWN ? m.size() : s);
      typedef typename M1::row_range_type Mr;
      typedef typename Mr::iterator IT;
      const int step = m.diagstep();
      IT it = m.get_row(0,0,N).begin();
      for(;N;--N) {
        MultXV_Helper<-1,UNKNOWN,ix,T,Mr>::call2(N,x,it);
        it.ShiftP(step);
      }
    }
  };

  // algo 3: Loop over columns
  template <int s, int ix, class T, class M1>
  struct MultXU_Helper<3,s,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      int N = (s == UNKNOWN ? m.size() : s);
      typedef typename M1::col_range_type Mc;
      typedef typename Mc::iterator IT;
      const int step = m.stepj();
      IT it = m.get_col(0,0,1).begin();
      int M=1;
      for(;N;--N) {
        MultXV_Helper<-1,UNKNOWN,ix,T,Mc>::call2(M++,x,it);
        it.ShiftP(step);
      }
    }
  };

  // algo 5: Fully unroll by rows
  template <int s, int ix, class T, class M1>
  struct MultXU_Helper<5,s,ix,T,M1>
  {
    template <int I, int M, int J, int N, bool iscomplex>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        Unroller<I,M/2,J,N,iscomplex>::unroll(x,m);
        Unroller<I+M/2,M-M/2,J+M/2,N-M/2,iscomplex>::unroll(x,m);
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
    { Unroller<0,s,0,s,M1::miscomplex>::unroll(x,m); }
  };

  // algo 6: Fully unroll by columns
  template <int s, int ix, class T, class M1>
  struct MultXU_Helper<6,s,ix,T,M1>
  {
    template <int I, int M, int J, int N, bool iscomplex>
    struct Unroller
    {
      static inline void unroll(const Scaling<ix,T>& x, M1& m)
      {
        Unroller<I,M-(N-N/2),J,N/2,iscomplex>::unroll(x,m);
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
    { Unroller<0,s,0,s,M1::miscomplex>::unroll(x,m); }
  };

  // algo 10: LowerTri, transpose:
  template <int s, int ix, class T, class M1>
  struct MultXU_Helper<10,s,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      typedef typename M1::transpose_type M1t;
      M1t mt = m.Transpose();
      MultXU_Helper<-1,s,ix,T,M1t>::call(x,mt);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int s, int ix, class T, class M1>
  struct MultXU_Helper<-1,s,ix,T,M1> 
  {
    static inline void call(const Scaling<ix,T>& x, M1& m)
    {
      TMVStaticAssert(!M1::munit || (ix == 1));
      typedef typename M1::value_type T1;
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
          (ix == 1) ? 0 :
          !ShapeTraits<M1::mshape>::upper ? 10 :
          unroll ? ( M1::mrowmajor ? 5 : 6 ) :
          M1::mrowmajor ? 2 : 3 ) };
      //std::cout<<"InlineMultXU: x = "<<T(x)<<std::endl;
      //std::cout<<"m = "<<TypeText(m)<<"  "<<m<<std::endl;
      //std::cout<<"algo = "<<algo<<std::endl;
      MultXU_Helper<algo,s,ix,T,M1>::call(x,m);
      //std::cout<<"m => "<<m<<std::endl;
    }
  };

  // m *= x
  template <int ix, class T, class M>
  inline void InlineMultXM(const Scaling<ix,T>& x, BaseMatrix_Tri_Mutable<M>& m)
  { MultXU_Helper<-1,M::msize,ix,T,M>::call(x,m.mat()); }

  // Defined in TMV_MultXU.cpp
  template <class T>
  void InstMultXM(const T x, UpperTriMatrixView<T,NonUnitDiag> m);
  template <class T>
  inline void InstMultXM(const T x, LowerTriMatrixView<T,NonUnitDiag> m)
  { InstMultXM(x,m.Transpose()); }


  //
  // Matrix * Scalar
  // -Matrix
  //

  // m2 = x * m1
  template <int ix, class T, class M1, class M2>
  inline void InlineMultXM(
      const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
      BaseMatrix_Tri_Mutable<M2>& m2)
  {
    TMVStaticAssert(ShapeTraits<M1::mshape>::upper == 
        int(ShapeTraits<M2::mshape>::upper));
    TMVStaticAssert(M1::munit || !M2::munit);
    TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
    TMVAssert(m1.size() == m2.size());
    enum { s = Sizes<M1::msize,M2::msize>::size };
    typedef typename M1::const_view_type M1v;
    typedef typename M2::view_type M2v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    AddUU_Helper<-1,s,false,ix,T,M1v,M2v>::call(x,m1v,m2v);
  }

  // Defined in TMV_MultXU.cpp
  template <class T1, DiagType D1, bool C1, class T2>
  void InstMultXM(const T2 x,
      const ConstUpperTriMatrixView<T1,D1,UNKNOWN,UNKNOWN,C1>& m1,
      UpperTriMatrixView<T2,NonUnitDiag> m2);

  template <class T1, DiagType D1, bool C1, class T2>
  void InstMultXM(const T2 x,
      const ConstLowerTriMatrixView<T1,D1,UNKNOWN,UNKNOWN,C1>& m1,
      LowerTriMatrixView<T2,NonUnitDiag> m2)
  { InstMultXM(x,m1.Transpose(),m2.Transpose()); }


} // namespace tmv

#endif
