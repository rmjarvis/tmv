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

#ifndef TMV_SwapM_H
#define TMV_SwapM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_MultXM.h"

namespace tmv {

  const int TMV_PERM_BLOCKSIZE = 32;

  //
  // Swap Matrices
  //

  template <int algo, int cs, int rs, class M1, class M2>
  struct SwapM_Helper;

  // algo 1: Linearize to vector version
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<1,cs,rs,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      typedef typename M1::linearview_type M1l;
      typedef typename M2::linearview_type M2l;
      M1l m1l = m1.LinearView();
      M2l m2l = m2.LinearView();
      enum { cs_rs = IntTraits2<cs,rs>::prod };
      SwapV_Helper<-1,cs_rs,M1l,M2l>::call(m1l,m2l);
    }
  };

  // algo 2: Loop over rows
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<2,cs,rs,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::row_type M1r;
      typedef typename M2::row_type M2r;
      typedef typename M1r::iterator IT1;
      typedef typename M2r::iterator IT2;
      const int step1 = m1.stepi();
      const int step2 = m2.stepi();
      IT1 it1 = m1.get_row(0).begin();
      IT2 it2 = m2.get_row(0).begin();
      for(;M;--M) {
        SwapV_Helper<-1,rs,M1r,M2r>::call2(N,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      } 
    }
  };

  // algo 3: Loop over columns
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<3,cs,rs,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
      const int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
      int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
      typedef typename M1::col_type M1c;
      typedef typename M2::col_type M2c;
      typedef typename M1c::iterator IT1;
      typedef typename M2c::iterator IT2;
      const int step1 = m1.stepj();
      const int step2 = m2.stepj();
      IT1 it1 = m1.get_col(0).begin();
      IT2 it2 = m2.get_col(0).begin();
      for(;N;--N) {
        SwapV_Helper<-1,cs,M1c,M2c>::call2(M,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      } 
    }
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<4,cs,rs,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
    {
#if TMV_OPT >= 2
      if (m1.CanLinearize() && m2.CanLinearize() &&
          m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
        SwapM_Helper<1,cs,rs,M1,M2>::call(m1,m2);
      else if (m1.isrm() && m2.isrm() || 
          ( !(m1.iscm() && m2.iscm()) &&
            (m1.colsize() > m1.rowsize()) ) )
        SwapM_Helper<2,cs,rs,M1,M2>::call(m1,m2);
      else 
        SwapM_Helper<3,cs,rs,M1,M2>::call(m1,m2);
#else
      enum { algo2 = (
#if TMV_OPT >= 1
          ( M1::mrowmajor && M2::mrowmajor ) ? 2 :
          ( M1::mcolmajor && M2::mcolmajor ) ? 3 :
          ( cs == UNKNOWN || rs == UNKNOWN ) ? 3 :
          ( cs > rs ) ? 2 : 
#endif
          3 ) };
      SwapM_Helper<algo2,cs,rs,M1,M2>::call(m1,m2);
#endif
    }
  };

  // algo 5: Fully unroll by rows
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<5,cs,rs,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(M1& m1, M2& m2)
      {
        Unroller<I,M/2,J,N>::unroll(m1,m2);
        Unroller<I+M/2,M-M/2,J,N>::unroll(m1,m2);
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
    { Unroller<0,cs,0,rs>::unroll(m1,m2); }
  };

  // algo 6: Fully unroll by columns
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<6,cs,rs,M1,M2>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(M1& m1, M2& m2)
      {
        Unroller<I,M,J,N/2>::unroll(m1,m2);
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
    { Unroller<0,cs,0,rs>::unroll(m1,m2); }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, class M1, class M2>
  struct SwapM_Helper<-1,cs,rs,M1,M2>
  {
    static inline void call(M1& m1, M2& m2)
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
      SwapM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
    }
  };

  template <class M1, class M2>
  inline void InlineSwap(
      BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
  {
    //std::cout<<"Inline Swap: m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same));
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    enum { cs = Sizes<M1::mcolsize,M2::mcolsize>::size };
    enum { rs = Sizes<M1::mrowsize,M2::mrowsize>::size };
    typedef typename M1::view_type M1v;
    typedef typename M2::view_type M2v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    SwapM_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    //std::cout<<"m1 => "<<m1.mat()<<std::endl;
    //std::cout<<"m2 => "<<m2.mat()<<std::endl;
  }

  // Defined in TMV_Matrix.cpp
  template <class T, bool C>
  void InstSwap(MatrixView<T,UNKNOWN,UNKNOWN,C> m1, MatrixView<T> m2); 

  template <bool checkalias, bool conj, bool inst, class M1, class M2>
  struct CallSwapm  // checkalias = true
  {
    static inline void call(M1& m1, M2& m2)
    {
      if (!SameStorage(m1,m2)) 
      {
        CallSwapm<false,conj,inst,M1,M2>::call(m1,m2);
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
  template <bool inst, class M1, class M2>
  struct CallSwapm<false,true,inst,M1,M2> // conj = true
  {
    static inline void call(M1& m1, M2& m2)
    { 
      typedef typename M1::conjugate_type M1c;
      typedef typename M2::conjugate_type M2c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      CallSwapm<false,false,inst,M1c,M2c>::call(m1c,m2c);
    }
  };
  template <class M1, class M2>
  struct CallSwapm<false,false,true,M1,M2> // inst = true
  {
    static inline void call(M1& m1, M2& m2)
    { InstSwap(m1.XView(),m2.XView()); }
  };
  template <class M1, class M2>
  struct CallSwapm<false,false,false,M1,M2> // inst = false
  {
    static inline void call(M1& m1, M2& m2)
    { InlineSwap(m1,m2); }
  };

  template <class M1, class M2> 
  inline void DoSwap(
      BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
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
        (M1::mcolmajor || M1::mrowmajor) &&
        (M2::mcolmajor || M2::mrowmajor) &&
        checkalias ) };
    CallSwapm<checkalias,M2::mconj,inst,M1,M2>::call(m1.mat(),m2.mat());
  }
  template <class M1, class M2> 
  inline void Swap(
      BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
  { DoSwap(m1,m2); }


  //
  // TransposeSelf
  //

  template <int algo, int size, class M1>
  struct TransposeSelf_Helper;

  // algo 1: Simple for loop
  template <int size, class M1>
  struct TransposeSelf_Helper<1,size,M1>
  {
    static inline void call(M1& m)
    {
      const size_t n = (size == UNKNOWN ? int(m.colsize()) : size);
      for(int i=1;i<n;++i) {
        typename M1::row_range_type v1 = m.get_row(i,0,i);
        typename M1::col_range_type v2 = m.get_col(i,0,i);
        Swap(v1,v2);
      }
    }
  };

  // algo 2: The same thing, but with iterators.
  template <int size, class M1>
  struct TransposeSelf_Helper<2,size,M1>
  {
    static inline void call(M1& m)
    {
      const size_t n = size == UNKNOWN ? int(m.colsize()) : size;
      if (n <= 1) return;
      typedef typename M1::row_type Mr;
      typedef typename M1::col_type Mc;
      typedef typename Mr::iterator IT1;
      typedef typename Mc::iterator IT2;
      IT1 it1 = m.row(1).begin();
      IT2 it2 = m.col(1).begin();
      const int step1 = m.stepi();
      const int step2 = m.stepj();
      for(int i=1;i<n;++i) {
        SwapV_Helper<-1,UNKNOWN,Mr,Mc>::call2(i,it1,it2);
        it1.ShiftP(step1);
        it2.ShiftP(step2);
      }
    }
  };

  // algo 3: The other way to do the loop.  
  // This way seems to be a little bit slower.
  template <int size, class M1>
  struct TransposeSelf_Helper<3,size,M1>
  {
    static inline void call(M1& m)
    {
      int n = size == UNKNOWN ? int(m.colsize()) : size;
      if (n <= 1) return;
      typedef typename M1::row_type Mr;
      typedef typename M1::col_type Mc;
      typedef typename Mr::iterator IT1;
      typedef typename Mc::iterator IT2;
      IT1 it1 = m.row(0).begin(); ++it1;
      IT2 it2 = m.col(0).begin(); ++it2;
      const int step = m.stepi() + m.stepj();
      for(--n;n;--n) {
        SwapV_Helper<-1,UNKNOWN,Mr,Mc>::call2(n,it1,it2);
        it1.ShiftP(step);
        it2.ShiftP(step);
      }
    }
  };

  // algo 5: Fully unroll
  template <int size, class M1>
  struct TransposeSelf_Helper<5,size,M1>
  {
    template <int I, int M, int J, int N>
    struct Unroller
    {
      static inline void unroll(M1& m)
      {
        Unroller<I,M/2,0,I>::unroll(m);
        Unroller<I+M/2,M-M/2,0,I+M/2>::unroll(m);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,1,J,N>
    {
      static inline void unroll(M1& m)
      {
        Unroller<I,1,J,N/2>::unroll(m);
        Unroller<I,1,J+N/2,N-N/2>::unroll(m);
      }
    };
    template <int I, int J, int N>
    struct Unroller<I,0,J,N>
    { static inline void unroll(M1& ) {} };
    template <int I, int J>
    struct Unroller<I,1,J,1>
    {
      static inline void unroll(M1& m)
      { TMV_SWAP(m.ref(I,J) , m.ref(J,I) ); }
    };
    template <int I, int J>
    struct Unroller<I,1,J,0>
    { static inline void unroll(M1& ) {} };

    static inline void call(M1& m)
    { Unroller<0,size,0,size>::unroll(m); }
  };

  // algo -1: Determine which algorithm to use
  template <int size, class M1>
  struct TransposeSelf_Helper<-1,size,M1>
  {
    static inline void call(M1& m)
    {
      enum { algo = (
#if TMV_OPT >= 1
          ( size != UNKNOWN && size < 8 ) ? 1 :
#endif
          2 ) };
      //std::cout<<"algo = "<<algo<<std::endl;
      TransposeSelf_Helper<algo,size,M1>::call(m);
    }
  };

  template <class M1>
  inline void InlineTransposeSelf(BaseMatrix_Rec_Mutable<M1>& m)
  {
    enum { size = Sizes<M1::mcolsize,M1::mrowsize>::size };
    typedef typename M1::view_type M1v;
    M1v m1v = m.View();
    TransposeSelf_Helper<-1,size,M1v>::call(m1v);
  }


  // Defined in TMV_Matrix.cpp
  template <class T>
  void InstTransposeSelf(MatrixView<T> m);

  template <class T>
  inline void InstTransposeSelf(MatrixView<T,UNKNOWN,UNKNOWN,true> m)
  { InstTransposeSelf(m.Conjugate()); }
  
  template <bool inst, class M>
  struct CallTransposeSelf // inst = true
  {
    static inline void call(M& m) 
    { InstTransposeSelf(m.XView()); } 
  };
  template <class M>
  struct CallTransposeSelf<false,M> // inst = false
  {
    static inline void call(M& m) 
    { InlineTransposeSelf(m); } 
  };

  template <class M>
  inline void TransposeSelf(BaseMatrix_Rec_Mutable<M>& m)
  {
    TMVStaticAssert((Sizes<M::mcolsize,M::mrowsize>::same)); 
    TMVAssert(m.colsize() == m.rowsize());
    typedef typename M::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        (M::mcolsize == UNKNOWN && M::mrowsize == UNKNOWN ) ) };
    CallTransposeSelf<inst,M>::call(m.mat());
  }


  //
  // PermuteRows
  //

  template <int algo, int cs, int rs, class M1> struct PermuteRows_Helper;

  // algo 1: Simple loop over rows
  template <int cs, int rs, class M1>
  struct PermuteRows_Helper<1,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      p += i1;
      for(int i=i1;i<i2;++i,++p) {
        TMVAssert(*p < int(m.colsize()));
        m.CSwapRows(i,*p);
      }
    }
  };
  
  // algo 2: Simple loop over columns
  template <int cs, int rs, class M1>
  struct PermuteRows_Helper<2,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      for (int j=0;j<N;++j) 
        m.get_col(j).Permute(p,i1,i2);
    }
  };

  // algo 3: Loop over rows with iterators
  template <int cs, int rs, class M1>
  struct PermuteRows_Helper<3,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      const int M = cs == UNKNOWN ? int(m.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      typedef typename M1::row_type M1r;
      typedef typename M1r::iterator IT;
      IT it1 = m.get_row(0).begin();
      const int stepi = m.stepi();
      it1.ShiftP(i1*stepi);
      p += i1;
      for(int i=i1;i<i2;++i,++p) {
        TMVAssert(*p < M);
        if (*p != i) {
          IT it2 = it1;
          it2.ShiftP((*p-i)*stepi);
          SwapV_Helper<-1,rs,M1r,M1r>::call2(N,it1,it2);
        }
        it1.ShiftP(stepi);
      }
    }
  };
  
  // algo 4: Loop over columns in blocks, then over rows with a block
  template <int cs, int rs, class M1>
  struct PermuteRows_Helper<4,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      const int M = cs == UNKNOWN ? int(m.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      typedef typename M1::row_type M1r;
      typedef typename M1r::iterator IT;

      int N_32 = (N>>5); // N_32 = N/32
      const int Nx = N - (N_32<<5); // Nx = N % 32
      enum { rsx = (rs == UNKNOWN ? UNKNOWN : (rs % 32) ) };
      IT it1 = m.get_row(0).begin();
      const int stepi = m.stepi();
      it1.ShiftP(i1*stepi);
      p += i1;
      if (N_32) do {
        const int* pi = p;
        IT it1i = it1;
        for(int i=i1;i<i2;++i,++pi) {
          TMVAssert(*pi < M);
          if (*pi != i) {
            IT it2 = it1i;
            it2.ShiftP((*pi-i)*stepi);
            SwapV_Helper<-1,32,M1r,M1r>::call2(32,it1i,it2);
          }
          it1i.ShiftP(stepi);
        }
        it1 += 32;
      } while (--N_32);
      if (Nx) {
        for(int i=i1;i<i2;++i,++p) {
          TMVAssert(*p < M);
          if (*p != i) {
            IT it2 = it1;
            it2.ShiftP((*p-i)*stepi);
            SwapV_Helper<-1,rsx,M1r,M1r>::call2(Nx,it1,it2);
          }
          it1.ShiftP(stepi);
        }
      }
    }
  };
  
  // algo -1: Determine which algorithm to use
  template <int cs, int rs, class M1>
  struct PermuteRows_Helper<-1,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int*const p, const int i1, const int i2)
    {
      enum { algo = (
#if TMV_OPT >= 1
          (M1::mcolmajor && (rs != UNKNOWN && rs <= 32)) ? 2 :
          M1::mcolmajor ? 4 : M1::mrowmajor ? 3 :
#endif
          1 ) };
      PermuteRows_Helper<algo,cs,rs,M1>::call(m,p,i1,i2);
    }
  };

  template <class M>
  inline void InlinePermuteRows(BaseMatrix_Rec_Mutable<M>& m, 
      const int*const p, const int i1, const int i2)
  {
    typedef typename M::view_type Mv;
    Mv mv = m.View();
    PermuteRows_Helper<-1,M::mcolsize,M::mrowsize,Mv>::call(mv,p,i1,i2); 
  }

  // Defined in TMV_Matrix.cpp
  template <class T>
  void InstPermuteRows(MatrixView<T> m, 
      const int*const p, const int i1, const int i2);

  template <class T>
  inline void InstPermuteRows(MatrixView<T,UNKNOWN,UNKNOWN,true> m, 
      const int*const p, const int i1, const int i2)
  { InstPermuteRows(m.Conjugate(),p,i1,i2); }

  template <bool inst, class M>
  struct CallPermuteRows // inst = true
  {
    static inline void call(M& m, 
        const int*const p, const int i1, const int i2) 
    { InstPermuteRows(m.XView(),p,i1,i2); } 
  };
  template <class M>
  struct CallPermuteRows<false,M> // inst = false
  {
    static inline void call(M& m, 
        const int*const p, const int i1, const int i2) 
    { InlinePermuteRows(m,p,i1,i2); } 
  };

  template <class M>
  inline void PermuteRows(BaseMatrix_Rec_Mutable<M>& m,
      const int*const p, const int i1, const int i2)
  {
    typedef typename M::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::mcolsize == UNKNOWN &&
        M::mrowsize == UNKNOWN) };
    CallPermuteRows<inst,M>::call(m.mat(),p,i1,i2);
  }

  template <int algo, int cs, int rs, class M1> struct ReversePermuteRows_Helper;

  // algo 1: Simple loop over rows
  template <int cs, int rs, class M1>
  struct ReversePermuteRows_Helper<1,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      p += i2-1;
      for(int i=i2-1;i>=i1;--i,--p) {
        TMVAssert(*p < int(m.colsize()));
        m.CSwapRows(i,*p);
      }
    }
  };
  
  // algo 2: Simple loop over columns
  template <int cs, int rs, class M1>
  struct ReversePermuteRows_Helper<2,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int*const p, const int i1, const int i2)
    {
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      for (int j=0;j<N;++j) 
        m.get_col(j).ReversePermute(p,i1,i2);
    }
  };

  // algo 3: Loop over rows with iterators
  template <int cs, int rs, class M1>
  struct ReversePermuteRows_Helper<3,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      const int M = cs == UNKNOWN ? int(m.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      typedef typename M1::row_type M1r;
      typedef typename M1r::iterator IT;
      IT it1 = m.get_row(0).begin();
      const int stepi = m.stepi();
      Prefetch_Write(it1.GetP()+i1*stepi);
      it1.ShiftP((i2-1)*stepi);
      p += i2-1;
      for(int i=i2-1;i>=i1;--i,--p) {
        TMVAssert(*p < M);
        if (*p != i) {
          IT it2 = it1;
          it2.ShiftP((*p-i)*stepi);
          SwapV_Helper<-1,rs,M1r,M1r>::call2(N,it1,it2);
        }
        it1.ShiftP(-stepi);
      }
    }
  };

  // algo 4: Loop over columns in blocks, then over rows with a block
  template <int cs, int rs, class M1>
  struct ReversePermuteRows_Helper<4,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int* p, const int i1, const int i2)
    {
      const int M = cs == UNKNOWN ? int(m.colsize()) : cs;
      const int N = rs == UNKNOWN ? int(m.rowsize()) : rs;
      typedef typename M1::value_type T1;
      typedef typename M1::row_type M1r;
      typedef typename M1r::iterator IT;

      int N_32 = (N>>5); // N_32 = N/32
      const int Nx = N - (N_32<<5); // Nx = N % 32
      enum { rsx = (rs == UNKNOWN ? UNKNOWN : (rs % 32) ) };
      IT it1 = m.get_row(0).begin();
      const int stepi = m.stepi();
      const int stepj = m.stepj();
      it1.ShiftP((i2-1)*stepi);
      p += i2-1;

      if (N_32) do {
        const int* pi = p;
        IT it1i = it1;
        for(int i=i2-1;i>=i1;--i,--pi) {
          TMVAssert(*pi < M);
          if (*pi != i) {
            IT it2 = it1i;
            it2.ShiftP((*pi-i)*stepi);
            SwapV_Helper<-1,32,M1r,M1r>::call2(32,it1i,it2);
          }
          it1i.ShiftP(-stepi);
        }
        it1 += 32;
      } while (--N_32);
      if (Nx) {
        for(int i=i2-1;i>=i1;--i,--p) {
          TMVAssert(*p < M);
          if (*p != i) {
            IT it2 = it1;
            it2.ShiftP((*p-i)*stepi);
            SwapV_Helper<-1,rsx,M1r,M1r>::call2(Nx,it1,it2);
          }
          it1.ShiftP(-stepi);
        }
      }
    }
  };
  
  // algo -1: Determine which algorithm to use
  template <int cs, int rs, class M1>
  struct ReversePermuteRows_Helper<-1,cs,rs,M1>
  {
    static inline void call(M1& m, 
        const int*const p, const int i1, const int i2)
    {
      enum { algo = (
#if TMV_OPT >= 1
          (M1::mrowmajor && (rs != UNKNOWN && rs <= 32)) ? 3 :
          M1::mcolmajor ? 2 : M1::mrowmajor ? 4 :
#endif
          1 ) };
      ReversePermuteRows_Helper<algo,cs,rs,M1>::call(m,p,i1,i2);
    }
  };

  template <class M>
  inline void InlineReversePermuteRows(BaseMatrix_Rec_Mutable<M>& m, 
      const int*const p, const int i1, const int i2)
  { 
    typedef typename M::view_type Mv;
    Mv mv = m.View();
    ReversePermuteRows_Helper<-1,M::mcolsize,M::mrowsize,Mv>::call(
      mv,p,i1,i2); 
  }

  // Defined in TMV_Matrix.cpp
  template <class T>
  void InstReversePermuteRows(MatrixView<T> m, 
      const int*const p, const int i1, const int i2);

  template <class T>
  inline void InstReversePermuteRows(MatrixView<T,UNKNOWN,UNKNOWN,true> m, 
      const int*const p, const int i1, const int i2)
  { InstReversePermuteRows(m.Conjugate(),p,i1,i2); }

  template <bool inst, class M>
  struct CallReversePermuteRows // inst = true
  {
    static inline void call(M& m, 
        const int*const p, const int i1, const int i2) 
    { InstReversePermuteRows(m.XView(),p,i1,i2); } 
  };
  template <class M>
  struct CallReversePermuteRows<false,M> // inst = false
  {
    static inline void call(M& m, 
        const int*const p, const int i1, const int i2) 
    { InlineReversePermuteRows(m,p,i1,i2); } 
  };

  template <class M>
  inline void ReversePermuteRows(BaseMatrix_Rec_Mutable<M>& m,
      const int*const p, const int i1, const int i2)
  {
    typedef typename M::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::mcolsize == UNKNOWN &&
        M::mrowsize == UNKNOWN) };
    CallReversePermuteRows<inst,M>::call(m.mat(),p,i1,i2);
  }

} // namespace tmv

#endif
