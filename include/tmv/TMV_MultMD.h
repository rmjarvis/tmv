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
// Boston, MA  02320-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_MultMD_H
#define TMV_MultMD_H

#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_AddVV.h"
#include "TMV_Vector.h"
#include "TMV_SmallVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Prefetch.h"
#include "TMV_DiagMatrix.h"
#include "TMV_CopyM.h"
#include "TMV_MultXM.h"

//#define PRINTALGO_MD

#ifdef PRINTALGO_MD
#include <iostream>
#endif

namespace tmv {

  //
  // Matrix * DiagMatrix
  //

  // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
  // It doesn't really seem to matter much either way.
#define TMV_ZeroIX (ix == 0)
//#define TMV_ZeroIX (ix != 1)

  template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper;

  // algo 0: cs or rs = 0, so nothing to do
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<0,cs,rs,add,ix,T,M1,M2,M3>
  { static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& ) {} };

  // algo 1: cs == 1, so simplifies to a MultDV function
  template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<1,1,rs,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M1::const_row_type M1r;
      typedef typename M2::const_diag_type M2d;
      typedef typename M3::row_type M3r;
      M1r m10 = m1.get_row(0);
      M2d m2d = m2.diag();
      M3r m30 = m3.get_row(0);
      ElemMultVV_Helper<-1,rs,add,ix,T,M1r,M2d,M3r>::call(x,m10,m2d,m30);
    }
  };

  // algo 2: rs == 1, so simplifies to an AddVV or MultXV function
  template <int cs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<2,cs,1,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M1::const_col_type M1c;
      typedef typename M3::col_type M3c;
      typedef typename Traits2<T,typename M2::value_type>::type PT2;
      M1c m10 = m1.get_col(0);
      PT2 xm2 = ZProd<false,false>::prod(x , m2.cref(0));
      M3c m30 = m3.get_col(0);
      AddVV_Helper<-1,cs,add,0,PT2,M1c,M3c>::call(xm2,m10,m30);
    }
  };

  // algo 3: column major loop
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<3,cs,rs,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      const int M = (cs == UNKNOWN ? m3.colsize() : cs);
      int N = (rs == UNKNOWN ? m3.rowsize() : rs);
#ifdef PRINTALGO_MD
      std::cout<<"algo 3: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

      typedef typename M1::const_col_type M1c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      const int Astepj = m1.stepj();
      IT1 A = m1.NonConj().get_col(0).begin();

      typedef typename M2::const_diag_type M2d;
      typedef typename M2d::const_nonconj_type::const_iterator IT2;
      typedef typename M2::value_type T2;
      typedef typename Traits2<T,T2>::type PT2;
      IT2 D = m2.NonConj().diag().begin();
      PT2 dj;

      typedef typename M3::col_type M3c;
      typedef typename M3c::iterator IT3;
      const int Bstepj = m3.stepj();
      IT3 B = m3.get_col(0).begin();
      enum { c2 = M2::mconj };

      if (N) do {
        dj = ZProd<false,c2>::prod(x , *D++);
        AddVV_Helper<-1,cs,add,0,PT2,M1c,M3c>::call2(M,dj,A,B);
        A.ShiftP(Astepj);
        B.ShiftP(Bstepj);
      } while (--N);
    }
  };

  // algo 4: The basic row major loop
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<4,cs,rs,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      int M = (cs == UNKNOWN ? m3.colsize() : cs);
      const int N = (rs == UNKNOWN ? m3.rowsize() : rs);
#ifdef PRINTALGO_MD
      std::cout<<"algo 4: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
      //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
#endif

      typedef typename M1::const_row_type M1r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      const int Astepi = m1.stepi();
      IT1 A = m1.NonConj().get_row(0).begin();

      typedef typename M2::const_diag_type M2d;
      typedef typename M2d::const_nonconj_type::const_iterator IT2;
      const IT2 D = m2.diag().NonConj().begin();

      typedef typename M3::row_type M3r;
      typedef typename M3r::iterator IT3;
      const int Bstepi = m3.stepi();
      IT3 B = m3.get_row(0).begin();

      if (M) do {
        ElemMultVV_Helper<-1,rs,add,ix,T,M1r,M2d,M3r>::call2(N,x,A,D,B);
        A.ShiftP(Astepi);
        B.ShiftP(Bstepi);
      } while (--M);
    }
  };

  // algo 5: rowmajor, ix==0 or m2.step != 1, so copy m2
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<5,cs,rs,add,ix,T,M1,M2,M3> // rs known
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX || M2::mstep != 1);
      const int M = (cs == UNKNOWN ? m3.colsize() : cs);
      const int N = (rs == UNKNOWN ? m3.rowsize() : rs);
      typedef typename M2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MD
      std::cout<<"algo 5: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
      //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
#endif
      typedef typename Traits2<T,T2>::type T2c;
      typedef SmallDiagMatrix<T2c,rs> M2c;
      typedef typename M2::const_diag_type M2d;
      typedef typename M2c::diag_type M2cd;
      typedef typename M2c::const_view_type M2cv;
      M2c m2c;
      M2d m2d = m2.diag();
      M2cd m2cd = m2c.diag();
      M2cv m2cv = m2c.View();
      AddVV_Helper<-1,rs,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
      MultMD_Helper<4,cs,rs,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
    }
  };
  template <int cs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<5,cs,UNKNOWN,add,ix,T,M1,M2,M3> // rs unknown
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX || M2::mstep != 1);
      const int M = (cs == UNKNOWN ? m3.colsize() : cs);
      const int rs = UNKNOWN;
      const int N = m3.rowsize();
      typedef typename M2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO_MD
      std::cout<<"algo 5: M,N,cs,rs,x = "<<M<<','<<N<<','<<cs<<','<<rs<<','<<T(x)<<std::endl;
      //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
#endif
      typedef typename Traits2<T,T2>::type T2c;
      typedef DiagMatrix<T2c> M2c;
      typedef typename M2::const_diag_type M2d;
      typedef typename M2c::diag_type M2cd;
      typedef typename M2c::const_view_type M2cv;
      M2c m2c(N);
      M2d m2d = m2.diag();
      M2cd m2cd = m2c.diag();
      M2cv m2cv = m2c.View();
      AddVV_Helper<-1,rs,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
#ifdef PRINTALGO_MD
      //std::cout<<"m2 = "<<m2.diag()<<std::endl;
      //std::cout<<"m2c => "<<m2c.diag()<<std::endl;
      //std::cout<<"x * m2 = "<<x*m2.diag().copy()<<std::endl;
#endif
      MultMD_Helper<4,cs,rs,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int cs, int rs, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultMD_Helper<-1,cs,rs,add,ix,T,M1,M2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      // Possible algorithms to choose from:
      //  0 = cs or rs == 0, so nothing to do
      //  1 = cs == 1: reduces to trivial MultDV function
      //  2 = rs == 1: reduces to trivial AddVV function
      //  3 = column major, simple for loop
      //  4 = row major
      //  5 = row major, ix==0 or m2.step() != 1, so copy m2

      enum { rm = M1::mrowmajor && M3::mrowmajor };
      enum { algo = (
          ( rs == 0 || cs == 0 ) ? 0 :
          ( cs == 1 ) ? 1 :
          ( rs == 1 ) ? 2 :
          rm ? ( 
#if TMV_OPT >= 1
            (TMV_ZeroIX || M2::mstep != 1) ? 5 : 
#endif
            4 ) :
          3 ) };
#ifdef PRINTALGO_MD
      std::cout<<"InlineMultMD: \n";
      std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
      std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
#endif
      MultMD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void InlineMultMM(const Scaling<ix,T>& x, 
      const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
      BaseMatrix_Rec_Mutable<M3>& m3)
  {
    //std::cout<<"InlineMultMD "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"add = "<<add<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3<<std::endl;
    TMVStaticAssert((Sizes<M3::mcolsize,M1::mcolsize>::same));
    TMVStaticAssert((Sizes<M3::mrowsize,M1::mrowsize>::same));
    TMVStaticAssert((Sizes<M3::mrowsize,M2::msize>::same));
    TMVAssert(m3.colsize() == m1.colsize());
    TMVAssert(m3.rowsize() == m1.rowsize());
    TMVAssert(m3.rowsize() == m2.size());
    const int cs = Sizes<M3::mcolsize,M1::mcolsize>::size;
    const int rs = Sizes<Sizes<M3::mrowsize,M1::mrowsize>::size,M2::msize>::size;
    typedef typename M1::const_view_type M1v;
    typedef typename M2::const_view_type M2v;
    typedef typename M3::view_type M3v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    M3v m3v = m3.View();
    MultMD_Helper<-1,cs,rs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    //std::cout<<"End InlineMultMD: m3 = "<<m3<<std::endl;
  }

  // D * M
  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void InlineMultMM(const Scaling<ix,T>& x,
      const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
      BaseMatrix_Mutable<M3>& m3)
  {
    typename M3::transpose_type m3t = m3.Transpose();
    InlineMultMM<add>(x,m2.Transpose(),m1.Transpose(),m3t);
  }


  // Defined in TMV_MultMD.cpp
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3);
  template <class T1, bool C1, class T2, bool C2, class T3>
  void InstAddMultMM(const T3 x,
      const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3);

  template <class T1, bool C1, class T2, bool C2, class T3>
  inline void InstMultMM(const T3 x,
      const ConstDiagMatrixView<T1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
  { InstMultMM(x,m2.Transpose(),m1.Transpose(),m3.Transpose()); }
  template <class T1, bool C1, class T2, bool C2, class T3>
  inline void InstAddMultMM(const T3 x,
      const ConstDiagMatrixView<T1,UNKNOWN,C1>& m1,
      const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3)
  { InstAddMultMM(x,m2.Transpose(),m1.Transpose(),m3.Transpose()); }

} // namespace tmv

#undef TMV_ZeroIX

#endif 
