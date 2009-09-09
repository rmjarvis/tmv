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


#ifndef TMV_MultUD_H
#define TMV_MultUD_H

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

//#define PRINTALGO

#ifdef PRINTALGO
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

  template <int algo, int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper;

  // algo 0: s = 0, so nothing to do
  template <bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<0,0,add,ix,T,M1,M2,M3>
  { static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& ) {} };

  // algo 1: s == 1, so simplifies to a scalar product
  template <bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<1,1,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      Maybe<add>::add(m3.ref(0), ZProd<false,false>::prod(x , 
            ZProd<false,false>::prod(m1.cref(0) , m2.cref(0)) ));
    }
  };

  // algo 11: UpperTri: column major loop
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<11,s,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO
      std::cout<<"algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      enum { unit = M1::munit };
      enum { c1 = M1::mconj };
      enum { c2 = M2::mconj };

      typedef typename M1::const_col_range_type M1c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      const int Astepj = m1.stepj();
      IT1 A = m1.get_col(0,0,1).NonConj().begin();

      typedef typename M2::const_diag_type M2d;
      typedef typename M2d::const_nonconj_type::const_iterator IT2;
      typedef typename M2::value_type T2;
      typedef typename Traits2<T,T2>::type PT2;
      IT2 D = m2.diag().NonConj().begin();
      PT2 dj;

      typedef typename M3::col_range_type M3c;
      typedef typename M3c::iterator IT3;
      const int Bstepj = m3.stepj();
      IT3 B = m3.get_col(0,0,1).begin();

      int len = 1;
      if (N) do {
        dj = ZProd<false,c2>::prod(x , *D++);
        AddVV_Helper<-1,UNKNOWN,add,0,PT2,M1c,M3c>::call2(len++,dj,A,B);
        A.ShiftP(Astepj);
        B.ShiftP(Bstepj);
      } while (--N);
    }
  };

  // algo 12: UpperTri: row major loop
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<12,s,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO
      std::cout<<"algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

      typedef typename M1::const_row_range_type M1r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      const int Adiagstep = m1.diagstep();
      IT1 A = m1.get_row(0,0,N).NonConj().begin();

      typedef typename M2::const_diag_type M2d;
      typedef typename M2d::const_nonconj_type::const_iterator IT2;
      IT2 D = m2.diag().NonConj().begin();

      typedef typename M3::row_range_type M3r;
      typedef typename M3r::iterator IT3;
      const int Bdiagstep = m3.diagstep();
      IT3 B = m3.get_row(0,0,N).begin();

      if (N) do {
        ElemMultVV_Helper<-1,UNKNOWN,add,ix,T,M1r,M2d,M3r>::call2(N,x,A,D++,B);
        A.ShiftP(Adiagstep);
        B.ShiftP(Bdiagstep);
      } while (--N);
    }
  };

  // TODO: Write an algo 15, unroll by cols

  // algo 17: UpperTri: m1 is unit diag, do calculation on offdiag part.
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<17,s,add,ix,T,M1,M2,M3> // s known
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      const int N = (s == UNKNOWN ? m3.size() : s);
      TMVStaticAssert(M1::munit);
      TMVStaticAssert(!M3::munit);
#ifdef PRINTALGO
      std::cout<<"algo 17: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef typename M1::const_offdiag_type M1o;
      typedef typename M2::const_subdiagmatrix_type M2s;
      typedef typename M3::offdiag_type M3o;
      enum { sm1 = IntTraits2<s,-1>::sum };
      M1o m1o = m1.OffDiag();
      M2s m2s = m2.SubDiagMatrix(1,N);
      M3o m3o = m3.OffDiag();
      MultUD_Helper<-1,sm1,add,ix,T,M1o,M2s,M3o>::call(x,m1o,m2s,m3o);
      typename M3::diag_type m3d = m3.diag();
      Maybe<add>::add(m3d , x * m2.diag());
    }
  };

  // algo 18: UpperTri: rowmajor, ix==0 or m2.step != 1, so copy m2
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<18,s,add,ix,T,M1,M2,M3> // s known
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX || M2::mstep != 1);
      const int N = (s == UNKNOWN ? m3.size() : s);
      typedef typename M2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 18: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef typename Traits2<T,T2>::type T2c;
      typedef SmallDiagMatrix<T2c,s> M2c;
      typedef typename M2::const_diag_type M2d;
      typedef typename M2c::diag_type M2cd;
      typedef typename M2c::const_view_type M2cv;
      M2c m2c;
      M2d m2d = m2.diag();
      M2cd m2cd = m2c.diag();
      M2cv m2cv = m2c.View();
      AddVV_Helper<-1,s,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
      MultUD_Helper<12,s,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
    }
  };
  template <bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<18,UNKNOWN,add,ix,T,M1,M2,M3> // s unknown
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX || M2::mstep != 1);
      const int s = UNKNOWN;
      const int N = m3.size();
      typedef typename M2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 18: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
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
      AddVV_Helper<-1,s,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
      MultUD_Helper<12,s,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
    }
  };


  // algo 21: LowerTri: column major loop
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<21,s,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO
      std::cout<<"algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      enum { unit = M1::munit };
      enum { c1 = M1::mconj };
      enum { c2 = M2::mconj };

      typedef typename M1::const_col_range_type M1c;
      typedef typename M1c::const_nonconj_type::const_iterator IT1;
      const int Adiagstep = m1.diagstep();
      IT1 A = m1.get_col(0,0,N).NonConj().begin();

      typedef typename M2::const_diag_type M2d;
      typedef typename M2d::const_nonconj_type::const_iterator IT2;
      typedef typename M2::value_type T2;
      typedef typename Traits2<T,T2>::type PT2;
      IT2 D = m2.diag().NonConj().begin();
      PT2 dj;

      typedef typename M3::col_range_type M3c;
      typedef typename M3c::iterator IT3;
      const int Bdiagstep = m3.diagstep();
      IT3 B = m3.get_col(0,0,N).begin();

      if (N) do {
        dj = ZProd<false,c2>::prod(x , *D++);
        AddVV_Helper<-1,UNKNOWN,add,0,PT2,M1c,M3c>::call2(N,dj,A,B);
        A.ShiftP(Adiagstep);
        B.ShiftP(Bdiagstep);
      } while (--N);
    }
  };

  // algo 22: LowerTri: row major loop
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<22,s,add,ix,T,M1,M2,M3>
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO
      std::cout<<"algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

      typedef typename M1::const_row_range_type M1r;
      typedef typename M1r::const_nonconj_type::const_iterator IT1;
      const int Astepi = m1.stepi();
      IT1 A = m1.get_row(0,0,1).NonConj().begin();

      typedef typename M2::const_diag_type M2d;
      typedef typename M2d::const_nonconj_type::const_iterator IT2;
      const IT2 D = m2.diag().NonConj().begin();

      typedef typename M3::row_range_type M3r;
      typedef typename M3r::iterator IT3;
      const int Bstepi = m3.stepi();
      IT3 B = m3.get_row(0,0,1).begin();

      int len = 1;
      if (N) do {
        ElemMultVV_Helper<-1,UNKNOWN,add,ix,T,M1r,M2d,M3r>::call2(
            len++,x,A,D,B);
        A.ShiftP(Astepi);
        B.ShiftP(Bstepi);
      } while (--N);
    }
  };

  // TODO: Write an algo 25, unroll by cols

  // algo 27: LowerTri: m1 is unit diag, do calculation on offdiag part.
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<27,s,add,ix,T,M1,M2,M3> // s known
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      const int N = (s == UNKNOWN ? m3.size() : s);
      TMVStaticAssert(M1::munit);
      TMVStaticAssert(!M3::munit);
#ifdef PRINTALGO
      std::cout<<"algo 27: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef typename M1::const_offdiag_type M1o;
      typedef typename M2::const_subdiagmatrix_type M2s;
      typedef typename M3::offdiag_type M3o;
      enum { sm1 = IntTraits2<s,-1>::sum };
      M1o m1o = m1.OffDiag();
      M2s m2s = m2.SubDiagMatrix(0,N-1);
      M3o m3o = m3.OffDiag();
      MultUD_Helper<-1,sm1,add,ix,T,M1o,M2s,M3o>::call(x,m1o,m2s,m3o);
      typename M3::diag_type m3d = m3.diag();
      Maybe<add>::add(m3d , x * m2.diag());
    }
  };

  // algo 28: LowerTri: rowmajor, ix==0 or m2.step != 1, so copy m2
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<28,s,add,ix,T,M1,M2,M3> // s known
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX || M2::mstep != 1);
      const int N = (s == UNKNOWN ? m3.size() : s);
      typedef typename M2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 28: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
      typedef typename Traits2<T,T2>::type T2c;
      typedef SmallDiagMatrix<T2c,s> M2c;
      typedef typename M2::const_diag_type M2d;
      typedef typename M2c::diag_type M2cd;
      typedef typename M2c::const_view_type M2cv;
      M2c m2c;
      M2d m2d = m2.diag();
      M2cd m2cd = m2c.diag();
      M2cv m2cv = m2c.View();
      AddVV_Helper<-1,s,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
      MultUD_Helper<22,s,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
    }
  };
  template <bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<28,UNKNOWN,add,ix,T,M1,M2,M3> // s unknown
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      TMVStaticAssert(TMV_ZeroIX || M2::mstep != 1);
      const int s = UNKNOWN;
      const int N = m3.size();
      typedef typename M2::value_type T2;
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
#ifdef PRINTALGO
      std::cout<<"algo 28: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
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
      AddVV_Helper<-1,s,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
      MultUD_Helper<22,s,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
    }
  };

  // algo -1: Determine which algorithm to use
  template <int s, bool add, int ix, class T, class M1, class M2, class M3>
  struct MultUD_Helper<-1,s,add,ix,T,M1,M2,M3>
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      // Possible algorithms to choose from:
      //
      // Trivial:
      //  0 = s == 0, so nothing to do
      //  1 = s == 1: reduces to trivial scalar product
      //
      // UpperTri:
      // 11 = column major
      // 12 = row major
      // 17 = m1 is unitdiag
      // 18 = row major, ix==0 or m2.step() != 1, so copy m2
      //
      // LowerTri:
      // 21 = column major
      // 22 = row major
      // 27 = m1 is unitdiag
      // 28 = row major, ix==0 or m2.step() != 1, so copy m2

      enum { rm = M1::mrowmajor && M3::mrowmajor };
      enum { algo = (
          ( s == 0 ) ? 0 :
          ( s == 1 ) ? 1 :
          ShapeTraits<M3::mshape>::upper ? ( // UpperTri
            M1::munit ? 17 :
            rm ? ( 
#if TMV_OPT >= 1
              (TMV_ZeroIX || M2::mstep != 1) ? 18 :
#endif
              12 ) :
            11 ) :
          ( // LowerTri
            M1::munit ? 27 :
            rm ? ( 
#if TMV_OPT >= 1
              (TMV_ZeroIX || M2::mstep != 1) ? 28 :
#endif
              22 ) :
            21 ) ) };
#ifdef PRINTALGO
      std::cout<<"InlineMultUD: \n";
      std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
      std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
      MultUD_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
    }
  };

  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void InlineMultMM(const Scaling<ix,T>& x, 
      const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
      BaseMatrix_Tri_Mutable<M3>& m3)
  {
    //std::cout<<"InlineMultMD "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"add = "<<add<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3<<std::endl;
    TMVStaticAssert(!M3::munit);
    TMVStaticAssert(ShapeTraits<M1::mshape>::upper == 
        int(ShapeTraits<M3::mshape>::upper));
    TMVStaticAssert((Sizes<M3::msize,M1::msize>::same));
    TMVStaticAssert((Sizes<M2::msize,M1::msize>::same));
    TMVStaticAssert((Sizes<M3::msize,M2::msize>::same));
    TMVAssert(m3.size() == m1.size());
    TMVAssert(m3.size() == m2.size());
    const int s = Sizes<Sizes<M3::msize,M1::msize>::size,M2::msize>::size;
    typedef typename M1::const_view_type M1v;
    typedef typename M2::const_view_type M2v;
    typedef typename M3::view_type M3v;
    M1v m1v = m1.View();
    M2v m2v = m2.View();
    M3v m3v = m3.View();
    MultUD_Helper<-1,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    //std::cout<<"End InlineMultMD: m3 = "<<m3<<std::endl;
  }

  // D * M
  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void InlineMultMM(const Scaling<ix,T>& x,
      const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Tri<M2>& m2,
      BaseMatrix_Mutable<M3>& m3)
  {
    typename M3::transpose_type m3t = m3.Transpose();
    InlineMultMM<add>(x,m2.Transpose(),m1.Transpose(),m3t);
  }


  // Defined in TMV_MultUD.cpp
  template <class T1, DiagType D2, bool C1, class T2, bool C2, class T3>
  void InstMultMM(const T3 x,
      const ConstUpperTriMatrixView<T1,D2,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, UpperTriMatrixView<T3> m3);
  template <class T1, DiagType D2, bool C1, class T2, bool C2, class T3>
  void InstAddMultMM(const T3 x,
      const ConstUpperTriMatrixView<T1,D2,UNKNOWN,UNKNOWN,C1>& m1,
      const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, UpperTriMatrixView<T3> m3);

  template <class T1, DiagType D2, bool C1, class T2, bool C2, class T3>
  inline void InstMultMM(const T3 x,
      const ConstDiagMatrixView<T1,UNKNOWN,C1>& m1,
      const ConstUpperTriMatrixView<T2,D2,UNKNOWN,UNKNOWN,C2>& m2,
      UpperTriMatrixView<T3> m3)
  { InstMultMM(x,m2.Transpose(),m1.Transpose(),m3.Transpose()); }
  template <class T1, DiagType D2, bool C1, class T2, bool C2, class T3>
  inline void InstAddMultMM(const T3 x,
      const ConstDiagMatrixView<T1,UNKNOWN,C1>& m1,
      const ConstUpperTriMatrixView<T2,D2,UNKNOWN,UNKNOWN,C2>& m2,
      UpperTriMatrixView<T3> m3)
  { InstAddMultMM(x,m2.Transpose(),m1.Transpose(),m3.Transpose()); }

  template <class M1, int ix, class T, class M2>
  inline void MultEqMM(BaseMatrix_Tri_Mutable<M1>& m1,
      const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
  { MultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }


} // namespace tmv

#undef TMV_ZeroIX

#endif 
