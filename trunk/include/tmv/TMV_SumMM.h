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


#ifndef TMV_SumMM_H
#define TMV_SumMM_H

#include "TMV_ProdXM.h"
#include "TMV_SumVV.h"

namespace tmv {

  //
  // Matrix += Matrix
  //

  template <bool checkalias, bool conj, bool rm, bool inst, int ix, class T, class M1, class M2>
  struct CallAddMM;
  template <bool conj, bool rm, bool inst, int ix, class T, class M1, class M2>
  struct CallAddMM<true,conj,rm,inst,ix,T,M1,M2> // checkalias = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { 
      if (!SameStorage(m1,m2)) 
      {
        CallAddMM<false,conj,rm,inst,ix,T,M1,M2>::call(x,m1,m2);
      }
      else 
      {
        if (ExactSameStorage(m1,m2)) {
          if (M1::mconj != int(M2::mconj)) {
            m2 += x * m1.copy();
          }
          else {  // They are initially equal.
            m2 *= T(x) + T(1);
          }
        } else {
          m2 += x * m1.copy();
        }
      }
    }
  };
  template <bool rm, bool inst, int ix, class T, class M1, class M2>
  struct CallAddMM<false,true,rm,inst,ix,T,M1,M2> // conj = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { 
      typedef typename M1::const_conjugate_type M1c;
      typedef typename M2::conjugate_type M2c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      CallAddMM<false,false,rm,inst,ix,T,M1c,M2c>::call(TMV_CONJ(x),m1c,m2c);
    }
  };
  template <bool inst, int ix, class T, class M1, class M2>
  struct CallAddMM<false,false,true,inst,ix,T,M1,M2> // rm = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { 
      typedef typename M1::const_transpose_type M1t;
      typedef typename M2::transpose_type M2t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      CallAddMM<false,false,false,inst,ix,T,M1t,M2t>::call(x,m1t,m2t);
    }
  };
  template <int ix, class T, class M1, class M2>
  struct CallAddMM<false,false,false,true,ix,T,M1,M2> // inst = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M2::value_type T2;
      InstAddMM(T2(x),m1.XView(),m2.XView()); 
    }
  };
  template <int ix, class T, class M1, class M2>
  struct CallAddMM<false,false,false,false,ix,T,M1,M2> // inst = false
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { InlineAddMM(x,m1,m2); }
  };

  template <int ix, class T, class M1, class M2>
  inline void AddMM(
      const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
      BaseMatrix_Mutable<M2>& m2)
  {
    typedef typename M1::value_type TM1;
    typedef typename M2::value_type TM2;
    enum { conj = M2::mconj };
    enum { rm = M2::mrowmajor && !M2::mcolmajor };
    enum { checkalias = (
        M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
        M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN) };
    enum { inst = (
        Traits<TM1>::isinst &&
        Traits<TM2>::isinst &&
        Traits2<T,TM2>::samebase &&
#ifdef TMV_INST_MIX
        Traits2<TM1,TM2>::samebase &&
#else
        Traits2<TM1,TM2>::sametype &&
#endif
        checkalias ) };
    CallAddMM<checkalias,conj,rm,inst,ix,T,M1,M2>::call(x,m1.mat(),m2.mat());
  }


  //
  // Matrix + Matrix
  //

  // Default unless overridden by specialized version
  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void InlineAddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1, 
      const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2, 
      BaseMatrix_Mutable<M3>& m3)
  { (m3 = x1*m1) += (x2*m2); }

  template <bool checkalias, bool conj, bool rm, bool inst, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct CallAddMM2;
  template <bool conj, bool rm, bool inst, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct CallAddMM2<true,conj,rm,inst,ix1,T1,M1,ix2,T2,M2,M3> // checkalias
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1,
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      if (!SameStorage(m1,m3) && !SameStorage(m2,m3)) 
      {
        CallAddMM2<false,conj,rm,inst,ix1,T1,M1,ix2,T2,M2,M3>::call(
            x1,m1,x2,m2,m3);
      }
      else if (!SameStorage(m2,m3)) // alias with m1 only
      {
        (m3 = x1 * m1) += x2 * m2;
      }
      else if (!SameStorage(m1,m3)) // alias with m2 only
      {
        (m3 = x2 * m2) += x1 * m1;
      }
      else // alias with both
      {
        typename M1::copy_type m1c = m1.copy();
        (m3 = x2 * m2) += x1 * m1c;
      }
    }
  };
  template <bool rm, bool inst, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct CallAddMM2<false,true,rm,inst,ix1,T1,M1,ix2,T2,M2,M3> // conj = true
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1,
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    { 
      typedef typename M1::const_conjugate_type M1c;
      typedef typename M2::const_conjugate_type M2c;
      typedef typename M3::conjugate_type M3c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      M3c m3c = m3.Conjugate();
      CallAddMM2<false,false,rm,inst,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
          TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
    }
  };
  template <bool inst, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct CallAddMM2<false,false,true,inst,ix1,T1,M1,ix2,T2,M2,M3> // rm = true
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
      CallAddMM2<false,false,false,inst,ix1,T1,M1t,ix2,T2,M2t,M3t>::call(
          x1,m1t,x2,m2t,m3t);
    }
  };
  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct CallAddMM2<false,false,false,true,ix1,T1,M1,ix2,T2,M2,M3> // inst = true
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1,
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    {
      // Rather than have an Inst version for this, just do it in two steps.
      (m3 = x1*m1) += x2*m2;
    }
  };
  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  struct CallAddMM2<false,false,false,false,ix1,T1,M1,ix2,T2,M2,M3> // inst = false
  {
    static inline void call(
        const Scaling<ix1,T1>& x1, const M1& m1,
        const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
    { InlineAddMM(x1,m1,x2,m2,m3); }
  };

  // m3 = m1 + m2
  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void AddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1, 
      const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2, 
      BaseMatrix_Mutable<M3>& m3)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same));
    TMVStaticAssert((Sizes<M2::mcolsize,M3::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M3::mrowsize>::same));
    TMVStaticAssert((Sizes<M2::mrowsize,M3::mrowsize>::same));
    TMVAssert(m1.colsize() == m3.colsize());
    TMVAssert(m2.colsize() == m3.colsize());
    TMVAssert(m1.rowsize() == m3.rowsize());
    TMVAssert(m2.rowsize() == m3.rowsize());

    typedef typename M1::value_type TM1;
    typedef typename M2::value_type TM2;
    typedef typename M3::value_type TM3;
    enum { conj = M3::mconj };
    enum { rm = M3::mrowmajor && !M3::mcolmajor };
    enum { checkalias = (
        M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
        M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN &&
        M3::mcolsize == UNKNOWN && M3::mrowsize == UNKNOWN) };
    enum { inst = (
        Traits<TM1>::isinst &&
        Traits<TM2>::isinst &&
        Traits<TM3>::isinst &&
        Traits2<T1,TM3>::samebase &&
        Traits2<T2,TM3>::samebase &&
#ifdef TMV_INST_MIX
        Traits2<TM1,TM3>::samebase &&
        Traits2<TM2,TM3>::samebase &&
#else
        Traits2<TM1,TM3>::sametype &&
        Traits2<TM2,TM3>::sametype &&
#endif
        checkalias ) };
    CallAddMM2<checkalias,conj,rm,inst,ix1,T1,M1,ix2,T2,M2,M3>::call(
        x1,m1.mat(),x2,m2.mat(),m3.mat());
  }


  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  class SumMM;

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  struct Traits<SumMM<ix1,T1,M1,ix2,T2,M2> >
  {
    typedef typename ProdXM<ix1,T1,M1>::value_type mtype1;
    typedef typename ProdXM<ix2,T2,M2>::value_type mtype2;
    typedef typename Traits2<mtype1,mtype2>::type value_type;

    enum { mcolsize = Sizes<M1::mcolsize,M2::mcolsize>::size };
    enum { mrowsize = Sizes<M1::mrowsize,M2::mrowsize>::size };
    enum { mshape = ShapeTraits2<M1::mshape,M2::mshape>::sum };
    enum { mfort = M1::mfort && M2::mfort };
    enum { mcalc = false };
    enum { rm1 = Traits<typename M1::calc_type>::mrowmajor };
    enum { rm2 = Traits<typename M2::calc_type>::mrowmajor };
    enum { cm1 = Traits<typename M1::calc_type>::mcolmajor };
    enum { mrowmajor = ( rm1 || (rm2 && !cm1) ) };

    typedef SumMM<ix1,T1,M1,ix2,T2,M2> type;
    typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,mrowmajor,mfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef typename TypeSelect<M1::mcalc&&M2::mcalc,const type,calc_type>::type eval_type;
    typedef InvalidType inverse_type;
  };

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  class SumMM :
    public BaseMatrix<SumMM<ix1,T1,M1,ix2,T2,M2> >
  {
  public:

    typedef SumMM<ix1,T1,M1,ix2,T2,M2> type;
    typedef typename Traits<type>::value_type value_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mshape = Traits<type>::mshape };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    inline SumMM(
        const T1& _x1, const BaseMatrix<M1>& _m1, 
        const T2& _x2, const BaseMatrix<M2>& _m2) :
      x1(_x1), m1(_m1.mat()), x2(_x2), m2(_m2.mat())
    {
      TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
      TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
      TMVAssert(m1.colsize() == m2.colsize());
      TMVAssert(m1.rowsize() == m2.rowsize());
    }

    inline const Scaling<ix1,T1>& GetX1() const { return x1; }
    inline const M1& GetM1() const { return m1; }
    inline const Scaling<ix2,T2>& GetX2() const { return x2; }
    inline const M2& GetM2() const { return m2; }

    inline size_t colsize() const { return m1.colsize(); }
    inline size_t rowsize() const { return m1.rowsize(); }

    inline value_type cref(int i, int j) const
    { return x1 * m1.cref(i,j) + x2 * m2.cref(i,j); }

    template <class M3>
    inline void AssignTo(BaseMatrix_Mutable<M3>& m3) const
    {
      TMVStaticAssert((ShapeTraits2<mshape,M3::mshape>::assignable)); 
      TMVStaticAssert((misreal || M3::miscomplex));
      TMVStaticAssert((Sizes<mcolsize,M3::mcolsize>::same)); 
      TMVStaticAssert((Sizes<mrowsize,M3::mrowsize>::same)); 
      TMVAssert(colsize() == m3.colsize());
      TMVAssert(rowsize() == m3.rowsize());
      typename M3::cview_type m3cv = m3.CView();
      AddMM(x1,m1.calc().CView(),x2,m2.calc().CView(),m3cv);
    }
  private:
    const Scaling<ix1,T1> x1;
    const M1& m1;
    const Scaling<ix2,T2> x2;
    const M2& m2;
  };

#define RT typename M2::real_type
  // m += m
  template <class M1, class M2>
  inline void AddEq(
      BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2) 
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVStaticAssert(M1::miscomplex || M2::misreal);
    typename M1::cview_type m1cv = m1.CView();
    AddMM(Scaling<1,RT>(),m2.calc().CView(),m1cv);
  }

  // m += xm
  template <class M1, int ix2, class T2, class M2>
  inline void AddEq(
      BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2) 
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVStaticAssert(M1::miscomplex || M2::misreal);
    typename M1::cview_type m1cv = m1.CView();
    AddMM(m2.GetX(),m2.GetM().calc().CView(),m1cv);
  }

  // m -= m
  template <class M1, class M2>
  inline void SubtractEq(
      BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVStaticAssert(M1::miscomplex || M2::misreal);
    typename M1::cview_type m1cv = m1.CView();
    AddMM(Scaling<-1,RT>(),m2.calc().CView(),m1cv);
  }

  // m -= xm
  template <class M1, int ix2, class T2, class M2>
  inline void SubtractEq(
      BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2) 
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVStaticAssert(M1::miscomplex || M2::misreal);
    typename M1::cview_type m1cv = m1.CView();
    AddMM(-m2.GetX(),m2.GetM().calc().CView(),m1cv);
  }

  // m + m
  template <class M1, class M2>
  inline SumMM<1,RT,M1,1,RT,M2> operator+(
      const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<1,RT,M1,1,RT,M2>(RT(1),m1,RT(1),m2); 
  }

  // xm + m
  template <int ix1, class T1, class M1, class M2>
  inline SumMM<ix1,T1,M1,1,RT,M2> operator+(
      const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<ix1,T1,M1,1,RT,M2>(T1(m1.GetX()),m1.GetM(),RT(1),m2); 
  }

  // m + xm
  template <class M1, int ix2, class T2, class M2>
  inline SumMM<1,RT,M1,ix2,T2,M2> operator+(
      const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<1,RT,M1,ix2,T2,M2>(RT(1),m1,T2(m2.GetX()),m2.GetM()); 
  }

  // xm + xm
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<ix1,T1,M1,ix2,T2,M2> operator+(
      const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<ix1,T1,M1,ix2,T2,M2>(
        T1(m1.GetX()),m1.GetM(),T2(m2.GetX()),m2.GetM()); 
  }

  // m - m
  template <class M1, class M2>
  inline SumMM<1,RT,M1,-1,RT,M2> operator-(
      const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<1,RT,M1,-1,RT,M2>(RT(1),m1,RT(-1),m2); 
  }

  // xm - m
  template <int ix1, class T1, class M1, class M2>
  inline SumMM<ix1,T1,M1,-1,RT,M2> operator-(
      const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<ix1,T1,M1,-1,RT,M2>(T1(m1.GetX()),m1.GetM(),RT(-1),m2); 
  }

  // m - xm
  template <class M1, int ix2, class T2, class M2>
  inline SumMM<1,RT,M1,-ix2,T2,M2> operator-(
      const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<1,RT,M1,-ix2,T2,M2>(RT(1),m1,-T2(m2.GetX()),m2.GetM()); 
  }

  // xm - xm
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<ix1,T1,M1,-ix2,T2,M2> operator-(
      const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return SumMM<ix1,T1,M1,-ix2,T2,M2>(
        T1(m1.GetX()),m1.GetM(),-T2(m2.GetX()),m2.GetM()); 
  }
#undef RT


  // Consolidate x*(xm+xm) type constructs:

#define RT typename SumMM<ix1,T1,M1,ix2,T2,M2>::real_type
#define CT typename SumMM<ix1,T1,M1,ix2,T2,M2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

  // -(xm+xm)
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<-ix1,T1,M1,-ix2,T2,M2> operator-(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    return SumMM<-ix1,T1,M1,-ix2,T2,M2>(
        -T1(smm.GetX1()),smm.GetM1(),-T2(smm.GetX2()),smm.GetM2());
  }

  // x * (xm+xm)
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,T1,M1,0,T2,M2> operator*(
      const int x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    return SumMM<0,T1,M1,0,T2,M2>(
        RT(x)*smm.GetX1(),smm.GetM1(), RT(x)*smm.GetX2(),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,T1,M1,0,T2,M2> operator*(
      const RT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    return SumMM<0,T1,M1,0,T2,M2>(
        x*smm.GetX1(),smm.GetM1(), x*smm.GetX2(),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,CT,M1,0,CT,M2> operator*(
      const CT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    return SumMM<0,CT,M1,0,CT,M2>(
        x*smm.GetX1(),smm.GetM1(), x*smm.GetX2(),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,CT,M1,0,CT,M2> operator*(
      const CCT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    return SumMM<0,CT,M1,0,CT,M2>(
        CT(x)*smm.GetX1(),smm.GetM1(), CT(x)*smm.GetX2(),smm.GetM2());
  }
  template <int ix, class T, int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2> operator*(
      const Scaling<ix,T>& x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    return SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2>(
        T(x)*smm.GetX1(),smm.GetM1(),T(x)*smm.GetX2(),smm.GetM2());
  }

  // (xm+xm)*x
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,T1,M1,0,T2,M2> operator*(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const int x)
  {
    return SumMM<0,T1,M1,0,T2,M2>(
        RT(x)*smm.GetX1(),smm.GetM1(), RT(x)*smm.GetX2(),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,T1,M1,0,T2,M2> operator*(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const RT x)
  {
    return SumMM<0,T1,M1,0,T2,M2>(
        x*smm.GetX1(),smm.GetM1(), x*smm.GetX2(),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,CT,M1,0,CT,M2> operator*(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CT x)
  {
    return SumMM<0,CT,M1,0,CT,M2>(
        x*smm.GetX1(),smm.GetM1(), x*smm.GetX2(),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,CT,M1,0,CT,M2> operator*(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CCT x)
  {
    return SumMM<0,CT,M1,0,CT,M2>(
        CT(x)*smm.GetX1(),smm.GetM1(), CT(x)*smm.GetX2(),smm.GetM2());
  }

  template <int ix, class T, int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2> operator*(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const Scaling<ix,T>& x)
  {
    return SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2>(
        T(x)*smm.GetX1(),smm.GetM1(),T(x)*smm.GetX2(),smm.GetM2());
  }

  // (xm+xm)/x
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,T1,M1,0,T2,M2> operator/(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const int x)
  {
    return SumMM<0,T1,M1,0,T2,M2>(
        smm.GetX1()/RT(x),smm.GetM1(), smm.GetX2()/RT(x),smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,T1,M1,0,T2,M2> operator/(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const RT x)
  {
    return SumMM<0,T1,M1,0,T2,M2>(
        smm.GetX1()/x,smm.GetM1(), smm.GetX2()/x,smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,CT,M1,0,CT,M2> operator/(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CT x)
  {
    return SumMM<0,CT,M1,0,CT,M2>(
        smm.GetX1()/x,smm.GetM1(), smm.GetX2()/x,smm.GetM2());
  }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<0,CT,M1,0,CT,M2> operator/(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CCT x)
  {
    return SumMM<0,CT,M1,0,CT,M2>(
        smm.GetX1()/CT(x),smm.GetM1(), smm.GetX2()/CT(x),smm.GetM2());
  }

  template <int ix, class T, int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2> operator/(
      const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const Scaling<ix,T>& x)
  {
    return SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2>(
        smm.GetX1()/T(x),smm.GetM1(),smm.GetX2()/T(x),smm.GetM2());
  }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



  // TypeText

  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline std::string TypeText(const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
  {
    std::ostringstream s;
    s << "SumMM< "<< ix1<<","<<TypeText(T1())<<" , "<<TypeText(smm.GetM1());
    s << " , "<<ix2<<","<<TypeText(T2())<<" , "<<TypeText(smm.GetM2())<<" >";
    return s.str();
  }



} // namespace tmv

#endif 
