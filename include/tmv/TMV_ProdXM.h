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

#ifndef TMV_ProdXM_H
#define TMV_ProdXM_H

#include "TMV_ProdXV.h"
#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

  //
  // Matrix *= Scalar 
  //

  template <bool conj, bool rm, bool inst, int ix, class T, class M>
  struct CallMultXM;

  template <bool conj, bool rm, bool inst, class T, class M>
  struct CallMultXM<conj,rm,inst,1,T,M> // ix == 1: don't do anything
  { static inline void call(const Scaling<1,T>& , M& ) {} };
  template <bool rm, bool inst, int ix, class T, class M>
  struct CallMultXM<true,rm,inst,ix,T,M> // conj = true
  {
    static inline void call(const Scaling<ix,T>& x, M& m) 
    {
      typedef typename M::conjugate_type Mc;
      Mc mc = m.Conjugate();
      CallMultXM<false,rm,inst,ix,T,Mc>::call(TMV_CONJ(x),mc);
    }
  };
  template <bool inst, int ix, class T, class M>
  struct CallMultXM<false,true,inst,ix,T,M> // rm = true
  {
    static inline void call(const Scaling<ix,T>& x, M& m) 
    {
      typedef typename M::transpose_type Mt;
      Mt mt = m.Transpose();
      CallMultXM<false,false,inst,ix,T,Mt>::call(x,mt);
    }
  };
  template <int ix, class T, class M>
  struct CallMultXM<false,false,true,ix,T,M> // inst = true
  {
    static inline void call(const Scaling<ix,T>& x, M& m)
    {
      typedef typename M::value_type T1;
      InstMultXM(T1(x),m.XView()); 
    }
  };
  template <int ix, class T, class M>
  struct CallMultXM<false,false,false,ix,T,M> // inst = false
  {
    static inline void call(const Scaling<ix,T>& x, M& m) 
    { InlineMultXM(x,m); } 
  };

  template <int ix, class T, class M>
  inline void MultXM(const Scaling<ix,T>& x, BaseMatrix_Mutable<M>& m)
  {
    typedef typename M::value_type T1;
    enum { conj = M::mconj };
    enum { rm = M::mrowmajor && !M::mcolmajor };
    enum { inst = (
        Traits<T>::isinst &&
        Traits<T1>::isinst &&
        Traits2<T,T1>::samebase &&
        M::mcolsize == UNKNOWN &&
        M::mrowsize == UNKNOWN) };
    return CallMultXM<conj,rm,inst,ix,T,M>::call(x,m.mat());
  }

  //
  // Matrix * Scalar
  // -Matrix
  //

  template <bool checkalias, bool conj, bool rm, bool inst, int ix, class T, class M1, class M2>
  struct CallMultXM2;
  template <bool conj, bool rm, bool inst, int ix, class T, class M1, class M2>
  struct CallMultXM2<true,conj,rm,inst,ix,T,M1,M2>  // checkalias = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      if (!SameStorage(m1,m2)) 
      {
        CallMultXM2<false,conj,rm,inst,ix,T,M1,M2>::call(x,m1,m2);
      }
      else 
      {
        (m2 = m1) *= x;
      }
    }
  };
  template <bool rm, bool inst, int ix, class T, class M1, class M2>
  struct CallMultXM2<false,true,rm,inst,ix,T,M1,M2> // conj = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M1::const_conjugate_type M1c;
      typedef typename M2::conjugate_type M2c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      CallMultXM2<false,false,rm,inst,ix,T,M1c,M2c>::call(TMV_CONJ(x),m1c,m2c);
    }
  };
  template <bool inst, int ix, class T, class M1, class M2>
  struct CallMultXM2<false,false,true,inst,ix,T,M1,M2> // rm = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M1::const_transpose_type M1t;
      typedef typename M2::transpose_type M2t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      CallMultXM2<false,false,false,inst,ix,T,M1t,M2t>::call(x,m1t,m2t);
    }
  };
  template <int ix, class T, class M1, class M2>
  struct CallMultXM2<false,false,false,true,ix,T,M1,M2> // inst = true
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      typedef typename M2::value_type T2;
      InstMultXM(T2(x),m1.XView(),m2.XView()); 
    }
  };
  template <int ix, class T, class M1, class M2>
  struct CallMultXM2<false,false,false,false,ix,T,M1,M2> // inst = false
  {
    static inline void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { InlineMultXM(x,m1,m2); }
  };

  template <int ix, class T, class M1, class M2>
  inline void MultXM(
      const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
      BaseMatrix_Mutable<M2>& m2)
  {
    typedef typename M1::value_type T1;
    typedef typename M2::value_type T2;
    enum { checkalias = (
        M1::mcolsize == UNKNOWN && M2::mcolsize == UNKNOWN &&
        M1::mrowsize == UNKNOWN && M2::mrowsize == UNKNOWN ) };
    enum { conj = M2::mconj };
    enum { rm = M2::mrowmajor && !M2::mcolmajor };
    enum { inst = (
        Traits<T>::isinst &&
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits2<T,T2>::samebase &&
        Traits2<T1,T2>::samebase &&
        checkalias ) };
    CallMultXM2<checkalias,conj,rm,inst,ix,T,M1,M2>::call(x,m1.mat(),m2.mat());
  }


  template <int ix, class T, class M>
  class ProdXM;

  template <int ix, class T, class M>
  struct Traits<ProdXM<ix,T,M> >
  {
    typedef typename Traits2<T,typename M::value_type>::type value_type;

    enum { mcolsize = M::mcolsize };
    enum { mrowsize = M::mrowsize };
    enum { mshape = ShapeTraits<M::mshape>::nonunit_shape };
    enum { mfort = M::mfort };
    enum { mcalc = false };
    enum { rm1 = Traits<typename M::calc_type>::mrowmajor };
    enum { mrowmajor = rm1 };

    typedef ProdXM<ix,T,M> type;
    typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,mrowmajor,mfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef type eval_type;
    typedef InvalidType inverse_type;
  };

  template <int ix, class T, class M>
  class ProdXM :
    public BaseMatrix<ProdXM<ix,T,M> >
  {
  public:

    typedef ProdXM<ix,T,M> type;
    typedef typename Traits<type>::value_type value_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mshape = Traits<type>::mshape };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    inline ProdXM(const T _x, const BaseMatrix<M>& _m) : 
      x(_x), m(_m.mat()) {}
    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const M& GetM() const { return m; }

    inline size_t colsize() const { return m.colsize(); }
    inline size_t rowsize() const { return m.rowsize(); }

    inline value_type cref(int i, int j) const
    { return x * m.cref(i,j); }

    template <class M2>
    inline void AssignTo(BaseMatrix_Mutable<M2>& m2) const
    {
      TMVStaticAssert((ShapeTraits2<mshape,M2::mshape>::assignable)); 
      TMVStaticAssert((Sizes<mcolsize,M2::mcolsize>::same));
      TMVStaticAssert((Sizes<mrowsize,M2::mrowsize>::same));
      TMVAssert(m2.colsize() == colsize());
      TMVAssert(m2.rowsize() == rowsize());
      TMVStaticAssert(misreal || M2::miscomplex);
      typename M2::cview_type m2cv = m2.CView();
      MultXM(x,m.calc().CView(),m2cv);
    }

  private:
    const Scaling<ix,T> x;
    const M& m;
  };

#define RT typename M::real_type
#define CT typename M::complex_type
#define CCT ConjRef<CT>

  // m *= x
  template <class M>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const int x)
  { 
    typename M::cview_type mcv = m.CView();
    MultXM(Scaling<0,RT>(RT(x)),mcv); 
  }

  template <class M>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const RT x)
  {
    typename M::cview_type mcv = m.CView();
    MultXM(Scaling<0,RT>(x),mcv); 
  }

  template <class M>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const CT x)
  {
    typename M::cview_type mcv = m.CView();
    MultXM(Scaling<0,CT>(x),mcv); 
  }

  template <class M>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const CCT x)
  {
    typename M::cview_type mcv = m.CView();
    MultXM(Scaling<0,CT>(CT(x)),mcv); 
  }

  template <class M, class T>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const Scaling<0,T>& x)
  {
    typename M::cview_type mcv = m.CView();
    MultXM(x,mcv); 
  }

  template <class M, class T>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const Scaling<1,T>& x)
  {}

  template <class M, class T>
  inline void MultEq(BaseMatrix_Mutable<M>& m, const Scaling<-1,T>& x)
  {
    typename M::cview_type mcv = m.CView();
    MultXM(x,mcv); 
  }

  // m /= x
  template <class M>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const int x)
  { MultEq(m,RT(1)/RT(x)); }

  template <class M>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const RT x)
  { MultEq(m,RT(1)/x); }

  template <class M>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const CT x)
  { MultEq(m,RT(1)/x); }

  template <class M>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const CCT x)
  { MultEq(m,RT(1)/CT(x)); }

  template <class M, class T>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const Scaling<0,T>& x)
  { MultEq(m,RT(1)/T(x)); }

  template <class M, class T>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const Scaling<1,T>& x)
  {}

  template <class M, class T>
  inline void DivEq(BaseMatrix_Mutable<M>& m, const Scaling<-1,T>& x)
  { MultEq(m,x); }

  // -m
  template <class M>
  inline ProdXM<-1,RT,M> operator-(const BaseMatrix<M>& m)
  { return ProdXM<-1,RT,M>(RT(-1),m); }

  // x * m
  template <class M>
  inline ProdXM<0,RT,M> operator*(const int x, const BaseMatrix<M>& m)
  { return ProdXM<0,RT,M>(RT(x),m); }

  template <class M>
  inline ProdXM<0,RT,M> operator*(const RT x, const BaseMatrix<M>& m)
  { return ProdXM<0,RT,M>(x,m); }

  template <class M>
  inline ProdXM<0,CT,M> operator*(const CT x, const BaseMatrix<M>& m)
  { return ProdXM<0,CT,M>(x,m); }

  template <class M>
  inline ProdXM<0,CT,M> operator*(const CCT x, const BaseMatrix<M>& m)
  { return CT(x)*m; }

  template <class M, int ix, class T>
  inline ProdXM<ix,T,M> operator*(const Scaling<ix,T>& x, 
      const BaseMatrix<M>& m)
  { return ProdXM<ix,T,M>(T(x),m); }

  // m * x
  template <class M>
  inline ProdXM<0,RT,M> operator*(const BaseMatrix<M>& m, const int x)
  { return RT(x)*m; }

  template <class M>
  inline ProdXM<0,RT,M> operator*(const BaseMatrix<M>& m, const RT x)
  { return x*m; }

  template <class M>
  inline ProdXM<0,CT,M> operator*(const BaseMatrix<M>& m, const CT x)
  { return x*m; }

  template <class M>
  inline ProdXM<0,CT,M> operator*(const BaseMatrix<M>& m, const CCT x)
  { return CT(x)*m; }

  template <class M, int ix, class T>
  inline ProdXM<ix,T,M> operator*(const BaseMatrix<M>& m,
      const Scaling<ix,T>& x)
  { return ProdXM<ix,T,M>(T(x),m); }

  // m / x
  template <class M>
  inline ProdXM<0,RT,M> operator/(const BaseMatrix<M>& m, const int x)
  { return (RT(1)/RT(x))*m; }

  template <class M>
  inline ProdXM<0,RT,M> operator/(const BaseMatrix<M>& m, const RT x)
  { return (RT(1)/x)*m; }

  template <class M>
  inline ProdXM<0,CT,M> operator/(const BaseMatrix<M>& m, const CT x)
  { return (RT(1)/x)*m; }

  template <class M>
  inline ProdXM<0,CT,M> operator/(const BaseMatrix<M>& m, const CCT x)
  { return (RT(1)/CT(x))*m; }

  template <class M, int ix, class T>
  inline ProdXM<ix,T,M> operator/(const BaseMatrix<M>& m,
      const Scaling<ix,T>& x)
  { return ProdXM<ix,T,M>(RT(1)/T(x),m); }

#undef RT
#undef CT
#undef CCT

  // Consolidate x*x*m type constructs:

#define RT typename ProdXM<ix,T,M>::real_type
#define CT typename ProdXM<ix,T,M>::complex_type
#define CCT ConjRef<CT>

  // -(x*m)
  template <int ix, class T, class M>
  inline ProdXM<-ix,T,M> operator-(const ProdXM<ix,T,M>& pxm)
  { return ProdXM<-ix,T,M>(-pxm.GetX(),pxm.GetM()); }

  // x*(x*m)
  template <int ix, class T, class M>
  inline ProdXM<0,T,M> operator*(const int x, const ProdXM<ix,T,M>& pxm)
  { return ProdXM<0,T,M>(RT(x)*pxm.GetX(),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,T,M> operator*(const RT x, const ProdXM<ix,T,M>& pxm)
  { return ProdXM<0,T,M>(x*pxm.GetX(),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,CT,M> operator*(const CT x, const ProdXM<ix,T,M>& pxm)
  { return ProdXM<0,CT,M>(x*pxm.GetX(),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,CT,M> operator*(const CCT x, const ProdXM<ix,T,M>& pxm)
  { return ProdXM<0,CT,M>(x*pxm.GetX(),pxm.GetM()); }

  template <int ix1, class T1, int ix, class T, class M>
  inline ProdXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
      const Scaling<ix1,T1>& x, const ProdXM<ix,T,M>& pxm)
  {
    return ProdXM<ix*ix1,typename Traits2<T1,T>::type,M>(
        T1(x)*pxm.GetX(),pxm.GetM()); 
  }

  // (x*m)*x
  template <int ix, class T, class M>
  inline ProdXM<0,T,M> operator*(const ProdXM<ix,T,M>& pxm, const int x)
  { return ProdXM<0,T,M>(RT(x)*pxm.GetX(),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,T,M> operator*(const ProdXM<ix,T,M>& pxm, const RT x)
  { return ProdXM<0,T,M>(x*pxm.GetX(),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,CT,M> operator*(const ProdXM<ix,T,M>& pxm, const CT x)
  { return ProdXM<0,CT,M>(x*pxm.GetX(),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,CT,M> operator*(const ProdXM<ix,T,M>& pxm, const CCT x)
  { return ProdXM<0,CT,M>(x*pxm.GetX(),pxm.GetM()); }

  template <int ix1, class T1, int ix, class T, class M>
  inline ProdXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
      const ProdXM<ix,T,M>& pxm, const Scaling<ix1,T1>& x)
  {
    return ProdXM<ix*ix1,typename Traits2<T1,T>::type,M>(
        T1(x)*pxm.GetX(),pxm.GetM()); 
  }

  // (x*m)/x
  template <int ix, class T, class M>
  inline ProdXM<0,RT,M> operator/(const ProdXM<ix,T,M>& pxm, const int x)
  { return ProdXM<0,RT,M>(pxm.GetX()/RT(x),pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,RT,M> operator/(const ProdXM<ix,T,M>& pxm, const RT x)
  { return ProdXM<0,RT,M>(pxm.GetX()/x,pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,CT,M> operator/(const ProdXM<ix,T,M>& pxm, const CT x)
  { return ProdXM<0,CT,M>(pxm.GetX()/x,pxm.GetM()); }

  template <int ix, class T, class M>
  inline ProdXM<0,CT,M> operator/(const ProdXM<ix,T,M>& pxm, const CCT x)
  { return ProdXM<0,CT,M>(pxm.GetX()/x,pxm.GetM()); }

  template <int ix1, class T1, int ix, class T, class M>
  inline ProdXM<ix*ix1,typename Traits2<T1,T>::type,M> operator/(
      const ProdXM<ix,T,M>& pxm, const Scaling<ix1,T1>& x)
  { 
    return ProdXM<ix*ix1,typename Traits2<T1,T>::type,M>(
        pxm.GetX()/T1(x),pxm.GetM()); 
  }

#undef RT
#undef CT
#undef CCT

#ifndef TMV_NO_ALIAS_CHECK
  // Have SameStorage look into a ProdXM object:
  template <int ix1, class T1, class M1, class M2>
  inline bool SameStorage(
      const ProdXM<ix1,T1,M1>& m1, const BaseMatrix_Calc<M2>& m2)
  { return SameStorage(m1.GetM(),m2); }
  template <class M1, int ix2, class T2, class M2>
  inline bool SameStorage(
      const BaseMatrix_Calc<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  { return SameStorage(m1,m2.GetM()); }
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline bool SameStorage(
      const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  { return SameStorage(m1.GetM(),m2.GetM()); }
#endif
  template <int ix, class T, class M>
  inline std::string TypeText(const ProdXM<ix,T,M>& pxm)
  {
    std::ostringstream s;
    s << "ProdXM< "<<ix<<","<<TypeText(T(pxm.GetX()));
    s << " , "<<TypeText(pxm.GetM())<<" >";
    return s.str();
  }

} // namespace tmv

#endif
