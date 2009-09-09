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


#ifndef TMV_ProdXV_H
#define TMV_ProdXV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

  //
  // Vector *= Scalar
  //

  template <bool inst, int ix, class T, class V>
  struct CallMultXV // inst = false, ix != 1
  {
    static inline void call(const Scaling<ix,T>& x, V& v)
    { InlineMultXV(x,v); }
  };
  template <int ix, class T, class V>
  struct CallMultXV<true,ix,T,V> // inst = true
  {
    static inline void call(const Scaling<ix,T>& x, V& v)
    {
      typedef typename V::value_type T1;
      InstMultXV(T1(x),v.XView()); 
    }
  };
  template <bool inst, class T, class V>
  struct CallMultXV<inst,1,T,V> // ix == 1: don't do anything
  { static inline void call(const Scaling<1,T>& , V& ) { } };

  template <int ix, class T, class V>
  inline void MultXV(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v)
  {
    typedef typename V::value_type T1;
    enum { inst = (
        Traits<T>::isinst &&
        Traits<T1>::isinst &&
        Traits2<T,T1>::samebase &&
        V::vsize == UNKNOWN) };
    return CallMultXV<inst,ix,T,V>::call(x,v.vec());
  }


  //
  // Vector * / Scalar
  // -Vector
  //

  // Define in TMV_MultXV.h
  template <int ix, class T, class V1, class V2>
  inline void InlineMultXV(const Scaling<ix,T>& x, 
      const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);

  // Defined in TMV_MultXV.cpp
  template <class T1, bool C1, class T2>
  void InstMultXV(const T2 x,
      const ConstVectorView<T1,UNKNOWN,C1>& v1, VectorView<T2> v2); 

  template <class V1, class T2>
  inline void InstMultXV(const T2 x, const V1& v1,
      VectorView<T2,UNKNOWN,true> v2)
  { InstMultXV(TMV_CONJ(x),v1.Conjugate(),v2.Conjugate()); }


  template <bool inst, int ix, class T, class V1, class V2>
  struct CallMultXV2 // inst = false
  {
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    { InlineMultXV(x,v1,v2); }
  };
  template <int ix, class T, class V1, class V2>
  struct CallMultXV2<true,ix,T,V1,V2> // inst = true
  {
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      typedef typename V2::value_type T2;
      InstMultXV(T2(x),v1.XView(),v2.XView()); 
    }
  };

  template <int ix, class T, class V1, class V2>
  inline void MultXV(const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
      BaseVector_Mutable<V2>& v2)
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    enum { inst = (
        Traits<T>::isinst &&
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits2<T,T1>::samebase &&
        Traits2<T1,T2>::sametype &&
        V1::vsize == UNKNOWN &&
        V2::vsize == UNKNOWN) };
    CallMultXV2<inst,ix,T,V1,V2>::call(x,v1.vec(),v2.vec());
  }

  template <int ix, class T, class V>
  class ProdXV;

  template <int ix, class T, class V>
  struct Traits<ProdXV<ix,T,V> >
  {
    typedef typename Traits2<T,typename V::value_type>::type value_type;
    enum { vsize = V::vsize };
    enum { vfort = V::vfort };
    enum { vcalc = false };
    typedef ProdXV<ix,T,V> type;
    typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef typename TypeSelect<V::vcalc,const type,calc_type>::type eval_type;
  };

  template <int ix, class T, class V>
  class ProdXV : 
    public BaseVector<ProdXV<ix,T,V> >
  {
  public:

    typedef ProdXV<ix,T,V> type;
    typedef typename Traits<type>::value_type value_type;
    enum { vsize = Traits<type>::vsize };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    inline ProdXV(const T _x, const BaseVector<V>& _v) : 
      x(_x), v(_v.vec()) {}
    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const V& GetV() const { return v; }

    inline size_t size() const { return v.size(); }
    inline value_type cref(int i) const 
    { return  x*v.cref(i); }

    template <class V2>
    inline void AssignTo(BaseVector_Mutable<V2>& v2) const
    {
      TMVStaticAssert((Sizes<V::vsize,V2::vsize>::same));
      TMVAssert(size() == v2.size());
      TMVStaticAssert(visreal || V2::viscomplex);
      typename V2::cview_type v2cv = v2.CView();
      MultXV(x,v.calc().CView(),v2cv);
    }
  private:
    const Scaling<ix,T> x;
    const V& v;
  };

#define RT typename V::real_type
#define CT typename V::complex_type
#define CCT ConjRef<CT>

  // v *= x
  template <class V>
  inline void MultEq(BaseVector_Mutable<V>& v, const int x)
  {
    typename V::cview_type vcv = v.CView();
    MultXV(Scaling<0,RT>(RT(x)),vcv); 
  }

  template <class V>
  inline void MultEq(BaseVector_Mutable<V>& v, const RT x)
  {
    typename V::cview_type vcv = v.CView();
    MultXV(Scaling<0,RT>(x),vcv); 
  }

  template <class V>
  inline void MultEq(BaseVector_Mutable<V>& v, const CT x)
  {
    typename V::cview_type vcv = v.CView();
    MultXV(Scaling<0,CT>(x),vcv); 
  }

  template <class V>
  inline void MultEq(BaseVector_Mutable<V>& v, const CCT x)
  {
    typename V::cview_type vcv = v.CView();
    MultXV(Scaling<0,CT>(CT(x)),vcv); 
  }
 
  template <class V, class T>
  inline void MultEq(BaseVector_Mutable<V>& v, const Scaling<0,T> x)
  {
    typename V::cview_type vcv = v.CView();
    MultXV(x,vcv); 
  }
 
  template <class V, class T>
  inline void MultEq(BaseVector_Mutable<V>& v, const Scaling<1,T> x)
  {}
 
  template <class V, class T>
  inline void MultEq(BaseVector_Mutable<V>& v, const Scaling<-1,T> x)
  {
    typename V::cview_type vcv = v.CView();
    MultXV(x,vcv); 
  }
 
  // v /= x
  template <class V>
  inline void DivEq(BaseVector_Mutable<V>& v, const int x)
  {
    typename V::cview_type vcv = v.CView();
    MultEq(vcv,RT(1)/RT(x)); 
  }

  template <class V>
  inline void DivEq(BaseVector_Mutable<V>& v, const RT x)
  {
    typename V::cview_type vcv = v.CView();
    MultEq(vcv,RT(1)/x); 
  }

  template <class V>
  inline void DivEq(BaseVector_Mutable<V>& v, const CT x)
  {
    typename V::cview_type vcv = v.CView();
    MultEq(vcv,RT(1)/x); 
  }

  template <class V>
  inline void DivEq(BaseVector_Mutable<V>& v, const CCT x)
  {
    typename V::cview_type vcv = v.CView();
    MultEq(vcv,RT(1)/CT(x)); 
  }

  template <class V, class T>
  inline void DivEq(BaseVector_Mutable<V>& v, const Scaling<0,T> x)
  {
    typename V::cview_type vcv = v.CView();
    MultEq(vcv,RT(1)/T(x)); 
  }
 
  template <class V, class T>
  inline void DivEq(BaseVector_Mutable<V>& v, const Scaling<1,T> x)
  {}
 
  template <class V, class T>
  inline void DivEq(BaseVector_Mutable<V>& v, const Scaling<-1,T> x)
  {
    typename V::cview_type vcv = v.CView();
    MultEq(vcv,x); 
  }
 
  // -v
  template <class V> 
  inline ProdXV<-1,RT,V> operator-(const BaseVector<V>& v)
  { return ProdXV<-1,RT,V>(RT(-1),v); }

  // x * v
  template <class V> 
  inline ProdXV<0,RT,V> operator*(const int x, const BaseVector<V>& v)
  { return ProdXV<0,RT,V>(RT(x),v); }

  template <class V> 
  inline ProdXV<0,RT,V> operator*(const RT x, const BaseVector<V>& v)
  { return ProdXV<0,RT,V>(x,v); }

  template <class V> 
  inline ProdXV<0,CT,V> operator*(const CT x, const BaseVector<V>& v)
  { return ProdXV<0,CT,V>(x,v); }

  template <class V> 
  inline ProdXV<0,CT,V> operator*(const CCT x, const BaseVector<V>& v)
  { return CT(x)*v; }

  template <class V, int ix, class T> 
  inline ProdXV<ix,T,V> operator*(const Scaling<ix,T> x, const BaseVector<V>& v)
  { return ProdXV<ix,T,V>(T(x),v); }

  // v * x
  template <class V> 
  inline ProdXV<0,RT,V> operator*(const BaseVector<V>& v, const int x)
  { return RT(x)*v; }

  template <class V> 
  inline ProdXV<0,RT,V> operator*(const BaseVector<V>& v, const RT x)
  { return x*v; }

  template <class V> 
  inline ProdXV<0,CT,V> operator*(const BaseVector<V>& v, const CT x)
  { return x*v; }

  template <class V> 
  inline ProdXV<0,CT,V> operator*(const BaseVector<V>& v, const CCT x)
  { return CT(x)*v; }

  template <class V, int ix, class T> 
  inline ProdXV<ix,T,V> operator*(const BaseVector<V>& v, const Scaling<ix,T> x)
  { return ProdXV<ix,T,V>(T(x),v); }

  // v / x
  template <class V> 
  inline ProdXV<0,RT,V> operator/(const BaseVector<V>& v, const int x)
  { return (RT(1)/RT(x))*v; }

  template <class V> 
  inline ProdXV<0,RT,V> operator/(const BaseVector<V>& v, const RT x)
  { return (RT(1)/x)*v; }

  template <class V> 
  inline ProdXV<0,CT,V> operator/(const BaseVector<V>& v, const CT x)
  { return (RT(1)/x)*v; }

  template <class V> 
  inline ProdXV<0,CT,V> operator/(const BaseVector<V>& v, const CCT x)
  { return (RT(1)/CT(x))*v; }

  template <class V, int ix, class T> 
  inline ProdXV<ix,T,V> operator/(const BaseVector<V>& v, const Scaling<ix,T> x)
  { return ProdXV<ix,T,V>(RT(1)/T(x),v); }

#undef RT
#undef CT
#undef CCT

  // Consolidate x*x*v type constructs:

#define RT typename ProdXV<ix,T,V>::real_type
#define CT typename ProdXV<ix,T,V>::complex_type
#define CCT ConjRef<CT>

  // -(x*v)
  template <int ix, class T, class V>
  inline ProdXV<-ix,T,V> operator-(const ProdXV<ix,T,V>& pxv)
  { return ProdXV<-ix,T,V>(-pxv.GetX(),pxv.GetV()); }

  // x * (x*v)
  template <int ix, class T, class V>
  inline ProdXV<0,T,V> operator*(const int x, const ProdXV<ix,T,V>& pxv)
  { return ProdXV<0,T,V>(RT(x)*pxv.GetX(),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,T,V> operator*(const RT x, const ProdXV<ix,T,V>& pxv)
  { return ProdXV<0,T,V>(x*pxv.GetX(),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,CT,V> operator*(const CT x, const ProdXV<ix,T,V>& pxv)
  { return ProdXV<0,CT,V>(x*pxv.GetX(),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,CT,V> operator*(const CCT x, const ProdXV<ix,T,V>& pxv)
  { return ProdXV<0,CT,V>(x*pxv.GetX(),pxv.GetV()); }

  template <int ix1, class T1, int ix, class T, class V>
  inline ProdXV<ix1*ix,typename Traits2<T,T1>::type,V> operator*(
      const Scaling<ix1,T1>& x, const ProdXV<ix,T,V>& pxv)
  { 
    return ProdXV<ix1*ix,typename Traits2<T,T1>::type,V>(
        T1(x)*pxv.GetX(),pxv.GetV()); 
  }

  // (x*v)*x
  template <int ix, class T, class V>
  inline ProdXV<0,T,V> operator*(const ProdXV<ix,T,V>& pxv, const int x)
  { return ProdXV<0,T,V>(RT(x)*pxv.GetX(),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,T,V> operator*(const ProdXV<ix,T,V>& pxv, const RT x)
  { return ProdXV<0,T,V>(x*pxv.GetX(),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,CT,V> operator*(const ProdXV<ix,T,V>& pxv, const CT x)
  { return ProdXV<0,CT,V>(x*pxv.GetX(),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,CT,V> operator*(const ProdXV<ix,T,V>& pxv, const CCT x)
  { return ProdXV<0,CT,V>(x*pxv.GetX(),pxv.GetV()); }

  template <int ix1, class T1, int ix, class T, class V>
  inline ProdXV<ix1*ix,typename Traits2<T,T1>::type,V> operator*(
      const ProdXV<ix,T,V>& pxv, const Scaling<ix1,T1>& x)
  {
    return ProdXV<ix1*ix,typename Traits2<T,T1>::type,V>(
        T1(x)*pxv.GetX(),pxv.GetV()); 
  }

  // (x*v)/x
  template <int ix, class T, class V>
  inline ProdXV<0,T,V> operator/(const ProdXV<ix,T,V>& pxv, const int x)
  { return ProdXV<0,T,V>(pxv.GetX()/RT(x),pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,T,V> operator/(const ProdXV<ix,T,V>& pxv, const RT x)
  { return ProdXV<0,T,V>(pxv.GetX()/x,pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,CT,V> operator/(const ProdXV<ix,T,V>& pxv, const CT x)
  { return ProdXV<0,CT,V>(pxv.GetX()/x,pxv.GetV()); }

  template <int ix, class T, class V>
  inline ProdXV<0,CT,V> operator/(const ProdXV<ix,T,V>& pxv, const CCT x)
  { return ProdXV<0,CT,V>(pxv.GetX()/x,pxv.GetV()); }

  template <int ix1, class T1, int ix, class T, class V>
  inline ProdXV<ix1*ix,typename Traits2<T,T1>::type,V> operator/(
      const ProdXV<ix,T,V>& pxv, const Scaling<ix1,T1>& x)
  { 
    return ProdXV<ix1*ix,typename Traits2<T,T1>::type,V>(
      pxv.GetX()/T1(x),pxv.GetV()); 
  }

#undef RT
#undef CT
#undef CCT

#ifndef TMV_NO_ALIAS_CHECK
  // Have SameStorage look into a ProdXV object:
  template <int ix1, class T1, class V1, class V2>
  inline bool SameStorage(
      const ProdXV<ix1,T1,V1>& v1, const BaseVector_Calc<V2>& v2)
  { return SameStorage(v1.GetV(),v2); }
  template <class V1, int ix2, class T2, class V2>
  inline bool SameStorage(
      const BaseVector_Calc<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { return SameStorage(v1,v2.GetV()); }
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline bool SameStorage(
      const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { return SameStorage(v1.GetV(),v2.GetV()); }
#endif


  template <int ix, class T, class V>
  inline std::string TypeText(const ProdXV<ix,T,V>& pxv)
  { 
    std::ostringstream s;
    s << "ProdXV< "<<ix<<","<<TypeText(T(pxv.GetX()));
    s << " , "<<TypeText(pxv.GetV())<<" >";
    return s.str();
  }


} // namespace tmv

#endif 
