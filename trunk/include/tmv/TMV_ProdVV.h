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


#ifndef TMV_ProdVV_H
#define TMV_ProdVV_H

#include "TMV_ProdXV.h"

namespace tmv {

  //
  // Vector * Vector
  //

#define PT typename Traits2<typename V1::value_type,typename V2::value_type>::type

  template <bool inst, class V1, class V2>
  struct CallMultVV // inst = false
  {
    static inline PT call(const V1& v1, const V2& v2)
    { return InlineMultVV(v1,v2); }
  };
  template <class V1, class V2>
  struct CallMultVV<true,V1,V2> // inst = true
  {
    static inline PT call(const V1& v1, const V2& v2)
    { return InstMultVV(v1.XView(),v2.XView()); }
  };

  template <class V1, class V2> 
  inline PT MultVV(
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
  {
    typedef typename V1::value_type T1;
    typedef typename V1::value_type T2;
    enum { inst = (
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
#ifdef TMV_INST_MIX
        Traits2<T1,T2>::samebase &&
#else
        Traits2<T1,T2>::sametype &&
#endif
        V1::vsize == UNKNOWN && V2::vsize == UNKNOWN) };
    return CallMultVV<inst,V1,V2>::call(v1.vec(),v2.vec());
  }


  // v * v
  template <class V1, class V2>
  inline PT operator*(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2) 
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    return MultVV(v1.calc().CView(),v2.calc().CView()); 
  }

#define PT2 typename Traits2<Tx,PT>::type
  // v * (x*v)
  template <class V1, int ix2, class Tx, class V2>
  inline PT2 operator*(const BaseVector<V1>& v1, const ProdXV<ix2,Tx,V2>& v2) 
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    return v2.GetX() * (v1 * v2.GetV());
  }

  // (x*v) * v
  template <int ix1, class Tx, class V1, class V2>
  inline PT2 operator*(const ProdXV<ix1,Tx,V1>& v1, const BaseVector<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    return v1.GetX() * (v1.GetV() * v2);
  }
#undef PT2

#define PT2 typename Traits2<Tx1,typename Traits2<Tx2,PT>::type>::type
  // (x*v) * (x*v)
  template <int ix1, class Tx1, class V1, int ix2, class Tx2, class V2>
  inline PT2 operator*(
      const ProdXV<ix1,Tx1,V1>& v1, const ProdXV<ix2,Tx2,V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    return v1.GetX() * v2.GetX() * (v1.GetV() * v2.GetV());
  }
#undef PT2

#undef PT


  //
  // ElementProd functions:
  //

  template <bool inst, bool add, int ix, class T, class V1, class V2, class V3>
  struct CallElemMultVV // inst = false
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
    { InlineElemMultVV<add>(x,v1,v2,v3); }
  };
  template <int ix, class T, class V1, class V2, class V3>
  struct CallElemMultVV<true,false,ix,T,V1,V2,V3> // inst = true, !add
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
    {
      typedef typename V3::value_type T3;
      InstElemMultVV(T3(x),v1.XView(),v2.XView(),v3.XView()); 
    }
  };
  template <int ix, class T, class V1, class V2, class V3>
  struct CallElemMultVV<true,true,ix,T,V1,V2,V3> // inst = true, add
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
    {
      typedef typename V3::value_type T3;
      InstAddElemMultVV(T3(x),v1.XView(),v2.XView(),v3.XView()); 
    }
  };

  template <bool add, int ix, class T, class V1, class V2, class V3>
  inline void ElemMultVV(const Scaling<ix,T>& x, 
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2,
      BaseVector_Mutable<V3>& v3)
  {
    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename V3::value_type T3;
    enum { inst = (
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits<T3>::isinst &&
        Traits2<T,T3>::samebase &&
#ifdef TMV_INST_MIX
        Traits2<T1,T3>::samebase &&
        Traits2<T2,T3>::samebase &&
#else
        Traits2<T1,T3>::sametype &&
        Traits2<T2,T3>::sametype &&
#endif
        V1::vconj == false &&
        V2::vconj == false &&
        V3::vconj == false &&
        V1::vsize == UNKNOWN &&
        V2::vsize == UNKNOWN &&
        V3::vsize == UNKNOWN) };
    return CallElemMultVV<inst,add,ix,T,V1,V2,V3>::call(
        x,v1.vec(),v2.vec(),v3.vec());
  }

  template <int ix, class T, class V1, class V2>
  class ElemProdVV;

  template <int ix, class T, class V1, class V2>
  struct Traits<ElemProdVV<ix,T,V1,V2> >
  {
    typedef typename ProdXV<ix,T,V1>::value_type vtype1;
    typedef typename V2::value_type vtype2;
    typedef typename Traits2<vtype1,vtype2>::type value_type;

    enum { vsize = Sizes<V1::vsize,V2::vsize>::size };
    enum { vfort = V1::vfort && V2::vfort };
    enum { vcalc = false };

    typedef ElemProdVV<ix,T,V1,V2> type;
    typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef typename TypeSelect<V1::vcalc&&V2::vcalc,const type,calc_type>::type eval_type;
  };

  template <int ix, class T, class V1, class V2>
  class ElemProdVV :
    public BaseVector<ElemProdVV<ix,T,V1,V2> >
  {
  public:

    typedef ElemProdVV<ix,T,V1,V2> type;
    typedef typename Traits<type>::value_type value_type;
    enum { vsize = Traits<type>::vsize };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    inline ElemProdVV(const T _x, const BaseVector<V1>& _v1,
        const BaseVector<V2>& _v2) : x(_x), v1(_v1.vec()), v2(_v2.vec())
    {
      TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
      TMVAssert(v1.size() == v2.size());
    }

    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const V1& GetV1() const { return v1; }
    inline const V2& GetV2() const { return v2; }

    inline size_t size() const { return v1.size(); }

    inline value_type cref(int i) const
    { return x * (v1.cref(i) * v2.cref(i)); }

    template <class V3>
    inline void AssignTo(BaseVector_Mutable<V3>& v3) const
    {
      TMVStaticAssert((Sizes<vsize,V3::vsize>::same));
      TMVAssert(size() == v3.size());
      TMVStaticAssert(visreal || V3::viscomplex);
      typename V3::cview_type v3cv = v3.CView();
      ElemMultVV<false>(x,v1.calc().CView(),v2.calc().CView(),v3cv);
    }
  private:
    const Scaling<ix,T> x;
    const V1& v1;
    const V2& v2;
  };

#define RT typename V2::real_type
  template <class V1, class V2> 
  inline ElemProdVV<1,RT,V1,V2> ElementProd(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVAssert(v1.size() == v2.size());
    return ElemProdVV<1,RT,V1,V2>(RT(1),v1,v2);
  }
#undef RT

  // v += [vv]
  template <class V3, int ix, class T, class V1, class V2>
  inline void AddEq(
      BaseVector_Mutable<V3>& v3, const ElemProdVV<ix,T,V1,V2>& vv)
  { 
    typename V3::cview_type v3cv = v3.CView();
    ElemMultVV<true>(vv.GetX(),vv.GetV1().calc().CView(),
        vv.GetV2().calc().CView(),v3cv); 
  }

  // v -= [vv]
  template <class V3, int ix, class T, class V1, class V2>
  inline void SubtractEq(
      BaseVector_Mutable<V3>& v3, const ElemProdVV<ix,T,V1,V2>& vv)
  {
    typename V3::cview_type v3cv = v3.CView();
    ElemMultVV<true>(-vv.GetX(),vv.GetV1().calc().CView(),
        vv.GetV2().calc().CView(),v3cv); 
  }

  // Consolidate x*(xvv) type constructs:

#define RT typename ElemProdVV<ix,T,V1,V2>::real_type
#define CT typename ElemProdVV<ix,T,V1,V2>::complex_type
#define CCT ConjRef<CT>

  // -[vv]
  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<-ix,T,V1,V2> operator-(const ElemProdVV<ix,T,V1,V2>& mv)
  { return ElemProdVV<-ix,T,V1,V2>(-mv.GetX(),mv.GetV1(),mv.GetV2()); }

  // x * [vv]
  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,T,V1,V2> operator*(
      const RT x, const ElemProdVV<ix,T,V1,V2>& mv)
  { return ElemProdVV<0,T,V1,V2>(x*mv.GetX(),mv.GetV1(),mv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,CT,V1,V2> operator*(
      const CT x, const ElemProdVV<ix,T,V1,V2>& mv)
  { return ElemProdVV<0,CT,V1,V2>(x*mv.GetX(),mv.GetV1(),mv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,CT,V1,V2> operator*(
      const CCT x, const ElemProdVV<ix,T,V1,V2>& mv)
  { return ElemProdVV<0,CT,V1,V2>(x*mv.GetX(),mv.GetV1(),mv.GetV2()); }

  template <int ix1, class T1, int ix, class T, class V1, class V2>
  inline ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2> operator*(
      const Scaling<ix1,T1>& x, const ElemProdVV<ix,T,V1,V2>& mv)
  {
    return ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2>(
        T1(x)*mv.GetX(),mv.GetV1(),mv.GetV2());
  }
  // [vv]*x
  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,T,V1,V2> operator*(
      const ElemProdVV<ix,T,V1,V2>& mv, const RT x)
  { return ElemProdVV<0,T,V1,V2>(x*mv.GetX(),mv.GetV1(),mv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,CT,V1,V2> operator*(
      const ElemProdVV<ix,T,V1,V2>& mv, const CT x)
  { return ElemProdVV<0,CT,V1,V2>(x*mv.GetX(),mv.GetV1(),mv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,CT,V1,V2> operator*(
      const ElemProdVV<ix,T,V1,V2>& mv, const CCT x)
  { return ElemProdVV<0,CT,V1,V2>(x*mv.GetX(),mv.GetV1(),mv.GetV2()); }

  template <int ix1, class T1, int ix, class T, class V1, class V2>
  inline ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2> operator*(
      const ElemProdVV<ix,T,V1,V2>& mv, const Scaling<ix1,T1>& x)
  {
    return ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2>(
        T1(x)*mv.GetX(),mv.GetV1(),mv.GetV2());
  }

  // [vv]/x
  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,T,V1,V2> operator/(
      const ElemProdVV<ix,T,V1,V2>& mv, const RT x)
  { return ElemProdVV<0,T,V1,V2>(mv.GetX()/x,mv.GetV1(),mv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,CT,V1,V2> operator/(
      const ElemProdVV<ix,T,V1,V2>& mv, const CT x)
  { return ElemProdVV<0,CT,V1,V2>(mv.GetX()/x,mv.GetV1(),mv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline ElemProdVV<0,CT,V1,V2> operator/(
      const ElemProdVV<ix,T,V1,V2>& mv, const CCT x)
  { return ElemProdVV<0,CT,V1,V2>(mv.GetX()/x,mv.GetV1(),mv.GetV2()); }

  template <int ix1, class T1, int ix, class T, class V1, class V2>
  inline ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2> operator/(
      const ElemProdVV<ix,T,V1,V2>& mv, const Scaling<ix1,T1>& x)
  {
    return ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2>(
        mv.GetX()/T1(x),mv.GetV1(),mv.GetV2());
  }

#undef RT
#undef CT
#undef CCT

  template <int ix, class T, class V1, class V2>
  inline std::string TypeText(const ElemProdVV<ix,T,V1,V2>& svv)
  {
    std::ostringstream s;
    s << "ElemProdVV< "<<ix<<","<<TypeText(T())<<" , ";
    s << TypeText(svv.GetV1())<<" , ";
    s << TypeText(svv.GetV2())<<" >";
    return s.str();
  }


} // namespace tmv

#endif 
