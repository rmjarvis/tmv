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


#ifndef TMV_SumVV_H
#define TMV_SumVV_H

#include "TMV_ProdXV.h"

namespace tmv {

  //
  // Vector += Vector
  //

  template <bool conj, bool inst, int ix, class T, class V1, class V2>
  struct CallAddVV;
  template <bool inst, int ix, class T, class V1, class V2>
  struct CallAddVV<true,inst,ix,T,V1,V2> // conj = true
  {
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    {
      typedef typename V1::const_conjugate_type V1c;
      typedef typename V2::conjugate_type V2c;
      V1c v1c = v1.Conjugate();
      V2c v2c = v2.Conjugate();
      CallAddVV<false,inst,ix,T,V1c,V2c>::call(TMV_CONJ(x),v1c,v2c);
    }
  };
  template <int ix, class T, class V1, class V2>
  struct CallAddVV<false,true,ix,T,V1,V2> // inst = true
  {
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    { 
      typedef typename V2::value_type T2;
      InstAddVV(T2(x),v1.XView(),v2.XView()); 
    }
  };
  template <int ix, class T, class V1, class V2>
  struct CallAddVV<false,false,ix,T,V1,V2> // inst = false
  {
    static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
    { InlineAddVV(x,v1,v2); }
  };

  template <int ix, class T, class V1, class V2>
  inline void AddVV(
      const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
      BaseVector_Mutable<V2>& v2)
  {
    typedef typename V1::value_type TV1;
    typedef typename V2::value_type TV2;
    enum { inst = (
        Traits<TV1>::isinst &&
        Traits<TV2>::isinst &&
        Traits2<T,TV2>::samebase &&
#ifdef TMV_INST_MIX
        Traits2<TV1,TV2>::samebase &&
#else
        Traits2<TV1,TV2>::sametype &&
#endif
        V1::vsize == UNKNOWN &&
        V2::vsize == UNKNOWN) };
    return CallAddVV<V2::vconj,inst,ix,T,V1,V2>::call(x,v1.vec(),v2.vec());
  }


  //
  // Vector + Vector
  //

  // Defined in TMV_AddVV.h
  template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  inline void InlineAddVV(
      const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
      const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
      BaseVector_Mutable<V3>& v3);

  template <bool inst, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct CallAddVV2 // inst = false
  {
    static inline void call(const Scaling<ix1,T1>& x1, const V1& v1, 
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    { InlineAddVV(x1,v1,x2,v2,v3); }
  };
  template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  struct CallAddVV2<true,ix1,T1,V1,ix2,T2,V2,V3> // inst = true
  {
    static inline void call(const Scaling<ix1,T1>& x1, const V1& v1, 
        const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
    {
      // There isn't much point in actually instantiating this, but 
      // if we have instantiated versions of AddVV and MultXV, then
      // we need to check for aliases. 
      // (The inline version doesn't need this check.)
      if (SameStorage(v1.vec(),v3.vec())) {
        if (SameStorage(v2.vec(),v3.vec())) {
          // Rather than using a temporary, just do the inline version.
          InlineAddVV(x1,v1,x2,v2,v3);
        } else {
          (v3 = x1 * v1) += x2 * v2;
        }
      } else {
        (v3 = x2 * v2) += x1 * v1;
      }
    }
  };

  // v3 = x1 * v1 + x2 * v2
  template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
  inline void AddVV(
      const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
      const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
      BaseVector_Mutable<V3>& v3)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
    TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
    TMVStaticAssert((Sizes<V2::vsize,V3::vsize>::same));
    TMVAssert(v1.size() == v3.size());
    TMVAssert(v2.size() == v3.size());

    typedef typename V1::value_type TV1;
    typedef typename V2::value_type TV2;
    typedef typename V3::value_type TV3;
    enum { inst = (
        Traits<TV1>::isinst &&
        Traits<TV2>::isinst &&
        Traits<TV3>::isinst &&
        Traits2<T1,TV3>::samebase &&
        Traits2<T2,TV3>::samebase &&
#ifdef TMV_INST_MIX
        Traits2<TV1,TV3>::samebase &&
        Traits2<TV2,TV3>::samebase &&
#else
        Traits2<TV1,TV3>::sametype &&
        Traits2<TV2,TV3>::sametype &&
#endif
        V1::vsize == UNKNOWN &&
        V2::vsize == UNKNOWN &&
        V3::vsize == UNKNOWN) };
    return CallAddVV2<inst,ix1,T1,V1,ix2,T2,V2,V3>::call(
        x1,v1.vec(),x2,v2.vec(),v3.vec());
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  class SumVV;

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  struct Traits<SumVV<ix1,T1,V1,ix2,T2,V2> >
  {
    typedef typename ProdXV<ix1,T1,V1>::value_type vtype1;
    typedef typename ProdXV<ix2,T2,V2>::value_type vtype2;
    typedef typename Traits2<vtype1,vtype2>::type value_type;

    enum { vsize = Sizes<V1::vsize,V2::vsize>::size };
    enum { vfort = V1::vfort && V2::vfort };
    enum { vcalc = false };

    typedef SumVV<ix1,T1,V1,ix2,T2,V2> type;
    typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef typename TypeSelect<V1::vcalc&&V2::vcalc,const type,calc_type>::type eval_type;
  };

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  class SumVV :
    public BaseVector<SumVV<ix1,T1,V1,ix2,T2,V2> >
  {
  public:

    typedef SumVV<ix1,T1,V1,ix2,T2,V2> type;
    typedef typename Traits<type>::value_type value_type;
    enum { vsize = Traits<type>::vsize };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    inline SumVV(const T1& _x1, const BaseVector<V1>& _v1,
        const T2& _x2, const BaseVector<V2>& _v2) :
      x1(_x1), v1(_v1.vec()), x2(_x2), v2(_v2.vec())
    { 
      TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
      TMVAssert(v1.size() == v2.size());
    }

    inline const Scaling<ix1,T1>& GetX1() const { return x1; }
    inline const V1& GetV1() const { return v1; }
    inline const Scaling<ix2,T2>& GetX2() const { return x2; }
    inline const V2& GetV2() const { return v2; }

    inline size_t size() const { return v1.size(); }
    inline value_type cref(int i) const
    { return x1 * v1.cref(i) + x2 * v2.cref(i); }

    template <class V3>
    inline void AssignTo(BaseVector_Mutable<V3>& v3) const
    {
      TMVStaticAssert((Sizes<vsize,V3::vsize>::same)); 
      TMVAssert(size() == v3.size());
      TMVStaticAssert(visreal || V3::viscomplex);
      typename V3::cview_type v3cv = v3.CView();
      AddVV(x1,v1.calc().CView(),x2,v2.calc().CView(),v3cv);
    }
  private:
    const Scaling<ix1,T1> x1;
    const V1& v1;
    const Scaling<ix2,T2> x2;
    const V2& v2;
  };

#define RT typename V2::real_type
  // v += v
  template <class V1, class V2>
  inline void AddEq(
      BaseVector_Mutable<V1>& v1, const BaseVector<V2>& v2) 
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    TMVStaticAssert(V1::viscomplex || V2::visreal);
    typename V1::cview_type v1cv = v1.CView();
    AddVV(Scaling<1,RT>(),v2.calc().CView(),v1cv);
  }

  // v += xv
  template <class V1, int ix2, class T2, class V2>
  inline void AddEq(
      BaseVector_Mutable<V1>& v1, const ProdXV<ix2,T2,V2>& v2) 
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    TMVStaticAssert(V1::viscomplex || V2::visreal);
    typename V1::cview_type v1cv = v1.CView();
    AddVV(v2.GetX(),v2.GetV().calc().CView(),v1cv);
  }

  // v -= v
  template <class V1, class V2>
  inline void SubtractEq(
      BaseVector_Mutable<V1>& v1, const BaseVector<V2>& v2)
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    TMVStaticAssert(V1::viscomplex || V2::visreal);
    typename V1::cview_type v1cv = v1.CView();
    AddVV(Scaling<-1,RT>(),v2.calc().CView(),v1cv);
  }

  // v -= xv
  template <class V1, int ix2, class T2, class V2>
  inline void SubtractEq(
      BaseVector_Mutable<V1>& v1, const ProdXV<ix2,T2,V2>& v2) 
  {
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    TMVStaticAssert(V1::viscomplex || V2::visreal);
    typename V1::cview_type v1cv = v1.CView();
    AddVV(-v2.GetX(),v2.GetV().calc().CView(),v1cv);
  }

  // v + v
  template <class V1, class V2>
  inline SumVV<1,RT,V1,1,RT,V2> operator+(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<1,RT,V1,1,RT,V2>(RT(1),v1,RT(1),v2); 
  }

  // xv + v
  template <int ix1, class T1, class V1, class V2>
  inline SumVV<ix1,T1,V1,1,RT,V2> operator+(
      const ProdXV<ix1,T1,V1>& v1, const BaseVector<V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<ix1,T1,V1,1,RT,V2>(T1(v1.GetX()),v1.GetV(),RT(1),v2); 
  }

  // v + xv
  template <class V1, int ix2, class T2, class V2>
  inline SumVV<1,RT,V1,ix2,T2,V2> operator+(
      const BaseVector<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<1,RT,V1,ix2,T2,V2>(RT(1),v1,T2(v2.GetX()),v2.GetV()); 
  }

  // xv + xv
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<ix1,T1,V1,ix2,T2,V2> operator+(
      const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<ix1,T1,V1,ix2,T2,V2>(
        T1(v1.GetX()),v1.GetV(),T2(v2.GetX()),v2.GetV()); 
  }

  // v - v
  template <class V1, class V2>
  inline SumVV<1,RT,V1,-1,RT,V2> operator-(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<1,RT,V1,-1,RT,V2>(RT(1),v1,RT(-1),v2); 
  }

  // xv - v
  template <int ix1, class T1, class V1, class V2>
  inline SumVV<ix1,T1,V1,-1,RT,V2> operator-(
      const ProdXV<ix1,T1,V1>& v1, const BaseVector<V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<ix1,T1,V1,-1,RT,V2>(T1(v1.GetX()),v1.GetV(),RT(-1),v2); 
  }

  // v - xv
  template <class V1, int ix2, class T2, class V2>
  inline SumVV<1,RT,V1,-ix2,T2,V2> operator-(
      const BaseVector<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<1,RT,V1,-ix2,T2,V2>(RT(1),v1,T2(-v2.GetX()),v2.GetV()); 
  }

  // xv - xv
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<ix1,T1,V1,-ix2,T2,V2> operator-(
      const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { 
    TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
    TMVAssert(v1.size() == v2.size());
    return SumVV<ix1,T1,V1,-ix2,T2,V2>(
        T1(v1.GetX()),v1.GetV(),T2(-v2.GetX()),v2.GetV()); 
  }
#undef RT


  // Consolidate x*(xv+xv) type constructs:

#define RT typename SumVV<ix1,T1,V1,ix2,T2,V2>::real_type
#define CT typename SumVV<ix1,T1,V1,ix2,T2,V2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

  // -(xv+xv)
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<-ix1,T1,V1,-ix2,T2,V2> operator-(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  {
    return SumVV<-ix1,T1,V1,-ix2,T2,V2>(
        -T1(svv.GetX1()),svv.GetV1(),-T2(svv.GetX2()),svv.GetV2()); 
  }

  // x * (xv+xv)
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,T1,V1,0,T2,V2> operator*(
      const int x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  {
    return SumVV<0,T1,V1,0,T2,V2>(
        RT(x)*svv.GetX1(),svv.GetV1(), RT(x)*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,T1,V1,0,T2,V2> operator*(
      const RT x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  {
    return SumVV<0,T1,V1,0,T2,V2>(
        x*svv.GetX1(),svv.GetV1(), x*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,CT,V1,0,CT,V2> operator*(
      const CT x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  {
    return SumVV<0,CT,V1,0,CT,V2>(
        x*svv.GetX1(),svv.GetV1(), x*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,CT,V1,0,CT,V2> operator*(
      const CCT x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  {
    return SumVV<0,CT,V1,0,CT,V2>(
        CT(x)*svv.GetX1(),svv.GetV1(), CT(x)*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix, class T, int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2> operator*(
      const Scaling<ix,T>& x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  {
    return SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2>(
        T(x)*svv.GetX1(),svv.GetV1(),T(x)*svv.GetX2(),svv.GetV2());
  }

  // (xv+xv)*x
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,T1,V1,0,T2,V2> operator*(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const int x)
  {
    return SumVV<0,T1,V1,0,T2,V2>(
        RT(x)*svv.GetX1(),svv.GetV1(), RT(x)*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,T1,V1,0,T2,V2> operator*(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const RT x)
  {
    return SumVV<0,T1,V1,0,T2,V2>(
        x*svv.GetX1(),svv.GetV1(), x*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,CT,V1,0,CT,V2> operator*(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CT x)
  { 
    return SumVV<0,CT,V1,0,CT,V2>(
        x*svv.GetX1(),svv.GetV1(), x*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,CT,V1,0,CT,V2> operator*(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CCT x)
  {
    return SumVV<0,CT,V1,0,CT,V2>(
        CT(x)*svv.GetX1(),svv.GetV1(), CT(x)*svv.GetX2(),svv.GetV2()); 
  }

  template <int ix, class T, int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2> operator*(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const Scaling<ix,T>& x)
  {
    return SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2>(
        T(x)*svv.GetX1(),svv.GetV1(),T(x)*svv.GetX2(),svv.GetV2());
  }

  // (xv+xv)/x
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,T1,V1,0,T2,V2> operator/(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const int x)
  {
    return SumVV<0,T1,V1,0,T2,V2>(
        svv.GetX1()/RT(x),svv.GetV1(), svv.GetX2()/RT(x),svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,T1,V1,0,T2,V2> operator/(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const RT x)
  {
    return SumVV<0,T1,V1,0,T2,V2>(
        svv.GetX1()/x,svv.GetV1(), svv.GetX2()/x,svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,CT,V1,0,CT,V2> operator/(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CT x)
  {
    return SumVV<0,CT,V1,0,CT,V2>(
        svv.GetX1()/x,svv.GetV1(), svv.GetX2()/x,svv.GetV2()); 
  }

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<0,CT,V1,0,CT,V2> operator/(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CCT x)
  {
    return SumVV<0,CT,V1,0,CT,V2>(
        svv.GetX1()/CT(x),svv.GetV1(), svv.GetX2()/CT(x),svv.GetV2()); 
  }

  template <int ix, class T, int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2> operator/(
      const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const Scaling<ix,T>& x)
  {
    return SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2>(
        svv.GetX1()/T(x),svv.GetV1(),svv.GetX2()/T(x),svv.GetV2());
  }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



  // TypeText

  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline std::string TypeText(const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
  { 
    std::ostringstream s;
    s << "SumVV< "<< ix1<<","<<TypeText(T1())<<" , "<<TypeText(svv.GetV1());
    s << " , "<<ix2<<","<<TypeText(T2())<<" , "<<TypeText(svv.GetV2())<<" >";
    return s.str();
  }

} // namespace tmv

#endif 
