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


#ifndef TMV_ProdMV_H
#define TMV_ProdMV_H

#include "TMV_ProdXM.h"
#include "TMV_ProdXV.h"
#include "TMV_Vector.h"
#include "TMV_SmallVector.h"
#include "TMV_MultMV.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
#endif

namespace tmv {

  //
  // Matrix * Vector
  //


  template <bool add, int ix, class T, class M1, class V2, class V3>
  inline void MultMV(const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
      const BaseVector_Calc<V2>& m2, BaseVector_Mutable<V3>& m3);

  template <bool checkalias, bool conj, bool inst, bool add, int ix, class T, class M1, class V2, class V3> struct CallMultMV;
  template <bool conj, bool inst, bool add, int ix, class T, class M1, class V2, class V3> 
  struct CallMultMV<true,conj,inst,add,ix,T,M1,V2,V3> // checkalias = true
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef RealType(T) RT;
      const Scaling<1,RT> one;

      if ( !SameStorage(m1.mat(),v3.vec()) &&
          // There are some cases where the vectors are allowed to have
          // the same storage, so we check for these first.
          ( (!ShapeTraits<M1::mshape>::lower && v2.step() >= v3.step()) ||
            (!ShapeTraits<M1::mshape>::upper && v2.step() <= v3.step()) ||
            !SameStorage(v2.vec(),v3.vec())) ) 
      {
        CallMultMV<false,conj,inst,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
      }
      else if (SameStorage(m1.mat(),v3.vec())) 
      {
        Maybe<add>::add( v3 , x * (m1*v2).calc() );
      } 
      else 
      {
        TMVAssert(SameStorage(v2.vec(),v3.vec()));
        Maybe<add>::add(v3 , m1 * (x*v2).calc() );
      } 
    }
  };
  template <bool inst, bool add, int ix, class T, class M1, class V2, class V3>
  struct CallMultMV<false,true,inst,add,ix,T,M1,V2,V3> // conj = true
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { 
      typedef typename M1::const_conjugate_type M1c;
      typedef typename V2::const_conjugate_type V2c;
      typedef typename V3::conjugate_type V3c;
      M1c m1c = m1.Conjugate();
      V2c v2c = v2.Conjugate();
      V3c v3c = v3.Conjugate();
      CallMultMV<false,false,inst,add,ix,T,M1c,V2c,V3c>::call(
          TMV_CONJ(x),m1c,v2c,v3c);
    }
  };
  template <bool add, int ix, class T, class M1, class V2, class V3>
  struct CallMultMV<false,false,false,add,ix,T,M1,V2,V3> // inst = false
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { InlineMultMV<add>(x,m1,v2,v3); }
  };
  template <int ix, class T, class M1, class V2, class V3>
  struct CallMultMV<false,false,true,true,ix,T,M1,V2,V3> // inst = true, add
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    { 
      typedef typename V3::value_type T3;
      InstAddMultMV(T3(x),m1.XView(),v2.XView(),v3.XView()); 
    }
  };
  template <int ix, class T, class M1, class V2, class V3>
  struct CallMultMV<false,false,true,false,ix,T,M1,V2,V3> // inst = true, !add
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
    {
      typedef typename V3::value_type T3;
      InstMultMV(T3(x),m1.XView(),v2.XView(),v3.XView()); 
    }
  };

  // v3 (+=) x * m1 * v2
  template <bool add, int ix, class T, class M1, class V2, class V3>
  inline void MultMV(const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
      const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,V3::vsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,V2::vsize>::same));
    TMVAssert(m1.colsize() == v3.size());
    TMVAssert(m1.rowsize() == v2.size());

    typedef typename M1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename V3::value_type T3;

    enum { checkalias = (
        M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
        V2::vsize == UNKNOWN && V3::vsize == UNKNOWN ) };
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
        checkalias ) };
#ifdef XDEBUG
    //std::cout<<"Start MultMV XDEBUG"<<std::endl;
    //std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    //std::cout<<"v3 = "<<TypeText(v3)<<"  "<<v3<<std::endl;
    Vector<typename V3::value_type> v3i = v3;
    //std::cout<<"v3i = "<<TypeText(v3i)<<"  "<<v3i<<std::endl;
    Vector<typename V3::value_type> v3c = v3;
    //std::cout<<"v3c = "<<TypeText(v3c)<<"  "<<v3c<<std::endl;
    if (!add) v3c.Zero();
    for(int i=0;i<v3.size();++i) {
      for(int j=0;j<v2.size();++j) {
        v3c.ref(i) += T(x) * m1.cref(i,j) * v2.cref(j);
      }
    }
    //std::cout<<"v3c => "<<v3c<<std::endl;
    //std::cout<<"checkalias, inst add = "<<checkalias<<" "<<inst<<" "<<add<<std::endl;
#endif
    CallMultMV<checkalias,V3::vconj,inst,add,ix,T,M1,V2,V3>::call(
        x,m1.mat(),v2.vec(),v3.vec());
#ifdef XDEBUG
    //std::cout<<"v3 => "<<v3<<std::endl;
    //std::cout<<"Norm(v3 - v3c) => "<<Norm(v3-v3c)<<std::endl;
    if (Norm(v3-v3c) > 1.e-6 * std::abs(T(x)) * Norm(m1) * Norm(v2)) {
      std::cout<<"MultMV:  add = "<<add<<std::endl;
      std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
      std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
      std::cout<<"v3 = "<<TypeText(v3)<<"  "<<v3i<<std::endl;
      std::cout<<"v3 -> "<<v3<<std::endl;
      std::cout<<"correct = "<<v3c<<std::endl;
      std::cout<<"diff = "<<v3-v3c<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(v3-v3c)<<std::endl;
      exit(1);
    }
#endif
  }

  // v3 (+=) x * v1 * m2
  template <bool add, int ix, class T, class V1, class M2, class V3>
  inline void MultVM(const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
      const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
  {
    TMVStaticAssert((Sizes<V1::vsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M2::mrowsize,V3::vsize>::same));
    TMVAssert(v1.size() == m2.colsize());
    TMVAssert(m2.rowsize() == v3.size());
    MultMV<add>(x,m2.Transpose(),v1,v3);
  }

  // Here we just need to check whether the shape of M2 requires 
  // that we copy v1 to new storage for a v *= m operation.
  template <bool copy, class V1, int ix, class T, class M2> struct CallMultEqVM;

  template <class V1, int ix, class T, class M2> 
  struct CallMultEqVM<true,V1,ix,T,M2> // copy = true
  {
    static inline void call(
        V1& v1, const Scaling<ix,T>& x, const M2& m2)
    { MultVM<false>(x,v1.copy().CView(),m2,v1); }
  };

  template <class V1, int ix, class T, class M2> 
  struct CallMultEqVM<false,V1,ix,T,M2> // copy = false
  {
    static inline void call(
        V1& v1, const Scaling<ix,T>& x, const M2& m2)
    { MultVM<false>(x,v1,m2,v1); }
  };

  template <class V1, int ix, class T, class M2>
  inline void MultEqVM(BaseVector_Mutable<V1>& v1,
      const Scaling<ix,T>& x, const BaseMatrix_Calc<M2>& m2)
  {
    enum { copy = (
        !ShapeTraits<M2::mshape>::upper && 
        !ShapeTraits<M2::mshape>::lower ) };
    CallMultEqVM<copy,V1,ix,T,M2>::call(v1.vec(),x,m2.mat());
  }

  template <int ix, class T, class M1, class V2>
  class ProdMV;

  template <int ix, class T, class M1, class V2>
  struct Traits<ProdMV<ix,T,M1,V2> >
  {
    typedef typename ProdXM<ix,T,M1>::value_type mtype1;
    typedef typename V2::value_type vtype2;
    typedef typename Traits2<mtype1,vtype2>::type value_type;

    enum { vsize = M1::mcolsize };
    enum { vfort = M1::mfort && V2::vfort };
    enum { vcalc = false };

    typedef ProdMV<ix,T,M1,V2> type;
    typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef const copy_type eval_type;
  };

  template <int ix, class T, class M1, class V2>
  class ProdMV :
    public BaseVector<ProdMV<ix,T,M1,V2> >
  {
  public:

    typedef ProdMV<ix,T,M1,V2> type;
    typedef typename Traits<type>::value_type value_type;

    enum { vsize = Traits<type>::vsize };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    inline ProdMV(const T _x, const BaseMatrix<M1>& _m1, 
        const BaseVector<V2>& _v2) : x(_x), m1(_m1.mat()), v2(_v2.vec())
    {
      TMVStaticAssert((Sizes<M1::mrowsize,V2::vsize>::same)); 
      TMVAssert(m1.rowsize() == v2.size());
    }

    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const M1& GetM() const { return m1; }
    inline const V2& GetV() const { return v2; }

    inline size_t size() const { return m1.colsize(); }

    inline value_type cref(int i) const
    { return x * (m1.get_row(i) * v2); }

    template <class V3>
    inline void AssignTo(BaseVector_Mutable<V3>& v3) const
    {
      TMVStaticAssert((visreal || V3::viscomplex));
      TMVStaticAssert((Sizes<vsize,V3::vsize>::same)); 
      TMVAssert(size() == v3.size());
      typename V3::cview_type v3cv = v3.CView();
      MultMV<false>(x,m1.calc().CView(),v2.calc().CView(),v3cv);
    }
  private:
    const Scaling<ix,T> x;
    const M1& m1;
    const V2& v2;
  };


  template <int ix, class T, class V1, class M2>
  class ProdVM;

  template <int ix, class T, class V1, class M2>
  struct Traits<ProdVM<ix,T,V1,M2> >
  {
    typedef typename V1::value_type vtype1;
    typedef typename ProdXM<ix,T,M2>::value_type mtype2;
    typedef typename Traits2<vtype1,mtype2>::type value_type;

    enum { vsize = M2::mrowsize };
    enum { vfort = V1::vfort && M2::mfort };
    enum { vcalc = false };

    typedef ProdVM<ix,T,V1,M2> type;
    typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef const copy_type eval_type;
  };

  template <int ix, class T, class V1, class M2>
  class ProdVM :
    public BaseVector<ProdVM<ix,T,V1,M2> >
  {
  public:

    typedef ProdVM<ix,T,V1,M2> type;
    typedef typename Traits<type>::value_type value_type;

    enum { vsize = Traits<type>::vsize };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { visreal = Traits<value_type>::isreal };
    enum { viscomplex = Traits<value_type>::iscomplex };

    inline ProdVM(const T _x, const BaseVector<V1>& _v1, 
        const BaseMatrix<M2>& _m2) : x(_x), v1(_v1.vec()), m2(_m2.mat())
    {
      TMVStaticAssert((Sizes<V1::vsize,M2::mcolsize>::same)); 
      TMVAssert(v1.size() == m2.colsize());
    }

    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const V1& GetV() const { return v1; }
    inline const M2& GetM() const { return m2; }

    inline size_t size() const { return m2.rowsize(); }

    inline value_type cref(int j) const
    { return x * (v1 * m2.get_col(j)); }

    template <class V3>
    inline void AssignTo(BaseVector_Mutable<V3>& v3) const
    {
      TMVStaticAssert((visreal || V3::viscomplex));
      TMVStaticAssert((Sizes<vsize,V3::vsize>::same)); 
      TMVAssert(size() == v3.size());
      typename V3::cview_type v3cv = v3.CView();
      MultVM<false>(x,v1.calc().CView(),m2.calc().CView(),v3cv);
    }
  private:
    const Scaling<ix,T> x;
    const V1& v1;
    const M2& m2;
  };


  // m * v
#define RT typename V::real_type
  template <class M, class V>
  inline ProdMV<1,RT,M,V> operator*(
      const BaseMatrix<M>& m, const BaseVector<V>& v)
  { return ProdMV<1,RT,M,V>(RT(1),m,v); }
#undef RT

  // m * xv
  template <class M, int ix, class T, class V>
  inline ProdMV<ix,T,M,V> operator*(
      const BaseMatrix<M>& m, const ProdXV<ix,T,V>& v)
  { return ProdMV<ix,T,M,V>(v.GetX(),m,v.GetV()); }

  // xm * v
  template <int ix, class T, class M, class V>
  inline ProdMV<ix,T,M,V> operator*(
      const ProdXM<ix,T,M>& m, const BaseVector<V>& v)
  { return ProdMV<ix,T,M,V>(m.GetX(),m.GetM(),v); }

  // xm * xv
#define PT typename Traits2<T1,T2>::type
  template <int ix1, class T1, class M, int ix2, class T2, class V>
  inline ProdMV<ix1*ix2,PT,M,V> operator*(
      const ProdXM<ix1,T1,M>& m, const ProdXV<ix2,T2,V>& v)
  { return ProdMV<ix1*ix2,PT,M,V>(m.GetX()*v.GetX(),m.GetM(),v.GetV()); }
#undef PT


  // v * m
#define RT typename V::real_type
  template <class M, class V>
  inline ProdVM<1,RT,V,M> operator*(
      const BaseVector<V>& v, const BaseMatrix<M>& m)
  { return ProdVM<1,RT,V,M>(RT(1),v,m); }
#undef RT

  // v * xm
  template <class M, int ix, class T, class V>
  inline ProdVM<ix,T,V,M> operator*(
      const BaseVector<V>& v, const ProdXM<ix,T,M>& m)
  { return ProdVM<ix,T,V,M>(m.GetX(),v,m.GetM()); }

  // xv * m
  template <int ix, class T, class M, class V>
  inline ProdVM<ix,T,V,M> operator*(
      const ProdXV<ix,T,V>& v, const BaseMatrix<M>& m)
  { return ProdVM<ix,T,V,M>(v.GetX(),v.GetV(),m); }

  // xv * xm
#define PT typename Traits2<T1,T2>::type
  template <int ix1, class T1, class M, int ix2, class T2, class V>
  inline ProdVM<ix1*ix2,PT,V,M> operator*(
      const ProdXV<ix1,T1,V>& v, const ProdXM<ix2,T2,M>& m)
  { return ProdVM<ix1*ix2,PT,V,M>(m.GetX()*v.GetX(),v.GetV(),m.GetM()); }
#undef PT



  // v *= m
#define RT typename V1::real_type
  template <class V1, class M2>
  inline void MultEq(
      BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
  {
    typename V1::cview_type v1cv = v1.CView();
    MultEqVM(v1cv,Scaling<1,RT>(),m2.calc().CView()); }
#undef RT

  // v *= xm
  template <class V1, int ix2, class T2, class M2>
  inline void MultEq(
      BaseVector_Mutable<V1>& v1, const ProdXM<ix2,T2,M2>& m2)
  {
    typename V1::cview_type v1cv = v1.CView();
    MultEqVM(v1cv,m2.GetX(),m2.GetM().calc().CView()); 
  }

  // v += mv
  template <class V3, int ix, class T, class M1, class V2>
  inline void AddEq(
      BaseVector_Mutable<V3>& v3, const ProdMV<ix,T,M1,V2>& mv)
  {
    typename V3::cview_type v3cv = v3.CView();
    MultMV<true>(mv.GetX(),mv.GetM().calc().CView(),
        mv.GetV().calc().CView(),v3.vec()); 
  }
  
  // v -= mv
  template <class V3, int ix, class T, class M1, class V2>
  inline void SubtractEq(
      BaseVector_Mutable<V3>& v3, const ProdMV<ix,T,M1,V2>& mv)
  {
    typename V3::cview_type v3cv = v3.CView();
    MultMV<true>(-mv.GetX(),mv.GetM().calc().CView(),
        mv.GetV().calc().CView(),v3.vec()); 
  }

  // v += vm
  template <class V3, int ix, class T, class M1, class V2>
  inline void AddEq(
      BaseVector_Mutable<V3>& v3, const ProdVM<ix,T,M1,V2>& vm)
  {
    typename V3::cview_type v3cv = v3.CView();
    MultVM<true>(vm.GetX(),vm.GetV().calc().CView(),
        vm.GetM().calc().CView(),v3.vec()); 
  }
  
  // v -= vm
  template <class V3, int ix, class T, class M1, class V2>
  inline void SubtractEq(
      BaseVector_Mutable<V3>& v3, const ProdVM<ix,T,M1,V2>& vm)
  {
    typename V3::cview_type v3cv = v3.CView();
    MultVM<true>(-vm.GetX(),vm.GetV().calc().CView(),
        vm.GetM().calc().CView(),v3.vec()); 
  }


  // Consolidate x*(xmv) type constructs:

#define RT typename ProdMV<ix,T,M1,V2>::real_type
#define CT typename ProdMV<ix,T,M1,V2>::complex_type
#define CCT ConjRef<CT>

  // -(mv)
  template <int ix, class T, class M1, class V2>
  inline ProdMV<-ix,T,M1,V2> operator-(const ProdMV<ix,T,M1,V2>& mv)
  { return ProdMV<-ix,T,M1,V2>(-mv.GetX(),mv.GetM(),mv.GetV()); }

  // x * (mv)
  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,T,M1,V2> operator*(const RT x, const ProdMV<ix,T,M1,V2>& mv)
  { return ProdMV<0,T,M1,V2>(x*mv.GetX(),mv.GetM(),mv.GetV()); }

  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,CT,M1,V2> operator*(const CT x, const ProdMV<ix,T,M1,V2>& mv)
  { return ProdMV<0,CT,M1,V2>(x*mv.GetX(),mv.GetM(),mv.GetV()); }

  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,CT,M1,V2> operator*(const CCT x, const ProdMV<ix,T,M1,V2>& mv)
  { return ProdMV<0,CT,M1,V2>(x*mv.GetX(),mv.GetM(),mv.GetV()); }

  template <int ix1, class T1, int ix, class T, class M1, class V2>
  inline ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
      const Scaling<ix1,T1>& x, const ProdMV<ix,T,M1,V2>& mv)
  { 
    return ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
      T1(x)*mv.GetX(),mv.GetM(),mv.GetV()); 
  }

  // (mv)*x
  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,T,M1,V2> operator*(const ProdMV<ix,T,M1,V2>& mv, const RT x)
  { return ProdMV<0,T,M1,V2>(x*mv.GetX(),mv.GetM(),mv.GetV()); }

  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,CT,M1,V2> operator*(const ProdMV<ix,T,M1,V2>& mv, const CT x)
  { return ProdMV<0,CT,M1,V2>(x*mv.GetX(),mv.GetM(),mv.GetV()); }

  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,CT,M1,V2> operator*(const ProdMV<ix,T,M1,V2>& mv, const CCT x)
  { return ProdMV<0,CT,M1,V2>(x*mv.GetX(),mv.GetM(),mv.GetV()); }

  template <int ix1, class T1, int ix, class T, class M1, class V2>
  inline ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
      const ProdMV<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
  { 
    return ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
      T1(x)*mv.GetX(),mv.GetM(),mv.GetV()); 
  }

  // (mv)/x
  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,T,M1,V2> operator/(const ProdMV<ix,T,M1,V2>& mv, const RT x)
  { return ProdMV<0,T,M1,V2>(mv.GetX()/x,mv.GetM(),mv.GetV()); }

  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,CT,M1,V2> operator/(const ProdMV<ix,T,M1,V2>& mv, const CT x)
  { return ProdMV<0,CT,M1,V2>(mv.GetX()/x,mv.GetM(),mv.GetV()); }

  template <int ix, class T, class M1, class V2>
  inline ProdMV<0,CT,M1,V2> operator/(const ProdMV<ix,T,M1,V2>& mv, const CCT x)
  { return ProdMV<0,CT,M1,V2>(mv.GetX()/x,mv.GetM(),mv.GetV()); }

  template <int ix1, class T1, int ix, class T, class M1, class V2>
  inline ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator/(
      const ProdMV<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
  { 
    return ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
      mv.GetX()/T1(x),mv.GetM(),mv.GetV()); 
  }

#undef RT
#undef CT
#undef CCT

#define RT typename ProdVM<ix,T,V1,M2>::real_type
#define CT typename ProdVM<ix,T,V1,M2>::complex_type
#define CCT ConjRef<CT>

  // -(vm)
  template <int ix, class T, class V1, class M2>
  inline ProdVM<-ix,T,V1,M2> operator-(const ProdVM<ix,T,V1,M2>& vm)
  { return ProdVM<-ix,T,V1,M2>(-vm.GetX(),vm.GetV(),vm.GetM()); }

  // x * (vm)
  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,T,V1,M2> operator*(const int x, const ProdVM<ix,T,V1,M2>& vm)
  { return ProdVM<0,T,V1,M2>(RT(x)*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,T,V1,M2> operator*(const RT x, const ProdVM<ix,T,V1,M2>& vm)
  { return ProdVM<0,T,V1,M2>(x*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,CT,V1,M2> operator*(const CT x, const ProdVM<ix,T,V1,M2>& vm)
  { return ProdVM<0,CT,V1,M2>(x*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,CT,V1,M2> operator*(const CCT x, const ProdVM<ix,T,V1,M2>& vm)
  { return ProdVM<0,CT,V1,M2>(x*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix1, class T1, int ix, class T, class V1, class M2>
  inline ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2> operator*(
      const Scaling<ix1,T1>& x, const ProdVM<ix,T,V1,M2>& vm)
  { 
    return ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2>(
      T1(x)*vm.GetX(),vm.GetV(),vm.GetM()); 
  }

  // (vm)*x
  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,T,V1,M2> operator*(const ProdVM<ix,T,V1,M2>& vm, const int x)
  { return ProdVM<0,T,V1,M2>(RT(x)*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,T,V1,M2> operator*(const ProdVM<ix,T,V1,M2>& vm, const RT x)
  { return ProdVM<0,T,V1,M2>(x*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,CT,V1,M2> operator*(const ProdVM<ix,T,V1,M2>& vm, const CT x)
  { return ProdVM<0,CT,V1,M2>(x*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,CT,V1,M2> operator*(const ProdVM<ix,T,V1,M2>& vm, const CCT x)
  { return ProdVM<0,CT,V1,M2>(x*vm.GetX(),vm.GetV(),vm.GetM()); }

  template <int ix1, class T1, int ix, class T, class V1, class M2>
  inline ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2> operator*(
      const ProdVM<ix,T,V1,M2>& vm, const Scaling<ix1,T1>& x)
  { 
    return ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2>(
      T1(x)*vm.GetX(),vm.GetV(),vm.GetM()); 
  }

  // (vm)/x
  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,T,V1,M2> operator/(const ProdVM<ix,T,V1,M2>& vm, const int x)
  { return ProdVM<0,T,V1,M2>(vm.GetX()/RT(x),vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,T,V1,M2> operator/(const ProdVM<ix,T,V1,M2>& vm, const RT x)
  { return ProdVM<0,T,V1,M2>(vm.GetX()/x,vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,CT,V1,M2> operator/(const ProdVM<ix,T,V1,M2>& vm, const CT x)
  { return ProdVM<0,CT,V1,M2>(vm.GetX()/x,vm.GetV(),vm.GetM()); }

  template <int ix, class T, class V1, class M2>
  inline ProdVM<0,CT,V1,M2> operator/(const ProdVM<ix,T,V1,M2>& vm, const CCT x)
  { return ProdVM<0,CT,V1,M2>(vm.GetX()/x,vm.GetV(),vm.GetM()); }

  template <int ix1, class T1, int ix, class T, class V1, class M2>
  inline ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2> operator/(
      const ProdVM<ix,T,V1,M2>& vm, const Scaling<ix1,T1>& x)
  { 
    return ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2>(
      vm.GetX()/T1(x),vm.GetV(),vm.GetM()); 
  }

#undef RT
#undef CT
#undef CCT

  // TypeText

  template <int ix, class T, class M1, class V2>
  inline std::string TypeText(const ProdMV<ix,T,M1,V2>& mv)
  {
    std::ostringstream s;
    s << "ProdMV< "<<ix<<","<<TypeText(T())<<" , ";
    s << TypeText(mv.GetM())<<" , "<<TypeText(mv.GetV())<<" >";
    return s.str();
  }

  template <int ix, class T, class V1, class M2>
  inline std::string TypeText(const ProdVM<ix,T,V1,M2>& vm)
  {
    std::ostringstream s;
    s << "ProdVM< "<<ix<<","<<TypeText(T())<<" , ";
    s << TypeText(vm.GetV())<<" , "<<TypeText(vm.GetM())<<" >";
    return s.str();
  }

} // namespace tmv

#ifdef XDEBUG
#undef XDEBUG
#endif

#endif 
