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


#ifndef TMV_OProdVV_H
#define TMV_OProdVV_H

#include "TMV_ProdXV.h"
#include "TMV_SumVV.h"

//#define XDEBUG_OPRODVV

#ifdef XDEBUG_OPRODVV
#include <iostream>
#endif

namespace tmv {

  //
  // Vector ^ Vector
  //

  template <bool add, int ix, class T, class V1, class V2, class M3>
  inline void Rank1Update(const Scaling<ix,T>& x,
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2,
      BaseMatrix_Mutable<M3>& m3);

  template <bool checkalias, bool conj, bool rm, bool inst, bool add, int ix, class T, class V1, class V2, class M3>
  struct CallRank1Update // checkalias = true
  { 
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef RealType(T) RT;
      const Scaling<1,RT> one;

      if (!SameStorage(v1.vec(),m3.mat()) &&
          !SameStorage(v1.vec(),m3.mat()) )
      {
        CallRank1Update<false,conj,rm,inst,add,ix,T,V1,V2,M3>::call(x,v1,v2,m3);
      }
      else if (SameStorage(v1.vec(),m3.mat())) 
      {
        Rank1Update<add>(one, (x*v1).calc() ,v2,m3);
      }
      else if (SameStorage(v2.vec(),m3.mat())) 
      {
        TMVAssert(SameStorage(v2.vec(),m3.mat()));
        Rank1Update<add>(one,v1, (x*v2).calc() ,m3);
      }
    }
  };
  template <bool rm, bool inst, bool add, int ix, class T, class V1, class V2, class M3>
  struct CallRank1Update<false,true,rm,inst,add,ix,T,V1,V2,M3> // conj = true
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename V1::const_conjugate_type V1c;
      typedef typename V2::const_conjugate_type V2c;
      typedef typename M3::conjugate_type M3c;
      V1c v1c = v1.Conjugate();
      V2c v2c = v2.Conjugate();
      M3c m3c = m3.Conjugate();
      CallRank1Update<false,false,rm,inst,add,ix,T,V1c,V2c,M3c>::call(
          TMV_CONJ(x),v1c,v2c,m3c); 
    }
  };
  template <bool inst, bool add, int ix, class T, class V1, class V2, class M3>
  struct CallRank1Update<false,false,true,inst,add,ix,T,V1,V2,M3> // rm = true
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename M3::transpose_type M3t;
      M3t m3t = m3.Transpose();
      CallRank1Update<false,false,false,inst,add,ix,T,V2,V1,M3t>::call(
          x,v2,v1,m3t);
    }
  };
  template <bool add, int ix, class T, class V1, class V2, class M3>
  struct CallRank1Update<false,false,false,false,add,ix,T,V1,V2,M3> // inst = false
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    { InlineRank1Update<add>(x,v1,v2,m3); }
  };
  template <int ix, class T, class V1, class V2, class M3>
  struct CallRank1Update<false,false,false,true,true,ix,T,V1,V2,M3> // inst = true, add
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    { 
      typedef typename M3::value_type T3;
      InstAddRank1Update(T3(x),v1.XView(),v2.XView(),m3.XView()); 
    }
  };
  template <int ix, class T, class V1, class V2, class M3>
  struct CallRank1Update<false,false,false,true,false,ix,T,V1,V2,M3> // inst = true, !add
  {
    static inline void call(
        const Scaling<ix,T>& x, const V1& v1, const V2& v2, M3& m3)
    {
      typedef typename M3::value_type T3;
      InstRank1Update(T3(x),v1.XView(),v2.XView(),m3.XView()); 
    }
  };

  // m3 (+=) x * v1 ^ v2
  template <bool add, int ix, class T, class V1, class V2, class M3>
  inline void Rank1Update(const Scaling<ix,T>& x,
      const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2,
      BaseMatrix_Mutable<M3>& m3)
  {
    TMVStaticAssert((Sizes<M3::mcolsize,V1::vsize>::same));
    TMVStaticAssert((Sizes<M3::mrowsize,V2::vsize>::same));
    TMVAssert(m3.colsize() == v1.size());
    TMVAssert(m3.rowsize() == v2.size());

    typedef typename V1::value_type T1;
    typedef typename V2::value_type T2;
    typedef typename M3::value_type T3;
    enum { checkalias = (
        M3::mcolsize == UNKNOWN &&
        M3::mrowsize == UNKNOWN &&
        V1::vsize == UNKNOWN &&
        V2::vsize == UNKNOWN ) };
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
#ifdef XDEBUG_OPRODVV
    std::cout<<"Start Rank1Update XDEBUG"<<std::endl;
    std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
    std::cout<<"v1 = "<<TypeText(v1)<<"  "<<v1<<std::endl;
    std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3.mat()<<std::endl;
    Matrix<typename M3::value_type> m3i = m3.mat();
    Matrix<typename M3::value_type> m3c = m3.mat();
    if (!add) m3c.Zero();
    for(int i=0;i<m3.colsize();++i) {
      for(int j=0;j<m3.rowsize();++j) {
        m3c.ref(i,j) += T(x) * v1.cref(i) * v2.cref(j);
      }
    }
    std::cout<<"m3c => "<<m3c<<std::endl;
#endif
    return CallRank1Update<checkalias,M3::mconj,M3::mrowmajor,inst,add,ix,T,V1,V2,M3>::call(
        x,v1.vec(),v2.vec(),m3.mat());
#ifdef XDEBUG_OPRODVV
    if (Norm(m3.mat()-m3c) > 1.e-6 * Norm(m3c)) {
      std::cout<<"Rank1Update:  add = "<<add<<std::endl;
      std::cout<<"checkalias = "<<checkalias<<"  inst = "<<inst<<std::endl;
      std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
      std::cout<<"v1 = "<<TypeText(v1)<<"  "<<v1<<std::endl;
      std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
      std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3i<<std::endl;
      std::cout<<"m3 -> "<<m3.mat()<<std::endl;
      std::cout<<"correct = "<<m3c<<std::endl;
      exit(1);
    }
#endif
  }


  template <int ix, class T, class V1, class V2>
  class OProdVV;

  template <int ix, class T, class V1, class V2>
  struct Traits<OProdVV<ix,T,V1,V2> >
  {
    typedef typename V1::value_type vtype1;
    typedef typename V2::value_type vtype2;
    typedef typename Traits2<vtype1,vtype2>::type value_type;

    typedef ProdXV<ix,T,V1> const_col_type;
    typedef ProdXV<ix,T,V2> const_row_type;

    enum { mcolsize = V1::vsize };
    enum { mrowsize = V2::vsize };
    enum { mshape = Rec };
    enum { mfort = V1::vfort && V2::vfort };
    enum { mcalc = false };
    enum { mrowmajor = false }; // arbitrary
    enum { mcolmajor = true };

    typedef OProdVV<ix,T,V1,V2> type;
    typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,mrowmajor,mfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef typename TypeSelect<V1::vcalc&&V1::vcalc,const type,calc_type>::type eval_type;
    typedef InvalidType inverse_type;
  };

  template <int ix, class T, class V1, class V2>
  class OProdVV :
    public BaseMatrix<OProdVV<ix,T,V1,V2> >
  {
  public:

    typedef OProdVV<ix,T,V1,V2> type;
    typedef typename Traits<type>::value_type value_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mshape = Traits<type>::mshape };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    inline OProdVV(const T _x, const BaseVector<V1>& _v1, 
        const BaseVector<V2>& _v2) : x(_x), v1(_v1.vec()), v2(_v2.vec()) {}

    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const V1& GetV1() const { return v1; }
    inline const V2& GetV2() const { return v2; }

    inline size_t colsize() const { return v1.size(); }
    inline size_t rowsize() const { return v2.size(); }

    inline value_type cref(int i, int j) const
    { return x * (v1.cref(i) * v2.cref(j)); }

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
      Rank1Update<false>(x,v1.calc().CView(),v2.calc().CView(),m3cv);
    }
  private:
    const Scaling<ix,T> x;
    const V1& v1;
    const V2& v2;
  };


  // v ^ v
#define RT typename V1::real_type
  template <class V1, class V2>
  inline OProdVV<1,RT,V1,V2> operator^(
      const BaseVector<V1>& v1, const BaseVector<V2>& v2)
  { return OProdVV<1,RT,V1,V2>(RT(1),v1,v2); }
#undef RT

  // v ^ xv
  template <class V1, int ix, class T, class V2>
  inline OProdVV<ix,T,V1,V2> operator^(
      const BaseVector<V1>& v1, const ProdXV<ix,T,V2>& v2)
  { return OProdVV<ix,T,V1,V2>(v2.GetX(),v1,v2.GetV()); }

  // xv ^ v
  template <int ix, class T, class V1, class V2>
  inline OProdVV<ix,T,V1,V2> operator^(
      const ProdXV<ix,T,V1>& v1, const BaseVector<V2>& v2)
  { return OProdVV<ix,T,V1,V2>(v1.GetX(),v1.GetV(),v2); }

  // xv ^ xv
#define PT typename Traits2<T1,T2>::type
  template <int ix1, class T1, class V1, int ix2, class T2, class V2>
  inline OProdVV<ix1*ix2,PT,V1,V2> operator^(
      const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
  { return OProdVV<ix1*ix2,PT,V1,V2>(v1.GetX()*v2.GetX(),v1.GetV(),v2.GetV()); }
#undef PT


  // m += vv
  template <class M3, int ix, class T, class V1, class V2>
  inline void AddEq(
      BaseMatrix_Mutable<M3>& m, const OProdVV<ix,T,V1,V2>& vv)
  {
    typename M3::cview_type mcv = m.CView();
    Rank1Update<true>(vv.GetX(),vv.GetV1().calc().CView(),
        vv.GetV2().calc().CView(),mcv); 
  }
  
  // m -= vv
  template <class M3, int ix, class T, class V1, class V2>
  inline void SubtractEq(
      BaseMatrix_Mutable<M3>& m, const OProdVV<ix,T,V1,V2>& vv)
  { 
    typename M3::cview_type mcv = m.CView();
    Rank1Update<true>(-vv.GetX(),vv.GetV1().calc().CView(),
      vv.GetV2().calc().CView(),mcv); 
  }


  // Consolidate x*(xmv) type constructs:

#define RT typename OProdVV<ix,T,V1,V2>::real_type
#define CT typename OProdVV<ix,T,V1,V2>::complex_type
#define CCT ConjRef<CT>

  // -(x*v)
  template <int ix, class T, class V1, class V2>
  inline OProdVV<-ix,T,V1,V2> operator-(const OProdVV<ix,T,V1,V2>& vv)
  { return OProdVV<-ix,T,V1,V2>(-vv.GetX(),vv.GetV1(),vv.GetV2()); }

  // x * (x*v)
  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,T,V1,V2> operator*(
      const RT x, const OProdVV<ix,T,V1,V2>& vv)
  { return OProdVV<0,T,V1,V2>(x*vv.GetX(),vv.GetV1(),vv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,CT,V1,V2> operator*(
      const CT x, const OProdVV<ix,T,V1,V2>& vv)
  { return OProdVV<0,CT,V1,V2>(x*vv.GetX(),vv.GetV1(),vv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,CT,V1,V2> operator*(
      const CCT x, const OProdVV<ix,T,V1,V2>& vv)
  { return OProdVV<0,CT,V1,V2>(x*vv.GetX(),vv.GetV1(),vv.GetV2()); }

  template <int ix1, class T1, int ix, class T, class V1, class V2>
  inline OProdVV<ix*ix1,typename Traits2<T,T1>::type,V1,V2> operator*(
      const Scaling<ix1,T1>& x, const OProdVV<ix,T,V1,V2>& vv)
  {
    return OProdVV<ix*ix1,typename Traits2<T,T1>::type,V1,V2>(
      T1(x)*vv.GetX(),vv.GetV1(),vv.GetV2()); 
  }

  // (x*v)*x
  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,T,V1,V2> operator*(
      const OProdVV<ix,T,V1,V2>& vv, const RT x)
  { return OProdVV<0,T,V1,V2>(x*vv.GetX(),vv.GetV1(),vv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,CT,V1,V2> operator*(
      const OProdVV<ix,T,V1,V2>& vv, const CT x)
  { return OProdVV<0,CT,V1,V2>(x*vv.GetX(),vv.GetV1(),vv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,CT,V1,V2> operator*(
      const OProdVV<ix,T,V1,V2>& vv, const CCT x)
  { return OProdVV<0,CT,V1,V2>(x*vv.GetX(),vv.GetV1(),vv.GetV2()); }

  template <int ix1, class T1, int ix, class T, class V1, class V2>
  inline OProdVV<ix*ix1,typename Traits2<T,T1>::type,V1,V2> operator*(
      const OProdVV<ix,T,V1,V2>& vv, const Scaling<ix1,T1>& x)
  {
    return OProdVV<ix*ix1,typename Traits2<T,T1>::type,V1,V2>(
      T1(x)*vv.GetX(),vv.GetV1(),vv.GetV2()); 
  }

  // (x*v)/x
  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,T,V1,V2> operator/(
      const OProdVV<ix,T,V1,V2>& vv, const RT x)
  { return OProdVV<0,T,V1,V2>(vv.GetX()/x,vv.GetV1(),vv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,CT,V1,V2> operator/(
      const OProdVV<ix,T,V1,V2>& vv, const CT x)
  { return OProdVV<0,CT,V1,V2>(vv.GetX()/x,vv.GetV1(),vv.GetV2()); }

  template <int ix, class T, class V1, class V2>
  inline OProdVV<0,CT,V1,V2> operator/(
      const OProdVV<ix,T,V1,V2>& vv, const CCT x)
  { return OProdVV<0,CT,V1,V2>(vv.GetX()/x,vv.GetV1(),vv.GetV2()); }

  template <int ix1, class T1, int ix, class T, class V1, class V2>
  inline OProdVV<ix*ix1,typename Traits2<T,T1>::type,V1,V2> operator/(
      const OProdVV<ix,T,V1,V2>& vv, const Scaling<ix1,T1>& x)
  {
    return OProdVV<ix*ix1,typename Traits2<T,T1>::type,V1,V2>(
      vv.GetX()/T1(x),vv.GetV1(),vv.GetV2()); 
  }

#undef RT
#undef CT
#undef CCT

  // TypeText

  template <int ix, class T, class V1, class V2>
  inline std::string TypeText(const OProdVV<ix,T,V1,V2>& svv)
  {
    std::ostringstream s;
    s << "OProdVV< "<<ix<<","<<TypeText(T())<<",";
    s << TypeText(svv.GetV1())<<" , "<<TypeText(svv.GetV2())<<" >";
    return s.str();
  }

} // namespace tmv

#endif 
