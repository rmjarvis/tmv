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


#ifndef TMV_ProdMM_H
#define TMV_ProdMM_H

#include "TMV_ProdXM.h"

//#define XDEBUG_PRODMM

#ifdef XDEBUG_PRODMM
#include <iostream>
#endif

namespace tmv {

  //
  // Matrix * Matrix
  //

  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void MultMM(const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
      const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);

  template <bool checkalias, bool conj, bool rm, bool inst, bool add, int ix, class T, class M1, class M2, class M3> struct CallMultMM;
  template <bool conj, bool rm, bool inst, bool add, int ix, class T, class M1, class M2, class M3> 
  struct CallMultMM<true,conj,rm,inst,add,ix,T,M1,M2,M3> // checkalias
  {
    static void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      //std::cout<<"CallMultMM checkalias:\n";
      //std::cout<<"m1 = "<<TypeText(m1)<<std::endl;
      //std::cout<<"m2 = "<<TypeText(m2)<<std::endl;
      //std::cout<<"m3 = "<<TypeText(m3)<<std::endl;
      //std::cout<<"m1.cptr = "<<m1.cptr()<<std::endl;
      //std::cout<<"m2.cptr = "<<m2.cptr()<<std::endl;
      //std::cout<<"m3.cptr = "<<m3.cptr()<<std::endl;
      //std::cout<<"SameStorage(m1,m3) = "<<SameStorage(m1,m3)<<std::endl;
      //std::cout<<"SameStorage(m2,m3) = "<<SameStorage(m2,m3)<<std::endl;
      // Check for m1 or m2 being aliases of m3.
      // Also check that all matrices are either rowmajor or colmajor here,
      // since they should use a temporary in that case too.
      // TODO: Make the temporary matrix in blocks if possible.
      typedef RealType(T) RT;
      const Scaling<1,RT> one;
      if (!SameStorage(m1,m3) && !SameStorage(m2,m3) )
      {
        //std::cout<<"not same storage"<<std::endl;
        CallMultMM<false,conj,rm,inst,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
      }
      else if (SameStorage(m1,m3))
      {
        //std::cout<<"m1,m2 same storage"<<std::endl;
        MultMM<add>(one, (x*m1).calc().CView(), m2, m3);
      }
      else  // SameStorage(m2,m3) 
      { 
        //std::cout<<"m2,m3 same storage"<<std::endl;
        MultMM<add>(one, m1, (x*m2).calc().CView(), m3); 
      }
    }
  };
  template <bool rm, bool inst, bool add, int ix, class T, class M1, class M2, class M3>
  struct CallMultMM<false,true,rm,inst,add,ix,T,M1,M2,M3> // conj = true
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M1::const_conjugate_type M1c;
      typedef typename M2::const_conjugate_type M2c;
      typedef typename M3::conjugate_type M3c;
      M1c m1c = m1.Conjugate();
      M2c m2c = m2.Conjugate();
      M3c m3c = m3.Conjugate();
      CallMultMM<false,false,rm,inst,add,ix,T,M1c,M2c,M3c>::call(
          TMV_CONJ(x),m1c,m2c,m3c);
    }
  };
  template <bool inst, bool add, int ix, class T, class M1, class M2, class M3>
  struct CallMultMM<false,false,true,inst,add,ix,T,M1,M2,M3> // rm = true
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M1::const_transpose_type M1t;
      typedef typename M2::const_transpose_type M2t;
      typedef typename M3::transpose_type M3t;
      M1t m1t = m1.Transpose();
      M2t m2t = m2.Transpose();
      M3t m3t = m3.Transpose();
      CallMultMM<false,false,false,inst,add,ix,T,M2t,M1t,M3t>::call(
          x,m2t,m1t,m3t);
    }
  };
  template <bool add, int ix, class T, class M1, class M2, class M3>
  struct CallMultMM<false,false,false,false,add,ix,T,M1,M2,M3> // inst = false
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { InlineMultMM<add>(x,m1,m2,m3); }
  };
  template <int ix, class T, class M1, class M2, class M3>
  struct CallMultMM<false,false,false,true,true,ix,T,M1,M2,M3> // add
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    { 
      typedef typename M3::value_type T3;
      typename M3::xview_type m3x = m3.XView();
      InstAddMultMM(T3(x),m1.XView(),m2.XView(),m3x); 
    }
  };
  template <int ix, class T, class M1, class M2, class M3>
  struct CallMultMM<false,false,false,true,false,ix,T,M1,M2,M3> // !add
  {
    static inline void call(
        const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
    {
      typedef typename M3::value_type T3;
      typename M3::xview_type m3x = m3.XView();
      InstMultMM(T3(x),m1.XView(),m2.XView(),m3x); 
    }
  };

  // m3 (+=) x * m1 * m2
  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void MultMM(
      const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
      const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same));
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mcolsize>::same));
    TMVStaticAssert((Sizes<M2::mrowsize,M3::mrowsize>::same));
    TMVAssert(m1.colsize() == m3.colsize());
    TMVAssert(m1.rowsize() == m2.colsize());
    TMVAssert(m2.rowsize() == m3.rowsize());

    typedef typename M1::value_type T1;
    typedef typename M2::value_type T2;
    typedef typename M3::value_type T3;
    enum { checkalias = (
        M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
        M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN &&
        M3::mcolsize == UNKNOWN && M3::mrowsize == UNKNOWN ) };
    enum { conj = M3::mconj };
    enum { rm = M3::mrowmajor && !M3::mcolmajor };
    enum { inst = (
        Traits<T1>::isinst &&
        Traits<T2>::isinst &&
        Traits<T3>::isinst &&
#ifdef TMV_INST_MIX
        Traits2<T1,T3>::samebase &&
        Traits2<T2,T3>::samebase &&
#else
        Traits2<T1,T3>::sametype &&
        Traits2<T2,T3>::sametype &&
#endif
        checkalias ) };
#ifdef XDEBUG_PRODMM
    //std::cout<<"Start MultMM XDEBUG"<<std::endl;
    //std::cout<<"add = "<<add<<std::endl;
    //std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
    //std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
    //std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
    //std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3.mat()<<std::endl;
    Matrix<typename M3::value_type> m3i = m3.mat();
    Matrix<typename M3::value_type> m3c = m3.mat();
    if (!add) m3c.Zero();
    for(int i=0;i<m3.colsize();++i) {
      for(int j=0;j<m3.rowsize();++j) {
        for(int k=0;k<m1.rowsize();++k) {
          m3c.ref(i,j) += T(x) * m1.cref(i,k) * m2.cref(k,j);
        }
      }
    }
    //std::cout<<"m3c => "<<m3c<<std::endl;
#endif
    CallMultMM<checkalias,conj,rm,inst,add,ix,T,M1,M2,M3>::call(
        x,m1.mat(),m2.mat(),m3.mat());
#ifdef XDEBUG_PRODMM
    if (Norm(m3.mat()-m3c) > 1.e-6 * Norm(m3c)) {
      std::cout<<"MultMM:  add = "<<add<<std::endl;
      std::cout<<"checkalias = "<<checkalias<<"  inst = "<<inst<<std::endl;
      std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
      std::cout<<"m1 = "<<TypeText(m1)<<"  "<<m1<<std::endl;
      std::cout<<"m2 = "<<TypeText(m2)<<"  "<<m2<<std::endl;
      std::cout<<"m3 = "<<TypeText(m3)<<"  "<<m3i<<std::endl;
      std::cout<<"m3 -> "<<m3.mat()<<std::endl;
      std::cout<<"correct = "<<m3c<<std::endl;
      std::cout<<"diff = "<<(m3c-m3.mat())<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(m3c-m3.mat())<<std::endl;
      exit(1);
    }
#endif
  }


  // Here we just need to check whether the shape of M2 requires 
  // that we copy m1 to new storage for a m *= m operation.
  template <bool copy, class M1, int ix, class T, class M2> struct CallMultEqMM;

  template <class M1, int ix, class T, class M2> 
  struct CallMultEqMM<true,M1,ix,T,M2> // copy = true
  {
    static inline void call(
        M1& m1, const Scaling<ix,T>& x, const M2& m2)
    { MultMM<false>(x,m1.copy().CView(),m2,m1); }
  };

  template <class M1, int ix, class T, class M2> 
  struct CallMultEqMM<false,M1,ix,T,M2> // copy = false
  {
    static inline void call(
        M1& m1, const Scaling<ix,T>& x, const M2& m2)
    { MultMM<false>(x,m1,m2,m1); }
  };

  template <class M1, int ix, class T, class M2>
  inline void MultEqMM(BaseMatrix_Mutable<M1>& m1,
      const Scaling<ix,T>& x, const BaseMatrix_Calc<M2>& m2)
  {
    enum { copy = (
        !ShapeTraits<M2::mshape>::upper && 
        !ShapeTraits<M2::mshape>::lower ) };
    CallMultEqMM<copy,M1,ix,T,M2>::call(m1.mat(),x,m2.mat());
  }

  template <int ix, class T, class M1, class M2>
  class ProdMM;

  template <int ix, class T, class M1, class M2>
  struct Traits<ProdMM<ix,T,M1,M2> >
  {
    typedef typename M1::value_type mtype1;
    typedef typename M2::value_type mtype2;
    typedef typename Traits2<mtype1,mtype2>::type value_type;

    enum { mcolsize = M1::mcolsize };
    enum { mrowsize = M2::mrowsize };
    enum { mshape = ShapeTraits2<M1::mshape,M2::mshape>::prod };
    enum { mfort = M1::mfort && M2::mfort };
    enum { mcalc = false };
    enum { rm1 = Traits<typename M1::calc_type>::mrowmajor };
    enum { rm2 = Traits<typename M2::calc_type>::mrowmajor };
    enum { cm1 = Traits<typename M1::calc_type>::mcolmajor };
    enum { cm2 = Traits<typename M2::calc_type>::mcolmajor };
    enum { mrowmajor = (
        ( (rm1 && rm2) || (!cm1 && !cm2 && mrowsize > int(mcolsize)) ) ) };

    typedef ProdMM<ix,T,M1,M2> type;
    typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,mrowmajor,mfort>::type copy_type;
    typedef const copy_type calc_type;
    typedef const copy_type eval_type;
    typedef InvalidType inverse_type;
  };

  template <int ix, class T, class M1, class M2>
  class ProdMM :
    public BaseMatrix<ProdMM<ix,T,M1,M2> >
  {
  public:

    typedef ProdMM<ix,T,M1,M2> type;
    typedef typename Traits<type>::value_type value_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mshape = Traits<type>::mshape };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    inline ProdMM(const T& _x, const BaseMatrix<M1>& _m1, 
        const BaseMatrix<M2>& _m2) : x(_x), m1(_m1.mat()), m2(_m2.mat())
    {
      TMVStaticAssert((Sizes<M1::mrowsize,M2::mcolsize>::same)); 
      TMVAssert(m1.rowsize() == m2.colsize());
    }

    inline const Scaling<ix,T>& GetX() const { return x; }
    inline const M1& GetM1() const { return m1; }
    inline const M2& GetM2() const { return m2; }

    inline size_t colsize() const { return m1.colsize(); }
    inline size_t rowsize() const { return m2.rowsize(); }

    inline value_type cref(int i, int j) const
    { return x * (m1.get_row(i) * m2.get_col(j)); }

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
      MultMM<false>(x,m1.calc().CView(),m2.calc().CView(),m3cv);
    }
  private:
    const Scaling<ix,T> x;
    const M1& m1;
    const M2& m2;
  };


  // m * m
#define RT typename M1::real_type
  template <class M1, class M2>
  inline ProdMM<1,RT,M1,M2> operator*(
      const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
  { return ProdMM<1,RT,M1,M2>(RT(1),m1,m2); }
#undef RT

  // m * xm
  template <class M1, int ix, class T, class M2>
  inline ProdMM<ix,T,M1,M2> operator*(
      const BaseMatrix<M1>& m1, const ProdXM<ix,T,M2>& m2)
  { return ProdMM<ix,T,M1,M2>(m2.GetX(),m1,m2.GetM()); }

  // xm * m
  template <int ix, class T, class M1, class M2>
  inline ProdMM<ix,T,M1,M2> operator*(
      const ProdXM<ix,T,M1>& m1, const BaseMatrix<M2>& m2)
  { return ProdMM<ix,T,M1,M2>(m1.GetX(),m1.GetM(),m2); }

  // xm * xm
#define PT typename Traits2<T1,T2>::type
  template <int ix1, class T1, class M1, int ix2, class T2, class M2>
  inline ProdMM<ix1*ix2,PT,M1,M2> operator*(
      const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  { return ProdMM<ix1*ix2,PT,M1,M2>(m1.GetX()*m2.GetX(),m1.GetM(),m2.GetM()); }
#undef PT


  // m *= m
#define RT typename M1::real_type
  template <class M1, class M2>
  inline void MultEq(
      BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
  {
    typename M1::cview_type m1cv = m1.CView();
    MultEqMM(m1cv,Scaling<1,RT>(),m2.calc().CView()); 
  }
#undef RT

  // m *= xm
  template <class M1, int ix2, class T2, class M2>
  inline void MultEq(
      BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
  {
    typename M1::cview_type m1cv = m1.CView();
    MultEqMM(m1cv,m2.GetX(),m2.GetM().calc().CView()); 
  }

  // m += mm
  template <class M3, int ix, class T, class M1, class M2>
  inline void AddEq(
      BaseMatrix_Mutable<M3>& m3, const ProdMM<ix,T,M1,M2>& mm)
  {
    typename M3::cview_type m3cv = m3.CView();
    MultMM<true>(mm.GetX(),mm.GetM1().calc().CView(),
        mm.GetM2().calc().CView(),m3cv); 
  }
  
  // m -= mm
  template <class M3, int ix, class T, class M1, class M2>
  inline void SubtractEq(
      BaseMatrix_Mutable<M3>& m3, const ProdMM<ix,T,M1,M2>& mm)
  { 
    typename M3::cview_type m3cv = m3.CView();
    MultMM<true>(-mm.GetX(),mm.GetM1().calc().CView(),
      mm.GetM2().calc().CView(),m3cv); 
  }


  // Consolidate x*(xmm) type constructs:

#define RT typename ProdMM<ix,T,M1,M2>::real_type
#define CT typename ProdMM<ix,T,M1,M2>::complex_type
#define CCT ConjRef<CT>

  // -(mm)
  template <int ix, class T, class M1, class M2>
  inline ProdMM<-ix,T,M1,M2> operator-(const ProdMM<ix,T,M1,M2>& mm)
  { return ProdMM<-ix,T,M1,M2>(-mm.GetX(),mm.GetM1(),mm.GetM2()); }

  // x * (mm)
  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,T,M1,M2> operator*(const int x, const ProdMM<ix,T,M1,M2>& mm)
  { return ProdMM<0,T,M1,M2>(RT(x)*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,T,M1,M2> operator*(const RT x, const ProdMM<ix,T,M1,M2>& mm)
  { return ProdMM<0,T,M1,M2>(x*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,CT,M1,M2> operator*(const CT x, const ProdMM<ix,T,M1,M2>& mm)
  { return ProdMM<0,CT,M1,M2>(x*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,CT,M1,M2> operator*(const CCT x, const ProdMM<ix,T,M1,M2>& mm)
  { return ProdMM<0,CT,M1,M2>(x*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix1, class T1, int ix, class T, class M1, class M2>
  inline ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
      const Scaling<ix1,T1>& x, const ProdMM<ix,T,M1,M2>& mm)
  { 
    return ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
      T1(x)*mm.GetX(),mm.GetM1(),mm.GetM2()); 
  }

  // (mm)*x
  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,T,M1,M2> operator*(const ProdMM<ix,T,M1,M2>& mm, const int x)
  { return ProdMM<0,T,M1,M2>(RT(x)*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,T,M1,M2> operator*(const ProdMM<ix,T,M1,M2>& mm, const RT x)
  { return ProdMM<0,T,M1,M2>(x*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,CT,M1,M2> operator*(const ProdMM<ix,T,M1,M2>& mm, const CT x)
  { return ProdMM<0,CT,M1,M2>(x*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,CT,M1,M2> operator*(const ProdMM<ix,T,M1,M2>& mm, const CCT x)
  { return ProdMM<0,CT,M1,M2>(x*mm.GetX(),mm.GetM1(),mm.GetM2()); }

  template <int ix1, class T1, int ix, class T, class M1, class M2>
  inline ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
      const ProdMM<ix,T,M1,M2>& mm, const Scaling<ix1,T1>& x)
  { 
    return ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
      T1(x)*mm.GetX(),mm.GetM1(),mm.GetM2()); 
  }

  // (mm)/x
  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,T,M1,M2> operator/(const ProdMM<ix,T,M1,M2>& mm, const int x)
  { return ProdMM<0,T,M1,M2>(mm.GetX()/RT(x),mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,T,M1,M2> operator/(const ProdMM<ix,T,M1,M2>& mm, const RT x)
  { return ProdMM<0,T,M1,M2>(mm.GetX()/x,mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,CT,M1,M2> operator/(const ProdMM<ix,T,M1,M2>& mm, const CT x)
  { return ProdMM<0,CT,M1,M2>(mm.GetX()/x,mm.GetM1(),mm.GetM2()); }

  template <int ix, class T, class M1, class M2>
  inline ProdMM<0,CT,M1,M2> operator/(const ProdMM<ix,T,M1,M2>& mm, const CCT x)
  { return ProdMM<0,CT,M1,M2>(mm.GetX()/x,mm.GetM1(),mm.GetM2()); }

  template <int ix1, class T1, int ix, class T, class M1, class M2>
  inline ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator/(
      const ProdMM<ix,T,M1,M2>& mm, const Scaling<ix1,T1>& x)
  { 
    return ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
      mm.GetX()/T1(x),mm.GetM1(),mm.GetM2()); 
  }

#undef RT
#undef CT
#undef CCT

  // TypeText

  template <int ix, class T, class M1, class M2>
  inline std::string TypeText(const ProdMM<ix,T,M1,M2>& mm)
  {
    std::ostringstream s;
    s << "ProdMM< "<<ix<<","<<TypeText(T())<<",";
    s << TypeText(mm.GetM1())<<" , "<<TypeText(mm.GetM2())<<" >";
    return s.str();
  }

} // namespace tmv

#endif 
