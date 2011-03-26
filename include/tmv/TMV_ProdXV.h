///////////////////////////////////////////////////////////////////////////////
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
#include "TMV_MultXV_Funcs.h"

namespace tmv {

    //
    // Scalar * Vector
    //

    // These first few are for when an argument is a composite vector
    // and needs to be calculated before running MultXV.
    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXV(
        const Scaling<ix,T>& x, const BaseVector<M1>& m1, 
        BaseVector_Mutable<M2>& m2)
    { MultXV(x,m1.calc(),m2.mat()); }
    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXV(
        const Scaling<ix,T>& x, const BaseVector<M1>& m1, 
        BaseVector_Mutable<M2>& m2)
    { NoAliasMultXV(x,m1.calc(),m2.mat()); }
    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXV(
        const Scaling<ix,T>& x, const BaseVector<M1>& m1, 
        BaseVector_Mutable<M2>& m2)
    { AliasMultXV(x,m1.calc(),m2.mat()); }

    // These are helpers to allow the caller to not use a Scaling object.
    template <bool add, class T, class V1, class V2>
    static inline void MultXV(
        const T& x, const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    { MultXV<add>(Scaling<0,T>(x),v1.vec(),v2.vec()); }
    template <bool add, class T, class V1, class V2>
    static inline void NoAliasMultXV(
        const T& x, const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    { NoAliasMultXV<add>(Scaling<0,T>(x),v1.vec(),v2.vec()); }
    template <bool add, class T, class V1, class V2>
    static inline void AliasMultXV(
        const T& x, const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    { AliasMultXV<add>(Scaling<0,T>(x),v1.vec(),v2.vec()); }

    template <bool add, class V1, class V2>
    static inline void MultXV(
        const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    { MultXV<add>(Scaling<1,typename V2::real_type>(),v1.vec(),v2.vec()); }
    template <bool add, class V1, class V2>
    static inline void NoAliasMultXV(
        const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        NoAliasMultXV<add>(
            Scaling<1,typename V2::real_type>(),v1.vec(),v2.vec()); 
    }
    template <bool add, class V1, class V2>
    static inline void AliasMultXV(
        const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        AliasMultXV<add>(
            Scaling<1,typename V2::real_type>(),v1.vec(),v2.vec()); 
    }

    template <class T, class V>
    static inline void Scale(const T& x, BaseVector_Mutable<V>& v)
    { Scale(Scaling<0,T>(x),v); }

    template <int ix, class T, class V>
    class ProdXV;

    template <int ix, class T, class V>
    struct Traits<ProdXV<ix,T,V> >
    {
        typedef typename Traits2<T,typename V::value_type>::type value_type;
        enum { _size = V::_size };
        enum { _fort = V::_fort };
        enum { _calc = false };
        typedef ProdXV<ix,T,V> type;
        typedef typename VCopyHelper<value_type,_size,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<V::_calc,const type,calc_type>::type 
            eval_type;
    };

    template <int ix, class T, class V>
    class ProdXV : 
        public BaseVector<ProdXV<ix,T,V> >
    {
    public:

        typedef ProdXV<ix,T,V> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ProdXV(const T _x, const BaseVector<V>& _v) : 
            x(_x), v(_v.vec()) {}
        const Scaling<ix,T>& getX() const { return x; }
        const V& getV() const { return v; }

        size_t size() const { return v.size(); }
        value_type cref(int i) const 
        { return  x*v.cref(i); }

        template <class V2>
        void assignTo(BaseVector_Mutable<V2>& v2) const
        {
            TMVStaticAssert((Sizes<type::_size,V2::_size>::same));
            TMVAssert(size() == v2.size());
            TMVStaticAssert(type::isreal || V2::iscomplex);
            MultXV<false>(x,v.vec(),v2.vec());
        }

        template <class V2>
        void newAssignTo(BaseVector_Mutable<V2>& v2) const
        {
            TMVStaticAssert((Sizes<type::_size,V2::_size>::same));
            TMVAssert(size() == v2.size());
            TMVStaticAssert(type::isreal || V2::iscomplex);
            NoAliasMultXV<false>(x,v.vec(),v2.vec());
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
    static inline void MultEq(BaseVector_Mutable<V>& v, const int x)
    { Scale(RT(x),v.vec()); }

    template <class V>
    static inline void MultEq(BaseVector_Mutable<V>& v, const RT x)
    { Scale(x,v.vec()); }

    template <class V>
    static inline void MultEq(BaseVector_Mutable<V>& v, const CT x)
    { Scale(x,v.vec()); }

    template <class V>
    static inline void MultEq(BaseVector_Mutable<V>& v, const CCT x)
    { Scale(CT(x),v.vec()); }

    template <class V, class T>
    static inline void MultEq(BaseVector_Mutable<V>& v, const Scaling<0,T> x)
    { Scale(x,v.vec()); }

    template <class V, class T>
    static inline void MultEq(BaseVector_Mutable<V>& v, const Scaling<1,T> x)
    {}

    template <class V, class T>
    static inline void MultEq(BaseVector_Mutable<V>& v, const Scaling<-1,T> x)
    { Scale(x,v.vec()); }

    // v /= x
    template <class V>
    static inline void DivEq(BaseVector_Mutable<V>& v, const int x)
    { Scale(RT(1)/RT(x),v.vec()); }

    template <class V>
    static inline void DivEq(BaseVector_Mutable<V>& v, const RT x)
    { Scale(RT(1)/x,v.vec()); }

    template <class V>
    static inline void DivEq(BaseVector_Mutable<V>& v, const CT x)
    { Scale(RT(1)/x,v.vec()); }

    template <class V>
    static inline void DivEq(BaseVector_Mutable<V>& v, const CCT x)
    { Scale(RT(1)/CT(x),v.vec()); }

    template <class V, class T>
    static inline void DivEq(BaseVector_Mutable<V>& v, const Scaling<0,T> x)
    { Scale(RT(1)/T(x),v.vec()); }

    template <class V, class T>
    static inline void DivEq(BaseVector_Mutable<V>& v, const Scaling<1,T> x)
    {}

    template <class V, class T>
    static inline void DivEq(BaseVector_Mutable<V>& v, const Scaling<-1,T> x)
    { Scale(x,v.vec()); }

    // -v
    template <class V>
    static inline ProdXV<-1,RT,V> operator-(const BaseVector<V>& v)
    { return ProdXV<-1,RT,V>(RT(-1),v); }

    // x * v
    template <class V>
    static inline ProdXV<0,RT,V> operator*(const int x, const BaseVector<V>& v)
    { return ProdXV<0,RT,V>(RT(x),v); }

    template <class V>
    static inline ProdXV<0,RT,V> operator*(const RT x, const BaseVector<V>& v)
    { return ProdXV<0,RT,V>(x,v); }

    template <class V>
    static inline ProdXV<0,CT,V> operator*(const CT x, const BaseVector<V>& v)
    { return ProdXV<0,CT,V>(x,v); }

    template <class V>
    static inline ProdXV<0,CT,V> operator*(const CCT x, const BaseVector<V>& v)
    { return CT(x)*v; }

    template <class V, int ix, class T>
    static inline ProdXV<ix,T,V> operator*(
        const Scaling<ix,T> x, const BaseVector<V>& v)
    { return ProdXV<ix,T,V>(T(x),v); }

    // v * x
    template <class V>
    static inline ProdXV<0,RT,V> operator*(const BaseVector<V>& v, const int x)
    { return RT(x)*v; }

    template <class V>
    static inline ProdXV<0,RT,V> operator*(const BaseVector<V>& v, const RT x)
    { return x*v; }

    template <class V>
    static inline ProdXV<0,CT,V> operator*(const BaseVector<V>& v, const CT x)
    { return x*v; }

    template <class V>
    static inline ProdXV<0,CT,V> operator*(const BaseVector<V>& v, const CCT x)
    { return CT(x)*v; }

    template <class V, int ix, class T>
    static inline ProdXV<ix,T,V> operator*(
        const BaseVector<V>& v, const Scaling<ix,T> x)
    { return ProdXV<ix,T,V>(T(x),v); }

    // v / x
    template <class V>
    static inline ProdXV<0,RT,V> operator/(const BaseVector<V>& v, const int x)
    { return (RT(1)/RT(x))*v; }

    template <class V>
    static inline ProdXV<0,RT,V> operator/(const BaseVector<V>& v, const RT x)
    { return (RT(1)/x)*v; }

    template <class V>
    static inline ProdXV<0,CT,V> operator/(const BaseVector<V>& v, const CT x)
    { return (RT(1)/x)*v; }

    template <class V>
    static inline ProdXV<0,CT,V> operator/(const BaseVector<V>& v, const CCT x)
    { return (RT(1)/CT(x))*v; }

    template <class V, int ix, class T>
    static inline ProdXV<ix,T,V> operator/(
        const BaseVector<V>& v, const Scaling<ix,T> x)
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
    static inline ProdXV<-ix,T,V> operator-(const ProdXV<ix,T,V>& pxv)
    { return ProdXV<-ix,T,V>(-pxv.getX(),pxv.getV()); }

    // x * (x*v)
    template <int ix, class T, class V>
    static inline ProdXV<0,T,V> operator*(
        const int x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,T,V>(RT(x)*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,T,V> operator*(
        const RT x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,T,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,CT,V> operator*(
        const CT x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,CT,V> operator*(
        const CCT x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix1, class T1, int ix, class T, class V>
    static inline ProdXV<ix1*ix,typename Traits2<T,T1>::type,V> operator*(
        const Scaling<ix1,T1>& x, const ProdXV<ix,T,V>& pxv)
    {
        return ProdXV<ix1*ix,typename Traits2<T,T1>::type,V>(
            T1(x)*pxv.getX(),pxv.getV()); 
    }

    // (x*v)*x
    template <int ix, class T, class V>
    static inline ProdXV<0,T,V> operator*(
        const ProdXV<ix,T,V>& pxv, const int x)
    { return ProdXV<0,T,V>(RT(x)*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,T,V> operator*(
        const ProdXV<ix,T,V>& pxv, const RT x)
    { return ProdXV<0,T,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,CT,V> operator*(
        const ProdXV<ix,T,V>& pxv, const CT x)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,CT,V> operator*(
        const ProdXV<ix,T,V>& pxv, const CCT x)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix1, class T1, int ix, class T, class V>
    static inline ProdXV<ix1*ix,typename Traits2<T,T1>::type,V> operator*(
        const ProdXV<ix,T,V>& pxv, const Scaling<ix1,T1>& x)
    {
        return ProdXV<ix1*ix,typename Traits2<T,T1>::type,V>(
            T1(x)*pxv.getX(),pxv.getV()); 
    }

    // (x*v)/x
    template <int ix, class T, class V>
    static inline ProdXV<0,T,V> operator/(
        const ProdXV<ix,T,V>& pxv, const int x)
    { return ProdXV<0,T,V>(pxv.getX()/RT(x),pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,T,V> operator/(
        const ProdXV<ix,T,V>& pxv, const RT x)
    { return ProdXV<0,T,V>(pxv.getX()/x,pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,CT,V> operator/(
        const ProdXV<ix,T,V>& pxv, const CT x)
    { return ProdXV<0,CT,V>(pxv.getX()/x,pxv.getV()); }

    template <int ix, class T, class V>
    static inline ProdXV<0,CT,V> operator/(
        const ProdXV<ix,T,V>& pxv, const CCT x)
    { return ProdXV<0,CT,V>(pxv.getX()/x,pxv.getV()); }

    template <int ix1, class T1, int ix, class T, class V>
    static inline ProdXV<ix1*ix,typename Traits2<T,T1>::type,V> operator/(
        const ProdXV<ix,T,V>& pxv, const Scaling<ix1,T1>& x)
    {
        return ProdXV<ix1*ix,typename Traits2<T,T1>::type,V>(
            pxv.getX()/T1(x),pxv.getV()); 
    }

#undef RT
#undef CT
#undef CCT

#ifndef TMV_NO_ALIAS_CHECK
    // Have SameStorage look into a ProdXV object:
    template <int ix1, class T1, class V1, class V2>
    static inline bool SameStorage(
        const ProdXV<ix1,T1,V1>& v1, const BaseVector_Calc<V2>& v2)
    { return SameStorage(v1.getV().vec(),v2.vec()); }
    template <class V1, int ix2, class T2, class V2>
    static inline bool SameStorage(
        const BaseVector_Calc<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    { return SameStorage(v1.vec(),v2.getV().vec()); }
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline bool SameStorage(
        const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    { return SameStorage(v1.getV().vec(),v2.getV().vec()); }
#endif


#ifdef TMV_DEBUG
    template <int ix, class T, class V>
    static inline std::string TMV_Text(const ProdXV<ix,T,V>& pxv)
    {
        std::ostringstream s;
        s << "ProdXV< "<<ix<<","<<TMV_Text(T(pxv.getX()));
        s << " , "<<TMV_Text(pxv.getV())<<" >";
        return s.str();
    }
#endif


} // namespace tmv

#endif 
