

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
    inline void MultXV(
        const Scaling<ix,T>& x, const BaseVector<M1>& m1, 
        BaseVector_Mutable<M2>& m2)
    { MultXV(x,m1.calc(),m2.mat()); }

    // Also allow x to be missing (taken to be 1) or a scalar.
    template <bool add, class V1, class V2>
    inline void MultXV(
        const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    { MultXV<add>(Scaling<1,typename V2::real_type>(),v1.calc(),v2.vec()); }
    template <bool add, class T, class V1, class V2>
    inline void MultXV(
        T x, const BaseVector<V1>& v1, BaseVector_Mutable<V2>& v2)
    { MultXV<add>(Scaling<0,T>(x),v1.calc(),v2.vec()); }
    template <class T, class V1>
    inline void Scale(T x, BaseVector_Mutable<V1>& v1)
    { Scale(Scaling<0,T>(x),v1.vec()); }


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
        enum { A = _fort ? FortranStyle : CStyle };
        typedef typename VCopyHelper<value_type,_size,A>::type copy_type;
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

        ProdXV(const Scaling<ix,T>& _x, const BaseVector<V>& _v) : 
            x(_x), v(_v.vec()) {}

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const V& getV() const { return v; }

        TMV_INLINE ptrdiff_t size() const { return v.size(); }
        value_type cref(ptrdiff_t i) const 
        { return  x*v.cref(i); }

        template <class V2>
        TMV_INLINE_ND void assignTo(BaseVector_Mutable<V2>& v2) const
        {
            TMVStaticAssert((Sizes<type::_size,V2::_size>::same));
            TMVAssert(size() == v2.size());
            TMVStaticAssert(type::isreal || V2::iscomplex);
            MultXV<false>(x,v.vec(),v2.vec());
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
    { Scale(RT(x),v.vec()); }

    template <class V>
    inline void MultEq(BaseVector_Mutable<V>& v, const RT x)
    { Scale(x,v.vec()); }

    template <class V>
    inline void MultEq(BaseVector_Mutable<V>& v, const CT x)
    { Scale(x,v.vec()); }

    template <class V>
    inline void MultEq(BaseVector_Mutable<V>& v, const CCT x)
    { Scale(CT(x),v.vec()); }

    // v /= x
    template <class V>
    inline void LDivEq(BaseVector_Mutable<V>& v, const int x)
    { Scale(RT(1)/RT(x),v.vec()); }

    template <class V>
    inline void LDivEq(BaseVector_Mutable<V>& v, const RT x)
    { Scale(RT(1)/x,v.vec()); }

    template <class V>
    inline void LDivEq(BaseVector_Mutable<V>& v, const CT x)
    { Scale(RT(1)/x,v.vec()); }

    template <class V>
    inline void LDivEq(BaseVector_Mutable<V>& v, const CCT x)
    { Scale(RT(1)/CT(x),v.vec()); }

    // -v
    template <class V>
    TMV_INLINE ProdXV<-1,RT,V> operator-(const BaseVector<V>& v)
    { return ProdXV<-1,RT,V>(RT(-1),v); }

    // x * v
    template <class V>
    TMV_INLINE ProdXV<0,RT,V> operator*(const int x, const BaseVector<V>& v)
    { return ProdXV<0,RT,V>(RT(x),v); }

    template <class V>
    TMV_INLINE ProdXV<0,RT,V> operator*(const RT x, const BaseVector<V>& v)
    { return ProdXV<0,RT,V>(x,v); }

    template <class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(const CT x, const BaseVector<V>& v)
    { return ProdXV<0,CT,V>(x,v); }

    template <class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(const CCT x, const BaseVector<V>& v)
    { return CT(x)*v; }

    // v * x
    template <class V>
    TMV_INLINE ProdXV<0,RT,V> operator*(const BaseVector<V>& v, const int x)
    { return RT(x)*v; }

    template <class V>
    TMV_INLINE ProdXV<0,RT,V> operator*(const BaseVector<V>& v, const RT x)
    { return x*v; }

    template <class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(const BaseVector<V>& v, const CT x)
    { return x*v; }

    template <class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(const BaseVector<V>& v, const CCT x)
    { return CT(x)*v; }

    // v / x
    template <class V>
    TMV_INLINE ProdXV<0,RT,V> operator/(const BaseVector<V>& v, const int x)
    { return (RT(1)/RT(x))*v; }

    template <class V>
    TMV_INLINE ProdXV<0,RT,V> operator/(const BaseVector<V>& v, const RT x)
    { return (RT(1)/x)*v; }

    template <class V>
    TMV_INLINE ProdXV<0,CT,V> operator/(const BaseVector<V>& v, const CT x)
    { return (RT(1)/x)*v; }

    template <class V>
    TMV_INLINE ProdXV<0,CT,V> operator/(const BaseVector<V>& v, const CCT x)
    { return (RT(1)/CT(x))*v; }

#undef RT
#undef CT
#undef CCT

    // Consolidate x*x*v type constructs:

#define RT typename ProdXV<ix,T,V>::real_type
#define CT typename ProdXV<ix,T,V>::complex_type
#define CCT ConjRef<CT>

    // -(x*v)
    template <int ix, class T, class V>
    TMV_INLINE ProdXV<-ix,T,V> operator-(const ProdXV<ix,T,V>& pxv)
    { return ProdXV<-ix,T,V>(-pxv.getX(),pxv.getV()); }

    // x * (x*v)
    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,T,V> operator*(
        const int x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,T,V>(RT(x)*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,T,V> operator*(
        const RT x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,T,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(
        const CT x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(
        const CCT x, const ProdXV<ix,T,V>& pxv)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    // (x*v)*x
    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,T,V> operator*(const ProdXV<ix,T,V>& pxv, const int x)
    { return ProdXV<0,T,V>(RT(x)*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,T,V> operator*(const ProdXV<ix,T,V>& pxv, const RT x)
    { return ProdXV<0,T,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(const ProdXV<ix,T,V>& pxv, const CT x)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,CT,V> operator*(const ProdXV<ix,T,V>& pxv, const CCT x)
    { return ProdXV<0,CT,V>(x*pxv.getX(),pxv.getV()); }

    // (x*v)/x
    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,T,V> operator/(const ProdXV<ix,T,V>& pxv, const int x)
    {
        return ProdXV<0,T,V>(
            ZProd<false,false>::quot(pxv.getX(),RT(x)),pxv.getV()); 
    }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,T,V> operator/(const ProdXV<ix,T,V>& pxv, const RT x)
    { 
        return ProdXV<0,T,V>(
            ZProd<false,false>::quot(pxv.getX(),x),pxv.getV()); 
    }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,CT,V> operator/(const ProdXV<ix,T,V>& pxv, const CT x)
    {
        return ProdXV<0,CT,V>(
            ZProd<false,false>::quot(pxv.getX(),x),pxv.getV()); 
    }

    template <int ix, class T, class V>
    TMV_INLINE ProdXV<0,CT,V> operator/(
        const ProdXV<ix,T,V>& pxv, const CCT x)
    { 
        return ProdXV<0,CT,V>(
            ZProd<false,false>::quot(pxv.getX(),CT(x)),pxv.getV()); 
    }

#undef RT
#undef CT
#undef CCT

    // Have SameStorage look into a ProdXV object:
    template <int ix1, class T1, class V1, class V2>
    TMV_INLINE bool SameStorage(
        const ProdXV<ix1,T1,V1>& v1, const BaseVector_Calc<V2>& v2)
    { return SameStorage(v1.getV().vec(),v2.vec()); }
    template <class V1, int ix2, class T2, class V2>
    TMV_INLINE bool SameStorage(
        const BaseVector_Calc<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    { return SameStorage(v1.vec(),v2.getV().vec()); }
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    TMV_INLINE bool SameStorage(
        const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    { return SameStorage(v1.getV().vec(),v2.getV().vec()); }


    template <int ix, class T, class V>
    inline std::string TMV_Text(const ProdXV<ix,T,V>& pxv)
    {
        std::ostringstream s;
        s << "ProdXV< "<<ix<<","<<TMV_Text(T(pxv.getX()));
        s << " , "<<TMV_Text(pxv.getV())<<" >";
        return s.str();
    }


} // namespace tmv

#endif 
