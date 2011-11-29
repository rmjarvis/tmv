

#ifndef TMV_ElemProdVV_H
#define TMV_ElemProdVV_H

#include "TMV_BaseVector.h"

namespace tmv {

    //
    // ElementProd functions:
    // (Element-wise Vector * Vector)
    //

    template <bool add, class T, class V1, class V2, class V3>
    static inline void ElemMultVV(
        const T& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    { ElemMultVV(Scaling<0,T>(x1),v1.vec(),v2.vec(),v3.vec()); }
    template <bool add, class T, class V1, class V2, class V3>
    static inline void NoAliasElemMultVV(
        const T& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    { NoAliasElemMultVV(Scaling<0,T>(x1),v1.vec(),v2.vec(),v3.vec()); }
    template <bool add, class T, class V1, class V2, class V3>
    static inline void InlineElemMultVV(
        const T& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    { InlineElemMultVV(Scaling<0,T>(x1),v1.vec(),v2.vec(),v3.vec()); }
    template <bool add, class T, class V1, class V2, class V3>
    static inline void AliasElemMultVV(
        const T& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    { AliasElemMultVV(Scaling<0,T>(x1),v1.vec(),v2.vec(),v3.vec()); }

    template <bool add, class V1, class V2, class V3>
    static inline void ElemMultVV(
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        ElemMultVV(
            Scaling<1,typename V3::real_type(),v1.vec(),v2.vec(),v3.vec()); 
    }
    template <bool add, class V1, class V2, class V3>
    static inline void NoAliasElemMultVV(
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        NoAliasElemMultVV(
            Scaling<1,typename V3::real_type(),v1.vec(),v2.vec(),v3.vec()); 
    }
    template <bool add, class V1, class V2, class V3>
    static inline void InlineElemMultVV(
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        InlineElemMultVV(
            Scaling<1,typename V3::real_type(),v1.vec(),v2.vec(),v3.vec()); 
    }
    template <bool add, class V1, class V2, class V3>
    static inline void AliasElemMultVV(
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        AliasElemMultVV(
            Scaling<1,typename V3::real_type(),v1.vec(),v2.vec(),v3.vec()); 
    }


    template <int ix, class T, class V1, class V2>
    class ElemProdVV;

    template <int ix, class T, class V1, class V2>
    struct Traits<ElemProdVV<ix,T,V1,V2> >
    {
        typedef typename ProdXV<ix,T,V1>::value_type vtype1;
        typedef typename V2::value_type vtype2;
        typedef typename Traits2<vtype1,vtype2>::type value_type;

        enum { _size = Sizes<V1::_size,V2::_size>::size };
        enum { _fort = V1::_fort && V2::_fort };
        enum { _calc = false };

        typedef ElemProdVV<ix,T,V1,V2> type;
        typedef typename VCopyHelper<value_type,_size,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<V1::_calc&&V2::_calc,
                const type,calc_type>::type eval_type;
    };

    template <int ix, class T, class V1, class V2>
    class ElemProdVV : 
        public BaseVector<ElemProdVV<ix,T,V1,V2> >
    {
    public:

        typedef ElemProdVV<ix,T,V1,V2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ElemProdVV(
            const T _x, const BaseVector<V1>& _v1,
            const BaseVector<V2>& _v2) : x(_x), v1(_v1.vec()), v2(_v2.vec())
        {
            TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
            TMVAssert(v1.size() == v2.size());
        }

        const Scaling<ix,T>& getX() const { return x; }
        const V1& getV1() const { return v1; }
        const V2& getV2() const { return v2; }

        size_t size() const { return v1.size(); }

        value_type cref(int i) const
        { return x * (v1.cref(i) * v2.cref(i)); }

        template <class V3>
        void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same));
            TMVAssert(size() == v3.size());
            TMVStaticAssert(type::isreal || V3::iscomplex);
            ElemMultVV<false>(x,v1.calc(),v2.calc(),v3.vec());
        }

        template <class V3>
        void newAssignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same));
            TMVAssert(size() == v3.size());
            TMVStaticAssert(type::isreal || V3::iscomplex);
            NoAliasElemMultVV<false>(x,v1.calc(),v2.calc(),v3.vec());
        }

    private:
        const Scaling<ix,T> x;
        const V1& v1;
        const V2& v2;
    };

#define RT typename V2::real_type
    template <class V1, class V2>
    static TMV_INLINE_ND ElemProdVV<1,RT,V1,V2> ElementProd(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        return ElemProdVV<1,RT,V1,V2>(RT(1),v1,v2);
    }
#undef RT

    // v += [vv]
    template <class V3, int ix, class T, class V1, class V2>
    static TMV_INLINE void AddEq(
        BaseVector_Mutable<V3>& v3, const ElemProdVV<ix,T,V1,V2>& vv)
    {
        ElemMultVV<true>(
            vv.getX(),vv.getV1().calc(),vv.getV2().calc(),v3.vec()); 
    }

    // v -= [vv]
    template <class V3, int ix, class T, class V1, class V2>
    static TMV_INLINE void SubtractEq(
        BaseVector_Mutable<V3>& v3, const ElemProdVV<ix,T,V1,V2>& vv)
    {
        ElemMultVV<true>(
            -vv.getX(),vv.getV1().calc(),vv.getV2().calc(),v3.vec()); 
    }

    // Consolidate x*(xvv) type constructs:

#define RT typename ElemProdVV<ix,T,V1,V2>::real_type
#define CT typename ElemProdVV<ix,T,V1,V2>::complex_type
#define CCT ConjRef<CT>

    // -[vv]
    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<-ix,T,V1,V2> operator-(
        const ElemProdVV<ix,T,V1,V2>& mv)
    { return ElemProdVV<-ix,T,V1,V2>(-mv.getX(),mv.getV1(),mv.getV2()); }

    // x * [vv]
    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,T,V1,V2> operator*(
        const RT x, const ElemProdVV<ix,T,V1,V2>& mv)
    { return ElemProdVV<0,T,V1,V2>(x*mv.getX(),mv.getV1(),mv.getV2()); }

    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,CT,V1,V2> operator*(
        const CT x, const ElemProdVV<ix,T,V1,V2>& mv)
    { return ElemProdVV<0,CT,V1,V2>(x*mv.getX(),mv.getV1(),mv.getV2()); }

    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,CT,V1,V2> operator*(
        const CCT x, const ElemProdVV<ix,T,V1,V2>& mv)
    { return ElemProdVV<0,CT,V1,V2>(x*mv.getX(),mv.getV1(),mv.getV2()); }

    template <int ix1, class T1, int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2> operator*(
        const Scaling<ix1,T1>& x, const ElemProdVV<ix,T,V1,V2>& mv)
    {
        return ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2>(
            T1(x)*mv.getX(),mv.getV1(),mv.getV2());
    }
    // [vv]*x
    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,T,V1,V2> operator*(
        const ElemProdVV<ix,T,V1,V2>& mv, const RT x)
    { return ElemProdVV<0,T,V1,V2>(x*mv.getX(),mv.getV1(),mv.getV2()); }

    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,CT,V1,V2> operator*(
        const ElemProdVV<ix,T,V1,V2>& mv, const CT x)
    { return ElemProdVV<0,CT,V1,V2>(x*mv.getX(),mv.getV1(),mv.getV2()); }

    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,CT,V1,V2> operator*(
        const ElemProdVV<ix,T,V1,V2>& mv, const CCT x)
    { return ElemProdVV<0,CT,V1,V2>(x*mv.getX(),mv.getV1(),mv.getV2()); }

    template <int ix1, class T1, int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2> operator*(
        const ElemProdVV<ix,T,V1,V2>& mv, const Scaling<ix1,T1>& x)
    {
        return ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2>(
            T1(x)*mv.getX(),mv.getV1(),mv.getV2());
    }

    // [vv]/x
    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,T,V1,V2> operator/(
        const ElemProdVV<ix,T,V1,V2>& mv, const RT x)
    { return ElemProdVV<0,T,V1,V2>(mv.getX()/x,mv.getV1(),mv.getV2()); }

    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,CT,V1,V2> operator/(
        const ElemProdVV<ix,T,V1,V2>& mv, const CT x)
    { return ElemProdVV<0,CT,V1,V2>(mv.getX()/x,mv.getV1(),mv.getV2()); }

    template <int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<0,CT,V1,V2> operator/(
        const ElemProdVV<ix,T,V1,V2>& mv, const CCT x)
    { return ElemProdVV<0,CT,V1,V2>(mv.getX()/x,mv.getV1(),mv.getV2()); }

    template <int ix1, class T1, int ix, class T, class V1, class V2>
    static TMV_INLINE ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2> operator/(
        const ElemProdVV<ix,T,V1,V2>& mv, const Scaling<ix1,T1>& x)
    {
        return ElemProdVV<ix1*ix,typename Traits2<T1,T>::type,V1,V2>(
            mv.getX()/T1(x),mv.getV1(),mv.getV2());
    }

#undef RT
#undef CT
#undef CCT

#ifdef TMV_TEXT
    template <int ix, class T, class V1, class V2>
    static inline std::string TMV_Text(const ElemProdVV<ix,T,V1,V2>& svv)
    {
        std::ostringstream s;
        s << "ElemProdVV< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(svv.getV1())<<" , ";
        s << TMV_Text(svv.getV2())<<" >";
        return s.str();
    }
#endif


} // namespace tmv

#endif 
