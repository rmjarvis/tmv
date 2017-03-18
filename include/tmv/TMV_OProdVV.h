

#ifndef TMV_OProdVV_H
#define TMV_OProdVV_H

#include "TMV_BaseVector.h"
#include "TMV_BaseMatrix.h"
#include "TMV_ProdXV.h"
#include "TMV_Rank1VVM_Funcs.h"

//#define XDEBUG_OPRODVV

#ifdef XDEBUG_OPRODVV
#include "TMV_Matrix.h"
#include <iostream>
#endif

namespace tmv {

    //
    // Vector ^ Vector
    //

    // These first few are for when an argument is a composite vector
    // and needs to be calculated before running Rank1Update
    template <bool add, int ix, class T, class V1, class V2, class M3>
    inline void Rank1Update(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1, 
        const BaseVector<V2>& v2, BaseMatrix_Mutable<M3>& m3)
    { Rank1Update<add>(x,v1.calc(),v2.calc(),m3.mat()); }

    // If everything is _Calc,  then there should be an overload
    // that says how to do the calculation.  This will give a 
    // compiler error on purpose.
    template <bool add, int ix, class T, class V1, class V2, class M3>
    inline void Rank1Update(
        const Scaling<ix,T>& , const BaseVector_Calc<V1>& , 
        const BaseVector_Calc<V2>& , BaseMatrix_Mutable<M3>& m3)
    { TMVStaticAssert(ix == 999); }

    // Also allow x to be missing (taken to be 1) or a scalar.
    template <bool add, class V1, class V2, class M3>
    inline void Rank1Update(
        const BaseVector<V1>& v1,
        const BaseVector<V2>& v2, BaseMatrix_Mutable<M3>& m3)
    {
        Rank1Update<add>(
            Scaling<1,typename M3::real_type>(),v1.vec(),v2.vec(),m3.mat()); 
    }
    template <bool add, class T, class V1, class V2, class M3>
    inline void Rank1Update(
        T x, const BaseVector<V1>& v1,
        const BaseVector<V2>& v2, BaseMatrix_Mutable<M3>& m3)
    { Rank1Update<add>(Scaling<0,T>(x),v1.vec(),v2.vec(),m3.mat()); }


#ifdef XDEBUG_OPRODVV
    template <bool add, int ix, class T, class V1, class V2, class M3>
    static void Rank1Update_Debug(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const BaseVector<V2>& v2, BaseMatrix_Mutable<M3>& m3)
    {
        //std::cout<<"Start Rank1Update XDEBUG"<<std::endl;
        //std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
        //std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1<<std::endl;
        //std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
        //std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3.mat()<<std::endl;
        Matrix<typename M3::value_type> m3i = m3.mat();
        Matrix<typename M3::value_type> m3c = m3.mat();
        if (!add) m3c.setZero();
        for(ptrdiff_t i=0;i<m3.colsize();++i) {
            for(ptrdiff_t j=0;j<m3.rowsize();++j) {
                m3c.ref(i,j) += T(x) * v1.cref(i) * v2.cref(j);
            }
        }
        //std::cout<<"m3c => "<<m3c<<std::endl;

        Rank1Update<add>(x,v1.vec(),v2.vec(),m3.mat());

        if (Norm(m3.mat()-m3c) > 1.e-6 * Norm(m3c)) {
            std::cout<<"Rank1Update:  add = "<<add<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3i<<std::endl;
            std::cout<<"m3 -> "<<m3.mat()<<std::endl;
            std::cout<<"correct = "<<m3c<<std::endl;
            exit(1);
        }
    }
#endif


    template <int ix, class T, class V1, class V2>
    class OProdVV;

    template <int ix, class T, class V1, class V2>
    struct Traits<OProdVV<ix,T,V1,V2> >
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        typedef typename Traits2<T1,T2>::type T12;
        typedef typename Traits2<T,T12>::type value_type;

        typedef ProdXV<ix,T,V1> const_col_type;
        typedef ProdXV<ix,T,V2> const_row_type;

        enum { _colsize = V1::_size };
        enum { _rowsize = V2::_size };
        enum { _nlo = IntTraits2<IntTraits<_colsize>::Sm1,0>::max };
        enum { _nhi = IntTraits2<IntTraits<_rowsize>::Sm1,0>::max };
        enum { _shape = Rec };
        enum { _fort = V1::_fort && V2::_fort };
        enum { _calc = false };
        enum { _rowmajor = false }; // arbitrary
        enum { _colmajor = true };

        typedef OProdVV<ix,T,V1,V2> type;
        enum { cs = _colsize };
        enum { rs = _rowsize };
        enum { A = ColMajor | (_fort ? FortranStyle : CStyle) };
        typedef typename MCopyHelper<value_type,Rec,cs,rs,A>:: type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (V1::_calc && V2::_calc),const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class V1, class V2>
    class OProdVV : 
        public BaseMatrix<OProdVV<ix,T,V1,V2> >
    {
    public:

        typedef OProdVV<ix,T,V1,V2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        OProdVV(const Scaling<ix,T>& _x, const BaseVector<V1>& _v1, 
                const BaseVector<V2>& _v2) :
            x(_x), v1(_v1.vec()), v2(_v2.vec()) {}

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const V1& getV1() const { return v1; }
        TMV_INLINE const V2& getV2() const { return v2; }

        TMV_INLINE ptrdiff_t colsize() const { return v1.size(); }
        TMV_INLINE ptrdiff_t rowsize() const { return v2.size(); }
        TMV_INLINE ptrdiff_t nlo() const { return TMV_MAX(colsize()-1,ptrdiff_t(0)); }
        TMV_INLINE ptrdiff_t nhi() const { return TMV_MAX(rowsize()-1,ptrdiff_t(0)); }

        value_type cref(ptrdiff_t i, ptrdiff_t j) const
        { return x * (v1.cref(i) * v2.cref(j)); }

        template <class M3>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M3::_shape>::assignable)); 
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_OPRODVV
            Rank1Update_Debug<false>(x,v1.vec(),v2.vec(),m3.mat());
#else
            Rank1Update<false>(x,v1.vec(),v2.vec(),m3.mat());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const V1& v1;
        const V2& v2;
    };


    // v ^ v
#define RT typename V1::real_type
    template <class V1, class V2>
    TMV_INLINE OProdVV<1,RT,V1,V2> operator^(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return OProdVV<1,RT,V1,V2>(RT(1),v1,v2); }
#undef RT

    // v ^ xv
    template <class V1, int ix, class T, class V2>
    TMV_INLINE OProdVV<ix,T,V1,V2> operator^(
        const BaseVector<V1>& v1, const ProdXV<ix,T,V2>& v2)
    { return OProdVV<ix,T,V1,V2>(v2.getX(),v1,v2.getV()); }

    // xv ^ v
    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<ix,T,V1,V2> operator^(
        const ProdXV<ix,T,V1>& v1, const BaseVector<V2>& v2)
    { return OProdVV<ix,T,V1,V2>(v1.getX(),v1.getV(),v2); }

    // xv ^ xv
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    TMV_INLINE OProdVV<ix1*ix2,PT,V1,V2> operator^(
        const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    {
        return OProdVV<ix1*ix2,PT,V1,V2>(
            v1.getX()*v2.getX(),v1.getV(),v2.getV()); 
    }
#undef PT


    // m += vv
    template <class M3, int ix, class T, class V1, class V2>
    inline void AddEq(
        BaseMatrix_Mutable<M3>& m, const OProdVV<ix,T,V1,V2>& vv)
    {
#ifdef XDEBUG_OPRODVV
        Rank1Update_Debug<true>(
            vv.getX(),vv.getV1().vec(),vv.getV2().vec(),m.mat()); 
#else
        Rank1Update<true>(
            vv.getX(),vv.getV1().vec(),vv.getV2().vec(),m.mat()); 
#endif
    }

    // m -= vv
    template <class M3, int ix, class T, class V1, class V2>
    inline void SubtractEq(
        BaseMatrix_Mutable<M3>& m, const OProdVV<ix,T,V1,V2>& vv)
    {
#ifdef XDEBUG_OPRODVV
        Rank1Update_Debug<true>(
            -vv.getX(),vv.getV1().vec(),vv.getV2().vec(),m.mat()); 
#else
        Rank1Update<true>(
            -vv.getX(),vv.getV1().vec(),vv.getV2().vec(),m.mat()); 
#endif
    }


    // Consolidate x*(xmv) type constructs:

#define RT typename OProdVV<ix,T,V1,V2>::real_type
#define CT typename OProdVV<ix,T,V1,V2>::complex_type
#define CCT ConjRef<CT>

    // -(x*v)
    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<-ix,T,V1,V2> operator-(const OProdVV<ix,T,V1,V2>& vv)
    { return OProdVV<-ix,T,V1,V2>(-vv.getX(),vv.getV1(),vv.getV2()); }

    // x * (x*v)
    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,T,V1,V2> operator*(
        const RT x, const OProdVV<ix,T,V1,V2>& vv)
    { return OProdVV<0,T,V1,V2>(x*vv.getX(),vv.getV1(),vv.getV2()); }

    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,CT,V1,V2> operator*(
        const CT x, const OProdVV<ix,T,V1,V2>& vv)
    { return OProdVV<0,CT,V1,V2>(x*vv.getX(),vv.getV1(),vv.getV2()); }

    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,CT,V1,V2> operator*(
        const CCT x, const OProdVV<ix,T,V1,V2>& vv)
    { return OProdVV<0,CT,V1,V2>(x*vv.getX(),vv.getV1(),vv.getV2()); }

    // (x*v)*x
    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,T,V1,V2> operator*(
        const OProdVV<ix,T,V1,V2>& vv, const RT x)
    { return OProdVV<0,T,V1,V2>(x*vv.getX(),vv.getV1(),vv.getV2()); }

    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,CT,V1,V2> operator*(
        const OProdVV<ix,T,V1,V2>& vv, const CT x)
    { return OProdVV<0,CT,V1,V2>(x*vv.getX(),vv.getV1(),vv.getV2()); }

    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,CT,V1,V2> operator*(
        const OProdVV<ix,T,V1,V2>& vv, const CCT x)
    { return OProdVV<0,CT,V1,V2>(x*vv.getX(),vv.getV1(),vv.getV2()); }

    // (x*v)/x
    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,T,V1,V2> operator/(
        const OProdVV<ix,T,V1,V2>& vv, const RT x)
    {
        return OProdVV<0,T,V1,V2>(
            ZProd<false,false>::quot(vv.getX(),x),vv.getV1(),vv.getV2()); }

    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,CT,V1,V2> operator/(
        const OProdVV<ix,T,V1,V2>& vv, const CT x)
    { 
        return OProdVV<0,CT,V1,V2>(
            ZProd<false,false>::quot(vv.getX(),x),vv.getV1(),vv.getV2()); }

    template <int ix, class T, class V1, class V2>
    TMV_INLINE OProdVV<0,CT,V1,V2> operator/(
        const OProdVV<ix,T,V1,V2>& vv, const CCT x)
    { 
        return OProdVV<0,CT,V1,V2>(
            ZProd<false,false>::quot(vv.getX(),CT(x)),vv.getV1(),vv.getV2()); }

#undef RT
#undef CT
#undef CCT

    // TMV_Text

    template <int ix, class T, class V1, class V2>
    inline std::string TMV_Text(const OProdVV<ix,T,V1,V2>& svv)
    {
        std::ostringstream s;
        s << "OProdVV< "<<ix<<","<<TMV_Text(T())<<",";
        s << TMV_Text(svv.getV1())<<" , "<<TMV_Text(svv.getV2())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
