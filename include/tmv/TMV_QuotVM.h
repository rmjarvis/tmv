

#ifndef TMV_QuotVM_H
#define TMV_QuotVM_H

#include "TMV_BaseVector.h"
#include "TMV_BaseMatrix.h"
#include "TMV_ProdXM.h"
#include "TMV_QuotXM.h"
#include "TMV_ProdXV.h"
#include "TMV_DivVM_Funcs.h"

//#define XDEBUG_QUOTVM

#ifdef XDEBUG_QUOTVM
#include "TMV_ProdMV.h"
#include "TMV_Vector.h"
#include <iostream>
#endif

namespace tmv {

    template <int ix, class T, class V1, class M2, class V3>
    inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1.calc(),m2.calc(),v3.vec()); }
    template <int ix, class T, class V1, class M2, class V3>
    inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(x,v1.calc(),m2.calc(),v3.vec()); }

    // Also allow x to be missing (taken to be 1) or a scalar.
    template <class V1, class M2, class V3>
    inline void LDiv(
        const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(Scaling<1,typename V3::real_type>(),v1.calc(),m2.calc(),v3.vec()); }
    template <class T, class V1, class M2, class V3>
    inline void LDiv(
        T x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(Scaling<0,T>(x),v1.calc(),m2.calc(),v3.vec()); }

    template <class V1, class M2, class V3>
    inline void RDiv(
        const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(Scaling<1,typename V3::real_type>(),v1.calc(),m2.calc(),v3.vec()); }
    template <class T, class V1, class M2, class V3>
    inline void RDiv(
        T x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    { RDiv(Scaling<0,T>(x),v1.calc(),m2.calc(),v3.vec()); }
    

    //
    // Vector (/%) Matrix
    //

#ifdef XDEBUG_QUOTVM
    template <int ix, class T, class V1, class M2, class V3>
    static void LDiv_Debug(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        // v3 = x*m2^-1*v1
        // --> m2*v3/x = v1
        Vector<typename V1::value_type> v1i = v1;
        Matrix<typename M2::value_type> m2i = m2;
        Vector<typename V3::value_type> v3i = v3;

        LDiv(x,v1.vec(),m2.mat(),v3.vec());

        Vector<typename V3::value_type> v1c = m2i*v3/x;
        const typename V3::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(v1i-v1c) > 1.e-6 * kappa * Norm(v1i)) {
            std::cout<<"LDivVM:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2i<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3i<<std::endl;
            std::cout<<"v3 -> "<<v3<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1i-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1i-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTVM
    template <int ix, class T, class V1, class M2, class V3>
    static void RDiv_Debug(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        // v3 = x*v1*m2^-1
        // --> v3*m2/x = v1
        Vector<typename V1::value_type> v1i = v1;
        Matrix<typename M2::value_type> m2i = m2;
        Vector<typename V3::value_type> v3i = v3;

        RDiv(x,v1.vec(),m2.mat(),v3.vec());

        Vector<typename V3::value_type> v1c = v3*m2i/x;
        const typename V3::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(v1i-v1c) > 1.e-6 * kappa * Norm(v1i)) {
            std::cout<<"RDivVM:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2i<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3i<<std::endl;
            std::cout<<"v3 -> "<<v3<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1i-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1i-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTVM
    template <class V1, class M2>
    static void LDivEq_Debug(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    {
        Vector<typename V1::value_type> v1i = v1;
        Matrix<typename M2::value_type> m2i = m2;

        LDivEq(v1.vec(),m2.mat());

        Vector<typename V1::value_type> v1c = m2i*v1;
        const typename V1::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(v1i-v1c) > 1.e-6 * kappa * Norm(v1i)) {
            std::cout<<"LDivEqVM:  \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2i<<std::endl;
            std::cout<<"v1 -> "<<v1<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1i-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1i-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTVM
    template <class V1, class M2>
    static void RDivEq_Debug(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    {
        Vector<typename V1::value_type> v1i = v1;
        Matrix<typename M2::value_type> m2i = m2;

        RDivEq(v1.vec(),m2.mat());

        Vector<typename V1::value_type> v1c = v1*m2i;
        const typename V1::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(v1i-v1c) > 1.e-6 * kappa * Norm(v1i)) {
            std::cout<<"RDivEqVM:  \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2i<<std::endl;
            std::cout<<"v1 -> "<<v1<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1i-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1i-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

    template <int ix, class T, class V1, class M2>
    class QuotVM;

    template <int ix, class T, class V1, class M2>
    struct Traits<QuotVM<ix,T,V1,M2> >
    {
        typedef typename ProdXV<ix,T,V1>::value_type vtype1;
        typedef typename M2::value_type mtype2;
        typedef typename Traits2<vtype1,mtype2>::type value_type;

        enum { _size = M2::_rowsize };
        enum { _fort = M2::_fort && V1::_fort };
        enum { _calc = false };

        typedef QuotVM<ix,T,V1,M2> type;
        enum { A = _fort ? FortranStyle : CStyle };
        typedef typename VCopyHelper<value_type,_size,A>::type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
    };

    template <int ix, class T, class V1, class M2>
    class QuotVM : 
        public BaseVector<QuotVM<ix,T,V1,M2> >
    {
    public:

        typedef QuotVM<ix,T,V1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        QuotVM(const Scaling<ix,T>& _x, const BaseVector<V1>& _v1,
               const BaseMatrix<M2>& _m2) : 
            x(_x), v1(_v1.vec()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M2::_colsize,V1::_size>::same)); 
            TMVAssert(m2.colsize() == v1.size());
        }

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const V1& getV() const { return v1; }
        TMV_INLINE const M2& getM() const { return m2; }

        TMV_INLINE int size() const { return m2.rowsize(); }

        template <class V3>
        TMV_INLINE_ND void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::isreal || V3::iscomplex));
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_QUOTVM
            LDiv_Debug(x,v1.vec(),m2.mat(),v3.vec());
#else
            LDiv(x,v1.vec(),m2.mat(),v3.vec());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const V1& v1;
        const M2& m2;
    };

    template <int ix, class T, class V1, class M2>
    class RQuotVM;

    template <int ix, class T, class V1, class M2>
    struct Traits<RQuotVM<ix,T,V1,M2> >
    {
        typedef typename ProdXV<ix,T,V1>::value_type vtype1;
        typedef typename M2::value_type mtype2;
        typedef typename Traits2<vtype1,mtype2>::type value_type;

        enum { _size = M2::_colsize };
        enum { _fort = M2::_fort && V1::_fort };
        enum { _calc = false };

        typedef RQuotVM<ix,T,V1,M2> type;
        enum { A = _fort ? FortranStyle : CStyle };
        typedef typename VCopyHelper<value_type,_size,A>::type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
    };

    template <int ix, class T, class V1, class M2>
    class RQuotVM : 
        public BaseVector<RQuotVM<ix,T,V1,M2> >
    {
    public:

        typedef RQuotVM<ix,T,V1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        RQuotVM(const Scaling<ix,T>& _x, const BaseVector<V1>& _v1,
                const BaseMatrix<M2>& _m2) : 
            x(_x), v1(_v1.vec()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M2::_rowsize,V1::_size>::same)); 
            TMVAssert(m2.rowsize() == v1.size());
        }

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const V1& getV() const { return v1; }
        TMV_INLINE const M2& getM() const { return m2; }

        TMV_INLINE int size() const { return m2.colsize(); }

        template <class V3>
        TMV_INLINE_ND void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::isreal || V3::iscomplex));
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_QUOTVM
            RDiv_Debug(x,v1.vec(),m2.mat(),v3.vec());
#else
            RDiv(x,v1.vec(),m2.mat(),v3.vec());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const V1& v1;
        const M2& m2;
    };



    // v / m
#define RT typename V::real_type
    template <class M, class V>
    TMV_INLINE QuotVM<1,RT,V,M> operator/(
        const BaseVector<V>& v, const BaseMatrix<M>& m)
    { return QuotVM<1,RT,V,M>(RT(1),v,m); }
#undef RT

    // v / xm
    template <class M, int ix, class T, class V>
    TMV_INLINE QuotVM<ix,T,V,M> operator/(
        const BaseVector<V>& v, const ProdXM<ix,T,M>& m)
    {
        return QuotVM<ix,T,V,M>(
            typename Traits<T>::real_type(1)/m.getX(),v,m.getM()); 
    }

    // xv / m
    template <int ix, class T, class M, class V>
    TMV_INLINE QuotVM<ix,T,V,M> operator/(
        const ProdXV<ix,T,V>& v, const BaseMatrix<M>& m)
    { return QuotVM<ix,T,V,M>(v.getX(),v.getV(),m); }

    // xv / xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M, int ix2, class T2, class V>
    TMV_INLINE QuotVM<ix1*ix2,PT,V,M> operator/(
        const ProdXV<ix1,T1,V>& v, const ProdXM<ix2,T2,M>& m)
    {
        return QuotVM<ix1*ix2,PT,V,M>(
            ZProd<false,false>::quot(v.getX(),m.getX()),v.getV(),m.getM()); 
    }
#undef PT

    // x/m * v
    template <int ix, class T, class M, class V>
    TMV_INLINE QuotVM<ix,T,V,M> operator*(
        const QuotXM<ix,T,M>& m, const BaseVector<V>& v)
    { return QuotVM<ix,T,V,M>(m.getX(),v,m.getM()); }

    // x/m * xv
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M, int ix2, class T2, class V>
    TMV_INLINE QuotVM<ix1*ix2,PT,V,M> operator*(
        const QuotXM<ix1,T1,M>& m, const ProdXV<ix2,T2,V>& v)
    { return QuotVM<ix1*ix2,PT,V,M>(m.getX()*v.getX(),v.getV(),m.getM()); }
#undef PT


    // v /= m
    template <class V1, class M2>
    inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    {
#ifdef XDEBUG_QUOTVM
        LDivEq_Debug(v1.vec(),m2.mat()); 
#else
        LDivEq(v1.vec(),m2.calc()); 
#endif
    }

    // v /= xm
    template <class V1, int ix2, class T2, class M2>
    inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const ProdXM<ix2,T2,M2>& m2)
    {
        typedef typename Traits<T2>::real_type RT;
#ifdef XDEBUG_QUOTVM
        LDivEq_Debug(v1.vec(),m2.getM().mat()); 
#else
        LDivEq(v1.vec(),m2.getM().mat()); 
#endif
        Scale(ZProd<false,false>::quot(RT(1),m2.getX()),v1.vec());
    }


    // v % m
#define RT typename V::real_type
    template <class M, class V>
    TMV_INLINE RQuotVM<1,RT,V,M> operator%(
        const BaseVector<V>& v, const BaseMatrix<M>& m)
    { return RQuotVM<1,RT,V,M>(RT(1),v,m); }
#undef RT

    // v % xm
    template <class V, class M, int ix, class T>
    TMV_INLINE RQuotVM<ix,T,V,M> operator%(
        const BaseVector<V>& v, const ProdXM<ix,T,M>& m)
    {
        typedef typename Traits<T>::real_type RT;
        return RQuotVM<ix,T,V,M>(
            ZProd<false,false>::quot(RT(1),m.getX()),v,m.getM()); 
    }

    // xv % m
    template <int ix, class T, class V, class M>
    TMV_INLINE RQuotVM<ix,T,V,M> operator%(
        const ProdXV<ix,T,V>& v, const BaseMatrix<M>& m)
    { return RQuotVM<ix,T,V,M>(v.getX(),v.getV(),m); }

    // xv % xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class V, int ix2, class T2, class M>
    TMV_INLINE RQuotVM<ix1*ix2,PT,V,M> operator%(
        const ProdXV<ix1,T1,V>& v, const ProdXM<ix2,T2,M>& m)
    { 
        return RQuotVM<ix1*ix2,PT,V,M>(
            ZProd<false,false>::quot(v.getX(),m.getX()),v.getV(),m.getM()); 
    }
#undef PT

    // v * x/m
    template <class V, int ix, class T, class M>
    TMV_INLINE RQuotVM<ix,T,V,M> operator*(
        const BaseVector<V>& v, const QuotXM<ix,T,M>& m)
    { return RQuotVM<ix,T,V,M>(m.getX(),v,m.getM()); }

    // xv * x/m
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class V, int ix2, class T2, class M>
    TMV_INLINE RQuotVM<ix1*ix2,PT,V,M> operator*(
        const ProdXV<ix1,T1,V>& v, const QuotXM<ix2,T2,M>& m)
    { return RQuotVM<ix1*ix2,PT,V,M>(m.getX()*v.getX(),v.getV(),m.getM()); }
#undef PT

    // v %= m
    template <class V1, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    {
#ifdef XDEBUG_QUOTVM
        RDivEq_Debug(v1.vec(),m2.mat()); 
#else
        RDivEq(v1.vec(),m2.calc()); 
#endif
    }

    // v %= xm
    template <class V1, int ix2, class T2, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const ProdXM<ix2,T2,M2>& m2)
    {
        typedef typename Traits<T2>::real_type RT;
#ifdef XDEBUG_QUOTVM
        RDivEq_Debug(v1.vec(),m2.getM().mat()); 
#else
        RDivEq(v1.vec(),m2.getM().calc()); 
#endif
        Scale(ZProd<false,false>::quot(RT(1),m2.getX()),v1.vec());
    }

    // v *= x/m
    template <class V1, int ix2, class T2, class M2>
    inline void MultEq(
        BaseVector_Mutable<V1>& v1, const QuotXM<ix2,T2,M2>& m2)
    {
#ifdef XDEBUG_QUOTVM
        RDivEq_Debug(v1.vec(),m2.getM().mat()); 
#else
        RDivEq(v1.vec(),m2.getM().calc()); 
#endif
        Scale(m2.getX(),v1.vec());
    }


    // Consolidate x*(xv/m) type constructs:

#define RT typename QuotVM<ix,T,V1,M2>::real_type
#define CT typename QuotVM<ix,T,V1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(v/m)
    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<-ix,T,V1,M2> operator-(
        const QuotVM<ix,T,V1,M2>& qvm)
    { return QuotVM<-ix,T,V1,M2>(-qvm.getX(),qvm.getV(),qvm.getM()); }

    // x * (v/m)
    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,T,V1,M2> operator*(
        const int x, const QuotVM<ix,T,V1,M2>& qvm)
    { return QuotVM<0,T,V1,M2>(RT(x)*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,T,V1,M2> operator*(
        const RT x, const QuotVM<ix,T,V1,M2>& qvm)
    { return QuotVM<0,T,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,CT,V1,M2> operator*(
        const CT x, const QuotVM<ix,T,V1,M2>& qvm)
    { return QuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,CT,V1,M2> operator*(
        const CCT x, const QuotVM<ix,T,V1,M2>& qvm)
    { return QuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    // (v/m)*x
    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,T,V1,M2> operator*(
        const QuotVM<ix,T,V1,M2>& qvm, const int x)
    { return QuotVM<0,T,V1,M2>(RT(x)*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,T,V1,M2> operator*(
        const QuotVM<ix,T,V1,M2>& qvm, const RT x)
    { return QuotVM<0,T,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,CT,V1,M2> operator*(
        const QuotVM<ix,T,V1,M2>& qvm, const CT x)
    { return QuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,CT,V1,M2> operator*(
        const QuotVM<ix,T,V1,M2>& qvm, const CCT x)
    { return QuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    // (v/m)/x
    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,T,V1,M2> operator/(
        const QuotVM<ix,T,V1,M2>& qvm, const int x)
    {
        return QuotVM<0,T,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),RT(x)),qvm.getV(),qvm.getM()); 
    }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,T,V1,M2> operator/(
        const QuotVM<ix,T,V1,M2>& qvm, const RT x)
    {
        return QuotVM<0,T,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),x),qvm.getV(),qvm.getM()); 
    }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,CT,V1,M2> operator/(
        const QuotVM<ix,T,V1,M2>& qvm, const CT x)
    {
        return QuotVM<0,CT,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),x),qvm.getV(),qvm.getM()); 
    }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE QuotVM<0,CT,V1,M2> operator/(
        const QuotVM<ix,T,V1,M2>& qvm, const CCT x)
    {
        return QuotVM<0,CT,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),CT(x)),qvm.getV(),qvm.getM()); 
    }

#undef RT
#undef CT
#undef CCT

#define RT typename RQuotVM<ix,T,V1,M2>::real_type
#define CT typename RQuotVM<ix,T,V1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(v/m)
    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<-ix,T,V1,M2> operator-(
        const RQuotVM<ix,T,V1,M2>& qvm)
    { return RQuotVM<-ix,T,V1,M2>(-qvm.getX(),qvm.getV(),qvm.getM()); }

    // x * (v/m)
    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,T,V1,M2> operator*(
        const int x, const RQuotVM<ix,T,V1,M2>& qvm)
    { return RQuotVM<0,T,V1,M2>(RT(x)*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,T,V1,M2> operator*(
        const RT x, const RQuotVM<ix,T,V1,M2>& qvm)
    { return RQuotVM<0,T,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,CT,V1,M2> operator*(
        const CT x, const RQuotVM<ix,T,V1,M2>& qvm)
    { return RQuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,CT,V1,M2> operator*(
        const CCT x, const RQuotVM<ix,T,V1,M2>& qvm)
    { return RQuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    // (v/m)*x
    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,T,V1,M2> operator*(
        const RQuotVM<ix,T,V1,M2>& qvm, const int x)
    { return RQuotVM<0,T,V1,M2>(RT(x)*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,T,V1,M2> operator*(
        const RQuotVM<ix,T,V1,M2>& qvm, const RT x)
    { return RQuotVM<0,T,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,CT,V1,M2> operator*(
        const RQuotVM<ix,T,V1,M2>& qvm, const CT x)
    { return RQuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,CT,V1,M2> operator*(
        const RQuotVM<ix,T,V1,M2>& qvm, const CCT x)
    { return RQuotVM<0,CT,V1,M2>(x*qvm.getX(),qvm.getV(),qvm.getM()); }

    // (v/m)/x
    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,T,V1,M2> operator/(
        const RQuotVM<ix,T,V1,M2>& qvm, const int x)
    {
        return RQuotVM<0,T,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),RT(x)),qvm.getV(),qvm.getM()); 
    }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,T,V1,M2> operator/(
        const RQuotVM<ix,T,V1,M2>& qvm, const RT x)
    {
        return RQuotVM<0,T,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),x),qvm.getV(),qvm.getM()); 
    }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,CT,V1,M2> operator/(
        const RQuotVM<ix,T,V1,M2>& qvm, const CT x)
    {
        return RQuotVM<0,CT,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),x),qvm.getV(),qvm.getM()); 
    }

    template <int ix, class T, class V1, class M2>
    TMV_INLINE RQuotVM<0,CT,V1,M2> operator/(
        const RQuotVM<ix,T,V1,M2>& qvm, const CCT x)
    {
        return RQuotVM<0,CT,V1,M2>(
            ZProd<false,false>::quot(qvm.getX(),x),qvm.getV(),qvm.getM()); 
    }

#undef RT
#undef CT
#undef CCT


    // TMV_Text

    template <int ix, class T, class V1, class M2>
    inline std::string TMV_Text(const QuotVM<ix,T,V1,M2>& qvm)
    {
        std::ostringstream s;
        s << "QuotVM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(qvm.getV())<<" , "<<TMV_Text(qvm.getM())<<" >";
        return s.str();
    }

    template <int ix, class T, class V1, class M2>
    inline std::string TMV_Text(const RQuotVM<ix,T,V1,M2>& qvm)
    {
        std::ostringstream s;
        s << "RQuotVM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(qvm.getV())<<" , "<<TMV_Text(qvm.getM())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
