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


#ifndef TMV_QuotVM_H
#define TMV_QuotVM_H

#include "TMV_ProdXM.h"
#include "TMV_QuotXM.h"
#include "TMV_ProdXV.h"

//#define XDEBUG_QUOTVM

#ifdef XDEBUG_QUOTVM
#include "TMV_ProdMV.h";
#include "TMV_Vector.h"
#include <iostream>
#endif

namespace tmv {

    //
    // Vector (/%) Matrix
    //

    // These are intentionally not defined to make sure we
    // get a compiler error if they are used.
    // All real calls should go through a more specific version than 
    // just the BaseMatrix_Calc's.
    template <class V1, class M2, class V3>
    inline void LDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, class M2, class V3>
    inline void NoAliasLDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, class M2, class V3>
    inline void InlineLDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, class M2, class V3>
    inline void AliasLDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <class V1, class M2, class V3>
    inline void RDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, class M2, class V3>
    inline void NoAliasRDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, class M2, class V3>
    inline void InlineRDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <class V1, class M2, class V3>
    inline void AliasRDivVM(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3);


#ifdef XDEBUG_QUOTVM
    template <class V1, class M2, class V3>
    inline void LDivVM_Debug(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        // v3 = m2^-1*v1
        // --> m2*v3 = v1
        Vector<typename V3::value_type> v3i = v3;

        LDivVM(v1.vec(),m2.mat(),v3.vec());

        Vector<typename V3::value_type> v1c = m2*v3;
        const typename V3::real_type kappa = Norm(m2.inverse())*Norm(m2);

        if (Norm(v1-v1c) > 1.e-6 * kappa * Norm(v1)) {
            std::cout<<"LDivVM:  \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3i<<std::endl;
            std::cout<<"v3 -> "<<v3<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTVM
    template <class V1, class M2, class V3>
    inline void RDivVM_Debug(
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Calc<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        // v3 = v1*m2^-1
        // --> v3*m2 = v1
        Vector<typename V3::value_type> v3i = v3;

        RDivVM(x,v1.vec(),m2.mat(),v3.vec());

        Vector<typename V3::value_type> v1c = v3*m2/x;
        const typename V3::real_type kappa = Norm(m2.inverse())*Norm(m2);

        if (Norm(v1-v1c) > 1.e-6 * kappa * Norm(v1)) {
            std::cout<<"RDivVM:  \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3i<<std::endl;
            std::cout<<"v3 -> "<<v3<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTVM
    template <class V1, class M2>
    inline void LDivEqVM_Debug(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        Vector<typename V1::value_type> v1i = v1;

        LDivEqVM(v1.vec(),m2.mat());

        Vector<typename V1::value_type> v1c = m2*v1i;
        const typename V1::real_type kappa = Norm(m2.inverse())*Norm(m2);

        if (Norm(v1-v1c) > 1.e-6 * kappa * Norm(v1i)) {
            std::cout<<"LDivEqVM:  \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"v1 -> "<<v1<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1-v1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTVM
    template <class V1, class M2>
    inline void RDivEqVM_Debug(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        Vector<typename V1::value_type> v1i = v1;

        RDivEqVM(v1.vec(),m2.mat());

        Vector<typename V1::value_type> v1c = v1i*m2;
        const typename V1::real_type kappa = Norm(m2.inverse())*Norm(m2);

        if (Norm(v1-v1c) > 1.e-6 * kappa * Norm(v1i)) {
            std::cout<<"RDivEqVM:  \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"v1 -> "<<v1<<std::endl;
            std::cout<<"v1c = "<<v1c<<std::endl;
            std::cout<<"diff = "<<v1-v1c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1-v1c)<<std::endl;
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

        enum { vsize = M2::mrowsize };
        enum { vfort = M2::mfort && V1::vfort };
        enum { vcalc = false };

        typedef QuotVM<ix,T,V1,M2> type;
        typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
    };

    template <int ix, class T, class V1, class M2>
    class QuotVM : public BaseVector<QuotVM<ix,T,V1,M2> >
    {
    public:

        typedef QuotVM<ix,T,V1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        inline QuotVM(
            const T _x, const BaseVector<V1>& _v1, const BaseMatrix<M2>& _m2
        ) : 
            x(_x), v1(_v1.vec()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M2::mcolsize,V1::vsize>::same)); 
            TMVAssert(m2.colsize() == v1.size());
        }

        inline const Scaling<ix,T>& getX() const { return x; }
        inline const V1& getV() const { return v1; }
        inline const M2& getM() const { return m2; }

        inline size_t size() const { return m2.rowsize(); }

        inline value_type cref(int i) const
        { return this->calc().cref(i); }

        template <class V3>
        inline void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::visreal || V3::viscomplex));
            TMVStaticAssert((Sizes<type::vsize,V3::vsize>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_QUOTVM
            LDivVM_Debug(v1.calc(),m2.calc(),v3.vec());
#else
            LDivVM(v1.calc(),m2.calc(),v3.vec());
#endif
            Scale(x,v3.vec());
        }

        template <class V3>
        inline void newAssignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::visreal || V3::viscomplex));
            TMVStaticAssert((Sizes<type::vsize,V3::vsize>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_QUOTVM
            LDivVM_Debug(v1.calc(),m2.calc(),v3.vec());
#else
            NoAliasLDivVM(v1.calc(),m2.calc(),v3.vec());
#endif
            Scale(x,v3.vec());
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

        enum { vsize = M2::mrowsize };
        enum { vfort = M2::mfort && V1::vfort };
        enum { vcalc = false };

        typedef RQuotVM<ix,T,V1,M2> type;
        typedef typename VCopyHelper<value_type,vsize,vfort>::type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
    };

    template <int ix, class T, class V1, class M2>
    class RQuotVM : public BaseVector<RQuotVM<ix,T,V1,M2> >
    {
    public:

        typedef RQuotVM<ix,T,V1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        inline RQuotVM(
            const T _x, const BaseVector<V1>& _v1, const BaseMatrix<M2>& _m2
        ) : 
            x(_x), v1(_v1.vec()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M2::mcolsize,V1::vsize>::same)); 
            TMVAssert(m2.colsize() == v1.size());
        }

        inline const Scaling<ix,T>& getX() const { return x; }
        inline const V1& getV() const { return v1; }
        inline const M2& getM() const { return m2; }

        inline size_t size() const { return m2.rowsize(); }

        inline value_type cref(int i) const
        { return this->calc().cref(i); }

        template <class V3>
        inline void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::visreal || V3::viscomplex));
            TMVStaticAssert((Sizes<type::vsize,V3::vsize>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_QUOTVM
            LDivVM_Debug(v1.calc(),m2.calc(),v3.vec());
#else
            LDivVM(v1.calc(),m2.calc(),v3.vec());
#endif
            Scale(x,v3.vec());
        }

        template <class V3>
        inline void newAssignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::visreal || V3::viscomplex));
            TMVStaticAssert((Sizes<type::vsize,V3::vsize>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_QUOTVM
            LDivVM_Debug(v1.calc(),m2.calc(),v3.vec());
#else
            NoAliasLDivVM(v1.calc(),m2.calc(),v3.vec());
#endif
            Scale(x,v3.vec());
        }

    private:
        const Scaling<ix,T> x;
        const V1& v1;
        const M2& m2;
    };



    // v / m
#define RT typename V::real_type
    template <class M, class V>
    inline QuotVM<1,RT,V,M> operator/(
        const BaseVector<V>& v, const BaseMatrix<M>& m)
    { return QuotVM<1,RT,V,M>(RT(1),v,m); }
#undef RT

    // v / xm
    template <class M, int ix, class T, class V>
    inline QuotVM<ix,T,V,M> operator/(
        const BaseVector<V>& v, const ProdXM<ix,T,M>& m)
    { return QuotVM<ix,T,V,M>(Traits<T>::real_type(1)/m.getX(),v,m.getM()); }

    // xv / m
    template <int ix, class T, class M, class V>
    inline QuotVM<ix,T,V,M> operator/(
        const ProdXV<ix,T,V>& v, const BaseMatrix<M>& m)
    { return QuotVM<ix,T,V,M>(v.getX(),v.getV(),m); }

    // xv / xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M, int ix2, class T2, class V>
    inline QuotVM<ix1*ix2,PT,V,M> operator/(
        const ProdXV<ix1,T1,V>& v, const ProdXM<ix2,T2,M>& m)
    { return QuotVM<ix1*ix2,PT,V,M>(m.getX()/v.getX(),v.getV(),m.getM()); }
#undef PT

    // x/m * v
    template <int ix, class T, class M, class V>
    inline QuotVM<ix,T,V,M> operator*(
        const QuotXM<ix,T,M>& m, const BaseVector<V>& v)
    { return QuotVM<ix,T,V,M>(m.getX(),v,m.getM()); }

    // x/m * xv
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M, int ix2, class T2, class V>
    inline QuotVM<ix1*ix2,PT,V,M> operator*(
        const QuotXM<ix1,T1,M>& m, const ProdXV<ix2,T2,V>& v)
    { return QuotVM<ix1*ix2,PT,V,M>(m.getX()*v.getX(),v.getV(),m.getM()); }
#undef PT


    // v /= m
    template <class V1, class M2>
    inline void DivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    { 
#ifdef XDEBUG_QUOTVM
        LDivEqVM_Debug(v1,m2.calc()); 
#else
        LDivEqVM(v1,m2.calc()); 
#endif
    }

    // v /= xm
    template <class V1, int ix2, class T2, class M2>
    inline void DivEq(
        BaseVector_Mutable<V1>& v1, const ProdXM<ix2,T2,M2>& m2)
    { 
#ifdef XDEBUG_QUOTVM
        LDivEqVM_Debug(v1,m2.getM().calc()); 
#else
        LDivEqVM(v1,m2.getM().calc()); 
#endif
        Scale(typename Traits<T2>::real_type(1)/m2.getX(),v1.vec());
    }


    // v % m
#define RT typename V::real_type
    template <class M, class V>
    inline RQuotVM<1,RT,V,M> operator%(
        const BaseVector<V>& v, const BaseMatrix<M>& m)
    { return RQuotVM<1,RT,V,M>(RT(1),v,m); }
#undef RT

    // v % xm
    template <class V, class M, int ix, class T>
    inline RQuotVM<ix,T,V,M> operator%(
        const BaseVector<V>& v, const ProdXM<ix,T,M>& m)
    { return RQuotVM<ix,T,V,M>(Traits<T>::real_type(1)/m.getX(),v,m.getM()); }

    // xv % m
    template <int ix, class T, class V, class M>
    inline RQuotVM<ix,T,V,M> operator%(
        const ProdXV<ix,T,V>& v, const BaseMatrix<M>& m)
    { return RQuotVM<ix,T,V,M>(v.getX(),v.getV(),m); }

    // xv % xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class V, int ix2, class T2, class M>
    inline RQuotVM<ix1*ix2,PT,V,M> operator%(
        const ProdXV<ix1,T1,V>& v, const ProdXM<ix2,T2,M>& m)
    { return RQuotVM<ix1*ix2,PT,V,M>(m.getX()/v.getX(),v.getV(),m.getM()); }
#undef PT

    // v * x/m
    template <class V, int ix, class T, class M>
    inline RQuotVM<ix,T,V,M> operator*(
        const BaseVector<V>& v, const QuotXM<ix,T,M>& m)
    { return RQuotVM<ix,T,V,M>(m.getX(),v,m.getM()); }

    // xv * x/m
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class V, int ix2, class T2, class M>
    inline RQuotVM<ix1*ix2,PT,V,M> operator*(
        const ProdXV<ix1,T1,V>& v, const QuotXM<ix2,T2,M>& m)
    { return RQuotVM<ix1*ix2,PT,V,M>(m.getX()*v.getX(),v.getV(),m.getM()); }
#undef PT

    // v %= m
    template <class V1, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    { 
#ifdef XDEBUG_QUOTVM
        RDivEqVM_Debug(v1,Scaling<1,typename V1::real_type>(),m2.calc()); 
#else
        RDivEqVM(v1,Scaling<1,typename V1::real_type>(),m2.calc()); 
#endif
    }

    // v %= xm
    template <class V1, int ix2, class T2, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const ProdXM<ix2,T2,M2>& m2)
    { 
#ifdef XDEBUG_QUOTVM
        RDivEqVM_Debug(v1,m2.getM().calc()); 
#else
        RDivEqVM(v1,m2.getM().calc()); 
#endif
        Scale(typename Traits<T2>::real_type(1)/m2.getX(),v1.vec());
    }

    // v *= x/m
    template <class V1, int ix2, class T2, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const QuotXM<ix2,T2,M2>& m2)
    { 
#ifdef XDEBUG_QUOTVM
        RDivEqVM_Debug(v1,m2.getX(),m2.getM().calc()); 
#else
        RDivEqVM(v1,m2.getX(),m2.getM().calc()); 
#endif
    }


    // Consolidate x*(xv/m) type constructs:

#define RT typename QuotVM<ix,T,M1,V2>::real_type
#define CT typename QuotVM<ix,T,M1,V2>::complex_type
#define CCT ConjRef<CT>

    // -(v/m)
    template <int ix, class T, class M1, class V2>
    inline QuotVM<-ix,T,M1,V2> operator-(const QuotVM<ix,T,M1,V2>& mv)
    { return QuotVM<-ix,T,M1,V2>(-mv.getX(),mv.getM(),mv.getV()); }

    // x * (v/m)
    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,T,M1,V2> operator*(
        const RT x, const QuotVM<ix,T,M1,V2>& mv)
    { return QuotVM<0,T,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,CT,M1,V2> operator*(
        const CT x, const QuotVM<ix,T,M1,V2>& mv)
    { return QuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,CT,M1,V2> operator*(
        const CCT x, const QuotVM<ix,T,M1,V2>& mv)
    { return QuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    inline QuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
        const Scaling<ix1,T1>& x, const QuotVM<ix,T,M1,V2>& mv)
    { 
        return QuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            T1(x)*mv.getX(),mv.getM(),mv.getV()); 
    }

    // (v/m)*x
    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,T,M1,V2> operator*(
        const QuotVM<ix,T,M1,V2>& mv, const RT x)
    { return QuotVM<0,T,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,CT,M1,V2> operator*(
        const QuotVM<ix,T,M1,V2>& mv, const CT x)
    { return QuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,CT,M1,V2> operator*(
        const QuotVM<ix,T,M1,V2>& mv, const CCT x)
    { return QuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    inline QuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
        const QuotVM<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
    { 
        return QuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            T1(x)*mv.getX(),mv.getM(),mv.getV()); 
    }

    // (v/m)/x
    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,T,M1,V2> operator/(
        const QuotVM<ix,T,M1,V2>& mv, const RT x)
    { return QuotVM<0,T,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,CT,M1,V2> operator/(
        const QuotVM<ix,T,M1,V2>& mv, const CT x)
    { return QuotVM<0,CT,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline QuotVM<0,CT,M1,V2> operator/(
        const QuotVM<ix,T,M1,V2>& mv, const CCT x)
    { return QuotVM<0,CT,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    inline QuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator/(
        const QuotVM<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
    { 
        return QuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            mv.getX()/T1(x),mv.getM(),mv.getV()); 
    }

#undef RT
#undef CT
#undef CCT

#define RT typename RQuotVM<ix,T,M1,V2>::real_type
#define CT typename RQuotVM<ix,T,M1,V2>::complex_type
#define CCT ConjRef<CT>

    // -(v/m)
    template <int ix, class T, class M1, class V2>
    inline RQuotVM<-ix,T,M1,V2> operator-(const RQuotVM<ix,T,M1,V2>& mv)
    { return RQuotVM<-ix,T,M1,V2>(-mv.getX(),mv.getM(),mv.getV()); }

    // x * (v/m)
    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,T,M1,V2> operator*(
        const RT x, const RQuotVM<ix,T,M1,V2>& mv)
    { return RQuotVM<0,T,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,CT,M1,V2> operator*(
        const CT x, const RQuotVM<ix,T,M1,V2>& mv)
    { return RQuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,CT,M1,V2> operator*(
        const CCT x, const RQuotVM<ix,T,M1,V2>& mv)
    { return RQuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    inline RQuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
        const Scaling<ix1,T1>& x, const RQuotVM<ix,T,M1,V2>& mv)
    { 
        return RQuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            T1(x)*mv.getX(),mv.getM(),mv.getV()); 
    }

    // (v/m)*x
    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,T,M1,V2> operator*(
        const RQuotVM<ix,T,M1,V2>& mv, const RT x)
    { return RQuotVM<0,T,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,CT,M1,V2> operator*(
        const RQuotVM<ix,T,M1,V2>& mv, const CT x)
    { return RQuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,CT,M1,V2> operator*(
        const RQuotVM<ix,T,M1,V2>& mv, const CCT x)
    { return RQuotVM<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    inline RQuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
        const RQuotVM<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
    { 
        return RQuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            T1(x)*mv.getX(),mv.getM(),mv.getV()); 
    }

    // (v/m)/x
    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,T,M1,V2> operator/(
        const RQuotVM<ix,T,M1,V2>& mv, const RT x)
    { return RQuotVM<0,T,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,CT,M1,V2> operator/(
        const RQuotVM<ix,T,M1,V2>& mv, const CT x)
    { return RQuotVM<0,CT,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    inline RQuotVM<0,CT,M1,V2> operator/(
        const RQuotVM<ix,T,M1,V2>& mv, const CCT x)
    { return RQuotVM<0,CT,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    inline RQuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator/(
        const RQuotVM<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
    { 
        return RQuotVM<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            mv.getX()/T1(x),mv.getM(),mv.getV()); 
    }

#undef RT
#undef CT
#undef CCT


    // TMV_Text

    template <int ix, class T, class M1, class V2>
    inline std::string TMV_Text(const QuotVM<ix,T,M1,V2>& qvm)
    {
        std::ostringstream s;
        s << "QuotVM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(qvm.getM())<<" , "<<TMV_Text(qvm.getV())<<" >";
        return s.str();
    }

    template <int ix, class T, class M1, class V2>
    inline std::string TMV_Text(const RQuotVM<ix,T,M1,V2>& qvm)
    {
        std::ostringstream s;
        s << "RQuotVM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(qvm.getM())<<" , "<<TMV_Text(qvm.getV())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
