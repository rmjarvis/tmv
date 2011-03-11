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


#ifndef TMV_ProdMV_H
#define TMV_ProdMV_H

#include "TMV_ProdXM.h"
#include "TMV_ProdXV.h"

//#define XDEBUG_PRODMV

#ifdef XDEBUG_PRODMV
#include "TMV_Vector.h"
#include <iostream>
#endif

namespace tmv {

    //
    // Matrix * Vector
    //

    // These first few are for when an argument is a composite matrix
    // and needs to be calculated before running MultMV.
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static void MultMV(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(x,m1.calc(),v2.calc(),v3.vec()); }
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static void NoAliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { NoAliasMultMV<add>(x,m1.calc(),v2.calc(),v3.vec()); }
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static void AliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1, 
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { AliasMultMV<add>(x,m1.calc(),v2.calc(),v3.vec()); }

    // These are helpers to allow the caller to not use a Scaling object.
    template <bool add, class T, class M1, class V2, class V3>
    static void MultMV(
        const T& x, const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { MultMV<add>(Scaling<0,T>(x),m1.mat(),v2.vec(),v3.vec()); }
    template <bool add, class T, class M1, class V2, class V3>
    static void NoAliasMultMV(
        const T& x, const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { NoAliasMultMV<add>(Scaling<0,T>(x),m1.mat(),v2.vec(),v3.vec()); }
    template <bool add, class T, class M1, class V2, class V3>
    static void AliasMultMV(
        const T& x, const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { AliasMultMV<add>(Scaling<0,T>(x),m1.mat(),v2.vec(),v3.vec()); }

    template <bool add, class M1, class V2, class V3>
    static void MultMV(
        const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        MultMV<add>(
            Scaling<1,typename V3::real_type>(),m1.mat(),v2.vec(),v3.vec()); 
    }
    template <bool add, class M1, class V2, class V3>
    static void NoAliasMultMV(
        const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { 
        NoAliasMultMV<add>(
            Scaling<1,typename V3::real_type>(),m1.mat(),v2.vec(),v3.vec()); 
    }
    template <bool add, class M1, class V2, class V3>
    static void AliasMultMV(
        const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    { 
        AliasMultMV<add>(
            Scaling<1,typename V3::real_type>(),m1.mat(),v2.vec(),v3.vec());
    }

#ifdef XDEBUG_PRODMV
    template <bool add, int ix, class T, class M1, class V2, class V3>
    static void MultMV_Debug(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseVector<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        Vector<typename V3::value_type> v3i = v3;
        Vector<typename V3::value_type> v3c = v3;
        if (!add) v3c.setZero();
        for(size_t i=0;i<v3.size();++i) {
            for(size_t j=0;j<v2.size();++j) {
                v3c.ref(i) += T(x) * m1.cref(i,j) * v2.cref(j);
            }
        }

        MultMV<add>(x,m1.mat(),v2.vec(),v3.vec());

        if (Norm(v3-v3c) > 1.e-6 * std::abs(T(x)) * Norm(m1) * Norm(v2)) {
            std::cout<<"MultMV:  add = "<<add<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<"  "<<v2<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3i<<std::endl;
            std::cout<<"v3 -> "<<v3<<std::endl;
            std::cout<<"correct = "<<v3c<<std::endl;
            std::cout<<"diff = "<<v3-v3c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v3-v3c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_PRODMV
    template <bool add, int ix, class T, class V1, class M2, class V3>
    static void MultVM_Debug(
        const Scaling<ix,T>& x, const BaseVector<V1>& v1,
        const BaseMatrix<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        Vector<typename V3::value_type> v3i = v3;
        Vector<typename V3::value_type> v3c = v3;
        if (!add) v3c.setZero();
        for(size_t i=0;i<v1.size();++i) {
            for(size_t j=0;j<v3.size();++j) {
                v3c.ref(j) += T(x) * v1.cref(i) *  m2.cref(i,j);
            }
        }

        MultVM<add>(x,v1.vec(),m2.mat(),v3.vec());

        if (Norm(v3-v3c) > 1.e-6 * std::abs(T(x)) * Norm(m2) * Norm(v1)) {
            std::cout<<"MultVM:  add = "<<add<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<"  "<<v3i<<std::endl;
            std::cout<<"v3 -> "<<v3<<std::endl;
            std::cout<<"correct = "<<v3c<<std::endl;
            std::cout<<"diff = "<<v3-v3c<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v3-v3c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_PRODMV
    template <class V1, int ix, class T, class M2>
    static void MultEqVM_Debug(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix<M2>& m2)
    {
        Vector<typename V1::value_type> v1i = v1;
        Vector<typename V1::value_type> v3 = v1;
        v3.setZero();
        for(size_t i=0;i<v1.size();++i) {
            for(size_t j=0;j<v1.size();++j) {
                v3.ref(j) += T(x) * v1.cref(i) *  m2.cref(i,j);
            }
        }

        MultEqVM(v1.vec(),x,m2.mat());

        if (Norm(v1-v3) > 1.e-6 * std::abs(T(x)) * Norm(m2) * Norm(v1i)) {
            std::cout<<"MultEqVM:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<"  "<<v1i<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
            std::cout<<"v1 -> "<<v1<<std::endl;
            std::cout<<"correct = "<<v3<<std::endl;
            std::cout<<"diff = "<<v1-v3<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(v1-v3)<<std::endl;
            exit(1);
        }
    }
#endif

    template <int ix, class T, class M1, class V2>
    class ProdMV;

    template <int ix, class T, class M1, class V2>
    struct Traits<ProdMV<ix,T,M1,V2> >
    {
        typedef typename ProdXM<ix,T,M1>::value_type mtype1;
        typedef typename V2::value_type vtype2;
        typedef typename Traits2<mtype1,vtype2>::type value_type;

        enum { _size = M1::_colsize };
        enum { _fort = M1::_fort && V2::_fort };
        enum { _calc = false };

        typedef ProdMV<ix,T,M1,V2> type;
        typedef typename VCopyHelper<value_type,_size,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (M1::_calc && V2::_calc),const type,calc_type>::type eval_type;
    };

    template <int ix, class T, class M1, class V2>
    class ProdMV : public BaseVector<ProdMV<ix,T,M1,V2> >
    {
    public:

        typedef ProdMV<ix,T,M1,V2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ProdMV(
            const T _x, const BaseMatrix<M1>& _m1, const BaseVector<V2>& _v2
        ) : 
            x(_x), m1(_m1.mat()), v2(_v2.vec())
        {
            TMVStaticAssert((Sizes<M1::_rowsize,V2::_size>::same)); 
            TMVAssert(m1.rowsize() == v2.size());
        }

        const Scaling<ix,T>& getX() const { return x; }
        const M1& getM() const { return m1; }
        const V2& getV() const { return v2; }

        size_t size() const { return m1.colsize(); }

        value_type cref(int i) const
        { return x * (m1.get_row(i) * v2); }

        template <class V3>
        void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::isreal || V3::iscomplex));
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_PRODMV
            MultMV_Debug<false>(x,m1.mat(),v2.vec(),v3.vec());
#else
            MultMV<false>(x,m1.mat(),v2.vec(),v3.vec());
#endif
        }

        template <class V3>
        void newAssignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::isreal || V3::iscomplex));
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_PRODMV
            MultMV_Debug<false>(x,m1.mat(),v2.vec(),v3.vec());
#else
            NoAliasMultMV<false>(x,m1.mat(),v2.vec(),v3.vec());
#endif
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

        enum { _size = M2::_rowsize };
        enum { _fort = V1::_fort && M2::_fort };
        enum { _calc = false };

        typedef ProdVM<ix,T,V1,M2> type;
        typedef typename VCopyHelper<value_type,_size,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (V1::_calc && M2::_calc),const type,calc_type>::type eval_type;
    };

    template <int ix, class T, class V1, class M2>
    class ProdVM : public BaseVector<ProdVM<ix,T,V1,M2> >
    {
    public:

        typedef ProdVM<ix,T,V1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        ProdVM(
            const T _x, const BaseVector<V1>& _v1, const BaseMatrix<M2>& _m2
        ) : 
            x(_x), v1(_v1.vec()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<V1::_size,M2::_colsize>::same)); 
            TMVAssert(v1.size() == m2.colsize());
        }

        const Scaling<ix,T>& getX() const { return x; }
        const V1& getV() const { return v1; }
        const M2& getM() const { return m2; }

        size_t size() const { return m2.rowsize(); }

        value_type cref(int j) const
        { return x * (v1 * m2.get_col(j)); }

        template <class V3>
        void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::isreal || V3::iscomplex));
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_PRODMV
            MultVM_Debug<false>(x,v1.vec(),m2.mat(),v3.vec());
#else
            MultVM<false>(x,v1.vec(),m2.mat(),v3.vec());
#endif
        }

        template <class V3>
        void newAssignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((type::isreal || V3::iscomplex));
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
#ifdef XDEBUG_PRODMV
            MultVM_Debug<false>(x,v1.vec(),m2.mat(),v3.vec());
#else
            NoAliasMultVM<false>(x,v1.vec(),m2.mat(),v3.vec());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const V1& v1;
        const M2& m2;
    };


    // m * v
#define RT typename V::real_type
    template <class M, class V>
    static ProdMV<1,RT,M,V> operator*(
        const BaseMatrix<M>& m, const BaseVector<V>& v)
    { return ProdMV<1,RT,M,V>(RT(1),m,v); }
#undef RT

    // m * xv
    template <class M, int ix, class T, class V>
    static ProdMV<ix,T,M,V> operator*(
        const BaseMatrix<M>& m, const ProdXV<ix,T,V>& v)
    { return ProdMV<ix,T,M,V>(v.getX(),m,v.getV()); }

    // xm * v
    template <int ix, class T, class M, class V>
    static ProdMV<ix,T,M,V> operator*(
        const ProdXM<ix,T,M>& m, const BaseVector<V>& v)
    { return ProdMV<ix,T,M,V>(m.getX(),m.getM(),v); }

    // xm * xv
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M, int ix2, class T2, class V>
    static ProdMV<ix1*ix2,PT,M,V> operator*(
        const ProdXM<ix1,T1,M>& m, const ProdXV<ix2,T2,V>& v)
    { return ProdMV<ix1*ix2,PT,M,V>(m.getX()*v.getX(),m.getM(),v.getV()); }
#undef PT


    // v * m
#define RT typename V::real_type
    template <class M, class V>
    static ProdVM<1,RT,V,M> operator*(
        const BaseVector<V>& v, const BaseMatrix<M>& m)
    { return ProdVM<1,RT,V,M>(RT(1),v,m); }
#undef RT

    // v * xm
    template <class M, int ix, class T, class V>
    static ProdVM<ix,T,V,M> operator*(
        const BaseVector<V>& v, const ProdXM<ix,T,M>& m)
    { return ProdVM<ix,T,V,M>(m.getX(),v,m.getM()); }

    // xv * m
    template <int ix, class T, class M, class V>
    static ProdVM<ix,T,V,M> operator*(
        const ProdXV<ix,T,V>& v, const BaseMatrix<M>& m)
    { return ProdVM<ix,T,V,M>(v.getX(),v.getV(),m); }

    // xv * xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M, int ix2, class T2, class V>
    static ProdVM<ix1*ix2,PT,V,M> operator*(
        const ProdXV<ix1,T1,V>& v, const ProdXM<ix2,T2,M>& m)
    { return ProdVM<ix1*ix2,PT,V,M>(m.getX()*v.getX(),v.getV(),m.getM()); }
#undef PT



    // v *= m
    template <class V1, class M2>
    static void MultEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix<M2>& m2)
    { 
#ifdef XDEBUG_PRODMV
        MultEqVM_Debug(v1,Scaling<1,typename V1::real_type>(),m2.mat()); 
#else
        MultEqVM(v1,Scaling<1,typename V1::real_type>(),m2.mat()); 
#endif
    }

    // v *= xm
    template <class V1, int ix2, class T2, class M2>
    static void MultEq(
        BaseVector_Mutable<V1>& v1, const ProdXM<ix2,T2,M2>& m2)
    { 
#ifdef XDEBUG_PRODMV
        MultEqVM_Debug(v1,m2.getX(),m2.getM().mat()); 
#else
        MultEqVM(v1,m2.getX(),m2.getM().mat()); 
#endif
    }

    // v += mv
    template <class V3, int ix, class T, class M1, class V2>
    static void AddEq(
        BaseVector_Mutable<V3>& v3, const ProdMV<ix,T,M1,V2>& mv)
    { 
#ifdef XDEBUG_PRODMV
        MultMV_Debug<true>(
            mv.getX(),mv.getM().mat(),mv.getV().vec(),v3.vec()); 
#else
        MultMV<true>(mv.getX(),mv.getM().mat(),mv.getV().vec(),v3.vec()); 
#endif
    }

    // v -= mv
    template <class V3, int ix, class T, class M1, class V2>
    static void SubtractEq(
        BaseVector_Mutable<V3>& v3, const ProdMV<ix,T,M1,V2>& mv)
    {
#ifdef XDEBUG_PRODMV
        MultMV_Debug<true>(
            -mv.getX(),mv.getM().mat(),mv.getV().vec(),v3.vec()); 
#else
        MultMV<true>(-mv.getX(),mv.getM().mat(),mv.getV().vec(),v3.vec()); 
#endif
    }

    // v += vm
    template <class V3, int ix, class T, class M1, class V2>
    static void AddEq(
        BaseVector_Mutable<V3>& v3, const ProdVM<ix,T,M1,V2>& vm)
    { 
#ifdef XDEBUG_PRODMV
        MultVM_Debug<true>(
            vm.getX(),vm.getV().vec(),vm.getM().mat(),v3.vec()); 
#else
        MultVM<true>(vm.getX(),vm.getV().vec(),vm.getM().mat(),v3.vec()); 
#endif
    }

    // v -= vm
    template <class V3, int ix, class T, class M1, class V2>
    static void SubtractEq(
        BaseVector_Mutable<V3>& v3, const ProdVM<ix,T,M1,V2>& vm)
    { 
#ifdef XDEBUG_PRODMV
        MultVM_Debug<true>(
            -vm.getX(),vm.getV().vec(),vm.getM().mat(),v3.vec()); 
#else
        MultVM<true>(-vm.getX(),vm.getV().vec(),vm.getM().mat(),v3.vec()); 
#endif
    }


    // Consolidate x*(xmv) type constructs:

#define RT typename ProdMV<ix,T,M1,V2>::real_type
#define CT typename ProdMV<ix,T,M1,V2>::complex_type
#define CCT ConjRef<CT>

    // -(mv)
    template <int ix, class T, class M1, class V2>
    static ProdMV<-ix,T,M1,V2> operator-(const ProdMV<ix,T,M1,V2>& mv)
    { return ProdMV<-ix,T,M1,V2>(-mv.getX(),mv.getM(),mv.getV()); }

    // x * (mv)
    template <int ix, class T, class M1, class V2>
    static ProdMV<0,T,M1,V2> operator*(
        const RT x, const ProdMV<ix,T,M1,V2>& mv)
    { return ProdMV<0,T,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    static ProdMV<0,CT,M1,V2> operator*(
        const CT x, const ProdMV<ix,T,M1,V2>& mv)
    { return ProdMV<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    static ProdMV<0,CT,M1,V2> operator*(
        const CCT x, const ProdMV<ix,T,M1,V2>& mv)
    { return ProdMV<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    static ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
        const Scaling<ix1,T1>& x, const ProdMV<ix,T,M1,V2>& mv)
    { 
        return ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            T1(x)*mv.getX(),mv.getM(),mv.getV()); 
    }

    // (mv)*x
    template <int ix, class T, class M1, class V2>
    static ProdMV<0,T,M1,V2> operator*(
        const ProdMV<ix,T,M1,V2>& mv, const RT x)
    { return ProdMV<0,T,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    static ProdMV<0,CT,M1,V2> operator*(
        const ProdMV<ix,T,M1,V2>& mv, const CT x)
    { return ProdMV<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    static ProdMV<0,CT,M1,V2> operator*(
        const ProdMV<ix,T,M1,V2>& mv, const CCT x)
    { return ProdMV<0,CT,M1,V2>(x*mv.getX(),mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    static ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator*(
        const ProdMV<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
    { 
        return ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            T1(x)*mv.getX(),mv.getM(),mv.getV()); 
    }

    // (mv)/x
    template <int ix, class T, class M1, class V2>
    static ProdMV<0,T,M1,V2> operator/(
        const ProdMV<ix,T,M1,V2>& mv, const RT x)
    { return ProdMV<0,T,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    static ProdMV<0,CT,M1,V2> operator/(
        const ProdMV<ix,T,M1,V2>& mv, const CT x)
    { return ProdMV<0,CT,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix, class T, class M1, class V2>
    static ProdMV<0,CT,M1,V2> operator/(
        const ProdMV<ix,T,M1,V2>& mv, const CCT x)
    { return ProdMV<0,CT,M1,V2>(mv.getX()/x,mv.getM(),mv.getV()); }

    template <int ix1, class T1, int ix, class T, class M1, class V2>
    static ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2> operator/(
        const ProdMV<ix,T,M1,V2>& mv, const Scaling<ix1,T1>& x)
    { 
        return ProdMV<ix1*ix,typename Traits2<T1,T>::type,M1,V2>(
            mv.getX()/T1(x),mv.getM(),mv.getV()); 
    }

#undef RT
#undef CT
#undef CCT

#define RT typename ProdVM<ix,T,V1,M2>::real_type
#define CT typename ProdVM<ix,T,V1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(vm)
    template <int ix, class T, class V1, class M2>
    static ProdVM<-ix,T,V1,M2> operator-(const ProdVM<ix,T,V1,M2>& vm)
    { return ProdVM<-ix,T,V1,M2>(-vm.getX(),vm.getV(),vm.getM()); }

    // x * (vm)
    template <int ix, class T, class V1, class M2>
    static ProdVM<0,T,V1,M2> operator*(
        const int x, const ProdVM<ix,T,V1,M2>& vm)
    { return ProdVM<0,T,V1,M2>(RT(x)*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,T,V1,M2> operator*(
        const RT x, const ProdVM<ix,T,V1,M2>& vm)
    { return ProdVM<0,T,V1,M2>(x*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,CT,V1,M2> operator*(
        const CT x, const ProdVM<ix,T,V1,M2>& vm)
    { return ProdVM<0,CT,V1,M2>(x*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,CT,V1,M2> operator*(
        const CCT x, const ProdVM<ix,T,V1,M2>& vm)
    { return ProdVM<0,CT,V1,M2>(x*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix1, class T1, int ix, class T, class V1, class M2>
    static ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2> operator*(
        const Scaling<ix1,T1>& x, const ProdVM<ix,T,V1,M2>& vm)
    { 
        return ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2>(
            T1(x)*vm.getX(),vm.getV(),vm.getM()); 
    }

    // (vm)*x
    template <int ix, class T, class V1, class M2>
    static ProdVM<0,T,V1,M2> operator*(
        const ProdVM<ix,T,V1,M2>& vm, const int x)
    { return ProdVM<0,T,V1,M2>(RT(x)*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,T,V1,M2> operator*(
        const ProdVM<ix,T,V1,M2>& vm, const RT x)
    { return ProdVM<0,T,V1,M2>(x*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,CT,V1,M2> operator*(
        const ProdVM<ix,T,V1,M2>& vm, const CT x)
    { return ProdVM<0,CT,V1,M2>(x*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,CT,V1,M2> operator*(
        const ProdVM<ix,T,V1,M2>& vm, const CCT x)
    { return ProdVM<0,CT,V1,M2>(x*vm.getX(),vm.getV(),vm.getM()); }

    template <int ix1, class T1, int ix, class T, class V1, class M2>
    static ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2> operator*(
        const ProdVM<ix,T,V1,M2>& vm, const Scaling<ix1,T1>& x)
    { 
        return ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2>(
            T1(x)*vm.getX(),vm.getV(),vm.getM()); 
    }

    // (vm)/x
    template <int ix, class T, class V1, class M2>
    static ProdVM<0,T,V1,M2> operator/(
        const ProdVM<ix,T,V1,M2>& vm, const int x)
    { return ProdVM<0,T,V1,M2>(vm.getX()/RT(x),vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,T,V1,M2> operator/(
        const ProdVM<ix,T,V1,M2>& vm, const RT x)
    { return ProdVM<0,T,V1,M2>(vm.getX()/x,vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,CT,V1,M2> operator/(
        const ProdVM<ix,T,V1,M2>& vm, const CT x)
    { return ProdVM<0,CT,V1,M2>(vm.getX()/x,vm.getV(),vm.getM()); }

    template <int ix, class T, class V1, class M2>
    static ProdVM<0,CT,V1,M2> operator/(
        const ProdVM<ix,T,V1,M2>& vm, const CCT x)
    { return ProdVM<0,CT,V1,M2>(vm.getX()/x,vm.getV(),vm.getM()); }

    template <int ix1, class T1, int ix, class T, class V1, class M2>
    static ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2> operator/(
        const ProdVM<ix,T,V1,M2>& vm, const Scaling<ix1,T1>& x)
    { 
        return ProdVM<ix1*ix,typename Traits2<T1,T>::type,V1,M2>(
            vm.getX()/T1(x),vm.getV(),vm.getM()); 
    }

#undef RT
#undef CT
#undef CCT

    // TMV_Text

    template <int ix, class T, class M1, class V2>
    static std::string TMV_Text(const ProdMV<ix,T,M1,V2>& mv)
    {
        std::ostringstream s;
        s << "ProdMV< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(mv.getM())<<" , "<<TMV_Text(mv.getV())<<" >";
        return s.str();
    }

    template <int ix, class T, class V1, class M2>
    static std::string TMV_Text(const ProdVM<ix,T,V1,M2>& vm)
    {
        std::ostringstream s;
        s << "ProdVM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(vm.getV())<<" , "<<TMV_Text(vm.getM())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
