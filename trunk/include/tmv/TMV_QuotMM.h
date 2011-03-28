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


#ifndef TMV_QuotMM_H
#define TMV_QuotMM_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXM.h"
#include "TMV_QuotXM.h"
#include "TMV_DivMM_Funcs.h"

//#define XDEBUG_QUOTMM

#ifdef XDEBUG_QUOTMM
#include "TMV_ProdMM.h"
#include "TMV_LUD.h"
#include <iostream>
#endif

namespace tmv {

    template <int ix, class T, class M1, class M2, class M3>
    static TMV_INLINE void LDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { LDiv(x,m1.calc(),m2.calc(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static TMV_INLINE void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { NoAliasLDiv(x,m1.calc(),m2.calc(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static TMV_INLINE void AliasLDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { AliasLDiv(x,m1.calc(),m2.calc(),m3.mat()); }

    template <class M1, class M2>
    static TMV_INLINE void LDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { LDivEq(m1.mat(),m2.calc()); }
    template <class M1, class M2>
    static TMV_INLINE void NoAliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { NoAliasLDivEq(m1.mat(),m2.calc()); }
    template <class M1, class M2>
    static TMV_INLINE void AliasLDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { AliasLDivEq(m1.mat(),m2.calc()); }

    template <int ix, class T, class M1, class M2, class M3>
    static TMV_INLINE void RDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { RDiv(x,m1.calc(),m2.calc(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static TMV_INLINE void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { NoAliasRDiv(x,m1.calc(),m2.calc(),m3.mat()); }
    template <int ix, class T, class M1, class M2, class M3>
    static TMV_INLINE void AliasRDiv(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { AliasRDiv(x,m1.calc(),m2.calc(),m3.mat()); }

#if 0
    // Defined below...
    template <class M1, class M2>
    static TMV_INLINE void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { RDivEq(m1.mat(),m2.calc()); }
#endif
    template <class M1, class M2>
    static TMV_INLINE void NoAliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { NoAliasRDivEq(m1.mat(),m2.calc()); }
    template <class M1, class M2>
    static TMV_INLINE void AliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    { AliasRDivEq(m1.mat(),m2.calc()); }


    //
    // Matrix (/%) Matrix
    //

#ifdef XDEBUG_QUOTMM
    template <int ix, class T, class M1, class M2, class M3>
    static void LDiv_Debug(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        // m3 = m2^-1*m1
        // --> m2*m3 = m1
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M2::value_type> m2i = m2.mat();
        Matrix<typename M3::value_type> m3i = m3.mat();

        LDiv(x,m1.mat(),m2.mat(),m3.mat());

        Matrix<typename M3::value_type> m3f = m3.mat();
        Matrix<typename M3::value_type> m1c = m2i*m3f/x;
        const typename M3::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(m1i-m1c) > 1.e-6 * kappa * Norm(m1i)) {
            std::cout<<"LDiv:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            if (m3.colsize() < 100 && m3.rowsize() < 100 && m1.colsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 -> "<<m3f<<std::endl;
                std::cout<<"m1c = "<<m1c<<std::endl;
                std::cout<<"diff = "<<m1i-m1c<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m1i-m1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTMM
    template <int ix, class T, class M1, class M2, class M3>
    static void RDiv_Debug(
        const Scaling<ix,T>& x, const BaseMatrix<M1>& m1,
        const BaseMatrix<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        // m3 = m1*m2^-1
        // --> m3*m2 = m1
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M2::value_type> m2i = m2.mat();
        Matrix<typename M3::value_type> m3i = m3.mat();

        RDiv(x,m1.mat(),m2.mat(),m3.mat());

        Matrix<typename M3::value_type> m3f = m3.mat();
        Matrix<typename M3::value_type> m1c = m3f*m2i/x;
        const typename M3::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(m1i-m1c) > 1.e-6 * kappa * Norm(m1i)) {
            std::cout<<"RDiv:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            if (m3.colsize() < 100 && m3.rowsize() < 100 && m1.colsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 -> "<<m3f<<std::endl;
                std::cout<<"m1c = "<<m1c<<std::endl;
                std::cout<<"diff = "<<m1i-m1c<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m1i-m1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTMM
    template <class M1, class M2>
    static void LDivEq_Debug(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M2::value_type> m2i = m2.mat();

        LDivEq(m1.mat(),m2.mat());

        Matrix<typename M1::value_type> m1f = m1.mat();
        Matrix<typename M1::value_type> m1c = m2i*m1f;
        const typename M1::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(m1i-m1c) > 1.e-6 * kappa * Norm(m1i)) {
            std::cout<<"LDivEq:  \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            if (m1.colsize() < 100 && m1.rowsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m1 -> "<<m1f<<std::endl;
                std::cout<<"m1c = "<<m1c<<std::endl;
                std::cout<<"diff = "<<m1i-m1c<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m1i-m1c)<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_QUOTMM
    template <class M1, class M2>
    static void RDivEq_Debug(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M2::value_type> m2i = m2.mat();

        RDivEq(m1.mat(),m2.calc());

        Matrix<typename M1::value_type> m1f = m1.mat();
        Matrix<typename M1::value_type> m1c = m1f*m2i;
        const typename M1::real_type kappa = Norm(m2i.inverse())*Norm(m2i);

        if (Norm(m1i-m1c) > 1.e-6 * kappa * Norm(m1i)) {
            std::cout<<"RDivEq:  \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            if (m1.colsize() < 100 && m1.rowsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m1 -> "<<m1f<<std::endl;
                std::cout<<"m1c = "<<m1c<<std::endl;
                std::cout<<"diff = "<<m1i-m1c<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m1i-m1c)<<std::endl;
            exit(1);
        }
    }
#endif

    template <int ix, class T, class M1, class M2>
    class QuotMM;

    template <int ix, class T, class M1, class M2>
    struct Traits<QuotMM<ix,T,M1,M2> >
    {
        typedef typename ProdXM<ix,T,M1>::value_type mtype1;
        typedef typename M2::value_type mtype2;
        typedef typename Traits2<mtype1,mtype2>::type value_type;

        enum { _colsize = M2::_rowsize };
        enum { _rowsize = M1::_rowsize };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { shape1 = Traits<QuotXM<ix,T,M2> >::_shape };
        enum { _shape = ShapeTraits2<shape1,M1::_shape>::prod };

        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { cm2 = Traits<typename M2::calc_type>::_colmajor };
        enum { _rowmajor = 
            (rm1 && rm2) || (!cm1 && !cm2 && _rowsize > int(_colsize)) };

        typedef QuotMM<ix,T,M1,M2> type;
        typedef typename MCopyHelper<value_type,_shape,_colsize,_rowsize,
                _rowmajor,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
        typedef typename Traits<calc_type>::inverse_type inverse_type;
    };

    template <int ix, class T, class M1, class M2>
    class QuotMM : 
        public BaseMatrix<QuotMM<ix,T,M1,M2> >
    {
    public:

        typedef QuotMM<ix,T,M1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        QuotMM(
            const T _x, const BaseMatrix<M1>& _m1, const BaseMatrix<M2>& _m2
        ) : 
            x(_x), m1(_m1.mat()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M2::_colsize,M1::_colsize>::same)); 
            TMVAssert(m2.colsize() == m1.colsize());
        }

        TMV_INLINE const Scaling<ix,T>& getX() const { return x; }
        TMV_INLINE const M1& getM1() const { return m1; }
        TMV_INLINE const M2& getM2() const { return m2; }

        TMV_INLINE size_t colsize() const { return m2.rowsize(); }
        TMV_INLINE size_t rowsize() const { return m1.rowsize(); }

        TMV_INLINE value_type cref(int i, int j) const
        { return this->calc().cref(i,j); }

        template <class M3>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_QUOTMM
            LDiv_Debug(x,m1.mat(),m2.mat(),m3.mat());
#else
            LDiv(x,m1.mat(),m2.mat(),m3.mat());
#endif
        }

        template <class M3>
        TMV_INLINE_ND void newAssignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_QUOTMM
            LDiv_Debug(x,m1.mat(),m2.mat(),m3.mat());
#else
            NoAliasLDiv(x,m1.mat(),m2.mat(),m3.mat());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const M1& m1;
        const M2& m2;
    };

    template <int ix, class T, class M1, class M2>
    class RQuotMM;

    template <int ix, class T, class M1, class M2>
    struct Traits<RQuotMM<ix,T,M1,M2> >
    {
        typedef typename ProdXM<ix,T,M1>::value_type mtype1;
        typedef typename M2::value_type mtype2;
        typedef typename Traits2<mtype1,mtype2>::type value_type;

        enum { _colsize = M1::_colsize };
        enum { _rowsize = M2::_colsize };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { shape1 = Traits<QuotXM<ix,T,M2> >::_shape };
        enum { _shape = ShapeTraits2<M1::_shape,shape1>::prod };

        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { cm2 = Traits<typename M2::calc_type>::_colmajor };
        enum { _rowmajor = 
            (rm1 && rm2) || (!cm1 && !cm2 && _rowsize > int(_colsize)) };

        typedef RQuotMM<ix,T,M1,M2> type;
        typedef typename MCopyHelper<value_type,_shape,_colsize,_rowsize,
                _rowmajor,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef const copy_type eval_type;
        typedef typename Traits<calc_type>::inverse_type inverse_type;
    };

    template <int ix, class T, class M1, class M2>
    class RQuotMM : 
        public BaseMatrix<RQuotMM<ix,T,M1,M2> >
    {
    public:

        typedef RQuotMM<ix,T,M1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        RQuotMM(
            const T _x, const BaseMatrix<M1>& _m1, const BaseMatrix<M2>& _m2
        ) : 
            x(_x), m1(_m1.mat()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
            TMVAssert(m1.rowsize() == m2.rowsize());
        }

        const Scaling<ix,T>& getX() const { return x; }
        const M1& getM1() const { return m1; }
        const M2& getM2() const { return m2; }

        size_t colsize() const { return m1.colsize(); }
        size_t rowsize() const { return m2.colsize(); }

        value_type cref(int i, int j) const
        { return this->calc().cref(i,j); }

        template <class M3>
        void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_QUOTMM
            RDiv_Debug(x,m1.mat(),m2.mat(),m3.mat());
#else
            RDiv(x,m1.mat(),m2.mat(),m3.mat());
#endif
        }

        template <class M3>
        void newAssignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_QUOTMM
            RDiv_Debug(x,m1.mat(),m2.mat(),m3.mat());
#else
            NoAliasRDiv(x,m1.mat(),m2.mat(),m3.mat());
#endif
        }

    private:
        const Scaling<ix,T> x;
        const M1& m1;
        const M2& m2;
    };



    // m / m
#define RT typename M2::real_type
    template <class M1, class M2>
    static TMV_INLINE QuotMM<1,RT,M1,M2> operator/(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return QuotMM<1,RT,M1,M2>(RT(1),m1,m2); }
#undef RT

    // m / xm
    template <class M1, int ix, class T, class M2>
    static TMV_INLINE QuotMM<ix,T,M1,M2> operator/(
        const BaseMatrix<M1>& m1, const ProdXM<ix,T,M2>& m2)
    {
        return QuotMM<ix,T,M1,M2>(
            typename Traits<T>::real_type(1)/m2.getX(),m1,m2.getM()); 
    }

    // xm / m
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<ix,T,M1,M2> operator/(
        const ProdXM<ix,T,M1>& m1, const BaseMatrix<M2>& m2)
    { return QuotMM<ix,T,M1,M2>(m1.getX(),m1.getM(),m2); }

    // xm / xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static TMV_INLINE QuotMM<ix1*ix2,PT,M1,M2> operator/(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return QuotMM<ix1*ix2,PT,M1,M2>(m1.getX()/m2.getX(),m1.getM(),m2.getM()); }
#undef PT

    // x/m * m
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<ix,T,M1,M2> operator*(
        const QuotXM<ix,T,M2>& m2, const BaseMatrix<M1>& m1)
    { return QuotMM<ix,T,M1,M2>(m2.getX(),m1,m2.getM()); }

    // x/m * xm
#define PT typename Traits2<T1,T2>::type
    template <int ix2, class T2, class M2, int ix1, class T1, class M1>
    static TMV_INLINE QuotMM<ix1*ix2,PT,M1,M2> operator*(
        const QuotXM<ix2,T2,M2>& m2, const ProdXM<ix1,T1,M1>& m1)
    {
        return QuotMM<ix1*ix2,PT,M1,M2>(
            m1.getX()*m2.getX(),m1.getM(),m2.getM()); 
    }
#undef PT


    // m /= m
    template <class M1, class M2>
    static TMV_INLINE void DivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
#ifdef XDEBUG_QUOTMM
        LDivEq_Debug(m1.mat(),m2.mat()); 
#else
        LDivEq(m1.mat(),m2.mat()); 
#endif
    }

    // m /= xm
    template <class M1, int ix2, class T2, class M2>
    static TMV_INLINE void DivEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
#ifdef XDEBUG_QUOTMM
        LDivEq_Debug(m1.mat(),m2.getM().mat()); 
#else
        LDivEq(m1.mat(),m2.getM().mat()); 
#endif
        Scale(typename Traits<T2>::real_type(1)/m2.getX(),m1.mat());
    }


    // m % m
#define RT typename M2::real_type
    template <class M1, class M2>
    static TMV_INLINE RQuotMM<1,RT,M1,M2> operator%(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return RQuotMM<1,RT,M1,M2>(RT(1),m1,m2); }
#undef RT

    // m % xm
    template <class M1, int ix, class T, class M2>
    static TMV_INLINE RQuotMM<ix,T,M1,M2> operator%(
        const BaseMatrix<M1>& m1, const ProdXM<ix,T,M2>& m2)
    {
        return RQuotMM<ix,T,M1,M2>(
            typename Traits<T>::real_type(1)/m2.getX(),m1,m2.getM()); 
    }

    // xm % m
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<ix,T,M1,M2> operator%(
        const ProdXM<ix,T,M1>& m1, const BaseMatrix<M2>& m2)
    { return RQuotMM<ix,T,M1,M2>(m1.getX(),m1.getM(),m2); }

    // xm % xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static TMV_INLINE RQuotMM<ix1*ix2,PT,M1,M2> operator%(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        return RQuotMM<ix1*ix2,PT,M1,M2>(
            m1.getX()/m2.getX(),m1.getM(),m2.getM()); 
    }
#undef PT

    // m * x/m
    template <class M1, int ix, class T, class M2>
    static TMV_INLINE RQuotMM<ix,T,M1,M2> operator*(
        const BaseMatrix<M1>& m1, const QuotXM<ix,T,M2>& m2)
    { return RQuotMM<ix,T,M1,M2>(m2.getX(),m1,m2.getM()); }

    // xm * x/m
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static TMV_INLINE RQuotMM<ix1*ix2,PT,M1,M2> operator*(
        const ProdXM<ix1,T1,M1>& m1, const QuotXM<ix2,T2,M2>& m2)
    {
        return RQuotMM<ix1*ix2,PT,M1,M2>(
            m1.getX()*m2.getX(),m1.getM(),m2.getM()); 
    }
#undef PT

    // m %= m
    template <class M1, class M2>
    static TMV_INLINE void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
#ifdef XDEBUG_QUOTMM
        RDivEq_Debug(m1.mat(),m2.mat()); 
#else
        RDivEq(m1.mat(),m2.calc()); 
#endif
    }

    // m %= xm
    template <class M1, int ix2, class T2, class M2>
    static TMV_INLINE void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
#ifdef XDEBUG_QUOTMM
        RDivEq_Debug(m1.mat(),m2.getM().mat()); 
#else
        RDivEq(m1.mat(),m2.getM().calc()); 
#endif
        Scale(typename Traits<T2>::real_type(1)/m2.getX(),m1.mat());
    }

    // m *= x/m
    template <class M1, int ix2, class T2, class M2>
    static TMV_INLINE void MultEq(
        BaseMatrix_Mutable<M1>& m1, const QuotXM<ix2,T2,M2>& m2)
    {
#ifdef XDEBUG_QUOTMM
        RDivEq_Debug(m1.mat(),m2.getM().mat()); 
#else
        RDivEq(m1.mat(),m2.getM().calc()); 
#endif
        Scale(m2.getX(),m1.mat());
    }


    // Consolidate x*(xm/m) type constructs:

#define RT typename QuotMM<ix,T,M1,M2>::real_type
#define CT typename QuotMM<ix,T,M1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(m/m)
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<-ix,T,M1,M2> operator-(const QuotMM<ix,T,M1,M2>& qmm)
    { return QuotMM<-ix,T,M1,M2>(-qmm.getX(),qmm.getM1(),qmm.getM2()); }

    // x * (m/m)
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,T,M1,M2> operator*(
        const RT x, const QuotMM<ix,T,M1,M2>& qmm)
    { return QuotMM<0,T,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,CT,M1,M2> operator*(
        const CT x, const QuotMM<ix,T,M1,M2>& qmm)
    { return QuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,CT,M1,M2> operator*(
        const CCT x, const QuotMM<ix,T,M1,M2>& qmm)
    { return QuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
        const Scaling<ix1,T1>& x, const QuotMM<ix,T,M1,M2>& qmm)
    {
        return QuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            T1(x)*qmm.getX(),qmm.getM1(),qmm.getM2()); 
    }

    // (m/m)*x
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,T,M1,M2> operator*(
        const QuotMM<ix,T,M1,M2>& qmm, const RT x)
    { return QuotMM<0,T,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,CT,M1,M2> operator*(
        const QuotMM<ix,T,M1,M2>& qmm, const CT x)
    { return QuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,CT,M1,M2> operator*(
        const QuotMM<ix,T,M1,M2>& qmm, const CCT x)
    { return QuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
        const QuotMM<ix,T,M1,M2>& qmm, const Scaling<ix1,T1>& x)
    {
        return QuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            T1(x)*qmm.getX(),qmm.getM1(),qmm.getM2()); 
    }

    // (m/m)/x
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,T,M1,M2> operator/(
        const QuotMM<ix,T,M1,M2>& qmm, const RT x)
    { return QuotMM<0,T,M1,M2>(qmm.getX()/x,qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,CT,M1,M2> operator/(
        const QuotMM<ix,T,M1,M2>& qmm, const CT x)
    { return QuotMM<0,CT,M1,M2>(qmm.getX()/x,qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<0,CT,M1,M2> operator/(
        const QuotMM<ix,T,M1,M2>& qmm, const CCT x)
    { return QuotMM<0,CT,M1,M2>(qmm.getX()/x,qmm.getM1(),qmm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    static TMV_INLINE QuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator/(
        const QuotMM<ix,T,M1,M2>& qmm, const Scaling<ix1,T1>& x)
    {
        return QuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            qmm.getX()/T1(x),qmm.getM1(),qmm.getM2()); 
    }

#undef RT
#undef CT
#undef CCT

#define RT typename RQuotMM<ix,T,M1,M2>::real_type
#define CT typename RQuotMM<ix,T,M1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(m/m)
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<-ix,T,M1,M2> operator-(const RQuotMM<ix,T,M1,M2>& qmm)
    { return RQuotMM<-ix,T,M1,M2>(-qmm.getX(),qmm.getM1(),qmm.getM2()); }

    // x * (m/m)
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,T,M1,M2> operator*(
        const RT x, const RQuotMM<ix,T,M1,M2>& qmm)
    { return RQuotMM<0,T,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,CT,M1,M2> operator*(
        const CT x, const RQuotMM<ix,T,M1,M2>& qmm)
    { return RQuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,CT,M1,M2> operator*(
        const CCT x, const RQuotMM<ix,T,M1,M2>& qmm)
    { return RQuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
        const Scaling<ix1,T1>& x, const RQuotMM<ix,T,M1,M2>& qmm)
    {
        return RQuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            T1(x)*qmm.getX(),qmm.getM1(),qmm.getM2()); 
    }

    // (m/m)*x
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,T,M1,M2> operator*(
        const RQuotMM<ix,T,M1,M2>& qmm, const RT x)
    { return RQuotMM<0,T,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,CT,M1,M2> operator*(
        const RQuotMM<ix,T,M1,M2>& qmm, const CT x)
    { return RQuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,CT,M1,M2> operator*(
        const RQuotMM<ix,T,M1,M2>& qmm, const CCT x)
    { return RQuotMM<0,CT,M1,M2>(x*qmm.getX(),qmm.getM1(),qmm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
        const RQuotMM<ix,T,M1,M2>& qmm, const Scaling<ix1,T1>& x)
    {
        return RQuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            T1(x)*qmm.getX(),qmm.getM1(),qmm.getM2()); 
    }

    // (m/m)/x
    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,T,M1,M2> operator/(
        const RQuotMM<ix,T,M1,M2>& qmm, const RT x)
    { return RQuotMM<0,T,M1,M2>(qmm.getX()/x,qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,CT,M1,M2> operator/(
        const RQuotMM<ix,T,M1,M2>& qmm, const CT x)
    { return RQuotMM<0,CT,M1,M2>(qmm.getX()/x,qmm.getM1(),qmm.getM2()); }

    template <int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<0,CT,M1,M2> operator/(
        const RQuotMM<ix,T,M1,M2>& qmm, const CCT x)
    { return RQuotMM<0,CT,M1,M2>(qmm.getX()/x,qmm.getM1(),qmm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    static TMV_INLINE RQuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator/(
        const RQuotMM<ix,T,M1,M2>& qmm, const Scaling<ix1,T1>& x)
    {
        return RQuotMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            qmm.getX()/T1(x),qmm.getM1(),qmm.getM2()); 
    }

#undef RT
#undef CT
#undef CCT


    // TMV_Text

#ifdef TMV_DEBUG
    template <int ix, class T, class M1, class M2>
    static inline std::string TMV_Text(const QuotMM<ix,T,M1,M2>& qmm)
    {
        std::ostringstream s;
        s << "QuotMM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(qmm.getM1())<<" , "<<TMV_Text(qmm.getM2())<<" >";
        return s.str();
    }

    template <int ix, class T, class M1, class M2>
    static inline std::string TMV_Text(const RQuotMM<ix,T,M1,M2>& qmm)
    {
        std::ostringstream s;
        s << "RQuotMM< "<<ix<<","<<TMV_Text(T())<<" , ";
        s << TMV_Text(qmm.getM1())<<" , "<<TMV_Text(qmm.getM2())<<" >";
        return s.str();
    }
#endif

} // namespace tmv

#endif 
