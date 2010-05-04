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


#ifndef TMV_ProdMM_H
#define TMV_ProdMM_H

#include "TMV_ProdXM.h"

//#define XDEBUG_PRODMM

#ifdef XDEBUG_PRODMM
#include "TMV_Matrix.h"
#include <iostream>
#endif

namespace tmv {

    //
    // Matrix * Matrix
    //

    // These first few are intentionally not defined to make sure we
    // get a compiler error if they are used.
    // All real calls should go through a more specific version than 
    // just the BaseMatrix_Calc's.
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3);

    // These are helpers to allow the caller to not use a Scaling object.
    template <bool add, class T, class M1, class M2, class M3>
    inline void MultMM(
        const T& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { MultMM<add>(Scaling<0,T>(x),m1.mat(),m2.mat(),m3.mat()); }
    template <bool add, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const T& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { NoAliasMultMM<add>(Scaling<0,T>(x),m1.mat(),m2.mat(),m3.mat()); }
    template <bool add, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const T& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { InlineMultMM<add>(Scaling<0,T>(x),m1.mat(),m2.mat(),m3.mat()); }
    template <bool add, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const T& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { AliasMultMM<add>(Scaling<0,T>(x),m1.mat(),m2.mat(),m3.mat()); }

    template <bool add, class M1, class M2, class M3>
    inline void MultMM(
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        MultMM<add>(
            Scaling<1,typename M3::real_type>(),m1.mat(),m2.mat(),m3.mat()); 
    }
    template <bool add, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { 
        NoAliasMultMM<add>(
            Scaling<1,typename M3::real_type>(),m1.mat(),m2.mat(),m3.mat()); 
    }
    template <bool add, class M1, class M2, class M3>
    inline void InlineMultMM(
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    { 
        InlineMultMM<add>(
            Scaling<1,typename M3::real_type>(),m1.mat(),m2.mat(),m3.mat()); 
    }
    template <bool add, class M1, class M2, class M3>
    inline void AliasMultMM(
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        AliasMultMM<add>(
            Scaling<1,typename M3::real_type>(),m1.mat(),m2.mat(),m3.mat()); 
    }

#ifdef XDEBUG_PRODMM
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM_Debug(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Calc<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M2::value_type> m2i = m2.mat();
        Matrix<typename M3::value_type> m3i = m3.mat();
        Matrix<typename M3::value_type> m3c = m3.mat();
        if (!add) m3c.setZero();
        for(size_t i=0;i<m3.colsize();++i) {
            for(size_t j=0;j<m3.rowsize();++j) {
                for(size_t k=0;k<m1.rowsize();++k) {
                    m3c.ref(i,j) += T(x) * m1.cref(i,k) * m2.cref(k,j);
                }
            }
        }
        //std::cout<<"Start MultMM XDEBUG"<<std::endl;
        //std::cout<<"add = "<<add<<std::endl;
        //std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
        //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        //std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
        //std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3.mat()<<std::endl;
        //std::cout<<"m3i = "<<TMV_Text(m3i)<<"  "<<m3i<<std::endl;
        //std::cout<<"m3c = "<<TMV_Text(m3c)<<"  "<<m3c<<std::endl;
        //std::cout<<"m3c => "<<m3c<<std::endl;

        MultMM<add>(x,m1.mat(),m2.mat(),m3.mat());

        //std::cout<<"diff = "<<m3.mat()-m3c<<std::endl;
        //std::cout<<"Norm(diff) = "<<Norm(m3.mat()-m3c)<<std::endl;
        if (Norm(m3.mat()-m3c) > 1.e-6 * std::abs(T(x))*Norm(m1)*Norm(m2)) {
            std::cout<<"MultMM:  add = "<<add<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            if (m3.colsize() < 100 && m3.rowsize() < 100 && m1.colsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 -> "<<m3.mat()<<std::endl;
                std::cout<<"correct = "<<m3c<<std::endl;
                std::cout<<"diff = "<<(m3c-m3.mat())<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m3c-m3.mat())<<std::endl;
            exit(1);
        }
    }
#endif

#ifdef XDEBUG_PRODMM
    template <class M1, int ix, class T, class M2>
    inline void MultEqMM_Debug(
        BaseMatrix_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M2>& m2)
    {
        Matrix<typename M1::value_type> m1i = m1.mat();
        Matrix<typename M1::value_type> m2i = m2.mat();
        Matrix<typename M1::value_type> m3 = m1.mat();
        m3.setZero();
        for(size_t i=0;i<m3.colsize();++i) {
            for(size_t j=0;j<m3.rowsize();++j) {
                for(size_t k=0;k<m1.rowsize();++k) {
                    m3.ref(i,j) += T(x) * m1.cref(i,k) * m2.cref(k,j);
                }
            }
        }

        MultEqMM(m1.mat(),x,m2.mat());

        if (Norm(m1.mat()-m3) > 1.e-6 * std::abs(T(x))*Norm(m1i)*Norm(m2)) {
            std::cout<<"MultEqMM:  \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            if (m1.colsize() < 100 && m1.rowsize() < 100) {
                std::cout<<"m1 = "<<m1i<<std::endl;
                std::cout<<"m2 = "<<m2i<<std::endl;
                std::cout<<"m1 -> "<<m1.mat()<<std::endl;
                std::cout<<"correct = "<<m3<<std::endl;
                std::cout<<"diff = "<<(m3-m1.mat())<<std::endl;
            }
            std::cout<<"Norm(diff) = "<<Norm(m3-m1.mat())<<std::endl;
            exit(1);
        }
    }
#endif

    template <int ix, class T, class M1, class M2>
    class ProdMM;

    template <int ix, class T, class M1, class M2>
    struct Traits<ProdMM<ix,T,M1,M2> >
    {
        typedef typename M1::value_type mtype1;
        typedef typename M2::value_type mtype2;
        typedef typename Traits2<mtype1,mtype2>::type value_type;

        enum { _colsize = M1::_colsize };
        enum { _rowsize = M2::_rowsize };
        enum { _shape = ShapeTraits2<M1::_shape,M2::_shape>::prod };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { cm2 = Traits<typename M2::calc_type>::_colmajor };
        enum { _rowmajor = 
            (rm1 && rm2) || (!cm1 && !cm2 && _rowsize > int(_colsize)) };

        typedef ProdMM<ix,T,M1,M2> type;
        typedef typename MCopyHelper<value_type,_shape,_colsize,_rowsize,
                _rowmajor,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (M1::_calc && M2::_calc),const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M1, class M2>
    class ProdMM : public BaseMatrix<ProdMM<ix,T,M1,M2> >
    {
    public:

        typedef ProdMM<ix,T,M1,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        inline ProdMM(
            const T& _x, const BaseMatrix<M1>& _m1, const BaseMatrix<M2>& _m2
        ) :
            x(_x), m1(_m1.mat()), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M1::_rowsize,M2::_colsize>::same)); 
            TMVAssert(m1.rowsize() == m2.colsize());
        }

        inline const Scaling<ix,T>& getX() const { return x; }
        inline const M1& getM1() const { return m1; }
        inline const M2& getM2() const { return m2; }

        inline size_t colsize() const { return m1.colsize(); }
        inline size_t rowsize() const { return m2.rowsize(); }

        inline value_type cref(int i, int j) const
        { return x * (m1.get_row(i) * m2.get_col(j)); }

        template <class M3>
        inline void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M3::_shape>::assignable)); 
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_PRODMM
            MultMM_Debug<false>(x,m1.calc(),m2.calc(),m3.mat());
#else
            MultMM<false>(x,m1.calc(),m2.calc(),m3.mat());
#endif
        }

        template <class M3>
        inline void newAssignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M3::_shape>::assignable)); 
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
#ifdef XDEBUG_PRODMM
            MultMM_Debug<false>(x,m1.calc(),m2.calc(),m3.mat());
#else
            NoAliasMultMM<false>(x,m1.calc(),m2.calc(),m3.mat());
#endif
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
    { return ProdMM<ix,T,M1,M2>(m2.getX(),m1,m2.getM()); }

    // xm * m
    template <int ix, class T, class M1, class M2>
    inline ProdMM<ix,T,M1,M2> operator*(
        const ProdXM<ix,T,M1>& m1, const BaseMatrix<M2>& m2)
    { return ProdMM<ix,T,M1,M2>(m1.getX(),m1.getM(),m2); }

    // xm * xm
#define PT typename Traits2<T1,T2>::type
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    inline ProdMM<ix1*ix2,PT,M1,M2> operator*(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        return ProdMM<ix1*ix2,PT,M1,M2>(
            m1.getX()*m2.getX(),m1.getM(),m2.getM()); 
    }
#undef PT


    // m *= m
#define RT typename M1::real_type
    template <class M1, class M2>
    inline void MultEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
#ifdef XDEBUG_PRODMM
        MultEqMM_Debug(m1.mat(),Scaling<1,RT>(),m2.calc()); 
#else
        MultEqMM(m1.mat(),Scaling<1,RT>(),m2.calc()); 
#endif
    }
#undef RT

    // m *= xm
    template <class M1, int ix2, class T2, class M2>
    inline void MultEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { 
#ifdef XDEBUG_PRODMM
        MultEqMM_Debug(m1.mat(),m2.getX(),m2.getM().calc()); 
#else
        MultEqMM(m1.mat(),m2.getX(),m2.getM().calc()); 
#endif
    }

    // m += mm
    template <class M3, int ix, class T, class M1, class M2>
    inline void AddEq(
        BaseMatrix_Mutable<M3>& m3, const ProdMM<ix,T,M1,M2>& mm)
    { 
#ifdef XDEBUG_PRODMM
        MultMM_Debug<true>(
            mm.getX(),mm.getM1().calc(),mm.getM2().calc(),m3.mat()); 
#else
        MultMM<true>(mm.getX(),mm.getM1().calc(),mm.getM2().calc(),m3.mat()); 
#endif
    }

    // m -= mm
    template <class M3, int ix, class T, class M1, class M2>
    inline void SubtractEq(
        BaseMatrix_Mutable<M3>& m3, const ProdMM<ix,T,M1,M2>& mm)
    { 
#ifdef XDEBUG_PRODMM
        MultMM_Debug<true>(
            -mm.getX(),mm.getM1().calc(),mm.getM2().calc(),m3.mat()); 
#else
        MultMM<true>(-mm.getX(),mm.getM1().calc(),mm.getM2().calc(),m3.mat()); 
#endif
    }


    // Consolidate x*(xmm) type constructs:

#define RT typename ProdMM<ix,T,M1,M2>::real_type
#define CT typename ProdMM<ix,T,M1,M2>::complex_type
#define CCT ConjRef<CT>

    // -(mm)
    template <int ix, class T, class M1, class M2>
    inline ProdMM<-ix,T,M1,M2> operator-(const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<-ix,T,M1,M2>(-mm.getX(),mm.getM1(),mm.getM2()); }

    // x * (mm)
    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,T,M1,M2> operator*(
        const int x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,T,M1,M2>(RT(x)*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,T,M1,M2> operator*(
        const RT x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,T,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,CT,M1,M2> operator*(
        const CT x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,CT,M1,M2> operator*(
        const CCT x, const ProdMM<ix,T,M1,M2>& mm)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    inline ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
        const Scaling<ix1,T1>& x, const ProdMM<ix,T,M1,M2>& mm)
    { 
        return ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            T1(x)*mm.getX(),mm.getM1(),mm.getM2()); 
    }

    // (mm)*x
    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,T,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const int x)
    { return ProdMM<0,T,M1,M2>(RT(x)*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,T,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const RT x)
    { return ProdMM<0,T,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,CT,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const CT x)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,CT,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const CCT x)
    { return ProdMM<0,CT,M1,M2>(x*mm.getX(),mm.getM1(),mm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    inline ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator*(
        const ProdMM<ix,T,M1,M2>& mm, const Scaling<ix1,T1>& x)
    { 
        return ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            T1(x)*mm.getX(),mm.getM1(),mm.getM2()); 
    }

    // (mm)/x
    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,T,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const int x)
    { return ProdMM<0,T,M1,M2>(mm.getX()/RT(x),mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,T,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const RT x)
    { return ProdMM<0,T,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,CT,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const CT x)
    { return ProdMM<0,CT,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix, class T, class M1, class M2>
    inline ProdMM<0,CT,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const CCT x)
    { return ProdMM<0,CT,M1,M2>(mm.getX()/x,mm.getM1(),mm.getM2()); }

    template <int ix1, class T1, int ix, class T, class M1, class M2>
    inline ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2> operator/(
        const ProdMM<ix,T,M1,M2>& mm, const Scaling<ix1,T1>& x)
    { 
        return ProdMM<ix1*ix,typename Traits2<T1,T>::type,M1,M2>(
            mm.getX()/T1(x),mm.getM1(),mm.getM2()); 
    }

#undef RT
#undef CT
#undef CCT

    // TMV_Text

    template <int ix, class T, class M1, class M2>
    inline std::string TMV_Text(const ProdMM<ix,T,M1,M2>& mm)
    {
        std::ostringstream s;
        s << "ProdMM< "<<ix<<","<<TMV_Text(T())<<",";
        s << TMV_Text(mm.getM1())<<" , "<<TMV_Text(mm.getM2())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
