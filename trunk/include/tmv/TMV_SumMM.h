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


#ifndef TMV_SumMM_H
#define TMV_SumMM_H

#include "TMV_BaseMatrix.h"
#include "TMV_ProdXM.h"
#include "TMV_AddMM_Funcs.h"

namespace tmv {

#if 1
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix<M2>& m2,
        BaseMatrix_Mutable<M3>& m3)
    { AddMM(x1,m1.calc(),x2,m2.calc(),m3.mat()); }
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix<M2>& m2,
        BaseMatrix_Mutable<M3>& m3)
    { NoAliasAddMM(x1,m1.calc(),x2,m2.calc(),m3.mat()); }
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    static inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix<M1>& m1,
        const Scaling<ix2,T2>& x2, const BaseMatrix<M2>& m2,
        BaseMatrix_Mutable<M3>& m3)
    { AliasAddMM(x1,m1.calc(),x2,m2.calc(),m3.mat()); }
#endif

    //
    // Matrix + Matrix
    //

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    class SumMM;

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    struct Traits<SumMM<ix1,T1,M1,ix2,T2,M2> >
    {
        typedef typename ProdXM<ix1,T1,M1>::value_type mtype1;
        typedef typename ProdXM<ix2,T2,M2>::value_type mtype2;
        typedef typename Traits2<mtype1,mtype2>::type value_type;

        enum { _colsize = Sizes<M1::_colsize,M2::_colsize>::size };
        enum { _rowsize = Sizes<M1::_rowsize,M2::_rowsize>::size };
        enum { _shape = ShapeTraits2<M1::_shape,M2::_shape>::sum };
        enum { _fort = M1::_fort && M2::_fort };
        enum { _calc = false };
        enum { rm1 = Traits<typename M1::calc_type>::_rowmajor };
        enum { rm2 = Traits<typename M2::calc_type>::_rowmajor };
        enum { cm1 = Traits<typename M1::calc_type>::_colmajor };
        enum { _rowmajor = ( rm1 || (rm2 && !cm1) ) };

        typedef SumMM<ix1,T1,M1,ix2,T2,M2> type;
        typedef typename MCopyHelper<value_type,_shape,_colsize,_rowsize,
                _rowmajor,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (M1::_calc && M2::_calc),const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    class SumMM : 
        public BaseMatrix<SumMM<ix1,T1,M1,ix2,T2,M2> >
    {
    public:

        typedef SumMM<ix1,T1,M1,ix2,T2,M2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        SumMM(
            const T1& _x1, const BaseMatrix<M1>& _m1, 
            const T2& _x2, const BaseMatrix<M2>& _m2) :
            x1(_x1), m1(_m1.mat()), x2(_x2), m2(_m2.mat())
        {
            TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
            TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
            TMVAssert(m1.colsize() == m2.colsize());
            TMVAssert(m1.rowsize() == m2.rowsize());
        }

        const Scaling<ix1,T1>& getX1() const { return x1; }
        const M1& getM1() const { return m1; }
        const Scaling<ix2,T2>& getX2() const { return x2; }
        const M2& getM2() const { return m2; }

        size_t colsize() const { return m1.colsize(); }
        size_t rowsize() const { return m1.rowsize(); }

        value_type cref(int i, int j) const
        { return x1 * m1.cref(i,j) + x2 * m2.cref(i,j); }

        template <class M3>
        void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M3::_shape>::assignable)); 
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
            AddMM(x1,m1.mat(),x2,m2.mat(),m3.mat());
        }

        template <class M3>
        void newAssignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::_shape,M3::_shape>::assignable)); 
            TMVStaticAssert((type::isreal || M3::iscomplex));
            TMVStaticAssert((Sizes<type::_colsize,M3::_colsize>::same)); 
            TMVStaticAssert((Sizes<type::_rowsize,M3::_rowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
            NoAliasAddMM(x1,m1.mat(),x2,m2.mat(),m3.mat());
        }

    private:
        const Scaling<ix1,T1> x1;
        const M1& m1;
        const Scaling<ix2,T2> x2;
        const M2& m2;
    };

#define RT typename M2::real_type
    // m += m
    template <class M1, class M2>
    static inline void AddEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2) 
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVStaticAssert(M1::iscomplex || M2::isreal);
        MultXM<true>(m2.mat(),m1.mat());
    }

    // m += xm
    template <class M1, int ix2, class T2, class M2>
    static inline void AddEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2) 
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVStaticAssert(M1::iscomplex || M2::isreal);
        MultXM<true>(m2.getX(),m2.getM().mat(),m1.mat());
    }

    // m -= m
    template <class M1, class M2>
    static inline void SubtractEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVStaticAssert(M1::iscomplex || M2::isreal);
        MultXM<true>(Scaling<-1,RT>(),m2.mat(),m1.mat());
    }

    // m -= xm
    template <class M1, int ix2, class T2, class M2>
    static inline void SubtractEq(
        BaseMatrix_Mutable<M1>& m1, const ProdXM<ix2,T2,M2>& m2) 
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVStaticAssert(M1::iscomplex || M2::isreal);
        MultXM<true>(-m2.getX(),m2.getM().mat(),m1.mat());
    }

    // m + m
    template <class M1, class M2>
    static inline SumMM<1,RT,M1,1,RT,M2> operator+(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<1,RT,M1,1,RT,M2>(RT(1),m1,RT(1),m2); 
    }

    // xm + m
    template <int ix1, class T1, class M1, class M2>
    static inline SumMM<ix1,T1,M1,1,RT,M2> operator+(
        const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<ix1,T1,M1,1,RT,M2>(T1(m1.getX()),m1.getM(),RT(1),m2); 
    }

    // m + xm
    template <class M1, int ix2, class T2, class M2>
    static inline SumMM<1,RT,M1,ix2,T2,M2> operator+(
        const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<1,RT,M1,ix2,T2,M2>(RT(1),m1,T2(m2.getX()),m2.getM()); 
    }

    // xm + xm
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<ix1,T1,M1,ix2,T2,M2> operator+(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<ix1,T1,M1,ix2,T2,M2>(
            T1(m1.getX()),m1.getM(),T2(m2.getX()),m2.getM()); 
    }

    // m - m
    template <class M1, class M2>
    static inline SumMM<1,RT,M1,-1,RT,M2> operator-(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<1,RT,M1,-1,RT,M2>(RT(1),m1,RT(-1),m2); 
    }

    // xm - m
    template <int ix1, class T1, class M1, class M2>
    static inline SumMM<ix1,T1,M1,-1,RT,M2> operator-(
        const ProdXM<ix1,T1,M1>& m1, const BaseMatrix<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<ix1,T1,M1,-1,RT,M2>(T1(m1.getX()),m1.getM(),RT(-1),m2); 
    }

    // m - xm
    template <class M1, int ix2, class T2, class M2>
    static inline SumMM<1,RT,M1,-ix2,T2,M2> operator-(
        const BaseMatrix<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<1,RT,M1,-ix2,T2,M2>(RT(1),m1,-T2(m2.getX()),m2.getM()); 
    }

    // xm - xm
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<ix1,T1,M1,-ix2,T2,M2> operator-(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        return SumMM<ix1,T1,M1,-ix2,T2,M2>(
            T1(m1.getX()),m1.getM(),-T2(m2.getX()),m2.getM()); 
    }
#undef RT


    // Consolidate x*(xm+xm) type constructs:

#define RT typename SumMM<ix1,T1,M1,ix2,T2,M2>::real_type
#define CT typename SumMM<ix1,T1,M1,ix2,T2,M2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

    // -(xm+xm)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<-ix1,T1,M1,-ix2,T2,M2> operator-(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<-ix1,T1,M1,-ix2,T2,M2>(
            -T1(smm.getX1()),smm.getM1(),-T2(smm.getX2()),smm.getM2());
    }

    // x * (xm+xm)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,T1,M1,0,T2,M2> operator*(
        const int x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            RT(x)*smm.getX1(),smm.getM1(), RT(x)*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,T1,M1,0,T2,M2> operator*(
        const RT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,CT,M1,0,CT,M2> operator*(
        const CT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,CT,M1,0,CT,M2> operator*(
        const CCT x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            CT(x)*smm.getX1(),smm.getM1(), CT(x)*smm.getX2(),smm.getM2());
    }
    template <int ix, class T, int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2> operator*(
        const Scaling<ix,T>& x, const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        return SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2>(
            T(x)*smm.getX1(),smm.getM1(),T(x)*smm.getX2(),smm.getM2());
    }

    // (xm+xm)*x
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,T1,M1,0,T2,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const int x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            RT(x)*smm.getX1(),smm.getM1(), RT(x)*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,T1,M1,0,T2,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const RT x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,CT,M1,0,CT,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            x*smm.getX1(),smm.getM1(), x*smm.getX2(),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,CT,M1,0,CT,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CCT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            CT(x)*smm.getX1(),smm.getM1(), CT(x)*smm.getX2(),smm.getM2());
    }

    template <int ix, class T, int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2> operator*(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const Scaling<ix,T>& x)
    {
        return SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2>(
            T(x)*smm.getX1(),smm.getM1(),T(x)*smm.getX2(),smm.getM2());
    }

    // (xm+xm)/x
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,T1,M1,0,T2,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const int x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            smm.getX1()/RT(x),smm.getM1(), smm.getX2()/RT(x),smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,T1,M1,0,T2,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const RT x)
    {
        return SumMM<0,T1,M1,0,T2,M2>(
            smm.getX1()/x,smm.getM1(), smm.getX2()/x,smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,CT,M1,0,CT,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            smm.getX1()/x,smm.getM1(), smm.getX2()/x,smm.getM2());
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<0,CT,M1,0,CT,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const CCT x)
    {
        return SumMM<0,CT,M1,0,CT,M2>(
            smm.getX1()/CT(x),smm.getM1(), smm.getX2()/CT(x),smm.getM2());
    }

    template <int ix, class T, int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2> operator/(
        const SumMM<ix1,T1,M1,ix2,T2,M2>& smm, const Scaling<ix,T>& x)
    {
        return SumMM<ix1*ix,TX1,M1,ix2*ix,TX2,M2>(
            smm.getX1()/T(x),smm.getM1(),smm.getX2()/T(x),smm.getM2());
    }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



    // TMV_Text

#ifdef TMV_DEBUG
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    static inline std::string TMV_Text(const SumMM<ix1,T1,M1,ix2,T2,M2>& smm)
    {
        std::ostringstream s;
        s << "SumMM< "<< ix1<<","<<TMV_Text(T1())
            <<" , "<<TMV_Text(smm.getM1())
            <<" , "<<ix2<<","<<TMV_Text(T2())
            <<" , "<<TMV_Text(smm.getM2())<<" >";
        return s.str();
    }
#endif



} // namespace tmv

#endif 
