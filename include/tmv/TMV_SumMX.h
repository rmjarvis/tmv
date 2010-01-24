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


#ifndef TMV_SumMX_H
#define TMV_SumMX_H

#include "TMV_ProdXM.h"

namespace tmv {

    //
    // Matrix + Scalar
    //

    template <class T, class M>
    inline void AddMX(const T& x, BaseMatrix_Mutable<M>& m)
    {
        TMVStaticAssert((Sizes<M::mcolsize,M::mrowsize>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        m.mat().diag().addToAll(x); 
    }

    template <int ix1, class T1, class M1, class T2, class M3>
    inline void AddMX(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const T2& x2, BaseMatrix_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same)); 
        TMVStaticAssert((Sizes<M1::mrowsize,M3::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        MultXM<false>(x1,m1.mat(),m3.mat());
        AddMX(x2,m3.mat());
    }

    template <int ix1, class T1, class M1, class T2, class M3>
    inline void NoAliasAddMX(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
        const T2& x2, BaseMatrix_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same)); 
        TMVStaticAssert((Sizes<M1::mrowsize,M3::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        NoAliasMultXM<false>(x1,m1.mat(),m3.mat());
        AddMX(x2,m3.mat());
    }

    template <int ix1, class T1, class M1, class T2>
    class SumMX;

    template <int ix1, class T1, class M1, class T2>
    struct Traits<SumMX<ix1,T1,M1,T2> >
    {
        typedef typename ProdXM<ix1,T1,M1>::value_type mtype1;
        typedef typename Traits2<mtype1,T2>::type value_type;

        enum { mcolsize = M1::mcolsize };
        enum { mrowsize = M1::mrowsize };
        enum { mshape = M1::mshape };
        enum { mfort = M1::mfort };
        enum { mcalc = false };
        enum { mrowmajor = M1::mrowmajor };
        enum { mcolmajor = !mrowmajor };

        typedef SumMX<ix1,T1,M1,T2> type;
        typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,
                mrowmajor,mfort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            M1::mcalc,const type,calc_type>::type eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix1, class T1, class M1, class T2>
    class SumMX : public BaseMatrix<SumMX<ix1,T1,M1,T2> >
    {
    public:

        typedef SumMX<ix1,T1,M1,T2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        inline SumMX(const T1& _x1, const BaseMatrix<M1>& _m1, const T2& _x2) :
            x1(_x1), m1(_m1.mat()), x2(_x2)
        {
            TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
            TMVAssert(m1.colsize() == m1.rowsize());
        }

        inline const Scaling<ix1,T1>& getX1() const { return x1; }
        inline const M1& getM1() const { return m1; }
        inline const T2& getX2() const { return x2; }

        inline size_t colsize() const { return m1.colsize(); }
        inline size_t rowsize() const { return m1.rowsize(); }
        inline size_t ls() const { return m1.ls(); }

        inline value_type cref(int i, int j) const
        { return (i == j) ? (x1 * m1.cref(i,j) + x2) : (x1 * m1.cref(i,j)); }

        template <class M3>
        inline void assignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::misreal || M3::miscomplex));
            TMVStaticAssert((Sizes<type::mcolsize,M3::mcolsize>::same)); 
            TMVStaticAssert((Sizes<type::mrowsize,M3::mrowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
            AddMX(x1,m1.calc(),x2,m3.mat());
        }

        template <class M3>
        inline void newAssignTo(BaseMatrix_Mutable<M3>& m3) const
        {
            TMVStaticAssert((type::misreal || M3::miscomplex));
            TMVStaticAssert((Sizes<type::mcolsize,M3::mcolsize>::same)); 
            TMVStaticAssert((Sizes<type::mrowsize,M3::mrowsize>::same)); 
            TMVAssert(colsize() == m3.colsize());
            TMVAssert(rowsize() == m3.rowsize());
            NoAliasAddMX(x1,m1.calc(),x2,m3.mat());
        }

    private:
        const Scaling<ix1,T1> x1;
        const M1& m1;
        const T2 x2;
    };

#define RT typename M1::real_type
#define CT typename M1::complex_type
#define CCT ConjRef<CT>
    // m += x
    template <class M1>
    inline void AddEq(BaseMatrix_Mutable<M1>& m1, const int x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        AddMX(RT(x2),m1.mat());
    }

    template <class M1>
    inline void AddEq(BaseMatrix_Mutable<M1>& m1, const RT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        AddMX(x2,m1.mat());
    }

    template <class M1>
    inline void AddEq(BaseMatrix_Mutable<M1>& m1, const CT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert(M1::miscomplex);
        AddMX(x2,m1.mat());
    }

    template <class M1>
    inline void AddEq(BaseMatrix_Mutable<M1>& m1, const CCT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert(M1::miscomplex);
        AddMX(CT(x2),m1.mat());
    }

    template <class M1, int ix2, class T2>
    inline void AddEq(BaseMatrix_Mutable<M1>& m1, const Scaling<ix2,T2>& x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert(M1::miscomplex || T2::misreal);
        AddMX(T2(x2),m1.mat());
    }

    // m -= x
    template <class M1>
    inline void SubtractEq(BaseMatrix_Mutable<M1>& m1, const int x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        AddMX(-RT(x2),m1.mat());
    }

    template <class M1>
    inline void SubtractEq(BaseMatrix_Mutable<M1>& m1, const RT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        AddMX(-x2,m1.mat());
    }

    template <class M1>
    inline void SubtractEq(BaseMatrix_Mutable<M1>& m1, const CT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert(M1::miscomplex);
        AddMX(-x2,m1.mat());
    }

    template <class M1>
    inline void SubtractEq(BaseMatrix_Mutable<M1>& m1, const CCT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert(M1::miscomplex);
        AddMX(-CT(x2),m1.mat());
    }

    template <class M1, int ix2, class T2>
    inline void SubtractEq(
        BaseMatrix_Mutable<M1>& m1, const Scaling<ix2,T2>& x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        TMVStaticAssert(M1::miscomplex || T2::misreal);
        AddMX(-T2(x2),m1.mat());
    }


    // m + x
    template <class M1>
    inline SumMX<1,RT,M1,RT> operator+(const BaseMatrix<M1>& m1, const int x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,RT>(RT(1),m1,RT(x2)); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,RT> operator+(const BaseMatrix<M1>& m1, const RT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,RT>(RT(1),m1,x2); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,CT> operator+(const BaseMatrix<M1>& m1, const CT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,CT>(RT(1),m1,x2); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,CT> operator+(const BaseMatrix<M1>& m1, const CCT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,CT>(RT(1),m1,CT(x2)); 
    }

    template <class M1, int ix2, class T2>
    inline SumMX<1,RT,M1,CT> operator+(
        const BaseMatrix<M1>& m1, const Scaling<ix2,T2> x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,T2>(RT(1),m1,T2(x2)); 
    }

    // x + m
    template <class M1>
    inline SumMX<1,RT,M1,RT> operator+(const int x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,RT>(RT(1),m1,RT(x2)); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,RT> operator+(const RT x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,RT>(RT(1),m1,x2); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,CT> operator+(const CT x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,CT>(RT(1),m1,x2); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,CT> operator+(const CCT x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,CT>(RT(1),m1,CT(x2)); 
    }

    template <class M1, int ix2, class T2>
    inline SumMX<1,RT,M1,CT> operator+(
        const Scaling<ix2,T2> x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,T2>(RT(1),m1,T2(x2)); 
    }

    // m - x
    template <class M1>
    inline SumMX<1,RT,M1,RT> operator-(const BaseMatrix<M1>& m1, const int x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,RT>(RT(1),m1,-RT(x2)); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,RT> operator-(const BaseMatrix<M1>& m1, const RT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,RT>(RT(1),m1,-x2); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,CT> operator-(const BaseMatrix<M1>& m1, const CT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,CT>(RT(1),m1,-x2); 
    }

    template <class M1>
    inline SumMX<1,RT,M1,CT> operator-(const BaseMatrix<M1>& m1, const CCT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,CT>(RT(1),m1,-CT(x2)); 
    }

    template <class M1, int ix2, class T2>
    inline SumMX<1,RT,M1,CT> operator-(
        const BaseMatrix<M1>& m1, const Scaling<ix2,T2> x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<1,RT,M1,T2>(RT(1),m1,-T2(x2)); 
    }

    // x - m
    template <class M1>
    inline SumMX<-1,RT,M1,RT> operator-(const int x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-1,RT,M1,RT>(RT(-1),m1,RT(x2)); 
    }

    template <class M1>
    inline SumMX<-1,RT,M1,RT> operator-(const RT x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-1,RT,M1,RT>(RT(-1),m1,x2); 
    }

    template <class M1>
    inline SumMX<-1,RT,M1,CT> operator-(const CT x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-1,RT,M1,CT>(RT(-1),m1,x2); 
    }

    template <class M1>
    inline SumMX<-1,RT,M1,CT> operator-(const CCT x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-1,RT,M1,CT>(RT(-1),m1,CT(x2)); 
    }

    template <class M1, int ix2, class T2>
    inline SumMX<-1,RT,M1,CT> operator-(
        const Scaling<ix2,T2> x2, const BaseMatrix<M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-1,RT,M1,T2>(RT(-1),m1,T2(x2)); 
    }

    // xm + x
    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,RT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const int x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),RT(x2)); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,RT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const RT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,CT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const CT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,CT> operator+(
        const ProdXM<ix1,T1,M1>& m1, const CCT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),CT(x2)); 
    }

    template <int ix1, class T1, class M1, int ix2, class T2>
    inline SumMX<ix1,T1,M1,T2> operator+(
        const ProdXM<ix1,T1,M1>& m1, const Scaling<ix2,T2> x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,T2>(T1(m1.getX()),m1.getM(),T2(x2)); 
    }

    // x + xm 
    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,RT> operator+(
        const int x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),int(x2)); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,RT> operator+(
        const RT x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,CT> operator+(
        const CT x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,CT> operator+(
        const CCT x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),CT(x2)); 
    }

    template <int ix1, class T1, class M1, int ix2, class T2>
    inline SumMX<ix1,T1,M1,T2> operator+(
        const Scaling<ix2,T2> x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,T2>(T1(m1.getX()),m1.getM(),T2(x2)); 
    }

    // xm - x
    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,RT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const int x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),-RT(x2)); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,RT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const RT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,RT>(T1(m1.getX()),m1.getM(),-x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,CT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const CT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),-x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<ix1,T1,M1,CT> operator-(
        const ProdXM<ix1,T1,M1>& m1, const CCT x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,CT>(T1(m1.getX()),m1.getM(),-CT(x2)); 
    }

    template <int ix1, class T1, class M1, int ix2, class T2>
    inline SumMX<ix1,T1,M1,T2> operator-(
        const ProdXM<ix1,T1,M1>& m1, const Scaling<ix2,T2> x2)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<ix1,T1,M1,T2>(T1(m1.getX()),m1.getM(),-T2(x2)); 
    }

    // x - xm 
    template <int ix1, class T1, class M1>
    inline SumMX<-ix1,T1,M1,RT> operator-(
        const int x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-ix1,T1,M1,RT>(-T1(m1.getX()),m1.getM(),RT(x2)); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<-ix1,T1,M1,RT> operator-(
        const RT x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-ix1,T1,M1,RT>(-T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<-ix1,T1,M1,CT> operator-(
        const CT x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-ix1,T1,M1,CT>(-T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1>
    inline SumMX<-ix1,T1,M1,CT> operator-(
        const CCT x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-ix1,T1,M1,CT>(-T1(m1.getX()),m1.getM(),x2); 
    }

    template <int ix1, class T1, class M1, int ix2, class T2>
    inline SumMX<-ix1,T1,M1,T2> operator-(
        const Scaling<ix2,T2> x2, const ProdXM<ix1,T1,M1>& m1)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same)); 
        TMVAssert(m1.colsize() == m1.rowsize());
        return SumMX<-ix1,T1,M1,T2>(-T1(m1.getX()),m1.getM(),T2(x2)); 
    }

#undef RT
#undef CT
#undef CCT


    // Consolidate x*(xm+x) type constructs:

#define RT typename SumMX<ix1,T1,M1,T2>::real_type
#define CT typename SumMX<ix1,T1,M1,T2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

    // -(xm+x)
    template <int ix1, class T1, class M1, class T2>
    inline SumMX<-ix1,T1,M1,T2> operator-(
        const SumMX<ix1,T1,M1,T2>& smx)
    { return SumMX<-ix1,T1,M1,T2>(-T1(smx.getX1()),smx.getM1(),-smx.getX2()); }

    // x * (xm+x)
    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,T1,M1,T2> operator*(
        const int x, const SumMX<ix1,T1,M1,T2>& smx)
    {
        return SumMX<0,T1,M1,T2>(
            RT(x)*smx.getX1(),smx.getM1(),RT(x)*smx.getX2());
    }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,T1,M1,T2> operator*(
        const RT x, const SumMX<ix1,T1,M1,T2>& smx)
    { return SumMX<0,T1,M1,T2>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,CT,M1,CT> operator*(
        const CT x, const SumMX<ix1,T1,M1,T2>& smx)
    { return SumMX<0,CT,M1,CT>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,CT,M1,CT> operator*(
        const CCT x, const SumMX<ix1,T1,M1,T2>& smx)
    {
        return SumMX<0,CT,M1,CT>(
            CT(x)*smx.getX1(),smx.getM1(),CT(x)*smx.getX2());
    }
    template <int ix, class T, int ix1, class T1, class M1, class T2>
    inline SumMX<ix1*ix,TX1,M1,TX2> operator*(
        const Scaling<ix,T>& x, const SumMX<ix1,T1,M1,T2>& smx)
    {
        return SumMX<ix1*ix,TX1,M1,TX2>(
            T(x)*smx.getX1(),smx.getM1(),x*smx.getX2());
    }

    // (xm+x)*x
    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,T1,M1,T2> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const int x)
    {
        return SumMX<0,T1,M1,T2>(
            RT(x)*smx.getX1(),smx.getM1(),RT(x)*smx.getX2());
    }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,T1,M1,T2> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const RT x)
    { return SumMX<0,T1,M1,T2>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,CT,M1,CT> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const CT x)
    { return SumMX<0,CT,M1,CT>(x*smx.getX1(),smx.getM1(),x*smx.getX2()); }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,CT,M1,CT> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const CCT x)
    {
        return SumMX<0,CT,M1,CT>(
            CT(x)*smx.getX1(),smx.getM1(),CT(x)*smx.getX2());
    }

    template <int ix, class T, int ix1, class T1, class M1, class T2>
    inline SumMX<ix1*ix,TX1,M1,TX2> operator*(
        const SumMX<ix1,T1,M1,T2>& smx, const Scaling<ix,T>& x)
    {
        return SumMX<ix1*ix,TX1,M1,TX2>(
            T(x)*smx.getX1(),smx.getM1(),x*smx.getX2());
    }

    // (xm+x)/x
    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,T1,M1,T2> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const int x)
    {
        return SumMX<0,T1,M1,T2>(
            smx.getX1()/RT(x),smx.getM1(),smx.getX2()/RT(x));
    }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,T1,M1,T2> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const RT x)
    {
        return SumMX<0,T1,M1,T2>(smx.getX1()/x,smx.getM1(),smx.getX2()/x);
    }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,CT,M1,CT> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const CT x)
    { return SumMX<0,CT,M1,CT>(smx.getX1()/x,smx.getM1(), smx.getX2()/x); }

    template <int ix1, class T1, class M1, class T2>
    inline SumMX<0,CT,M1,CT> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const CCT x)
    {
        return SumMX<0,CT,M1,CT>(
            smx.getX1()/CT(x),smx.getM1(),smx.getX2()/CT(x));
    }

    template <int ix, class T, int ix1, class T1, class M1, class T2>
    inline SumMX<ix1*ix,TX1,M1,TX2> operator/(
        const SumMX<ix1,T1,M1,T2>& smx, const Scaling<ix,T>& x)
    {
        return SumMX<ix1*ix,TX1,M1,TX2>(
            smx.getX1()/T(x),smx.getM1(),smx.getX2()/T(x));
    }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



    // TMV_Text

    template <int ix1, class T1, class M1, class T2>
    inline std::string TMV_Text(const SumMX<ix1,T1,M1,T2>& smx)
    {
        std::ostringstream s;
        s << "SumMX< "<< ix1<<","<<TMV_Text(T1())
            <<" , "<<TMV_Text(smx.getM1())
            <<" , "<<TMV_Text(T2())<<" >";
        return s.str();
    }



} // namespace tmv

#endif 
