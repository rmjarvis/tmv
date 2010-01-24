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

#ifndef TMV_ProdXM_H
#define TMV_ProdXM_H

#include "TMV_ProdXV.h"
#include "TMV_BaseMatrix.h"

namespace tmv {
    
    //
    // Scalar * Matrix
    //

    // These first few are intentionally not defined to make sure we
    // get a compiler error if they are used.
    // All real calls should go through a more specific version than 
    // just the BaseMatrix_Calc's.
    template <bool add, int ix, class T, class M1, class M2>
    inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        BaseMatrix_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        BaseMatrix_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    inline void InlineMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        BaseMatrix_Mutable<M2>& m2);
    template <bool add, int ix, class T, class M1, class M2>
    inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1, 
        BaseMatrix_Mutable<M2>& m2);

    // These are helpers to allow the caller to not use a Scaling object.
    template <bool add, class T, class M1, class M2>
    inline void MultXM(
        const T& x, const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MultXM<add>(Scaling<0,T>(x),m1.mat(),m2.mat()); }
    template <bool add, class T, class M1, class M2>
    inline void NoAliasMultXM(
        const T& x, const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { NoAliasMultXM<add>(Scaling<0,T>(x),m1.mat(),m2.mat()); }
    template <bool add, class T, class M1, class M2>
    inline void InlineMultXM(
        const T& x, const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { InlineMultXM<add>(Scaling<0,T>(x),m1.mat(),m2.mat()); }
    template <bool add, class T, class M1, class M2>
    inline void AliasMultXM(
        const T& x, const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { AliasMultXM<add>(Scaling<0,T>(x),m1.mat(),m2.mat()); }

    template <bool add, class M1, class M2>
    inline void MultXM(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { MultXM<add>(Scaling<1,typename M2::real_type>(),m1.mat(),m2.mat()); }
    template <bool add, class M1, class M2>
    inline void NoAliasMultXM(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        NoAliasMultXM<add>(
            Scaling<1,typename M2::real_type>(),m1.mat(),m2.mat()); 
    }
    template <bool add, class M1, class M2>
    inline void InlineMultXM(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    {
        InlineMultXM<add>(
            Scaling<1,typename M2::real_type>(),m1.mat(),m2.mat()); 
    }
    template <bool add, class M1, class M2>
    inline void AliasMultXM(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2)
    { 
        AliasMultXM<add>(
            Scaling<1,typename M2::real_type>(),m1.mat(),m2.mat()); 
    }

    template <class T, class M>
    inline void Scale(const T& x, BaseMatrix_Mutable<M>& m)
    { Scale(Scaling<0,T>(x),m.mat()); }

    template <class T, class M>
    inline void InlineScale(const T& x, BaseMatrix_Mutable<M>& m)
    { InlineScale(Scaling<0,T>(x),m.mat()); }


    template <int ix, class T, class M>
    class ProdXM;

    template <int ix, class T, class M>
    struct Traits<ProdXM<ix,T,M> >
    {
        typedef typename Traits2<T,typename M::value_type>::type value_type;

        enum { mcolsize = M::mcolsize };
        enum { mrowsize = M::mrowsize };
        enum { mshape = ShapeTraits<M::mshape>::nonunit_shape };
        enum { mfort = M::mfort };
        enum { mcalc = false };
        enum { rm1 = Traits<typename M::calc_type>::mrowmajor };
        enum { mrowmajor = rm1 };

        typedef ProdXM<ix,T,M> type;
        typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,
                mrowmajor,mfort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<M::mcalc,const type,calc_type>::type 
            eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M>
    class ProdXM : public BaseMatrix<ProdXM<ix,T,M> >
    {
    public:

        typedef ProdXM<ix,T,M> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        inline ProdXM(const T _x, const BaseMatrix<M>& _m) : 
            x(_x), m(_m.mat()) {}
        inline const Scaling<ix,T>& getX() const { return x; }
        inline const M& getM() const { return m; }

        inline size_t colsize() const { return m.colsize(); }
        inline size_t rowsize() const { return m.rowsize(); }

        inline value_type cref(int i, int j) const
        { return x * m.cref(i,j); }

        template <class M2>
        inline void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::mshape,M2::mshape>::assignable)); 
            TMVStaticAssert((Sizes<type::mcolsize,M2::mcolsize>::same));
            TMVStaticAssert((Sizes<type::mrowsize,M2::mrowsize>::same));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVStaticAssert(type::misreal || M2::miscomplex);
            MultXM<false>(x,m.calc(),m2.mat());
        }

        template <class M2>
        inline void newAssignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((
                    ShapeTraits2<type::mshape,M2::mshape>::assignable)); 
            TMVStaticAssert((Sizes<type::mcolsize,M2::mcolsize>::same));
            TMVStaticAssert((Sizes<type::mrowsize,M2::mrowsize>::same));
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVStaticAssert(type::misreal || M2::miscomplex);
            NoAliasMultXM<false>(x,m.calc(),m2.mat());
        }

    private:
        const Scaling<ix,T> x;
        const M& m;
    };

#define RT typename M::real_type
#define CT typename M::complex_type
#define CCT ConjRef<CT>

    // m *= x
    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const int x)
    { Scale(RT(x),m.mat()); }

    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const RT x)
    { Scale(x,m.mat()); }

    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const CT x)
    { Scale(x,m.mat()); }

    template <class M>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const CCT x)
    { Scale(CT(x),m.mat()); }

    template <class M, class T>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const Scaling<0,T>& x)
    { Scale(x,m.mat()); }

    template <class M, class T>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const Scaling<1,T>& x)
    {}

    template <class M, class T>
    inline void MultEq(BaseMatrix_Mutable<M>& m, const Scaling<-1,T>& x)
    { Scale(x,m.mat()); }

    // m /= x
    template <class M>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const int x)
    { Scale(RT(1)/RT(x),m.mat()); }

    template <class M>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const RT x)
    { Scale(RT(1)/x,m.mat()); }

    template <class M>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const CT x)
    { Scale(RT(1)/x,m.mat()); }

    template <class M>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const CCT x)
    { Scale(RT(1)/CT(x),m.mat()); }

    template <class M, class T>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const Scaling<0,T>& x)
    { Scale(RT(1)/T(x),m.mat()); }

    template <class M, class T>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const Scaling<1,T>& x)
    {}

    template <class M, class T>
    inline void DivEq(BaseMatrix_Mutable<M>& m, const Scaling<-1,T>& x)
    { Scale(x,m.mat()); }

    // -m
    template <class M>
    inline ProdXM<-1,RT,M> operator-(const BaseMatrix<M>& m)
    { return ProdXM<-1,RT,M>(RT(-1),m); }

    // x * m
    template <class M>
    inline ProdXM<0,RT,M> operator*(const int x, const BaseMatrix<M>& m)
    { return ProdXM<0,RT,M>(RT(x),m); }

    template <class M>
    inline ProdXM<0,RT,M> operator*(const RT x, const BaseMatrix<M>& m)
    { return ProdXM<0,RT,M>(x,m); }

    template <class M>
    inline ProdXM<0,CT,M> operator*(const CT x, const BaseMatrix<M>& m)
    { return ProdXM<0,CT,M>(x,m); }

    template <class M>
    inline ProdXM<0,CT,M> operator*(const CCT x, const BaseMatrix<M>& m)
    { return CT(x)*m; }

    template <class M, int ix, class T>
    inline ProdXM<ix,T,M> operator*(const Scaling<ix,T>& x, 
                                    const BaseMatrix<M>& m)
    { return ProdXM<ix,T,M>(T(x),m); }

    // m * x
    template <class M>
    inline ProdXM<0,RT,M> operator*(const BaseMatrix<M>& m, const int x)
    { return RT(x)*m; }

    template <class M>
    inline ProdXM<0,RT,M> operator*(const BaseMatrix<M>& m, const RT x)
    { return x*m; }

    template <class M>
    inline ProdXM<0,CT,M> operator*(const BaseMatrix<M>& m, const CT x)
    { return x*m; }

    template <class M>
    inline ProdXM<0,CT,M> operator*(const BaseMatrix<M>& m, const CCT x)
    { return CT(x)*m; }

    template <class M, int ix, class T>
    inline ProdXM<ix,T,M> operator*(const BaseMatrix<M>& m,
                                    const Scaling<ix,T>& x)
    { return ProdXM<ix,T,M>(T(x),m); }

    // m / x
    template <class M>
    inline ProdXM<0,RT,M> operator/(const BaseMatrix<M>& m, const int x)
    { return (RT(1)/RT(x))*m; }

    template <class M>
    inline ProdXM<0,RT,M> operator/(const BaseMatrix<M>& m, const RT x)
    { return (RT(1)/x)*m; }

    template <class M>
    inline ProdXM<0,CT,M> operator/(const BaseMatrix<M>& m, const CT x)
    { return (RT(1)/x)*m; }

    template <class M>
    inline ProdXM<0,CT,M> operator/(const BaseMatrix<M>& m, const CCT x)
    { return (RT(1)/CT(x))*m; }

    template <class M, int ix, class T>
    inline ProdXM<ix,T,M> operator/(const BaseMatrix<M>& m,
                                    const Scaling<ix,T>& x)
    { return ProdXM<ix,T,M>(RT(1)/T(x),m); }

#undef RT
#undef CT
#undef CCT

    // Consolidate x*x*m type constructs:

#define RT typename ProdXM<ix,T,M>::real_type
#define CT typename ProdXM<ix,T,M>::complex_type
#define CCT ConjRef<CT>

    // -(x*m)
    template <int ix, class T, class M>
    inline ProdXM<-ix,T,M> operator-(const ProdXM<ix,T,M>& pxm)
    { return ProdXM<-ix,T,M>(-pxm.getX(),pxm.getM()); }

    // x*(x*m)
    template <int ix, class T, class M>
    inline ProdXM<0,T,M> operator*(const int x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,T,M>(RT(x)*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,T,M> operator*(const RT x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,T,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,CT,M> operator*(const CT x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,CT,M> operator*(const CCT x, const ProdXM<ix,T,M>& pxm)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    inline ProdXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
        const Scaling<ix1,T1>& x, const ProdXM<ix,T,M>& pxm)
    {
        return ProdXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            T1(x)*pxm.getX(),pxm.getM()); 
    }

    // (x*m)*x
    template <int ix, class T, class M>
    inline ProdXM<0,T,M> operator*(const ProdXM<ix,T,M>& pxm, const int x)
    { return ProdXM<0,T,M>(RT(x)*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,T,M> operator*(const ProdXM<ix,T,M>& pxm, const RT x)
    { return ProdXM<0,T,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,CT,M> operator*(const ProdXM<ix,T,M>& pxm, const CT x)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,CT,M> operator*(const ProdXM<ix,T,M>& pxm, const CCT x)
    { return ProdXM<0,CT,M>(x*pxm.getX(),pxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    inline ProdXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
        const ProdXM<ix,T,M>& pxm, const Scaling<ix1,T1>& x)
    {
        return ProdXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            T1(x)*pxm.getX(),pxm.getM()); 
    }

    // (x*m)/x
    template <int ix, class T, class M>
    inline ProdXM<0,RT,M> operator/(const ProdXM<ix,T,M>& pxm, const int x)
    { return ProdXM<0,RT,M>(pxm.getX()/RT(x),pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,RT,M> operator/(const ProdXM<ix,T,M>& pxm, const RT x)
    { return ProdXM<0,RT,M>(pxm.getX()/x,pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,CT,M> operator/(const ProdXM<ix,T,M>& pxm, const CT x)
    { return ProdXM<0,CT,M>(pxm.getX()/x,pxm.getM()); }

    template <int ix, class T, class M>
    inline ProdXM<0,CT,M> operator/(const ProdXM<ix,T,M>& pxm, const CCT x)
    { return ProdXM<0,CT,M>(pxm.getX()/x,pxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    inline ProdXM<ix*ix1,typename Traits2<T1,T>::type,M> operator/(
        const ProdXM<ix,T,M>& pxm, const Scaling<ix1,T1>& x)
    { 
        return ProdXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            pxm.getX()/T1(x),pxm.getM()); 
    }

#undef RT
#undef CT
#undef CCT

#ifndef TMV_NO_ALIAS_CHECK
    // Have SameStorage look into a ProdXM object:
    template <int ix1, class T1, class M1, class M2>
    inline bool SameStorage(
        const ProdXM<ix1,T1,M1>& m1, const BaseMatrix_Calc<M2>& m2)
    { return SameStorage(m1.getM(),m2); }
    template <class M1, int ix2, class T2, class M2>
    inline bool SameStorage(
        const BaseMatrix_Calc<M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return SameStorage(m1,m2.getM()); }
    template <int ix1, class T1, class M1, int ix2, class T2, class M2>
    inline bool SameStorage(
        const ProdXM<ix1,T1,M1>& m1, const ProdXM<ix2,T2,M2>& m2)
    { return SameStorage(m1.getM(),m2.getM()); }
#endif
    template <int ix, class T, class M>
    inline std::string TMV_Text(const ProdXM<ix,T,M>& pxm)
    {
        std::ostringstream s;
        s << "ProdXM< "<<ix<<","<<TMV_Text(T(pxm.getX()));
        s << " , "<<TMV_Text(pxm.getM())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
