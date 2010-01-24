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

#ifndef TMV_QuotXM_H
#define TMV_QuotXM_H

#include "TMV_ProdXV.h"
#include "TMV_ProdXM.h"
#include "TMV_BaseMatrix.h"

namespace tmv {
    
    //
    // Scalar * Matrix
    //

    // These are intentionally not defined to make sure we
    // get a compiler error if they are used.
    // All real calls should go through a more specific version than 
    // just the BaseMatrix_Calc's.
    template <int ix, class T, class M1, class M2>
    inline void Invert(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void NoAliasInvert(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void InlineInvert(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void AliasInvert(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& m2);


    template <int ix, class T, class M>
    class QuotXM;

    template <int ix, class T, class M>
    struct Traits<QuotXM<ix,T,M> >
    {
        typedef typename Traits2<T,typename M::value_type>::type value_type;

        enum { mcolsize = M::mcolsize };
        enum { mrowsize = M::mrowsize };
        enum { mshape = ShapeTraits<M::mshape>::nonunit_shape };
        enum { mfort = M::mfort };
        enum { mcalc = false };
        enum { rm1 = Traits<typename M::calc_type>::mrowmajor };
        enum { mrowmajor = rm1 };

        typedef QuotXM<ix,T,M> type;
        typedef typename MCopyHelper<value_type,mshape,mcolsize,mrowsize,
                mrowmajor,mfort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<M::mcalc,const type,calc_type>::type 
            eval_type;
        typedef InvalidType inverse_type;
    };

    template <int ix, class T, class M>
    class QuotXM : public BaseMatrix<QuotXM<ix,T,M> >
    {
    public:

        typedef QuotXM<ix,T,M> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        inline QuotXM(const T _x, const BaseMatrix<M>& _m) : 
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
            Invert(x,m.calc(),m2.mat());
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
            NoAliasInvert(x,m.calc(),m2.mat());
        }

    private:
        const Scaling<ix,T> x;
        const M& m;
    };

#define RT typename M::real_type
#define CT typename M::complex_type
#define CCT ConjRef<CT>

    // x / m
    template <class M>
    inline QuotXM<0,RT,M> operator/(const int x, const BaseMatrix<M>& m)
    { return QuotXM<0,RT,M>(RT(x),m); }

    template <class M>
    inline QuotXM<0,RT,M> operator/(const RT x, const BaseMatrix<M>& m)
    { return QuotXM<0,RT,M>(x,m); }

    template <class M>
    inline QuotXM<0,CT,M> operator/(const CT x, const BaseMatrix<M>& m)
    { return QuotXM<0,CT,M>(x,m); }

    template <class M>
    inline QuotXM<0,CT,M> operator/(const CCT x, const BaseMatrix<M>& m)
    { return CT(x)/m; }

    template <int ix, class T, class M>
    inline QuotXM<ix,T,M> operator/(
        const Scaling<ix,T>& x, const BaseMatrix<M>& m)
    { return QuotXM<ix,T,M>(T(x),m); }

    // x / xm
    template <int ix, class T, class M>
    inline QuotXM<0,T,M> operator/(const int x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<0,T,M>(RT(x)/pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,T,M> operator/(const RT x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<0,T,M>(x/pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator/(const CT x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<0,CT,M>(x/pxm.getX(),pxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator/(const CCT x, const ProdXM<ix,T,M>& pxm)
    { return CT(x)/pxm; }

    template <int ix1, class T1, int ix, class T, class M>
    inline QuotXM<ix1*ix,T,M> operator/(
        const Scaling<ix,T>& x, const ProdXM<ix,T,M>& pxm)
    { return QuotXM<ix1*ix,T,M>(T(x)/pxm.getX(),pxm.getM()); }

#undef RT
#undef CT
#undef CCT

    // Consolidate x*x/m type constructs:

#define RT typename QuotXM<ix,T,M>::real_type
#define CT typename QuotXM<ix,T,M>::complex_type
#define CCT ConjRef<CT>

    // -(x*m)
    template <int ix, class T, class M>
    inline QuotXM<-ix,T,M> operator-(const QuotXM<ix,T,M>& qxm)
    { return QuotXM<-ix,T,M>(-qxm.getX(),qxm.getM()); }

    // x*(x*m)
    template <int ix, class T, class M>
    inline QuotXM<0,T,M> operator*(const int x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,T,M>(RT(x)*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,T,M> operator*(const RT x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,T,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator*(const CT x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator*(const CCT x, const QuotXM<ix,T,M>& qxm)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    inline QuotXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
        const Scaling<ix1,T1>& x, const QuotXM<ix,T,M>& qxm)
    {
        return QuotXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            T1(x)*qxm.getX(),qxm.getM()); 
    }

    // (x*m)*x
    template <int ix, class T, class M>
    inline QuotXM<0,T,M> operator*(const QuotXM<ix,T,M>& qxm, const int x)
    { return QuotXM<0,T,M>(RT(x)*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,T,M> operator*(const QuotXM<ix,T,M>& qxm, const RT x)
    { return QuotXM<0,T,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator*(const QuotXM<ix,T,M>& qxm, const CT x)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator*(const QuotXM<ix,T,M>& qxm, const CCT x)
    { return QuotXM<0,CT,M>(x*qxm.getX(),qxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    inline QuotXM<ix*ix1,typename Traits2<T1,T>::type,M> operator*(
        const QuotXM<ix,T,M>& qxm, const Scaling<ix1,T1>& x)
    {
        return QuotXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            T1(x)*qxm.getX(),qxm.getM()); 
    }

    // (x*m)/x
    template <int ix, class T, class M>
    inline QuotXM<0,RT,M> operator/(const QuotXM<ix,T,M>& qxm, const int x)
    { return QuotXM<0,RT,M>(qxm.getX()/RT(x),qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,RT,M> operator/(const QuotXM<ix,T,M>& qxm, const RT x)
    { return QuotXM<0,RT,M>(qxm.getX()/x,qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator/(const QuotXM<ix,T,M>& qxm, const CT x)
    { return QuotXM<0,CT,M>(qxm.getX()/x,qxm.getM()); }

    template <int ix, class T, class M>
    inline QuotXM<0,CT,M> operator/(const QuotXM<ix,T,M>& qxm, const CCT x)
    { return QuotXM<0,CT,M>(qxm.getX()/x,qxm.getM()); }

    template <int ix1, class T1, int ix, class T, class M>
    inline QuotXM<ix*ix1,typename Traits2<T1,T>::type,M> operator/(
        const QuotXM<ix,T,M>& qxm, const Scaling<ix1,T1>& x)
    { 
        return QuotXM<ix*ix1,typename Traits2<T1,T>::type,M>(
            qxm.getX()/T1(x),qxm.getM()); 
    }

#undef RT
#undef CT
#undef CCT

    template <int ix, class T, class M>
    inline std::string TMV_Text(const QuotXM<ix,T,M>& qxm)
    {
        std::ostringstream s;
        s << "QuotXM< "<<ix<<","<<TMV_Text(T(qxm.getX()));
        s << " , "<<TMV_Text(qxm.getM())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
