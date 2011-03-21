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


#ifndef TMV_SumVV_H
#define TMV_SumVV_H

#include "TMV_ProdXV.h"
#include "TMV_AddVV_Funcs.h"

namespace tmv {

    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    static inline void AddVV(
        const Scaling<ix1,T1>& x1, const BaseVector<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    { AddVV(x1,v1.calc(),x2,v2.calc(),v3.vec()); }
    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    static inline void NoAliasAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    { NoAliasAddVV(x1,v1.calc(),x2,v2.calc(),v3.vec()); }
    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    static inline void AliasAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    { AliasAddVV(x1,v1.calc(),x2,v2.calc(),v3.vec()); }

    //
    // Vector + Vector
    //

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    class SumVV;

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    struct Traits<SumVV<ix1,T1,V1,ix2,T2,V2> >
    {
        typedef typename ProdXV<ix1,T1,V1>::value_type vtype1;
        typedef typename ProdXV<ix2,T2,V2>::value_type vtype2;
        typedef typename Traits2<vtype1,vtype2>::type value_type;

        enum { _size = Sizes<V1::_size,V2::_size>::size };
        enum { _fort = V1::_fort && V2::_fort };
        enum { _calc = false };

        typedef SumVV<ix1,T1,V1,ix2,T2,V2> type;
        typedef typename VCopyHelper<value_type,_size,_fort>::type copy_type;
        typedef const copy_type calc_type;
        typedef typename TypeSelect<
            (V1::_calc && V2::_calc),const type,calc_type>::type eval_type;
    };

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    class SumVV : 
        public BaseVector<SumVV<ix1,T1,V1,ix2,T2,V2> >
    {
    public:

        typedef SumVV<ix1,T1,V1,ix2,T2,V2> type;
        typedef typename Traits<type>::value_type value_type;

        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        SumVV(const T1& _x1, const BaseVector<V1>& _v1,
                     const T2& _x2, const BaseVector<V2>& _v2) :
            x1(_x1), v1(_v1.vec()), x2(_x2), v2(_v2.vec())
        {
            TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
            TMVAssert(v1.size() == v2.size());
        }

        const Scaling<ix1,T1>& getX1() const { return x1; }
        const V1& getV1() const { return v1; }
        const Scaling<ix2,T2>& getX2() const { return x2; }
        const V2& getV2() const { return v2; }

        size_t size() const { return v1.size(); }
        value_type cref(int i) const
        { return x1 * v1.cref(i) + x2 * v2.cref(i); }

        template <class V3>
        void assignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
            TMVStaticAssert(type::isreal || V3::iscomplex);
            AddVV(x1,v1.vec(),x2,v2.vec(),v3.vec());
        }

        template <class V3>
        void newAssignTo(BaseVector_Mutable<V3>& v3) const
        {
            TMVStaticAssert((Sizes<type::_size,V3::_size>::same)); 
            TMVAssert(size() == v3.size());
            TMVStaticAssert(type::isreal || V3::iscomplex);
            NoAliasAddVV(x1,v1.vec(),x2,v2.vec(),v3.vec());
        }

    private:
        const Scaling<ix1,T1> x1;
        const V1& v1;
        const Scaling<ix2,T2> x2;
        const V2& v2;
    };

#define RT typename V2::real_type
    // v += v
    template <class V1, class V2>
    static inline void AddEq(
        BaseVector_Mutable<V1>& v1, const BaseVector<V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        TMVStaticAssert(V1::iscomplex || V2::isreal);
        MultXV<true>(v2.vec(),v1);
    }

    // v += xv
    template <class V1, int ix2, class T2, class V2>
    static inline void AddEq(
        BaseVector_Mutable<V1>& v1, const ProdXV<ix2,T2,V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        TMVStaticAssert(V1::iscomplex || V2::isreal);
        MultXV<true>(v2.getX(),v2.getV().vec(),v1);
    }

    // v -= v
    template <class V1, class V2>
    static inline void SubtractEq(
        BaseVector_Mutable<V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        TMVStaticAssert(V1::iscomplex || V2::isreal);
        MultXV<true>(Scaling<-1,RT>(),v2.vec(),v1);
    }

    // v -= xv
    template <class V1, int ix2, class T2, class V2>
    static inline void SubtractEq(
        BaseVector_Mutable<V1>& v1, const ProdXV<ix2,T2,V2>& v2) 
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        TMVStaticAssert(V1::iscomplex || V2::isreal);
        MultXV<true>(-v2.getX(),v2.getV().vec(),v1);
    }

    // v + v
    template <class V1, class V2>
    static inline SumVV<1,RT,V1,1,RT,V2> operator+(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<1,RT,V1,1,RT,V2>(RT(1),v1,RT(1),v2); 
    }

    // xv + v
    template <int ix1, class T1, class V1, class V2>
    static inline SumVV<ix1,T1,V1,1,RT,V2> operator+(
        const ProdXV<ix1,T1,V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<ix1,T1,V1,1,RT,V2>(T1(v1.getX()),v1.getV(),RT(1),v2); 
    }

    // v + xv
    template <class V1, int ix2, class T2, class V2>
    static inline SumVV<1,RT,V1,ix2,T2,V2> operator+(
        const BaseVector<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<1,RT,V1,ix2,T2,V2>(RT(1),v1,T2(v2.getX()),v2.getV()); 
    }

    // xv + xv
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<ix1,T1,V1,ix2,T2,V2> operator+(
        const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<ix1,T1,V1,ix2,T2,V2>(
            T1(v1.getX()),v1.getV(),T2(v2.getX()),v2.getV()); 
    }

    // v - v
    template <class V1, class V2>
    static inline SumVV<1,RT,V1,-1,RT,V2> operator-(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<1,RT,V1,-1,RT,V2>(RT(1),v1,RT(-1),v2); 
    }

    // xv - v
    template <int ix1, class T1, class V1, class V2>
    static inline SumVV<ix1,T1,V1,-1,RT,V2> operator-(
        const ProdXV<ix1,T1,V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<ix1,T1,V1,-1,RT,V2>(T1(v1.getX()),v1.getV(),RT(-1),v2); 
    }

    // v - xv
    template <class V1, int ix2, class T2, class V2>
    static inline SumVV<1,RT,V1,-ix2,T2,V2> operator-(
        const BaseVector<V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<1,RT,V1,-ix2,T2,V2>(RT(1),v1,T2(-v2.getX()),v2.getV()); 
    }

    // xv - xv
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<ix1,T1,V1,-ix2,T2,V2> operator-(
        const ProdXV<ix1,T1,V1>& v1, const ProdXV<ix2,T2,V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return SumVV<ix1,T1,V1,-ix2,T2,V2>(
            T1(v1.getX()),v1.getV(),T2(-v2.getX()),v2.getV()); 
    }
#undef RT


    // Consolidate x*(xv+xv) type constructs:

#define RT typename SumVV<ix1,T1,V1,ix2,T2,V2>::real_type
#define CT typename SumVV<ix1,T1,V1,ix2,T2,V2>::complex_type
#define CCT ConjRef<CT>
#define TX1 typename Traits2<T,T1>::type
#define TX2 typename Traits2<T,T2>::type

    // -(xv+xv)
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<-ix1,T1,V1,-ix2,T2,V2> operator-(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        return SumVV<-ix1,T1,V1,-ix2,T2,V2>(
            -T1(svv.getX1()),svv.getV1(),-T2(svv.getX2()),svv.getV2()); 
    }

    // x * (xv+xv)
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,T1,V1,0,T2,V2> operator*(
        const int x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        return SumVV<0,T1,V1,0,T2,V2>(
            RT(x)*svv.getX1(),svv.getV1(), RT(x)*svv.getX2(),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,T1,V1,0,T2,V2> operator*(
        const RT x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        return SumVV<0,T1,V1,0,T2,V2>(
            x*svv.getX1(),svv.getV1(), x*svv.getX2(),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,CT,V1,0,CT,V2> operator*(
        const CT x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        return SumVV<0,CT,V1,0,CT,V2>(
            x*svv.getX1(),svv.getV1(), x*svv.getX2(),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,CT,V1,0,CT,V2> operator*(
        const CCT x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        return SumVV<0,CT,V1,0,CT,V2>(
            CT(x)*svv.getX1(),svv.getV1(), CT(x)*svv.getX2(),svv.getV2()); 
    }

    template <int ix, class T, int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2> operator*(
        const Scaling<ix,T>& x, const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        return SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2>(
            T(x)*svv.getX1(),svv.getV1(),T(x)*svv.getX2(),svv.getV2());
    }

    // (xv+xv)*x
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,T1,V1,0,T2,V2> operator*(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const int x)
    {
        return SumVV<0,T1,V1,0,T2,V2>(
            RT(x)*svv.getX1(),svv.getV1(), RT(x)*svv.getX2(),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,T1,V1,0,T2,V2> operator*(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const RT x)
    {
        return SumVV<0,T1,V1,0,T2,V2>(
            x*svv.getX1(),svv.getV1(), x*svv.getX2(),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,CT,V1,0,CT,V2> operator*(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CT x)
    {
        return SumVV<0,CT,V1,0,CT,V2>(
            x*svv.getX1(),svv.getV1(), x*svv.getX2(),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,CT,V1,0,CT,V2> operator*(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CCT x)
    {
        return SumVV<0,CT,V1,0,CT,V2>(
            CT(x)*svv.getX1(),svv.getV1(), CT(x)*svv.getX2(),svv.getV2()); 
    }

    template <int ix, class T, int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2> operator*(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const Scaling<ix,T>& x)
    {
        return SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2>(
            T(x)*svv.getX1(),svv.getV1(),T(x)*svv.getX2(),svv.getV2());
    }

    // (xv+xv)/x
    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,T1,V1,0,T2,V2> operator/(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const int x)
    {
        return SumVV<0,T1,V1,0,T2,V2>(
            svv.getX1()/RT(x),svv.getV1(), svv.getX2()/RT(x),svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,T1,V1,0,T2,V2> operator/(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const RT x)
    {
        return SumVV<0,T1,V1,0,T2,V2>(
            svv.getX1()/x,svv.getV1(), svv.getX2()/x,svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,CT,V1,0,CT,V2> operator/(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CT x)
    {
        return SumVV<0,CT,V1,0,CT,V2>(
            svv.getX1()/x,svv.getV1(), svv.getX2()/x,svv.getV2()); 
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<0,CT,V1,0,CT,V2> operator/(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const CCT x)
    {
        return SumVV<0,CT,V1,0,CT,V2>(
            svv.getX1()/CT(x),svv.getV1(), svv.getX2()/CT(x),svv.getV2()); 
    }

    template <int ix, class T, int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2> operator/(
        const SumVV<ix1,T1,V1,ix2,T2,V2>& svv, const Scaling<ix,T>& x)
    {
        return SumVV<ix1*ix,TX1,V1,ix2*ix,TX2,V2>(
            svv.getX1()/T(x),svv.getV1(),svv.getX2()/T(x),svv.getV2());
    }

#undef RT
#undef CT
#undef CCT
#undef TX1
#undef TX2



    // TMV_Text

    template <int ix1, class T1, class V1, int ix2, class T2, class V2>
    static inline std::string TMV_Text(const SumVV<ix1,T1,V1,ix2,T2,V2>& svv)
    {
        std::ostringstream s;
        s << "SumVV< "<< ix1<<","<<TMV_Text(T1())
            <<" , "<<TMV_Text(svv.getV1())
            <<" , "<<ix2<<","<<TMV_Text(T2())
            <<" , "<<TMV_Text(svv.getV2())<<" >";
        return s.str();
    }

} // namespace tmv

#endif 
