

#ifndef TMV_AddVV_H
#define TMV_AddVV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"
#include "TMV_MultXV_Funcs.h"

namespace tmv {

    // Defined in TMV_AddVV.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddVV(
        const T3 x1, const ConstVectorView<T1,C1>& v1,
        const T3 x2, const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddVV(
        const T3 x1, const ConstVectorView<T1,C1>& v1,
        const T3 x2, const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    //
    // Vector = x1 * Vector + x2 * Vector
    //

    template <int algo, ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper;

    // algo 0: size = 0, nothing to do
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<0,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& , const V1& ,
            const Scaling<ix2,T2>& , const V2& , V3& ) 
        {}
        static TMV_INLINE void call2(
            ptrdiff_t , const Scaling<ix1,T1>& , IT1 ,
            const Scaling<ix2,T2>& , IT2 , IT3 )
        {}
    };  

    // algo 1: complex vectors with unit step, convert to real version
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<1,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        typedef typename V1::const_flatten_type V1f;
        typedef typename V2::const_flatten_type V2f;
        typedef typename V3::flatten_type V3f;
        typedef typename V1f::const_nonconj_type::const_iterator IT1f;
        typedef typename V2f::const_nonconj_type::const_iterator IT2f;
        typedef typename V3f::iterator IT3f;
        enum { size2 = size == Unknown ? Unknown : (size<<1) };
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            V1f v1f = v1.flatten();
            V2f v2f = v2.flatten();
            V3f v3f = v3.flatten();
            AddVV_Helper<-4,size2,ix1,T1,V1f,ix2,T2,V2f,V3f>::call(
                x1,v1f,x2,v2f,v3f);
        }
        static inline void call2(
            ptrdiff_t n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            IT1f it1f = it1.flatten();
            IT2f it2f = it2.flatten();
            IT3f it3f = it3.flatten();
            const ptrdiff_t n2 = n<<1;
            AddVV_Helper<-4,size2,ix1,T1,V1f,ix2,T2,V2f,V3f>::call2(
                n2,x1,it1f,x2,it2f,it3f);
        }
    };

    // algo 11: simple for loop
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<11,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const ptrdiff_t n = size == Unknown ? v3.size() : size;
            for(ptrdiff_t i=0;i<n;++i) 
                v3.ref(i) = ZSum::sum(
                    ZProd<false,false>::prod(x1 , v1.cref(i)),
                    ZProd<false,false>::prod(x2 , v2.cref(i)));
        }
        static void call2(
            ptrdiff_t n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
            for(;n;--n) 
                *it3++ = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , *it1++), 
                    ZProd<false,c2>::prod(x2 , *it2++)); 
        }
    };  

    // algo 12: 2 at a time
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<12,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const ptrdiff_t n = size == Unknown ? v3.size() : size;
            call2(n,x1,v1.begin().nonConj(),x2,v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            const ptrdiff_t n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            ptrdiff_t n_2 = (n>>1);
            const ptrdiff_t nb = n-(n_2<<1);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;

            if (n_2) do {
                it3[0] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[0]),
                    ZProd<false,c2>::prod(x2 , it2[0]));
                it3[1] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[1]),
                    ZProd<false,c2>::prod(x2 , it2[1]));
                it1 += 2; it2 += 2; it3 += 2;
            } while (--n_2);
            if (nb) {
                *it3 = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , *it1),
                    ZProd<false,c2>::prod(x2 , *it2)); 
            }
        }
    };

    // algo 13: 4 at a time
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<13,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const ptrdiff_t n = size == Unknown ? v3.size() : size;
            call2(n,x1,v1.begin().nonConj(),x2,v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            const ptrdiff_t n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            ptrdiff_t n_4 = (n>>2);
            ptrdiff_t nb = n-(n_4<<2);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;

            if (n_4) do {
                it3[0] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[0]),
                    ZProd<false,c2>::prod(x2 , it2[0]));
                it3[1] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[1]),
                    ZProd<false,c2>::prod(x2 , it2[1]));
                it3[2] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[2]),
                    ZProd<false,c2>::prod(x2 , it2[2]));
                it3[3] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[3]),
                    ZProd<false,c2>::prod(x2 , it2[3]));
                it1 += 4; it2 += 4; it3 += 4;
            } while (--n_4);
            if (nb) do {
                *it3++ = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , *it1++),
                    ZProd<false,c2>::prod(x2 , *it2++)); 
            } while (--nb);
        }
    };

    // algo 14: 8 at a time
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<14,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const ptrdiff_t n = size == Unknown ? v3.size() : size;
            call2(n,x1,v1.begin().nonConj(),x2,v2.begin().nonConj(),v3.begin());
        }
        static void call2(
            const ptrdiff_t n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            ptrdiff_t n_8 = (n>>3);
            ptrdiff_t nb = n-(n_8<<3);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;

            if (n_8) do {
                it3[0] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[0]),
                    ZProd<false,c2>::prod(x2 , it2[0]));
                it3[1] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[1]),
                    ZProd<false,c2>::prod(x2 , it2[1]));
                it3[2] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[2]),
                    ZProd<false,c2>::prod(x2 , it2[2]));
                it3[3] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[3]),
                    ZProd<false,c2>::prod(x2 , it2[3]));
                it3[4] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[4]),
                    ZProd<false,c2>::prod(x2 , it2[4]));
                it3[5] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[5]),
                    ZProd<false,c2>::prod(x2 , it2[5]));
                it3[6] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[6]),
                    ZProd<false,c2>::prod(x2 , it2[6]));
                it3[7] = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , it1[7]),
                    ZProd<false,c2>::prod(x2 , it2[7]));
                it1 += 8; it2 += 8; it3 += 8;
            } while (--n_8);
            if (nb) do {
                *it3++ = ZSum::sum(
                    ZProd<false,c1>::prod(x1 , *it1++),
                    ZProd<false,c2>::prod(x2 , *it2++)); 
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<15,size,ix1,T1,V1,ix2,T2,V2,V3> // known size, unroll
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const V1& v1,
                const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
            {
                Unroller<I,N/2>::unroll(x1,v1,x2,v2,v3);
                Unroller<I+N/2,N-N/2>::unroll(x1,v1,x2,v2,v3);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const V1& v1,
                const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
            {
                v3.ref(I) = ZSum::sum(
                    ZProd<false,false>::prod(x1 , v1.cref(I)),
                    ZProd<false,false>::prod(x2 , v2.cref(I))); 
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const V1& v1,
                const Scaling<ix2,T2>& x2, const V2& v2, V3& v3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        { Unroller<0,size>::unroll(x1,v1,x2,v2,v3); }
    };

    // algo 90: call inst 
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<90,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAddVV(xx1,v1.xView(),xx2,v2.xView(),v3.xView());
        }
    };

    // algo 91: call inst alias
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<91,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAliasAddVV(xx1,v1.xView(),xx2,v2.xView(),v3.xView());
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<97,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            AddVV_Helper<-2,s,ix1,T1,V1c,ix2,T2,V2c,V3c>::call(
                TMV_CONJ(x1),v1c,TMV_CONJ(x2),v2c,v3c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<197,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            AddVV_Helper<99,s,ix1,T1,V1c,ix2,T2,V2c,V3c>::call(
                TMV_CONJ(x1),v1c,TMV_CONJ(x2),v2c,v3c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<98,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const bool s1 = SameStorage(v1,v3);
            const bool s2 = SameStorage(v2,v3);

            if (!s1 && !s2) {
                // No aliasing
                AddVV_Helper<-2,s,ix1,T1,V1,ix2,T2,V2,V3>::call(
                    x1,v1,x2,v2,v3);
            } else if (!s2) {
                // Alias with v1 only, do v1 first
                MultXV<false>(x1,v1,v3);
                typename V3::noalias_type v3na = v3.noAlias();
                MultXV<true>(x2,v2,v3na);
            } else if (!s1) {
                // Alias with v2 only, do v2 first
                MultXV<false>(x2,v2,v3);
                typename V3::noalias_type v3na = v3.noAlias();
                MultXV<true>(x1,v1,v3na);
            } else {
                // Need a temporary
                typename V1::copy_type v1c = v1;
                MultXV<false>(x2,v2,v3);
                typename V3::noalias_type v3na = v3.noAlias();
                MultXV<true>(x1,v1c,v3na);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<99,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type TV1;
            typedef typename V2::value_type TV2;
            typedef typename V3::value_type TV3;
            const bool inst =
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TV1,TV2>::samebase &&
                Traits2<TV1,TV3>::samebase &&
#else
                Traits2<TV1,TV2>::sametype &&
                Traits2<TV1,TV3>::sametype &&
#endif
                Traits<TV3>::isinst;
            const bool conj = V3::_conj;
            const int algo = 
                conj ? 197 :
                inst ? 91 :
                98;
            AddVV_Helper<algo,s,ix1,T1,V1,ix2,T2,V2,V3>::call(x1,v1,x2,v2,v3);
        }
    };

    // algo -4: No branches or copies
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-4,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        typedef typename V3::value_type T3;
        typedef typename V3::real_type RT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        enum { allreal = V1::isreal && V2::isreal && V3::isreal };
        enum { allcomplex = (
                V1::iscomplex && V2::iscomplex && V3::iscomplex ) };
        enum { allunit = V1::_step == 1 && V2::_step == 1 && V3::_step == 1 };
        enum { flatten = (
                allunit && allcomplex &&
                Traits<T1>::isreal && Traits<T2>::isreal &&
                V1::_conj == int(V3::_conj) && V2::_conj == int(V3::_conj) ) };
        enum { algo =  (
                size == 0 ? 0 : 
                TMV_OPT == 0 ? 11 :
                flatten ? 1 :
                (sizeof(RT) == 8 && allunit && allreal) ? 12 :
                (sizeof(RT) == 4 && allunit && allreal) ? 13 :
                11 ) };
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            const int algo1 = 
                (size != Unknown && size <= int(128/sizeof(T3))) ? 15 :
                algo;
            AddVV_Helper<algo1,size,ix1,T1,V1,ix2,T2,V2,V3>::call(
                x1,v1,x2,v2,v3);
        }
        static TMV_INLINE void call2(
            const ptrdiff_t n, const Scaling<ix1,T1>& x1, const IT1& it1,
            const Scaling<ix2,T2>& x2, const IT2& it2, const IT3& it3)
        {
            TMVStaticAssert(!V3::_conj);
            AddVV_Helper<algo,size,ix1,T1,V1,ix2,T2,V2,V3>::call2(
                n,x1,it1,x2,it2,it3);
        } 
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t size, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-3,size,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            AddVV_Helper<-4,size,ix1,T1,V1,ix2,T2,V2,V3>::call(
                x1,v1,x2,v2,v3);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-2,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type TV1;
            typedef typename V2::value_type TV2;
            typedef typename V3::value_type TV3;
            const bool inst =
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TV1,TV2>::samebase &&
                Traits2<TV1,TV3>::samebase &&
#else
                Traits2<TV1,TV2>::sametype &&
                Traits2<TV1,TV3>::sametype &&
#endif
                Traits<TV3>::isinst;
            const bool conj = V3::_conj;
            const int algo = 
                conj ? 97 :
                inst ? 90 :
                -4;
            AddVV_Helper<algo,s,ix1,T1,V1,ix2,T2,V2,V3>::call(x1,v1,x2,v2,v3);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t s, int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-1,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const int algo = 
                V3::_checkalias ? 99 : 
                -2;
            AddVV_Helper<algo,s,ix1,T1,V1,ix2,T2,V2,V3>::call(x1,v1,x2,v2,v3);
        }
    };

    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    inline void AddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const ptrdiff_t s = Sizes<Sizes<V1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        AddVV_Helper<-1,s,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(x1,v1v,x2,v2v,v3v);
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    inline void InlineAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const ptrdiff_t s = Sizes<Sizes<V1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        AddVV_Helper<-4,s,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(x1,v1v,x2,v2v,v3v);
    }

    template <int ix1, class T1, class V1, int ix2, class T2, class V2, class V3>
    inline void InlineAliasAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const ptrdiff_t s = Sizes<Sizes<V1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        AddVV_Helper<98,s,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(x1,v1v,x2,v2v,v3v);
    }

} // namespace tmv

#endif 
