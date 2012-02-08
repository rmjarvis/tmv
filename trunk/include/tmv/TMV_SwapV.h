
#ifndef TMV_SwapV_H
#define TMV_SwapV_H

#include "TMV_BaseVector.h"

namespace tmv {

    // Defined in TMV_Vector.cpp
    template <class T, int C>
    void InstSwap(VectorView<T,C> v1, VectorView<T> v2);
    template <class T, int C>
    void InstAliasSwap(VectorView<T,C> v1, VectorView<T> v2);

    //
    // Swap Vectors
    //

    template <int algo, ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper;

    // algo 2: complex vectors with unit step, convert to real version
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<2,s,V1,V2>
    {
        static inline void call(V1& v1, V2& v2)
        {
            typedef typename V1::flatten_type V1f;
            typedef typename V2::flatten_type V2f;
            const ptrdiff_t s2 = IntTraits<s>::twoS;
            V1f v1f = v1.flatten();
            V2f v2f = v2.flatten();
            SwapV_Helper<-3,s2,V1f,V2f>::call(v1f,v2f);
        }
    };

    // algo 11: simple for loop
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<11,s,V1,V2>
    {
        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            for(ptrdiff_t i=0;i<n;++i) TMV_SWAP(v1.ref(i),v2.ref(i));
        }
        static void call2(ptrdiff_t n, IT1 it1, IT2 it2)
        { for(;n;--n) TMV_SWAP(*it1++,*it2++); }
    };

    // algo 12: 2 at a time
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<12,s,V1,V2>
    {
        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            call2(n,v1.begin(),v2.begin());
        }
        static void call2(const ptrdiff_t n, IT1 it1, IT2 it2)
        {
            typedef typename V1::value_type T1;
            T1 t0, t1;

            ptrdiff_t n_2 = (n>>1);
            const ptrdiff_t nb = n-(n_2<<1);
            if (n_2) do {
                t0 = it1[0];
                t1 = it1[1];
                it1[0] = it2[0];
                it1[1] = it2[1]; it1 += 2;
                it2[0] = t0;
                it2[1] = t1; it2 += 2;
            } while (--n_2);
            if (nb) {
                t0 = *it1;
                *it1 = *it2;
                *it2 = t0;
            }
        }
    };

    // algo 13: 4 at a time
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<13,s,V1,V2>
    {
        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            call2(n,v1.begin(),v2.begin());
        }
        static void call2(const ptrdiff_t n, IT1 it1, IT2 it2)
        {
            typedef typename V1::value_type T1;
            T1 t0, t1, t2, t3;

            ptrdiff_t n_4 = (n>>2);
            ptrdiff_t nb = n-(n_4<<2);

            if (n_4) do {
                t0 = it1[0];
                t1 = it1[1];
                t2 = it1[2];
                t3 = it1[3];
                it1[0] = it2[0];
                it1[1] = it2[1];
                it1[2] = it2[2];
                it1[3] = it2[3]; it1 += 4;
                it2[0] = t0;
                it2[1] = t1;
                it2[2] = t2;
                it2[3] = t3; it2 += 4;
            } while (--n_4);
            if (nb) do {
                t0 = *it1;
                *it1++ = *it2;
                *it2++ = t0;
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<15,s,V1,V2>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(V1& v1, V2& v2)
            {
                Unroller<I,N/2>::unroll(v1,v2);
                Unroller<I+N/2,N-N/2>::unroll(v1,v2);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void unroll(V1& v1, V2& v2)
            { TMV_SWAP(v1.ref(I),v2.ref(I)); }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        { static TMV_INLINE void unroll(V1& v1, V2& v2) {} };
        static inline void call(V1& v1, V2& v2)
        { Unroller<0,s>::unroll(v1,v2); }
    };

    // algo 21: complex vectors, but not unit step
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<21,s,V1,V2>
    {
        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            call2(n,v1.begin(),v2.begin());
        }
        static void call2(ptrdiff_t n, IT1 it1, IT2 it2)
        {
            typedef typename V1::real_type RT;
            RT t0, t1;
            if (n) do {
                t0 = real(*it1);
                t1 = imag(*it1);
                real(*it1) = real(*it2);
                imag(*it1) = imag(*it2); ++it1;
                real(*it2) = t0;
                imag(*it2) = t1; ++it2;
            } while (--n);
        }
    };

    // algo 22: complex vectors, but v1 is conjugate
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<22,s,V1,V2>
    {
        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            call2(n,v1.begin(),v2.begin());
        }
        static void call2(ptrdiff_t n, IT1 it1, IT2 it2)
        {
            typedef typename IT1::nonconj_type IT1n;
            typedef typename V1::real_type RT;
            RT t0, t1;
            IT1n it1n = it1.nonConj();
            if (n) do {
                t0 = real(*it1n);
                t1 = imag(*it1n);
                real(*it1n) = real(*it2);
                imag(*it1n) = -imag(*it2); ++it1n;
                real(*it2) = t0;
                imag(*it2) = -t1; ++it2;
            } while (--n);
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<90,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        { InstSwap(v1.xView(),v2.xView()); }
    };

    // algo 91: Call inst alias
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<91,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        { InstAliasSwap(v1.xView(),v2.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<97,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        {
            typedef typename V1::conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            SwapV_Helper<-2,s,V1c,V2c>::call(v1c,v2c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<197,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        {
            typedef typename V1::conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            SwapV_Helper<99,s,V1c,V2c>::call(v1c,v2c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<98,s,V1,V2>
    {
        static void call(V1& v1, V2& v2)
        {
            if ( !SameStorage(v1,v2) ||
                 v1.step()*v2.step() < 0) {
                // No aliasing (or no clobbering)
                SwapV_Helper<-2,s,V1,V2>::call(v1,v2);
            } else if (ExactSameStorage(v1,v2)) {
                // They are already equal modulo a conjugation.
                Maybe<V1::_conj != int(V2::_conj)>::conjself(v2);
            } else {
                // Need a temporary
                typename V1::copy_type v1c = v1;
                v1.noAlias() = v2;
                v2.noAlias() = v1c;
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<99,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                (s == Unknown || s > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo = 
                V2::_conj ? 197 :
                inst ? 91 :
                98;
            SwapV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<-3,s,V1,V2>
    {
        typedef typename V1::iterator IT1;
        typedef typename V2::iterator IT2;
        typedef typename V1::real_type RT;

        enum { allunit = V1::_step == 1 && V2::_step == 1 };
        enum { algo = (
                TMV_OPT == 0 ? 11 :
                // Strangely, algo 2 doesn't seem to be faster.
                //(V1::iscomplex && allunit && !V1::_conj) ? 2 :
                (V1::iscomplex) ? (V1::_conj ? 22 : 21) :
                (sizeof(RT) == 8 && allunit) ? 12 :
                (sizeof(RT) == 4 && allunit) ? 13 :
                11 ) };

        static TMV_INLINE void call(V1& v1, V2& v2)
        {
            TMVStaticAssert(!V2::_conj);
            const int algo1 = 
                s != Unknown && s <= ptrdiff_t(128/sizeof(RT)) ? 15 :
                algo;
            SwapV_Helper<algo1,s,V1,V2>::call(v1,v2); 
        }
        static TMV_INLINE void call2(ptrdiff_t n, IT1 it1, IT2 it2)
        {
            TMVStaticAssert(!V2::_conj);
            SwapV_Helper<algo,s,V1,V2>::call2(n,it1,it2); 
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<-2,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                (s == Unknown || s > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo = 
                V2::_conj ? 97 :
                inst ? 90 :
                -3;
            SwapV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t s, class V1, class V2>
    struct SwapV_Helper<-1,s,V1,V2>
    {
        static TMV_INLINE void call(V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            // This seems impossible at first glance, but it is if
            // the steps have opposite sign:
            const bool noclobber = 
                VStepHelper<V1,V2>::noclobber &&
                VStepHelper<V2,V1>::noclobber &&
                !VStepHelper<V1,V2>::same;
            const bool checkalias =
                (V1::_checkalias || V2::_checkalias) && !noclobber;
            const int algo = 
                checkalias ? 99 : 
                -2;
            SwapV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
    };

    template <class V1, class V2>
    inline void DoSwap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const ptrdiff_t s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        SwapV_Helper<-1,s,V1v,V2v>::call(v1v,v2v);
    }

    template <class V1, class V2>
    inline void InlineSwap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const ptrdiff_t s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        SwapV_Helper<-3,s,V1v,V2v>::call(v1v,v2v);
    }

    template <class V1, class V2>
    inline void InlineAliasSwap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const ptrdiff_t s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        SwapV_Helper<98,s,V1v,V2v>::call(v1v,v2v);
    }

    template <class V1, class V2>
    TMV_INLINE void Swap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    { DoSwap(v1,v2); }

} // namespace tmv

#endif
