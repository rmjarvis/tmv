
#ifndef TMV_CopyV_H
#define TMV_CopyV_H

#include "TMV_BaseVector.h"
#include <cstring>  // for memmove

namespace tmv {

    // Defined in TMV_Vector.cpp
    template <class T1, int C1, class T2>
    void InstCopy(const ConstVectorView<T1,C1>& v1, VectorView<T2> v2); 
    template <class T1, int C1, class T2>
    void InstAliasCopy(const ConstVectorView<T1,C1>& v1, VectorView<T2> v2); 

    //
    // Copy Vectors
    //

    template <int algo, ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper;

    // algo 0: s = 0, nothing to do
    template <class V1, class V2>
    struct CopyV_Helper<0,0,V1,V2>
    {
        typedef typename V1::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static TMV_INLINE void call(const V1&, V2& ) {}
        static TMV_INLINE void call2(ptrdiff_t , IT1 , IT2 ) {}
    };

    // algo 11: simple for loop
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<11,s,V1,V2>
    {
        typedef typename V1::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            for(ptrdiff_t i=0;i<n;++i) v2.ref(i) = v1.cref(i); 
        }
        static void call2(ptrdiff_t n, IT1 it1, IT2 it2)
        { for(;n;--n) *it2++ = *it1++; }
    };

    // algo 15: fully unroll
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<15,s,V1,V2>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void dounroll(const V1& v1, V2& v2)
            {
                Unroller<I,N/2>::dounroll(v1,v2);
                Unroller<I+N/2,N-N/2>::dounroll(v1,v2);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void dounroll(const V1& v1, V2& v2)
            { v2.ref(I) = v1.cref(I); }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        { static TMV_INLINE void dounroll(const V1& v1, V2& v2) {} };
        static inline void call(const V1& v1, V2& v2)
        { Unroller<0,s>::dounroll(v1,v2); }
    };

    // algo 21: memmove
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<21,s,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        {
            const ptrdiff_t n = s == Unknown ? v1.size() : s;
            memmove(v2.ptr(),v1.cptr(),n*sizeof(typename V2::value_type));
        }
        static inline void call2(
            const ptrdiff_t n, 
            typename V1::const_iterator it1, typename V2::iterator it2)
        { memmove(it2.get(),it1.get(),n*sizeof(typename V2::value_type)); }
    };

    // algo 22: std::copy
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<22,s,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        { std::copy(v1.begin(),v1.end(),v2.begin()); }
        static inline void call2(
            const ptrdiff_t n, 
            typename V1::const_iterator it1, typename V2::iterator it2)
        { std::copy(it1,it1+n,it2); }
    };

    // algo 90: Call inst
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<90,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        { InstCopy(v1.xView(),v2.xView()); }
    };

    // algo 91: Call inst alias
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<91,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        { InstAliasCopy(v1.xView(),v2.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<97,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            CopyV_Helper<-2,s,V1c,V2c>::call(v1c,v2c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<197,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            CopyV_Helper<99,s,V1c,V2c>::call(v1c,v2c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<98,s,V1,V2>
    {
        static void call(const V1& v1, V2& v2)
        {
            if ( !SameStorage(v1,v2) ||
                 v1.step()*v2.step() < 0 || 
                 std::abs(v2.step()) < std::abs(v1.step())) {
                // No aliasing (or no clobbering)
                CopyV_Helper<-2,s,V1,V2>::call(v1,v2);
            } else if (ExactSameStorage(v1,v2)) {
                // They are already equal modulo a conjugation.
                Maybe<V1::_conj != int(V2::_conj)>::conjself(v2); 
            } else if (v2.ptr()) { 
                // Need a temporary
                v2.noAlias() = v1.copy();
            } else {
                // else v2.ptr == 0, so don't need to do anything.
                // v1 and v2 are degenerate
                TMVAssert(v1.cptr() == 0);
                TMVAssert(v2.cptr() == 0);
                TMVAssert(v2.size() == 0);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<99,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                V2::_conj ? 197 :
                inst ? 91 :
                98;
            CopyV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<-3,s,V1,V2>
    {
        typedef typename V1::const_iterator IT1;
        typedef typename V2::iterator IT2;

        static void call(const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const int algo = 
                s == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                ( s != Unknown && s <= 8 ) ? 15 :
                ( Traits2<T1,T2>::sametype && 
                  V1::_conj == int(V2::_conj) &&
                  V1::_step == 1 && V2::_step == 1 ) ? 21 :
                11;
            CopyV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
        static void call2(const ptrdiff_t n, IT1 it1, IT2 it2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const int algo = 
                s == 0 ? 0 :
                ( TMV_OPT == 0 ) ? 11 :
                ( Traits2<T1,T2>::sametype && 
                  V1::_conj == int(V2::_conj) &&
                  V1::_step == 1 && V2::_step == 1 ) ? 21 :
                11;
            CopyV_Helper<algo,s,V1,V2>::call2(n,it1,it2);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<-2,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                V2::_conj ? 97 :
                inst ? 90 :
                -3;
            CopyV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t s, class V1, class V2>
    struct CopyV_Helper<-1,s,V1,V2>
    {
        static TMV_INLINE void call(const V1& v1, V2& v2)
        {
            const bool samestep = VStepHelper<V1,V2>::same;
            const bool noclobber = VStepHelper<V1,V2>::noclobber;
            const bool checkalias =
                V2::_checkalias && (samestep || !noclobber);
            const int algo = 
                checkalias ? 99 :
                -2;
            CopyV_Helper<algo,s,V1,V2>::call(v1,v2);
        }
    };

    template <class V1, class V2>
    inline void Copy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const ptrdiff_t s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        CopyV_Helper<-1,s,V1v,V2v>::call(v1v,v2v); 
    }

    template <class V1, class V2>
    inline void InlineCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const ptrdiff_t s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        CopyV_Helper<-3,s,V1v,V2v>::call(v1v,v2v); 
    }

    template <class V1, class V2>
    inline void InlineAliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const ptrdiff_t s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        CopyV_Helper<98,s,V1v,V2v>::call(v1v,v2v); 
    }


} // namespace tmv

#endif
