
#ifndef TMV_ScaleV_H
#define TMV_ScaleV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

    // Defined in TMV_ScaleV.cpp
    template <class T>
    void InstScale(const T x, VectorView<T> v);  

    //
    // Vector *= x
    //

    template <int algo, int s, int ix, class T, class V>
    struct ScaleV_Helper;

    // algo 0: trivial: s == 0 or ix == 1, so nothing to do
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<0,s,ix,T,V>
    { static TMV_INLINE void call(const Scaling<ix,T>& , V& ) {} };

    // algo 1: complex vector with unit step, convert to real version.
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<1,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        typedef typename V::flatten_type Vf;
        typedef typename Vf::iterator ITf;
        typedef typename V::real_type RT;
        enum { s2 = IntTraits<s>::twoS };
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            Vf vf = v.flatten();
            ScaleV_Helper<-3,s2,ix,T,Vf>::call(x,vf);
        }
        static void call2(int n, const Scaling<ix,T>& x, IT it)
        {
            ITf itf = it.flatten();
            const int n2 = n<<1;
            ScaleV_Helper<-3,s2,ix,T,Vf>::call2(n2,x,itf);
        }
    };

    // algo 11: simple loop
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<11,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT it)
        {
            if (n) do {
                *it = ZProd<false,false>::prod(x,*it); ++it;
            } while (--n);
        }
    };

    // algo 12: 2 at a time
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<12,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(const int n, const Scaling<ix,T>& x, IT it)
        {
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) do {
                it[0] = ZProd<false,false>::prod(x,it[0]);
                it[1] = ZProd<false,false>::prod(x,it[1]); it += 2;
            } while (--n_2);
            if (nb) *it = ZProd<false,false>::prod(x,*it);
        }
    };

    // algo 13: 4 at a time
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<13,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(const int n, const Scaling<ix,T>& x, IT it)
        {
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) do {
                it[0] = ZProd<false,false>::prod(x,it[0]);
                it[1] = ZProd<false,false>::prod(x,it[1]);
                it[2] = ZProd<false,false>::prod(x,it[2]);
                it[3] = ZProd<false,false>::prod(x,it[3]); it += 4;
            } while (--n_4);
            if (nb) do {
                *it = ZProd<false,false>::prod(x,*it); ++it;
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<15,s,ix,T,V>
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(const Scaling<ix,T>& x, V& v)
            {
                Unroller<I,N/2>::unroll(x,v);
                Unroller<I+N/2,N-N/2>::unroll(x,v);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(const Scaling<ix,T>& x, V& v)
            { v.ref(I) = ZProd<false,false>::prod(x,v.cref(I)); }
        };
        template <int I>
        struct Unroller<I,0>
        { static inline void unroll(const Scaling<ix,T>& , V& ) {} };

        static inline void call(const Scaling<ix,T>& x, V& v)
        { Unroller<0,s>::unroll(x,v); }
    };

#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<31,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::_step == 1;
            if (unit) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A++ *= x; 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT A1 = A+1;
                double x1 = x;
                __m128d xx = _mm_set1_pd(x1);
                __m128d xA;
                do {
                    Maybe<unit>::sse_load(xA,A.get(),A1.get());
                    xA = _mm_mul_pd(xA,xx);
                    Maybe<unit>::sse_store(A.get(),A1.get(),xA);
                    A+=2; A1+=2;
                } while (--n_2);
            }
            if (nb) *A++ *= x;
        }
    };

    // algo 32: double precision SSE2: x real, v complex
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<32,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            if (n) {
                double x1 = x;
                __m128d xx = _mm_set1_pd(x1);
                __m128d xA;
                do {
                    Maybe<true>::sse_load(xA,A.get());
                    xA = _mm_mul_pd(xA,xx);
                    Maybe<true>::sse_store(A.get(),xA);
                    ++A;
                } while (--n);
            }
        }
    };

    // algo 33: double precision SSE2: all complex
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<33,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            if (n) {
                std::complex<double> xx(x);
                double xr = real(xx);
                double xi = imag(xx);
                __m128d xxr = _mm_set1_pd(xr);
                __m128d xxi = _mm_set_pd(xi , -xi);
                __m128d xA;
                __m128d x0, x1, x2; // temp vars
                if (TMV_Aligned(A.get()) ) {
                    do {
                        // r = xr * Ar - xi * Ai
                        // i = xr * Ai + xi * Ar
                        Maybe<true>::sse_load(xA,A.get());
                        x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                        x1 = _mm_mul_pd(xxr,xA);
                        x2 = _mm_mul_pd(xxi,x0);
                        xA = _mm_add_pd(x1,x2);
                        Maybe<true>::sse_store(A.get(),xA);
                        ++A;
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.get());
                        x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                        x1 = _mm_mul_pd(xxr,xA);
                        x2 = _mm_mul_pd(xxi,x0);
                        xA = _mm_add_pd(x1,x2);
                        Maybe<true>::sse_storeu(A.get(),xA);
                        ++A;
                    } while (--n);
                }
            }
        }
    };
#endif

#ifdef __SSE__
    // algo 41: single precision SSE: all real
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<41,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::_step == 1;
            if (unit ) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A++ *= x; 
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) {
                IT A1 = A+1;
                IT A2 = A+2;
                IT A3 = A+3;
                float x1 = x;
                __m128 xx = _mm_set1_ps(x1);
                __m128 xA;
                do {
                    Maybe<unit>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    xA = _mm_mul_ps(xA,xx);
                    Maybe<unit>::sse_store(
                        A.get(),A1.get(),A2.get(),A3.get(),xA);
                    A+=4; A1+=4; A2+=4; A3+=4;
                } while (--n_4);
            }
            if (nb) do { *A++ *= x; } while (--nb);
        }
    };

    // algo 42: single precision SSE: x real, v complex
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<42,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::_step == 1;
            if (unit ) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A++ *= x; 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT A1 = A+1;
                float x1 = x;
                __m128 xx = _mm_set1_ps(x1);
                __m128 xA;
                do {
                    Maybe<unit>::sse_load(xA,A.get(),A1.get());
                    xA = _mm_mul_ps(xA,xx);
                    Maybe<unit>::sse_store(A.get(),A1.get(),xA);
                    A+=2; A1+=2;
                } while (--n_2);
            }
            if (nb) *A++ *= x; 
        }
    };

    // algo 43: single precision SSE: all complex
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<43,s,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = s == TMV_UNKNOWN ? int(v.size()) : s;
            call2(n,x,v.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::_step == 1;

            if (unit) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A = ZProd<false,false>::prod(x,*A); ++A;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT A1 = A+1;
                std::complex<float> xx(x);
                float xr = real(xx);
                float xi = imag(xx);
                __m128 xxr = _mm_set1_ps(xr);
                __m128 xxi = _mm_set_ps(xi , -xi , xi , -xi);
                __m128 xA;
                __m128 x0, x1, x2; // temp vars
                do {
                    // r = xr * Ar - xi * Ai
                    // i = xr * Ai + xi * Ar
                    Maybe<unit>::sse_load(xA,A.get(),A1.get());
                    x0 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(2,3,0,1));
                    x1 = _mm_mul_ps(xxr,xA); // xr*Ar, xr*Ai
                    x2 = _mm_mul_ps(xxi,x0); // -xi*Ai, xi*Ar
                    xA = _mm_add_ps(x1,x2);
                    Maybe<unit>::sse_store(A.get(),A1.get(),xA);
                    A += 2; A1 += 2;
                } while (--n_2);
            }
            if (nb) *A = ZProd<false,false>::prod(x,*A);
        }
    };
#endif

    // algo 90: Call inst
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<90,s,ix,T,V>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, V& v)
        {
            typedef typename V::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstScale(xx,v.xView());
        }
    };

    // algo 97: Conjugate
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<97,s,ix,T,V>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ScaleV_Helper<-2,s,ix,T,Vc>::call(TMV_CONJ(x),vc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<-3,s,ix,T,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::iterator IT;
        enum { unit = V::_step == 1 };
        enum { vfloat = Traits2<RT,float>::sametype };
        enum { vdouble = Traits2<RT,double>::sametype };
        enum { xreal = Traits<T>::isreal };
        enum { xcomplex = Traits<T>::iscomplex };
        enum { vreal = V::isreal };
        enum { vcomplex = V::iscomplex };
#if 0
        enum { algo = 11 };
#else
        enum { algo = (
                (s == 0 || ix == 1) ? 0 :
                (unit && xreal && vcomplex) ? 1 :
                TMV_OPT == 0 ? 11 :
#ifdef __SSE__
                (vfloat && xreal && vreal) ? 41 :
                (vfloat && xreal && vcomplex) ? 42 :
                (vfloat && xcomplex && vcomplex) ? 43 :
#endif
#ifdef __SSE2__
                (vdouble && xreal && vreal) ? 31 :
                (vdouble && xreal && vcomplex) ? 32 :
                (vdouble && xcomplex && vcomplex) ? 33 :
#endif
                (vreal && sizeof(RT) == 4) ? 13 :
                (vreal && sizeof(RT) == 8) ? 12 :
                11 ) };
#endif
        static TMV_INLINE void call(const Scaling<ix,T>& x, V& v)
        {
            TMVStaticAssert(!V::_conj);
            ScaleV_Helper<algo,s,ix,T,V>::call(x,v); 
        }
        static TMV_INLINE void call2(
            const int n, const Scaling<ix,T>& x, const IT& it)
        {
            TMVStaticAssert(!V::_conj);
            ScaleV_Helper<algo,s,ix,T,V>::call2(n,x,it); 
        }
    };

    // algo -2: Check for inst
    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<-2,s,ix,T,V>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, V& v)
        {
            typedef typename V::value_type T1;
            const bool inst =
                (s == TMV_UNKNOWN || s > 16) &&
                Traits<T1>::isinst;
            const bool conj = V::_conj;
            const int algo = 
                ix == 1 ? 0 :
                conj ? 97 :
                inst ? 90 : 
                -3;
            ScaleV_Helper<algo,s,ix,T,V>::call(x,v);
        }
    };

    template <int s, int ix, class T, class V>
    struct ScaleV_Helper<-1,s,ix,T,V>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, V& v)
        { ScaleV_Helper<-2,s,ix,T,V>::call(x,v); }
    };

    template <int ix, class T, class V>
    static inline void Scale(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ScaleV_Helper<-2,V::_size,ix,T,Vv>::call(x,vv);
    }

    template <int ix, class T, class V>
    static inline void InlineScale(
        const Scaling<ix,T>& x, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ScaleV_Helper<-3,V::_size,ix,T,Vv>::call(x,vv);
    }

} // namespace tmv

#endif 
