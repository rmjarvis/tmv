

#ifndef TMV_MultXV_H
#define TMV_MultXV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"
#include "TMV_MultXV_Funcs.h"
#include "TMV_CopyV.h"
#include "TMV_ScaleV.h"
#include "TMV_Array.h"

#ifdef PRINTALGO_XV
#include <iostream>
#include "TMV_VectorIO.h"
#endif

namespace tmv {

    // Defined in TMV_MultXV.cpp
    template <class T1, int C1, class T2>
    void InstMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2);
    template <class T1, int C1, class T2>
    void InstAddMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2);

    template <class T1, int C1, class T2>
    void InstAliasMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2);
    template <class T1, int C1, class T2>
    void InstAliasAddMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2);

    //
    // Vector = x * Vector
    //

    template <int algo, int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper;

    // algo 0: s == 0, nothing to do
    template <bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<0,0,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static TMV_INLINE void call(const Scaling<ix,T>& , const V1& , V2& ) {}
        static TMV_INLINE void call2(
            int n, const Scaling<ix,T>& , IT1 , IT2 ) {}
    };

    // algo 1: trivial: ix == 1, !add, so call Copy
    template <int s, class T, class V1, class V2>
    struct MultXV_Helper<1,s,false,1,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static TMV_INLINE void call(const Scaling<1,T>& , const V1& v1, V2& v2)
        { CopyV_Helper<-3,s,V1,V2>::call(v1,v2); }
        static TMV_INLINE void call2(int n, const Scaling<1,T>& , IT1 A, IT2 B)
        { CopyV_Helper<-3,s,V1,V2>::call2(A,B); }
    };

    // algo 101: same as 1, but use -1 algo
    template <int s, class T, class V1, class V2>
    struct MultXV_Helper<101,s,false,1,T,V1,V2>
    {
        static TMV_INLINE void call(const Scaling<1,T>& , const V1& v1, V2& v2)
        { CopyV_Helper<-1,s,V1,V2>::call(v1,v2); }
    };

    // algo 201: same as 1, but use -2 algo
    template <int s, class T, class V1, class V2>
    struct MultXV_Helper<201,s,false,1,T,V1,V2>
    {
        static TMV_INLINE void call(const Scaling<1,T>& , const V1& v1, V2& v2)
        { CopyV_Helper<-2,s,V1,V2>::call(v1,v2); }
    };

    // algo 2: complex vectors with unit step, convert to real version
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<2,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V1::const_flatten_type V1f;
            typedef typename V2::flatten_type V2f;
            const int s2 = s == TMV_UNKNOWN ? TMV_UNKNOWN : (s<<1);
            V1f v1f = v1.flatten();
            V2f v2f = v2.flatten();
            MultXV_Helper<-4,s2,add,ix,T,V1f,V2f>::call(x,v1f,v2f);
        }
        static inline void call2(
            const int n, const Scaling<ix,T>& x, 
            const IT1& A, const IT2& B)
        {
            typedef typename V1::const_flatten_type V1f;
            typedef typename V2::flatten_type V2f;
            typedef typename V1f::const_iterator IT1f;
            typedef typename V2f::iterator IT2f;
            const int s2 = s == TMV_UNKNOWN ? TMV_UNKNOWN : (s<<1);
            const int n2 = s == TMV_UNKNOWN ? 2*n : s2;
            IT1f Af = A.flatten();
            IT2f Bf = B.flatten();
            MultXV_Helper<-4,s2,add,ix,T,V1f,V2f>::call2(n2,x,Af,Bf);
        }
    };

    // algo 3: copy v1 to new storage
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<3,s,true,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        { 
            typename V2::noalias_type v2na = v2.noAlias();
            MultXV<true>(x,v1.copy(),v2); 
        }
    };

    // algo 4: copy v1 to v2 and turn into a MultEq op
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<4,s,false,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            CopyV_Helper<-3,s,V1,V2>::call(v1,v2);
            ScaleV_Helper<-3,s,ix,T,V2>::call(x,v2);
        }
    };

    // algo 11: simple for loop
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<11,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool c1 = V1::_conj;
            if (n) do {
                Maybe<add>::add(*B++ , ZProd<false,c1>::prod(x , *A++)); 
            } while (--n);
        }
    };  

    // algo 12: 2 at a time
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<12,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(const int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);
            const bool c1 = V1::_conj;

            if (n_2) do {
                Maybe<add>::add(B[0] , ZProd<false,c1>::prod(x , A[0]));
                Maybe<add>::add(B[1] , ZProd<false,c1>::prod(x , A[1]));
                A += 2; B += 2;
            } while (--n_2);
            if (nb) {
                Maybe<add>::add(*B , ZProd<false,c1>::prod(x , *A)); 
            }
        }
    };

    // algo 13: 4 at a time
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<13,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(const int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            const bool c1 = V1::_conj;

            if (n_4) do {
                Maybe<add>::add(B[0] , ZProd<false,c1>::prod(x , A[0]));
                Maybe<add>::add(B[1] , ZProd<false,c1>::prod(x , A[1]));
                Maybe<add>::add(B[2] , ZProd<false,c1>::prod(x , A[2]));
                Maybe<add>::add(B[3] , ZProd<false,c1>::prod(x , A[3]));
                A += 4; B += 4;
            } while (--n_4);
            if (nb) do {
                Maybe<add>::add(*B++ , ZProd<false,c1>::prod(x , *A++)); 
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<15,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const V1& v1, V2& v2)
            {
                Unroller<I,N/2>::unroll(x,v1,v2);
                Unroller<I+N/2,N-N/2>::unroll(x,v1,v2);
            }
            static inline void unroll2(
                const Scaling<ix,T>& x, const IT1& A, const IT2& B)
            {
                Unroller<I,N/2>::unroll2(x,A,B);
                Unroller<I+N/2,N-N/2>::unroll2(x,A,B);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const V1& v1, V2& v2)
            {
                Maybe<add>::add(
                    v2.ref(I), ZProd<false,false>::prod(x , v1.cref(I))); 
            }
            static inline void unroll2(
                const Scaling<ix,T>& x, const IT1& A, const IT2& B)
            {
                const bool c1 = V1::_conj;
                Maybe<add>::add(B[I] , ZProd<false,c1>::prod(x , A[I])); 
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(const Scaling<ix,T>& , const V1& , V2& ) 
            {}
            static inline void unroll2(
                const Scaling<ix,T>& , const IT1& , const IT2& ) 
            {}
        };
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        { Unroller<0,s>::unroll(x,v1,v2); }
        static inline void call2(
            const int , const Scaling<ix,T>& x, const IT1& A, const IT2& B)
        { Unroller<0,s>::unroll2(x,A,B); }
    };

#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<31,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    Maybe<add>::add(*B++ , x * *A++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , x * *A++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xx = _mm_set1_pd(double(x));
                __m128d xA,xB;
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    xB = _mm_mul_pd(xx,xA);
                    Maybe2<add,unit2>::sse_add(B.get(),B1.get(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) Maybe<add>::add(*B , x * *A);
        }
    };

    // algo 32: double precision SSE2: x real v1 real v2 complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<32,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , x * *A++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xx = _mm_set1_pd(double(x));
                __m128d xzero = _mm_set1_pd(0.);
                __m128d xA,xB;
                do {
                    Maybe<unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    xB = _mm_mul_pd(xx,xA);
                    Maybe<unit2>::sse_addu(B.get(),_mm_unpacklo_pd(xB,xzero));
                    Maybe<unit2>::sse_addu(B1.get(),_mm_unpackhi_pd(xB,xzero));
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) Maybe<add>::add(*B , x * *A);
        }
    };

    // algo 33: double precision SSE2: x real v1 complex v2 complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<33,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            // TODO: This must be wrong.  c1 isn't used!
            // But I'm not using this function right now, so I won't 
            // bother trying to fix it.
            //const bool c1 = V1::_conj;
            if (n) {
                __m128d xx = _mm_set1_pd(double(x));
                __m128d xA,xB;
                if (TMV_Aligned(A.get()) ) {
                    do {
                        Maybe<true>::sse_load(xA,A.get()); ++A;
                        xB = _mm_mul_pd(xx,xA); 
                        Maybe2<add,true>::sse_add(B.get(),xB); ++B;
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.get()); ++A;
                        xB = _mm_mul_pd(xx,xA); 
                        Maybe2<add,true>::sse_addu(B.get(),xB); ++B;
                    } while (--n);
                }
            }
        }
    };

    // algo 34: double precision SSE2: x complex v1 real v2 complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<34,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<double> >::sametype));
            const bool unit1 = V1::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , ZProd<false,false>::prod(x , *A++)); 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xr = _mm_set1_pd(real(x.x));
                __m128d xi = _mm_set1_pd(imag(x.x));
                __m128d xA,xBr,xBi;
                do {
                    Maybe<unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    xBr = _mm_mul_pd(xr,xA);
                    xBi = _mm_mul_pd(xi,xA);
                    Maybe2<add,true>::sse_add(B.get(),_mm_unpacklo_pd(xBr,xBi));
                    Maybe2<add,true>::sse_add(B1.get(),_mm_unpackhi_pd(xBr,xBi));
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) Maybe<add>::add(*B , ZProd<false,false>::prod(x , *A)); 
        }
    };

    // algo 35: double precision SSE2: all complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<35,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<double> >::sametype));
            const bool c1 = V1::_conj;
            if (n) {
                double xr = real(x.x);
                double mxr = Maybe<c1>::select(-xr,xr);
                double xi = imag(x.x);
                double mxi = Maybe<c1>::select(xi,-xi);
                // B = x * A
                // Br = xr * Ar - xi * Ai 
                // Bi = xr * Ai + xi * Ar
                __m128d xxr = _mm_set_pd(mxr, xr);
                __m128d xxi = _mm_set_pd(xi, mxi);
                __m128d xA,xB;
                __m128d x0, x1, x2; // temp values
                do {
                    Maybe<true>::sse_load(xA,A.get()); ++A;
                    x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                    x1 = _mm_mul_pd(xxr,xA);
                    x2 = _mm_mul_pd(xxi,x0);
                    xB = _mm_add_pd(x1,x2);
                    Maybe2<add,true>::sse_add(B.get(),xB); ++B;
                } while (--n);
            }
        }
    };
#endif

#ifdef __SSE__
    // algo 41: single precision SSE: all real
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<41,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    Maybe<add>::add(*B++ , x * *A++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , x * *A++);
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT1 A1 = A+1;
                IT1 A2 = A+2;
                IT1 A3 = A+3;

                IT2 B1 = B+1;
                IT2 B2 = B+2;
                IT2 B3 = B+3;

                __m128 xx = _mm_set1_ps(float(x));
                __m128 xA,xB;
                do {
                    // If !unit2, then we know that A is aligned, so we 
                    // can use load, otherwise use loadu.
                    Maybe2<!unit2,unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xB = _mm_mul_ps(xx,xA);
                    Maybe2<add,unit2>::sse_add(
                        B.get(),B1.get(),B2.get(),B3.get(),xB);
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do {
                Maybe<add>::add(*B++ , x * *A++);
            } while (--nb);
        }
    };

    // algo 42: single precision SSE: x real v1 real v2 complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<42,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , x * *A++);
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT1 A1 = A+1;
                IT1 A2 = A+2;
                IT1 A3 = A+3;

                IT2 B1 = B+1;
                IT2 B2 = B+2;
                IT2 B3 = B+3;

                __m128 xx = _mm_set1_ps(float(x));
                __m128 xzero = _mm_set1_ps(0.F);
                __m128 xA,xB;
                do {
                    Maybe<unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xB = _mm_mul_ps(xx,xA);
                    Maybe2<add,unit2>::sse_addu(
                        B.get(),B1.get(),_mm_unpacklo_ps(xB,xzero));
                    Maybe2<add,unit2>::sse_addu(
                        B2.get(),B3.get(),_mm_unpackhi_ps(xB,xzero));
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do {
                Maybe<add>::add(*B++ , x * *A++);
            } while (--nb);
        }
    };

    // algo 43: single precision SSE: x real v1 complex v2 complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<43,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;
            const bool c1 = V1::_conj;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    Maybe<add>::add(*B++ , ZProd<false,c1>::prod(x , *A++)); 
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , ZProd<false,c1>::prod(x , *A++)); 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                const float mx = Maybe<c1>::select( -float(x) , float(x) );
                // These look backwards, but order is from hi to lo values.
                __m128 xx = _mm_set_ps(mx, float(x), mx, float(x));
                __m128 xA,xB;
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    xB = _mm_mul_ps(xx,xA); 
                    Maybe2<add,unit2>::sse_add(B.get(),B1.get(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) Maybe<add>::add(*B , ZProd<false,c1>::prod(x , *A)); 
        }
    };

    // algo 44: single precision SSE: x complex v1 real v2 complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<44,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<float> >::sametype));
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , ZProd<false,false>::prod(x , *A++)); 
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT1 A1 = A+1;
                IT1 A2 = A+2;
                IT1 A3 = A+3;

                IT2 B1 = B+1;
                IT2 B2 = B+2;
                IT2 B3 = B+3;

                __m128 xr = _mm_set1_ps(real(x.x));
                __m128 xi = _mm_set1_ps(imag(x.x));
                __m128 xA,xBr,xBi;
                do {
                    Maybe<unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xBr = _mm_mul_ps(xr,xA);
                    xBi = _mm_mul_ps(xi,xA);
                    Maybe2<add,unit2>::sse_addu(
                        B.get(),B1.get(),_mm_unpacklo_ps(xBr,xBi));
                    Maybe2<add,unit2>::sse_addu(
                        B2.get(),B3.get(),_mm_unpackhi_ps(xBr,xBi));
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do {
                Maybe<add>::add(*B++ , ZProd<false,false>::prod(x , *A++)); 
            } while (--nb);
        }
    };

    // algo 45: single precision SSE: all complex
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<45,s,add,ix,T,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = s == TMV_UNKNOWN ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<float> >::sametype));
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;
            const bool c1 = V1::_conj;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    Maybe<add>::add(*B++ , ZProd<false,c1>::prod(x , *A++)); 
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    Maybe<add>::add(*B++ , ZProd<false,c1>::prod(x , *A++)); 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                float xr = real(x.x);
                float mxr = Maybe<c1>::select(-xr,xr);
                float xi = imag(x.x);
                float mxi = Maybe<c1>::select(xi,-xi);
                // B = x * A
                // Br = xr * Ar - xi * Ai 
                // Bi = xr * Ai + xi * Ar 
                __m128 xxr = _mm_set_ps(mxr, xr, mxr, xr);
                __m128 xxi = _mm_set_ps(xi, mxi, xi, mxi);
                __m128 xA,xB;
                __m128 x0, x1, x2; // temp values
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    x0 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(2,3,0,1));
                    x1 = _mm_mul_ps(xxr,xA);
                    x2 = _mm_mul_ps(xxi,x0);
                    xB = _mm_add_ps(x1,x2);
                    Maybe2<add,unit2>::sse_add(B.get(),B1.get(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) Maybe<add>::add(*B , ZProd<false,c1>::prod(x , *A)); 
        }
    };
#endif

    // algo 90: Call inst
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<90,s,true,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultXV(xx,v1.xView(),v2.xView()); 
        }
    };
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<90,s,false,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultXV(xx,v1.xView(),v2.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<91,s,true,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultXV(xx,v1.xView(),v2.xView()); 
        }
    };
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<91,s,false,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultXV(xx,v1.xView(),v2.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<97,s,add,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            MultXV_Helper<-2,s,add,ix,T,V1c,V2c>::call(TMV_CONJ(x),v1c,v2c);
        }
    };

    // algo 197: Conjugate
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<197,s,add,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            MultXV_Helper<99,s,add,ix,T,V1c,V2c>::call(TMV_CONJ(x),v1c,v2c);
        }
    };

    // algo 98: Inline check for aliases
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<98,s,false,ix,T,V1,V2>
    {
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            if (!SameStorage(v1,v2)) {
                // No aliasing 
                MultXV_Helper<-2,s,false,ix,T,V1,V2>::call(x,v1,v2);
            } else {
                // Let Copy handle the aliasing
                Copy(v1,v2);
                Scale(x,v2);
            }
        }
    };
    template <int s, int ix, class T, class V1, class V2>
    struct MultXV_Helper<98,s,true,ix,T,V1,V2>
    {
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            if ( !SameStorage(v1,v2) || 
                 ExactSameStorage(v1,v2) || 
                 v1.step()*v2.step() < 0 || 
                 std::abs(v2.step()) < std::abs(v1.step()) ) {
                // No aliasing (or no clobering)
                MultXV_Helper<-2,s,true,ix,T,V1,V2>::call(x,v1,v2);
            } else {
                // Need a temporary
                typename V2::noalias_type v2na = v2.noAlias();
                MultXV<true>(x,v1.copy(),v2na);
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<99,s,add,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                (s == TMV_UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                s == 0 ? 0 : 
                V2::_conj ? 197 :
                inst ? 91 :
                98;
            MultXV_Helper<algo,s,add,ix,T,V1,V2>::call(x,v1,v2);
        }
    };

    // algo -4: No branches or copies.
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<-4,s,add,ix,T,V1,V2>
    {
        typedef typename V1::value_type VT1;
        typedef typename V1::real_type RT1;
        typedef typename V2::value_type VT2;
        typedef typename V2::real_type RT2;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        enum { unit = V1::_step == 1 || V2::_step == 1 };
        enum { allunit = V1::_step == 1 && V2::_step == 1 };
        enum { allfloat = 
            Traits2<RT1,float>::sametype && Traits2<RT2,float>::sametype };
        enum { alldouble = 
            Traits2<RT1,double>::sametype && Traits2<RT2,double>::sametype };
        enum { xreal = Traits<T>::isreal };
        enum { xcomplex = Traits<T>::iscomplex };
        enum { v1real = V1::isreal };
        enum { v2real = V2::isreal };
        enum { v1complex = V1::iscomplex };
        enum { v2complex = V2::iscomplex };
        enum { allcomplex = v1complex && v2complex };
        enum { flatten = 
                xreal && allunit && allcomplex && V1::_conj==int(V2::_conj) };
        enum { algo = (
                s == 0 ? 0 : 
                ( ix == 1 && !add ) ? 1 :
                TMV_OPT == 0 ? 11 :
                flatten ? 2 :
                ( s != TMV_UNKNOWN && s <= int(128/sizeof(VT2)) ) ? 15 :
#if 0
                // The SSE routines don't seem to be faster.
#ifdef __SSE__
                ( allfloat && xreal && v1real && v2real ) ? 41 :
                ( allfloat && xreal && v1real && v2complex ) ? 42 :
                ( allfloat && xreal && v1complex && v2complex ) ? 43 :
                ( allfloat && xcomplex && v1real && v2complex ) ? 44 :
                ( allfloat && xcomplex && v1complex && v2complex ) ? 45 :
#endif
#ifdef __SSE2__
                ( alldouble && xreal && v1real && v2real ) ? 31 :
                ( alldouble && xreal && v1real && v2complex ) ? 32 :
                ( alldouble && xreal && v1complex && v2complex ) ? 33 :
                ( alldouble && xcomplex && v1real && v2complex ) ? 34 :
                ( alldouble && xcomplex && v1complex && v2complex ) ? 35 :
#endif
#endif
                ( allunit && xreal && v1real && sizeof(RT2) == 4 ) ? 13 :
                ( allunit && xreal && v1real && sizeof(RT2) == 8 ) ? 12 :
                11 ) };
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            TMVStaticAssert(!V2::_conj);
#ifdef PRINTALGO_XV
            std::cout<<"InlineMultXV: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"s = "<<s<<" = "<<v2.size()<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
            MultXV_Helper<algo,s,add,ix,T,V1,V2>::call(x,v1,v2); 
        }
        static TMV_INLINE void call2(
            int n, const Scaling<ix,T>& x, const IT1& it1, const IT2& it2)
        {
            TMVStaticAssert(!V2::_conj);
#ifdef PRINTALGO_XV
            std::cout<<"InlineMultXV: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"s = "<<s<<" = "<<n<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
            MultXV_Helper<algo,s,add,ix,T,V1,V2>::call2(n,x,it1,it2); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<-3,s,add,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        { MultXV_Helper<-4,s,add,ix,T,V1,V2>::call(x,v1,v2); }
    };

    // algo -2: Check for inst
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<-2,s,add,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                (s == TMV_UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                s == 0 ? 0 : 
                ( ix == 1 && !add ) ? 201 :
                V2::_conj ? 97 :
                inst ? 90 :
                -4;
            MultXV_Helper<algo,s,add,ix,T,V1,V2>::call(x,v1,v2);
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, int ix, class T, class V1, class V2>
    struct MultXV_Helper<-1,s,add,ix,T,V1,V2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const bool noclobber = VStepHelper<V1,V2>::noclobber;
            const bool checkalias =
                V2::_checkalias && !noclobber;
            const int algo = 
                s == 0 ? 0 : 
                ( ix == 1 && !add ) ? 101 :
                checkalias ? 99 : 
                -2;
            MultXV_Helper<algo,s,add,ix,T,V1,V2>::call(x,v1,v2);
        }
    };

    template <bool add, int ix, class T, class V1, class V2>
    inline void MultXV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        MultXV_Helper<-1,s,add,ix,T,V1v,V2v>::call(x,v1v,v2v);
    }

    template <bool add, int ix, class T, class V1, class V2>
    inline void InlineMultXV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        MultXV_Helper<-4,s,add,ix,T,V1v,V2v>::call(x,v1v,v2v);
    }

    template <bool add, int ix, class T, class V1, class V2>
    inline void InlineAliasMultXV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int s = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        MultXV_Helper<98,s,add,ix,T,V1v,V2v>::call(x,v1v,v2v);
    }

} // namespace tmv

#endif 
