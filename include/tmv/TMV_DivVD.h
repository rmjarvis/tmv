

#ifndef TMV_DivVD_H
#define TMV_DivVD_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_BaseVector.h"
#include "TMV_Array.h"

namespace tmv {

    // Defined in TMV_DivVD.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstElemDivVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasElemDivVV(
        const T3 x, const ConstVectorView<T1,C1>& v1,
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    //
    // Element-wise division:
    // v3(i) = x * v1(i) / v2(i)
    //

    template <int algo, int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper;

    // algo 0: s == 0, nothing to do
    template <int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<0,0,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static TMV_INLINE void call(
            const Scaling<ix,T>&, const V1&, const V2&, V3&) {}
        static TMV_INLINE void call2(
            int, const Scaling<ix,T>&, IT1, IT2, IT3) {}
    };

    // algo 11: simple for loop
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<11,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v3.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 11: N,s = "<<n<<','<<s<<std::endl;
            std::cout<<"c1,c2 = "<<c1<<','<<c2<<std::endl;
#endif
            if (n) do {
                *C++ = ZProd<false,false>::prod(
                    x, ZProd<c1,c2>::quot(*A++,*B++));
            } while (--n);
        }
    };  

#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<31,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 31: N,s = "<<n<<','<<s<<std::endl;
#endif
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    *B++ = x / *A++;
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = x / *A++;
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
                    xB = _mm_div_pd(xx,xA);
                    Maybe<unit2>::sse_store(B.get(),B1.get(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = x / *A;
        }
    };

    // algo 32: double precision SSE2: x real v1 real v2 complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<32,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 32: N,s = "<<n<<','<<s<<std::endl;
#endif
            const bool unit1 = V1::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = x / *A++;
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
                    xB = _mm_div_pd(xx,xA);
                    Maybe<true>::sse_storeu(
                        B.get(),_mm_unpacklo_pd(xB,xzero));
                    Maybe<true>::sse_storeu(
                        B1.get(),_mm_unpackhi_pd(xB,xzero));
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = x / *A;
        }
    };

    // algo 33: double precision SSE2: x real v1 complex v2 complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<33,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 33: N,s = "<<n<<','<<s<<std::endl;
#endif
            const bool c1 = V1::_conj;
            if (n) {
                const double mone = Maybe<c1>::select( double(x) , -double(x) );
                // These look backwards, but order is from hi to lo values.
                __m128d xconj = _mm_set_pd(mone, double(x));
                __m128d xA,xB;
                __m128d xAc, xnorm, x1, x2; // temp values
                if ( TMV_Aligned(A.get()) && TMV_Aligned(B.get()) ) {
                    do {
                        Maybe<true>::sse_load(xA,A.get()); ++A;
                        xAc = _mm_mul_pd(xconj,xA); // x*conj(xA)
                        x1 = _mm_mul_pd(xA,xA);
                        x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                        xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                        xB = _mm_div_pd(xAc,xnorm);  // = x/xA
                        Maybe<true>::sse_store(B.get(),xB); ++B;
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.get()); ++A;
                        xAc = _mm_mul_pd(xconj,xA); // x*conj(xA)
                        x1 = _mm_mul_pd(xA,xA);
                        x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                        xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                        xB = _mm_div_pd(xAc,xnorm);  // = x/xA
                        Maybe<true>::sse_storeu(B.get(),xB); ++B;
                    } while (--n);
                }
            }
        }
    };

    // algo 34: double precision SSE2: x complex v1 real v2 complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<34,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 34: N,s = "<<n<<','<<s<<std::endl;
#endif
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<double> >::sametype));
            const bool unit1 = V1::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xone = _mm_set1_pd(1.);
                __m128d xr = _mm_set1_pd(real(x.x));
                __m128d xi = _mm_set1_pd(imag(x.x));
                __m128d xA,xBr,xBi,xAinv;
                do {
                    Maybe<unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    xAinv = _mm_div_pd(xone,xA);
                    xBr = _mm_mul_pd(xr,xAinv);
                    xBi = _mm_mul_pd(xi,xAinv);
                    Maybe<true>::sse_store(B.get(),_mm_unpacklo_pd(xBr,xBi));
                    Maybe<true>::sse_store(B1.get(),_mm_unpackhi_pd(xBr,xBi));
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = x / *A;
        }
    };

    // algo 35: double precision SSE2: all complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<35,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 35: N,s = "<<n<<','<<s<<std::endl;
#endif
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<double> >::sametype));
            const bool c1 = V1::_conj;
            if (n) {
                double xr = real(x.x);
                double mxr = Maybe<c1>::select(xr,-xr);
                double xi = imag(x.x);
                double mxi = Maybe<c1>::select(-xi,xi);
                // B = x * conj(A) / norm(A)
                // Br = xr * Ar + xi * Ai / norm
                // Bi = -xr * Ai + xi * Ar / norm
                __m128d xxr = _mm_set_pd(mxr, xr);
                __m128d xxi = _mm_set_pd(xi, mxi);
                __m128d xA,xB;
                __m128d xnorm, x0, x1, x2, x3, x4, x5; // temp values
                do {
                    Maybe<true>::sse_load(xA,A.get()); ++A;
                    x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                    x1 = _mm_mul_pd(xA,xA);
                    x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                    xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                    x3 = _mm_mul_pd(xxr,xA);
                    x4 = _mm_mul_pd(xxi,x0);
                    x5 = _mm_add_pd(x3,x4);
                    xB = _mm_div_pd(x5,xnorm);  // = x/xA
                    Maybe<true>::sse_store(B.get(),xB); ++B;
                } while (--n);
            }
        }
    };
#endif

#ifdef __SSE__
    // algo 41: single precision SSE: all real
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<41,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 41: N,s = "<<n<<','<<s<<std::endl;
#endif
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;
            const bool unit3 = V2::_step == 1;

            if (unit3) {
                while (n && !TMV_Aligned(C.get()) ) {
                    *C++ = x * *A++ / *B++;
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *C++ = x * *A++ / *B++;
                    --n;
                }
            } else if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    *C++ = x * *A++ / *B++;
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

                IT3 C1 = C+1;
                IT3 C2 = C+2;
                IT3 C3 = C+3;

                __m128 xx = _mm_set1_ps(float(x));
                __m128 xA,xB,xC;
                do {
                    Maybe2<!unit3,unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    Maybe2<!(unit3||unit1),unit2>::sse_load(
                        xB,B.get(),B1.get(),B2.get(),B3.get());
                    B+=4; B1+=4; B2+=4; B3+=4;
                    xC = _mm_div_ps(xA,xB);
                    Maybe<ix!=1>::sse_mult(xx,xC);
                    Maybe<unit2>::sse_store(
                        C.get(),C1.get(),C2.get(),C3.get(),xC);
                    C+=4; C1+=4; C2+=4; C3+=4;
                } while (--n_4);
            }

            if (nb) do { *C++ = x * *A++ / *B++; } while (--nb);
        }
    };

    // algo 42: single precision SSE: x real v1 real v2 complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<42,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 42: N,s = "<<n<<','<<s<<std::endl;
#endif
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = x / *A++;
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
                    xB = _mm_div_ps(xx,xA);
                    Maybe<unit2>::sse_storeu(
                        B.get(),B1.get(),_mm_unpacklo_ps(xB,xzero));
                    Maybe<unit2>::sse_storeu(
                        B2.get(),B3.get(),_mm_unpackhi_ps(xB,xzero));
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do { *B++ = x / *A++; } while (--nb);
        }
    };

    // algo 43: single precision SSE: x real v1 complex v2 complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<43,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 43: N,s = "<<n<<','<<s<<std::endl;
#endif
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;
            const bool c1 = V1::_conj;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                const float mx = Maybe<c1>::select( float(x) , -float(x) );
                // These look backwards, but order is from hi to lo values.
                __m128 xconj = _mm_set_ps(mx, float(x), mx, float(x));
                __m128 xA,xB;
                __m128 xAc, xnorm, x1, x2; // temp values
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    xAc = _mm_mul_ps(xconj,xA); // conj(xA)
                    x1 = _mm_mul_ps(xA,xA);
                    x2 = _mm_shuffle_ps(x1,x1,_MM_SHUFFLE(2,3,0,1));
                    xnorm = _mm_add_ps(x1,x2); // = norm(xA)
                    xB = _mm_div_ps(xAc,xnorm);  // = x/xA
                    Maybe<unit2>::sse_store(B.get(),B1.get(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = ZProd<false,c1>::quot(x , *A);
        }
    };

    // algo 44: single precision SSE: x complex v1 real v2 complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<44,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 44: N,s = "<<n<<','<<s<<std::endl;
#endif
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<float> >::sametype));
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = x / *A++;
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

                __m128 xone = _mm_set1_ps(1.F);
                __m128 xr = _mm_set1_ps(real(x.x));
                __m128 xi = _mm_set1_ps(imag(x.x));
                __m128 xA,xBr,xBi,xAinv;
                do {
                    Maybe<unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xAinv = _mm_div_ps(xone,xA);
                    xBr = _mm_mul_ps(xr,xAinv);
                    xBi = _mm_mul_ps(xi,xAinv);
                    Maybe<unit2>::sse_storeu(
                        B.get(),B1.get(),_mm_unpacklo_ps(xBr,xBi));
                    Maybe<unit2>::sse_storeu(
                        B2.get(),B3.get(),_mm_unpackhi_ps(xBr,xBi));
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do { *B++ = x / *A++; } while (--nb);
        }
    };

    // algo 45: single precision SSE: all complex
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<45,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = s == Unknown ? v2.size() : s;
            call2(n,x,v1.begin().nonConj(),v2.begin().nonConj(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
#ifdef PRINTALGO_DivD
            std::cout<<"ElemDivVV algo 45: N,s = "<<n<<','<<s<<std::endl;
#endif
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<float> >::sametype));
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;
            const bool c1 = V1::_conj;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                float xr = real(x.x);
                float mxr = Maybe<c1>::select(xr,-xr);
                float xi = imag(x.x);
                float mxi = Maybe<c1>::select(-xi,xi);
                // B = x * conj(A) / norm(A)
                // Br = xr * Ar + xi * Ai / norm
                // Bi = -xr * Ai + xi * Ar / norm
                __m128 xxr = _mm_set_ps(mxr, xr, mxr, xr);
                __m128 xxi = _mm_set_ps(xi, mxi, xi, mxi);
                __m128 xA,xB;
                __m128 xnorm, x0, x1, x2, x3, x4, x5; // temp values
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    x0 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(2,3,0,1));
                    x1 = _mm_mul_ps(xA,xA);
                    // It doesn't seem to matter which way we calculate x2.
                    // They take almost identical time.
                    //x2 = _mm_mul_ps(x0,x0);
                    x2 = _mm_shuffle_ps(x1,x1,_MM_SHUFFLE(2,3,0,1));
                    xnorm = _mm_add_ps(x1,x2); // = norm(xA)
                    x3 = _mm_mul_ps(xxr,xA);
                    x4 = _mm_mul_ps(xxi,x0);
                    x5 = _mm_add_ps(x3,x4);
                    xB = _mm_div_ps(x5,xnorm);  // = x/xA
                    Maybe<unit2>::sse_store(B.get(),B1.get(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = ZProd<false,c1>::quot(x , *A);
        }
    };
#endif

    // algo 90: Call inst
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<90,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstElemDivVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<91,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasElemDivVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<97,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            ElemDivVV_Helper<-2,s,ix,T,V1c,V2c,V3c>::call(
                TMV_CONJ(x),v1c,v2c,v3c);
        }
    };

    // algo 197: Conjugate
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<197,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            ElemDivVV_Helper<99,s,ix,T,V1c,V2c,V3c>::call(
                TMV_CONJ(x),v1c,v2c,v3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<98,s,ix,T,V1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber1 =
                !SameStorage(v1,v3) ||
                ExactSameStorage(v1,v3) || 
                v1.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v1.step());
            const bool noclobber2 =
                !SameStorage(v2,v3) ||
                ExactSameStorage(v2,v3) || 
                v2.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v2.step());
            if (noclobber1) {
                if (noclobber2) {
                    // No aliasing (or no clobering)
                    ElemDivVV_Helper<-2,s,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
                } else {
                    // Need a temporary for v2
                    typename V3::noalias_type v3na = v3.noAlias();
                    ElemDivVV(x,v1,v2.copy(),v3na);
                }
            } else {
                if (noclobber2) {
                    // Need a temporary for v1
                    typename V3::noalias_type v3na = v3.noAlias();
                    ElemDivVV(x,v1.copy(),v2,v3na);
                } else {
                    // Need a temporary for v3
                    typename V3::copy_type v3c(v3.size());
                    ElemDivVV(x,v1,v2,v3c);
                    typename V3::noalias_type v3na = v3.noAlias();
                    Copy(v3c,v3na);
                }
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<99,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst =
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                V3::_conj ? 197 : 
                inst ? 91 : 
                98;
            ElemDivVV_Helper<algo,s,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    // algo -4: No branches or copies
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-4,s,ix,T,V1,V2,V3>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        enum { algo = (
                // SSE stuff doesn't seem to be faster.  Huh.
                11 ) };
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
#ifdef PRINTALGO_DivD
            std::cout<<"Inline ElemDivVV: s = "<<s<<std::endl;
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<std::endl;
#endif
            ElemDivVV_Helper<algo,s,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
        static TMV_INLINE void call2(
            int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            TMVStaticAssert(!V3::_conj);
            ElemDivVV_Helper<algo,s,ix,T,V1,V2,V3>::call2(n,x,A,B,C);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-3,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            ElemDivVV_Helper<-4,s,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    // algo -2: Check for inst
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-2,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst =
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                V3::_conj ? 97 : 
                inst ? 90 : 
                -3;
            ElemDivVV_Helper<algo,s,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    // algo -1: Check for aliases?
    template <int s, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-1,s,ix,T,V1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber =
                VStepHelper<V1,V3>::noclobber &&
                VStepHelper<V2,V3>::noclobber;
            const bool checkalias =
                V3::_checkalias && !noclobber;
            const int algo =
                checkalias ? 99 : 
                -2;
            ElemDivVV_Helper<algo,s,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    template <int ix, class T, class V1, class V2, class V3>
    inline void ElemDivVV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemDivVV_Helper<-1,s,ix,T,V1v,V2v,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void InlineElemDivVV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemDivVV_Helper<-3,s,ix,T,V1v,V2v,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void InlineAliasElemDivVV(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::_size,V2::_size>::size;
        const int s = Sizes<s12,V3::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        ElemDivVV_Helper<98,s,ix,T,V1v,V2v,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    static void LDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        if (m2.isSingular()) ThrowSingular("DiagMatrix");
        ElemDivVV(x,v1,m2.diag(),v3);
    }

} // namespace tmv

#endif 
