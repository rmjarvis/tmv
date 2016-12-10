

#ifndef TMV_MultUV_H
#define TMV_MultUV_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseVector.h"
#include "TMV_MultVV.h"
#include "TMV_MultXV.h"
#include "TMV_Prefetch.h"
#include "TMV_MultMV_Funcs.h"

#ifdef PRINTALGO_UV
#include <iostream>
#endif

namespace tmv {

    // Defined in TMV_MultUV.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMV(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMV(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMV(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMV(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMV(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMV(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMV(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    //
    // Matrix * Vector
    //

    // UNROLL is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_UV_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_UV_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_UV_UNROLL 9
#else
#define TMV_UV_UNROLL 0
#endif

    // The minimum size to copy a vector if its step == Unknown.
#define TMV_UV_COPYSIZE 4

    // The crossover memory size to start using prefetch commands.
    // This is undoubtedly a function of the L1 (and L2?) cache size,
    // but 2KBytes is probably not too bad for most machines.
    // (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_UV_PREFETCH 2048


    // Note: all algorithms here are designed to work if v2 and v3 are 
    // the same storage.  So it works for things like v *= U without 
    // any need for a temporary copy.
    template <int algo, ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper;

    // algo 0: s = 0, so nothing to do
    template <bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<0,0,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const M1& , const V2& , V3& ) 
        {} 
    };

    // algo 1: s == 1, so simplifies to a scalar product
    template <bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<1,1,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 1: N,s,x = "<<1<<','<<1<<','<<T(x)<<std::endl;
#endif
            Maybe<add>::add( 
                v3.ref(0) , 
                ZProd<false,false>::prod(
                    x, Maybe<!M1::_unit>::prod(m1.cref(0,0) , v2.cref(0)))); 
        }
    };

    // algo 11: The basic column major loop for UpperTri
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<11,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            ptrdiff_t N = (s == Unknown ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1::const_diag_type M1d;
            PT2 Xj;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M1d::const_nonconj_type::const_iterator IT1d;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;
            IT1 A0j = m1.get_col(0,0,1).begin().nonConj();
            IT2 X = v2.begin().nonConj();
            IT3 Y0 = v3.begin();
            IT3 Y = Y0;
            const ptrdiff_t Astepj = m1.stepj();

            const bool dopref = N * sizeof(T1) >= TMV_UV_PREFETCH;

            Prefetch_Read(m1.cptr());
            Prefetch_Read(v2.cptr());
            Prefetch_MultiWrite(v3.ptr());

            for(ptrdiff_t j=0;N--;++j) {
                // loop from j = 0 .. N-1
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X++);
                    // y.subVector(0,j) += x(j) * A.col(j,0,j);
                    MultXV_Helper<-4,Unknown,true,0,PT2,M1c,V3>::call2(
                        j,Scaling<0,PT2>(Xj),A0j,Y0);
                    A0j.shiftP(Astepj);
                    if (dopref) Prefetch_Read(A0j.get());
                    Maybe<add>::add(*Y++,Maybe<!unit>::prod(m1.cref(j,j),Xj));
                } else {
                    A0j.shiftP(Astepj);
                    if (dopref) Prefetch_Read(A0j.get());
                    ++X; 
                    Maybe<!add>::set(*Y++,T3(0));
                }
            }
        }
    };

    // algo 12: The basic row major loop for UpperTri
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<12,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            ptrdiff_t N = (s == Unknown ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT;
            typedef typename M1::const_row_sub_type M1r;
            PT Yi;

            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            // Actually Aii is the address of A(i,i+1)
            IT1 Aii = m1.get_row(0,1,N).begin().nonConj();
            IT2 X = v2.begin().nonConj();
            IT3 Y = v3.begin();
            const ptrdiff_t Adiagstep = m1.diagstep();
            const ptrdiff_t Astepi = m1.stepi();

            const bool dopref = N * sizeof(T1) >= TMV_UV_PREFETCH;

            Prefetch_Read(m1.cptr());
            Prefetch_MultiRead(v2.cptr());
            Prefetch_Write(v3.ptr());

            // [ A00 A01 A02 ] [  0 ]   [ A01 B2 ]
            // [  0  A11 A12 ] [ B2 ] = [ A11 B2 ]
            // [  0   0  A22 ] [  0 ]   [    0   ]
            while (N > 0 && v2(N-1) == T2(0)) Maybe<!add>::set(v3(--N),T3(0));
            ptrdiff_t i1=0;
            while (N > 0 && *X == T2(0)) { ++i1; --N; ++X; ++Aii; }
            for(ptrdiff_t i=0;i<i1;++i) {
                Yi = MultVV_Helper<-4,Unknown,M1r,V2>::call2(N,Aii-1,X);
                Aii.shiftP(Astepi);
                if (dopref) Prefetch_Read(Aii.get()-1);
                Maybe<add>::add(*Y++, x * Yi);
            }

            for(ptrdiff_t i=i1;N--;++i) {
                // loop from i = 0 .. N-1
                // Yi = A.row(i,i,N) * x.subVector(i,N)
                //    = A(i,i) * x(i)
                //      + A.row(i,i+1,N) * x.subVector(i+1,N)
                Yi = Maybe<!unit>::template zprod<false,c2>(m1.cref(i,i),*X++);
                Yi += MultVV_Helper<-4,Unknown,M1r,V2>::call2(N,Aii,X);
                Aii.shiftP(Adiagstep);
                if (dopref) Prefetch_Read(Aii.get()-1);
                Maybe<add>::add(*Y++, x * Yi);
            }
        }
    };

    // algo 15: UpperTri: unroll by rows
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<15,s,add,ix,T,M1,V2,V3>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                Unroller<I,N/2>::unroll(x,m1,v2,v3);
                Unroller<I+N/2,N-N/2>::unroll(x,m1,v2,v3);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                typedef typename M1::const_row_sub_type M1r;
                typedef typename V2::const_subvector_type V2s;
                const bool unit = M1::_unit;
                Maybe<add>::add(
                    v3.ref(I) , 
                    ZProd<false,false>::prod(
                        x , Maybe<!unit>::prod(m1.cref(I,I) , v2.cref(I)) +
                        MultVV_Helper<-4,s-I-1,M1r,V2s>::call(
                            m1.get_row(I,I+1,s),v2.cSubVector(I+1,s))));
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& , const M1& , const V2& , V3& ) {}
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 15: N,s,x = "<<s<<','<<s<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,v2,v3); 
        }
    };

    // algo 21: The basic column major loop for LowerTri
    // This is designed to work if v2 and v3 are the same storage.
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<21,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const ptrdiff_t N = (s == Unknown ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T,typename V2::value_type>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1::const_diag_type M1d;
            PT2 Xj;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M1d::const_nonconj_type::const_iterator IT1d;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;

            // Actually Ajj is the address of A(j,j+1)
            IT1 Ajj = m1.get_col(N-1,N,N).begin().nonConj();
            IT2 X = v2.begin().nonConj() + N-1;
            IT3 Y = v3.begin() + N;
            const ptrdiff_t Adiagstep = m1.diagstep();

            const bool dopref = N * Adiagstep * sizeof(T1) >= TMV_UV_PREFETCH;

            Prefetch_Read(v2.cptr());
            Prefetch_MultiWrite(v3.ptr());
            if (dopref) Prefetch_Read(Ajj.get());
            else Prefetch_Read(m1.cptr());

            for(ptrdiff_t j=N,len=0;j--;++len) {
                // loop from j = N-1 .. 0
                if (*X != T2(0)) {
                    Xj = ZProd<false,c2>::prod(x , *X--);
                    // y.subVector(j+1,N) += x(j) * A.col(j,j+1,N);
                    MultXV_Helper<-4,Unknown,true,0,PT2,M1c,V3>::call2(
                        len,Scaling<0,PT2>(Xj),Ajj,Y--);
                    Ajj.shiftP(-Adiagstep);
                    if (dopref) Prefetch_Read(Ajj.get()-1);
                    // y(j) (+)= x(j) * A(j,j);
                    Maybe<add>::add(*Y , Maybe<!unit>::prod(m1.cref(j,j),Xj));
                } else {
                    Ajj.shiftP(-Adiagstep);
                    if (dopref) Prefetch_Read(Ajj.get()-1);
                    --X;
                    Maybe<!add>::set(*--Y,T3(0));
                }
            }
        }
    };

    // algo 22: The basic row major loop for LowerTri
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<22,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            ptrdiff_t N = (s == Unknown ? m1.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            typedef typename Traits2<T1,T2>::type PT;
            typedef typename M1::const_row_sub_type M1r;
            PT Yi;

            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename V2::const_nonconj_type::const_iterator IT2;
            typedef typename V3::iterator IT3;

            const bool unit = M1::_unit;
            const bool c2 = V2::_conj;

            IT1 Ai0 = m1.get_row(N-1,0,N-1).begin().nonConj();
            IT2 X0 = v2.begin().nonConj();
            IT2 X = X0 + N-1;
            IT3 Y = v3.begin() + N-1;
            const ptrdiff_t Astepi = m1.stepi();

            const bool dopref = 
                N * m1.diagstep() * sizeof(T1) >= TMV_UV_PREFETCH;

            Prefetch_MultiRead(v2.cptr());
            Prefetch_Write(v3.ptr());
            if (dopref) Prefetch_Read(Ai0.get());
            else Prefetch_Read(m1.cptr());

            // [ A00  0   0  ] [  0 ]   [    0   ]
            // [ A10 A11  0  ] [ B2 ] = [ A11 B2 ]
            // [ A20 A21 A22 ] [  0 ]   [ A21 B2 ]
            ptrdiff_t i=N-1;
            ptrdiff_t i1=0;
            while (N > 0 && *X0 == T2(0)) {
                Maybe<!add>::set(v3(i1++),T3(0));
                --N; ++Ai0; ++X0;
            }
            ptrdiff_t i2=i;
            while (N > 0 && *X == T2(0)) { --i2; --N; --X; }
            for(;i>i2;--i) {
                Yi = MultVV_Helper<-4,Unknown,M1r,V2>::call2(N,Ai0,X0);
                Ai0.shiftP(-Astepi);
                if (dopref) Prefetch_Read(Ai0.get());
                Maybe<add>::add(*Y--, x * Yi);
            }

            for(;N--;--i) {
                // loop from i = N-1 .. 0
                // Yi = A.row(i,0,i+1) * x.subVector(0,i+1)
                //    = A(i,i) * x(i)
                //      + A.row(i,0,i) * x.subVector(0,i)
                Yi = MultVV_Helper<-4,Unknown,M1r,V2>::call2(N,Ai0,X0);
                Ai0.shiftP(-Astepi);
                if (dopref) Prefetch_Read(Ai0.get());
                Yi += Maybe<!unit>::template zprod<false,c2>(m1.cref(i,i),*X--);
                Maybe<add>::add(*Y--, x * Yi);
            } 
        }
    };

    // algo 25: LowerTri: unroll by rows
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<25,s,add,ix,T,M1,V2,V3>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                Unroller<I+N/2,N-N/2>::unroll(x,m1,v2,v3);
                Unroller<I,N/2>::unroll(x,m1,v2,v3);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
            {
                typedef typename M1::const_row_sub_type M1r;
                typedef typename V2::const_subvector_type V2s;
                const bool unit = M1::_unit;
                Maybe<add>::add(
                    v3.ref(I) , 
                    ZProd<false,false>::prod(
                        x ,  Maybe<!unit>::prod(m1.cref(I,I) , v2.cref(I)) +
                        MultVV_Helper<-4,I,M1r,V2s>::call(
                            m1.get_row(I,0,I),v2.cSubVector(0,I))));
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix,T>& , const M1& , const V2& , V3& ) {}
        };

        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 25: N,s,x = "<<s<<','<<s<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,v2,v3); 
        }
    };

    // algo 43: colmajor, v3.step == Unknown, so maybe copy v3
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<43,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const ptrdiff_t N = s == Unknown ? m1.size() : s;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 43: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            if (N > TMV_UV_COPYSIZE) {
                MultUV_Helper<85,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else 
                MultUV_Helper<-4,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 53: rowmajor, v2.step == Unknown, so maybe copy v2
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<53,s,add,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const ptrdiff_t N = s == Unknown ? m1.size() : s;
#ifdef PRINTALGO_UV
            std::cout<<"UV algo 53: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            if (N > TMV_UV_COPYSIZE) {
                MultUV_Helper<81,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else 
                MultUV_Helper<-4,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo 81: copy v2
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<81,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            const ptrdiff_t N = s == Unknown ? m1.size() : s;
            std::cout<<"UV algo 81: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typename V3::noalias_type v3na = v3.noAlias();
            MultMV<add>(x,m1,v2.copy(),v3na);
        }
    };

    template <int ix, class T, class M, class V> class ProdMV;

    // algo 85: v3c = x*m1*v2, v3 (+)= v3c
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<85,s,add,ix,T,M1,V2,V3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
#ifdef PRINTALGO_UV
            const ptrdiff_t N = s == Unknown ? m1.size() : s;
            std::cout<<"UV algo 85: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typename V3::noalias_type v3na = v3.noAlias();
            MultXV<add>(one,ProdMV<ix,T,M1,V2>(x,m1,v2).calc(),v3na);
        }
    };

    // algo 90: call inst
    template <ptrdiff_t s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<90,s,false,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };
    template <ptrdiff_t s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<90,s,true,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };

    // algo 91: call inst alias
    template <ptrdiff_t s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<91,s,false,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };
    template <ptrdiff_t s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<91,s,true,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename V3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMV(xx,m1.xView(),v2.xView(),v3.xView());
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<97,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultUV_Helper<-2,s,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<197,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            M1c m1c = m1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            MultUV_Helper<99,s,add,ix,T,M1c,V2c,V3c>::call(
                TMV_CONJ(x),m1c,v2c,v3c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<98,s,true,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            if ( !SameStorage(m1,v3) &&
                 ( (M1::_upper && v2.step()>=v3.step()) ||
                   (M1::_lower && v2.step()<=v3.step()) ||
                   !SameStorage(v2.vec(),v3.vec())) ) {
                // No aliasing (or no clobbering)
                MultUV_Helper<-2,s,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else if (SameStorage(m1,v3)) {
                // Use temporary for v3
                MultUV_Helper<85,s,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                // SameStorage(v2,v3)
                // Use temporary for v2
                MultUV_Helper<81,s,true,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            }
        }
    };
    // if !add, then don't need the temporary if v2,v3 are aliased.
    template <ptrdiff_t s, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<98,s,false,ix,T,M1,V2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            if ( !SameStorage(m1,v3) &&
                 ( (M1::_upper && v2.step()>=v3.step()) ||
                   (M1::_lower && v2.step()<=v3.step()) ||
                   !SameStorage(v2.vec(),v3.vec())) ) {
                // No aliasing (or no clobbering)
                MultUV_Helper<-2,s,false,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else if (SameStorage(m1,v3)) {
                // Use temporary for v3
                MultUV_Helper<85,s,false,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
            } else {
                Copy(v2,v3);
                typename V3::noalias_type v3na = v3.noAlias();
                MultMV<false>(x,m1,v3,v3na);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<99,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
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
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                V3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -4: No branches or copies
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-4,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            const ptrdiff_t s2 = s > 20 ? Unknown : s;
            const ptrdiff_t s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const ptrdiff_t nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_UV_UNROLL;
            const int algo = 
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                M1::_upper ? (
                    unroll ? 15 :
                    M1::_colmajor ? 11 :
                    M1::_rowmajor ? 12 :
                    V2::_step == 1 ? 12 : V3::_step == 1 ? 11 : 12 ) :
                ( // lowertri
                    unroll ? 25 :
                    M1::_colmajor ? 21 :
                    M1::_rowmajor ? 22 :
                    V2::_step == 1 ? 22 : V3::_step == 1 ? 21 : 22 );
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-3,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::_conj);
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar product
            //
            // UpperTri:
            // 11 = column major, simple for loop
            // 12 = row major, simple for loop
            // 15 = fully unroll by rows
            //
            // LowerTri:
            // 21 = column major, simple for loop
            // 22 = row major, simple for loop
            // 25 = fully unroll by rows
            //
            // Copy a vector to new storage:
            // 81 = copy v2
            // 85 = temp v3 = x*m1*v2

            const ptrdiff_t s2 = s > 20 ? Unknown : s;
            const ptrdiff_t s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const ptrdiff_t nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_UV_UNROLL;
            const int algo = 
                ( s == 0 ) ? 0 : // trivial - nothing to do
                ( s == 1 ) ? 1 : // trivial - s = 1

                M1::_upper ? 
                TMV_OPT == 0 ? ( M1::_colmajor ? 11 : 12 ) :
                unroll ? 15 :
                M1::_colmajor ? (
                    V3::_step == Unknown ? (
#ifdef TMV_OPT_SCALE
                        s == Unknown ? 43 :
#endif
                        s == Unknown ? 85 :
                        s > TMV_UV_COPYSIZE ? 85 : 11 ) :
                    11 ) :
                M1::_rowmajor ? (
                    V2::_step == Unknown ? (
#ifdef TMV_OPT_SCALE
                        s == Unknown ? 53 :
#endif
                        s == Unknown ? 81 :
                        s > TMV_UV_COPYSIZE ? 81 : 12 ) :
                    12 ) :
                V2::_step == 1 ? 12 : V3::_step == 1 ? 11 : 12 :

                // lowertri
                TMV_OPT == 0 ? ( M1::_colmajor ? 21 : 22 ) :
                unroll ? 25 :
                M1::_colmajor ? (
                    V3::_step == Unknown ? (
#ifdef TMV_OPT_SCALE
                        s == Unknown ? 43 :
#endif
                        s == Unknown ? 85 :
                        s > TMV_UV_COPYSIZE ? 85 : 21 ) :
                    21 ) :
                M1::_rowmajor ? (
                    V2::_step == Unknown ? (
#ifdef TMV_OPT_SCALE
                        s == Unknown ? 53 :
#endif
                        s == Unknown ? 81 :
                        s > TMV_UV_COPYSIZE ? 81 : 22 ) :
                    22 ) :
                V2::_step == 1 ? 22 : V3::_step == 1 ? 21 : 22;
#ifdef PRINTALGO_UV
            std::cout<<"InlineMultUV: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"v3 = "<<TMV_Text(v3)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
#ifdef XDEBUG_MV
            typedef typename V3::real_type RT;
            typedef typename V3::value_type T3;
            Matrix<T3> m1c = m1;
            Vector<T3> v2c = v2;
            Vector<T3> v3i = v3;
            Vector<T3> v3c = v3;
            MultMV<add>(x,m1c,v2c,v3c);
#endif
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
#ifdef XDEBUG_MV
            if (Norm(v3-v3c) > 1.e-3*(Norm(m1c)*Norm(v2c)+(add?Norm(v3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"v2 = "<<v2c<<std::endl;
                std::cout<<"v3 = "<<v3i<<std::endl;
                std::cout<<"v3 => "<<v3<<std::endl;
                std::cout<<"Correct v3 = "<<v3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-2,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            typedef typename M1::value_type T1;
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
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                V3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t s, bool add, int ix, class T, class M1, class V2, class V3>
    struct MultUV_Helper<-1,s,add,ix,T,M1,V2,V3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const V2& v2, V3& v3)
        {
            const int algo = 
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                V3::_checkalias ? 99 : 
                -2;
            MultUV_Helper<algo,s,add,ix,T,M1,V2,V3>::call(x,m1,v2,v3);
        }
    };

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void MultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V2::_size,V3::_size>::same));
        TMVAssert(m1.size() == v3.size());
        TMVAssert(m1.size() == v2.size());
        TMVAssert(v2.size() == v3.size());
        const ptrdiff_t s = Sizes<Sizes<M1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultUV_Helper<-1,s,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V2::_size,V3::_size>::same));
        TMVAssert(m1.size() == v3.size());
        TMVAssert(m1.size() == v2.size());
        TMVAssert(v2.size() == v3.size());
        const ptrdiff_t s = Sizes<Sizes<M1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultUV_Helper<-3,s,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <bool add, int ix, class T, class M1, class V2, class V3>
    inline void InlineAliasMultMV(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1::_size,V3::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,V2::_size>::same));
        TMVStaticAssert((Sizes<V2::_size,V3::_size>::same));
        TMVAssert(m1.size() == v3.size());
        TMVAssert(m1.size() == v2.size());
        TMVAssert(v2.size() == v3.size());
        const ptrdiff_t s = Sizes<Sizes<M1::_size,V2::_size>::size,V3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        MultUV_Helper<98,s,add,ix,T,M1v,V2v,V3v>::call(x,m1v,v2v,v3v);
    }

    template <class V1, int ix, class T, class M2>
    TMV_INLINE void MultEqVM(
        BaseVector_Mutable<V1>& v1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { InlineMultMV<false>(x,m2.transpose(),v1.vec(),v1.vec()); }

} // namespace tmv

#undef TMV_UV_UNROLL
#undef TMV_UV_COPYSIZE
#undef TMV_UV_PREFETCH

#endif 
