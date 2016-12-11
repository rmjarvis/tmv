///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_SmallVectorArithFunc_H
#define TMV_SmallVectorArithFunc_H

#define CT std::complex<T>

namespace tmv {

    template <typename T1, typename T2>
    struct OKTypeHelper { enum { ok=true }; };
    template <typename T>
    struct OKTypeHelper<CT,T> { enum { ok=false }; };

    template <typename T1, typename T2, typename T3>
    struct OKTypeHelper2 { enum { ok=true }; };
    template <typename T>
    struct OKTypeHelper2<T,CT,T> { enum { ok=false }; };
    template <typename T>
    struct OKTypeHelper2<CT,T,T> { enum { ok=false }; };
    template <typename T>
    struct OKTypeHelper2<CT,CT,T> { enum { ok=false }; };

#define OKTypes(T1,T2) (OKTypeHelper<T1,T2>::ok)
#define OKTypes2(T1,T2,T3) (OKTypeHelper2<T1,T2,T3>::ok)


    // v *= x
    template <ptrdiff_t N, typename T, typename Tv, bool ok>
    struct MultXVHelper
    { MultXVHelper(const T, Tv*) {} };
    template <ptrdiff_t N, typename T, typename Tv>
    struct MultXVHelper<N,T,Tv,true>
    {
        MultXVHelper(const T x, Tv* v)
        { for(ptrdiff_t i=0;i<N;++i) v[i] *= x; }
    };
    template <ptrdiff_t N, typename T, typename Tv>
    inline void MultXV(const T x, Tv* v)
    { MultXVHelper<N,T,Tv,OKTypes(T,Tv)>(x,v); }

    // v2 = x * v1
    template <ptrdiff_t N, typename T, typename T1, typename T2, bool ok>
    struct MultXVHelper2
    { MultXVHelper2(const T, const T1*, T2*) {} };
    template <ptrdiff_t N, typename T, typename T1, typename T2>
    struct MultXVHelper2<N,T,T1,T2,true>
    {
        MultXVHelper2(const T x, const T1* v1, T2* v2)
        { for(ptrdiff_t i=0;i<N;++i) v2[i] = x * v1[i]; }
    };
    template <ptrdiff_t N, typename T, typename T1, typename T2>
    inline void MultXV( const T x, const T1* v1, T2* v2)
    { MultXVHelper2<N,T,T1,T2,OKTypes2(T,T1,T2)>(x,v1,v2); }

    // v2 += x * v1
    template <ptrdiff_t N, typename T1, typename T2, bool ok>
    struct AddVV_1Helper
    { AddVV_1Helper(const T1*, T2*) {} };
    template <ptrdiff_t N, typename T1, typename T2>
    struct AddVV_1Helper<N,T1,T2,true>
    {
        inline AddVV_1Helper(const T1* v1, T2* v2)
        { for(ptrdiff_t i=0;i<N;++i) v2[i] += v1[i]; }
    };
    template <ptrdiff_t N, typename T1, typename T2>
    inline void AddVV_1(const T1* v1, T2* v2)
    { AddVV_1Helper<N,T1,T2,OKTypes(T1,T2)>(v1,v2); }

    template <ptrdiff_t N, typename T1, typename T2, bool ok>
    struct AddVV_m1Helper
    { AddVV_m1Helper(const T1*, T2*) {} };
    template <ptrdiff_t N, typename T1, typename T2>
    struct AddVV_m1Helper<N,T1,T2,true>
    {
        inline AddVV_m1Helper(const T1* v1, T2* v2)
        { for(ptrdiff_t i=0;i<N;++i) v2[i] -= v1[i]; }
    };
    template <ptrdiff_t N, typename T1, typename T2>
    inline void AddVV_m1(const T1* v1, T2* v2)
    { AddVV_m1Helper<N,T1,T2,OKTypes(T1,T2)>(v1,v2); }

    template <ptrdiff_t N, typename T, typename T1, typename T2, bool ok>
    struct AddVVHelper
    { AddVVHelper(const T, const T1*, T2*) {} };
    template <ptrdiff_t N, typename T, typename T1, typename T2>
    struct AddVVHelper<N,T,T1,T2,true>
    {
        inline AddVVHelper(const T x, const T1* v1, T2* v2)
        { for(ptrdiff_t i=0;i<N;++i) v2[i] += x*v1[i]; }
    };
    template <ptrdiff_t N, typename T, typename T1, typename T2>
    inline void AddVV(const T x, const T1* v1, T2* v2)
    { AddVVHelper<N,T,T1,T2,OKTypes2(T,T1,T2)>(x,v1,v2); }

    // v3 = x1 * v1 + x2 * v2
    template <ptrdiff_t N, typename T1, typename T2, typename T3>
    inline void AddVV_1_1(const T1* v1, const T2* v2, T3* v3)
    { SmallVectorCopy<N>(v1,v3); AddVV_1<N>(v2,v3); }
    template <ptrdiff_t N, typename T1, typename T2, typename T3>
    inline void AddVV_1_m1(const T1* v1, const T2* v2, T3* v3)
    { SmallVectorCopy<N>(v1,v3); AddVV_m1<N>(v2,v3); }
    template <ptrdiff_t N, typename T, typename T1, typename T2, typename T3>
    inline void AddVV_1_x(const T1* v1, const T x2, const T2* v2, T3* v3)
    { MultXV<N>(x2,v2,v3); AddVV_1<N>(v1,v3); }
    template <ptrdiff_t N, typename T, typename T1, typename T2, typename T3>
    inline void AddVV_x_1(const T x1, const T1* v1, const T2* v2, T3* v3)
    { MultXV<N>(x1,v1,v3); AddVV_1<N>(v2,v3); }
    template <ptrdiff_t N, typename T, typename T1, typename T2, typename T3>
    inline void AddVV_x_m1(const T x1, const T1* v1, const T2* v2, T3* v3)
    { MultXV<N>(x1,v1,v3); AddVV_m1<N>(v2,v3); }
    template <ptrdiff_t N, typename T, typename T1, typename T2, typename T3>
    inline void AddVV(const T x1, const T1* v1, const T x2, const T2* v2, T3* v3)
    { MultXV<N>(x1,v1,v3); AddVV<N>(x2,v2,v3); }

    // v1 * v2 (dot product)
    template <typename T1, typename T2>
    struct ProdType { typedef T1 Tprod; };
    template <typename T>
    struct ProdType<T,std::complex<T> > { typedef std::complex<T> Tprod; };
#define ProductType(T1,T2) typename ProdType<T1,T2>::Tprod

    template <ptrdiff_t N, typename T, typename T1, typename T2>
    struct MultVVHelper
    {
        inline MultVVHelper(T& res, const T1* v1, const T2* v2)
        { for(ptrdiff_t i=0;i<N;++i) res += v1[i]*v2[i]; }
    };
    template <ptrdiff_t N, typename T1, typename T2>
    inline ProductType(T1,T2) MultVV(const T1* v1, const T2* v2)
    {
        ProductType(T1,T2) res(0);
        MultVVHelper<N,ProductType(T1,T2),T1,T2>(res,v1,v2);
        return res;
    }
#undef ProductType

    template <typename T, typename T1, typename T2, typename T3, ptrdiff_t N, int A1, int A2, int A3>
    inline void ElemMultVV(
        const T alpha, const SmallVector<T1,N,A1>& v1,
        const SmallVector<T2,N,A2>& v2, SmallVector<T3,N,A3>& v3)
    {
        if (alpha == T(1))
            for(ptrdiff_t i=0;i<N;++i) v3.ref(i) = v1.cref(i) * v2.cref(i);
        else
            for(ptrdiff_t i=0;i<N;++i) v3.ref(i) = alpha * v1.cref(i) * v2.cref(i);
    }
    template <typename T, typename T1, typename T2, ptrdiff_t N, int A1, int A2, int A3>
    inline void ElemMultVV(
        const CT , const SmallVector<T1,N,A1>& ,
        const SmallVector<T2,N,A2>& , SmallVector<T,N,A3>& )
    { TMVAssert(TMV_FALSE); }


    template <typename T, typename T1, typename T2, typename T3, ptrdiff_t N, int A1, int A2, int A3>
    inline void AddElemMultVV(
        const T alpha, const SmallVector<T1,N,A1>& v1,
        const SmallVector<T2,N,A2>& v2, SmallVector<T3,N,A3>& v3)
    {
        if (alpha == T(1))
            for(ptrdiff_t i=0;i<N;++i) v3.ref(i) += v1.cref(i) * v2.cref(i);
        else if (alpha == T(-1))
            for(ptrdiff_t i=0;i<N;++i) v3.ref(i) -= v1.cref(i) * v2.cref(i);
        else
            for(ptrdiff_t i=0;i<N;++i) v3.ref(i) += alpha * v1.cref(i) * v2.cref(i);
    }
    template <typename T, typename T1, typename T2, ptrdiff_t N, int A1, int A2, int A3>
    inline void AddElemMultVV(
        const CT , const SmallVector<T1,N,A1>& ,
        const SmallVector<T2,N,A2>& , SmallVector<T,N,A3>& )
    { TMVAssert(TMV_FALSE); }

#undef OKTypes
#undef OKTypes2

} // namespace tmv

#undef CT

#endif
