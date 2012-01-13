
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

#define EXPLICIT_ALIAS

#include "TMV_TestVectorArith.h"

template <int N, class T> static void DoTestSmallVectorArith_2a()
{
    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+2);

    tmv::SmallVector<std::complex<T>,N> ca = a*std::complex<T>(2,-1);;
    tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(-5,1);

    TestVectorArith2(a,ca,b,cb,"SmallVector CC");

#if (XTEST & 32)
    tmv::SmallVector<T,N,tmv::FortranStyle> af = a;
    tmv::SmallVector<T,N,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> caf = ca;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> cbf = cb;

    TestVectorArith2(af,caf,bf,cbf,"SmallVector FF");
    TestVectorArith2(a,ca,bf,cbf,"SmallVector CF");
    TestVectorArith2(af,caf,b,cb,"SmallVector FC");
#endif
}

template <class T> void TestSmallVectorArith_2a()
{
    DoTestSmallVectorArith_2a<2,T>();
    DoTestSmallVectorArith_2a<5,T>();
#if (XTEST & 2)
    DoTestSmallVectorArith_2a<1,T>();
    DoTestSmallVectorArith_2a<3,T>();
    DoTestSmallVectorArith_2a<4,T>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallVectorArith_2a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallVectorArith_2a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallVectorArith_2a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallVectorArith_2a<int>();
#endif
