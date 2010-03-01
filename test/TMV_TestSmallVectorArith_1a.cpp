
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_Vec.h"
#include <fstream>
#include <cstdio>

#include "TMV_TestVectorArith.h"

template <int N, class T> static void DoTestSmallVectorArith_1a()
{
    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+2);

    tmv::SmallVector<std::complex<T>,N> ca = a*std::complex<T>(2,-1);;
    tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(-5,1);

    TestVectorArith1<T>(a,ca,"SmallVector C");

#if (XTEST & 32)
    tmv::SmallVector<T,N,tmv::FortranStyle> af = a;
    tmv::SmallVector<T,N,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> caf = ca;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> cbf = cb;

    TestVectorArith1<T>(af,caf,"SmallVector F");
#endif
}

template <class T> void TestSmallVectorArith_1a()
{
    DoTestSmallVectorArith_1a<2,T>();
    DoTestSmallVectorArith_1a<5,T>();
#if (XTEST & 2)
    DoTestSmallVectorArith_1a<1,T>();
    DoTestSmallVectorArith_1a<3,T>();
    DoTestSmallVectorArith_1a<4,T>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallVectorArith_1a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallVectorArith_1a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallVectorArith_1a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallVectorArith_1a<int>();
#endif
