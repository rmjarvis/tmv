
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#include <fstream>
#include <fstream>
#include <cstdio>

#include "TMV_TestVectorArith.h"

template <int N, class T> static void DoTestSmallVectorArith_2c()
{
    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+2);

    tmv::SmallVector<std::complex<T>,N> ca = a*std::complex<T>(2,-1);;
    tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(-5,1);

    tmv::VectorView<T> av = a.view();
    tmv::VectorView<std::complex<T> > cav = ca.view();
    tmv::VectorView<T> bv = b.view();
    tmv::VectorView<std::complex<T> > cbv = cb.view();

    TestVectorArith2(a,ca,bv,cbv,"Vector/SmallVector");
}

template <class T> void TestSmallVectorArith_2c()
{
    DoTestSmallVectorArith_2c<2,T>();
    DoTestSmallVectorArith_2c<5,T>();
#if (XTEST & 2)
    DoTestSmallVectorArith_2c<1,T>();
    DoTestSmallVectorArith_2c<3,T>();
    DoTestSmallVectorArith_2c<4,T>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallVectorArith_2c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallVectorArith_2c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallVectorArith_2c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallVectorArith_2c<int>();
#endif
