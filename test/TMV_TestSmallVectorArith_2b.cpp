// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_Vec.h"
#include <fstream>
#include <fstream>
#include <cstdio>

#include "TMV_TestVectorArith.h"

template <int N, class T> static void DoTestSmallVectorArith_2b()
{
  tmv::SmallVector<T,N> a;
  for(int i=0;i<N;++i) a(i) = T(i+10);
  tmv::SmallVector<T,N> b;
  for(int i=0;i<N;++i) b(i) = T(-3*i+2);

  tmv::SmallVector<std::complex<T>,N> ca = a*std::complex<T>(2,-1);;
  tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(-5,1);

  tmv::VectorView<T> av = a.View();
  tmv::VectorView<std::complex<T> > cav = ca.View();
  tmv::VectorView<T> bv = b.View();
  tmv::VectorView<std::complex<T> > cbv = cb.View();

  TestVectorArith2<T>(av,cav,b,cb,"SmallVector/Vector");
}

template <class T> void TestSmallVectorArith_2b()
{
  DoTestSmallVectorArith_2b<2,T>();
  DoTestSmallVectorArith_2b<5,T>();
#ifdef XTEST
  DoTestSmallVectorArith_2b<1,T>();
  DoTestSmallVectorArith_2b<3,T>();
  DoTestSmallVectorArith_2b<4,T>();
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallVectorArith_2b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallVectorArith_2b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallVectorArith_2b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallVectorArith_2b<int>();
#endif
