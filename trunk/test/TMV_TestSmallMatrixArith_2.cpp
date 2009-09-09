// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2()
{
  TestSmallMatrixArith_2a<T>();
  TestSmallMatrixArith_2b<T>();
  TestSmallMatrixArith_2c<T>();
  TestSmallMatrixArith_2d<T>();
  TestSmallMatrixArith_2e<T>();
  TestSmallMatrixArith_2f<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2<int>();
#endif
