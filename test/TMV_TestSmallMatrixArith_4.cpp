// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4()
{
  TestSmallMatrixArith_4a<T>();
  TestSmallMatrixArith_4b<T>();
  TestSmallMatrixArith_4c<T>();
  TestSmallMatrixArith_4d<T>();
  TestSmallMatrixArith_4e<T>();
  TestSmallMatrixArith_4f<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4<int>();
#endif
