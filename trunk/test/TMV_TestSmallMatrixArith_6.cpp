// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6()
{
  TestSmallMatrixArith_6a<T>();
  TestSmallMatrixArith_6b<T>();
  TestSmallMatrixArith_6c<T>();
  TestSmallMatrixArith_6d<T>();
  TestSmallMatrixArith_6e<T>();
  TestSmallMatrixArith_6f<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6<int>();
#endif
