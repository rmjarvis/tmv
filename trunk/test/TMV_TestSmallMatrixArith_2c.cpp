// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2c()
{
  TestSmallMatrixArith_2<T,4,4>("4 4");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2c<int>();
#endif
