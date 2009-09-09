// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1a()
{
  TestSmallMatrixArith_1<T,2,2>("2 2");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1a<int>();
#endif
