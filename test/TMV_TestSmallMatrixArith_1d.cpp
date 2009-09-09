// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define NONSQUARE
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1d()
{
  TestSmallMatrixArith_1<T,3,6>("3 6");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1d<int>();
#endif
