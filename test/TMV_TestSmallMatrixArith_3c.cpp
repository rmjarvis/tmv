// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3c()
{
  TestSmallMatrixArith_3<T,4,4>("4 4");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3c<int>();
#endif
