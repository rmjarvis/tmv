// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define NONSQUARE
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3f()
{
#ifdef XTEST
  TestSmallMatrixArith_3<T,339,607>("339 607");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3f<int>();
#endif
