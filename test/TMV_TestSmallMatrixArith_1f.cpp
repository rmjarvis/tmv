
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1f()
{
#if (XTEST & 2)
    TestSmallMatrixArith_1<T,339,607>("339 607");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1f<int>();
#endif
