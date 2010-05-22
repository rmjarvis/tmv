
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5f()
{
#if (XTEST & 2)
    TestSmallMatrixArith_5<T,2,1,12>("2 1 12");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5f<int>();
#endif
