
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5f()
{
    TestSmallMatrixArith_5<T,2,4,12>("2 4 12");
#if XTEST & 2
    TestSmallMatrixArith_5<T,4,2,12>("4 2 12");
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
