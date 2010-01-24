
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1e()
{
#ifdef XTEST
#if (XTEST & 2)
    TestSmallMatrixArith_1<T,1,10>("1 10");
#endif
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1e<int>();
#endif
