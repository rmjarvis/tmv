
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_7.h"

template <class T> void TestSmallMatrixArith_7e()
{
#if (XTEST & 2)
    TestSmallMatrixArith_7<T,1,10>("1 10");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_7e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_7e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_7e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_7e<int>();
#endif
