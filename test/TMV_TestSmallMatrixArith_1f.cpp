
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1f()
{
    TestSmallMatrixArith_1<T,3,6>("3 6");
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
