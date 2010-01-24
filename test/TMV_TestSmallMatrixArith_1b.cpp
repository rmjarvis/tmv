
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1b()
{
    TestSmallMatrixArith_1<T,3,3>("3 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1b<int>();
#endif
