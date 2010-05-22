
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6b()
{
    TestSmallMatrixArith_6<T,3,3,3>("3 3 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6b<int>();
#endif
