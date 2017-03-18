
#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2c()
{
    TestSmallMatrixArith_2<T,3,3>("3 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2c<int>();
#endif
