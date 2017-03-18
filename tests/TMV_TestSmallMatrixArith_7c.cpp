
#include "TMV_TestSmallMatrixArith_7.h"

template <class T> void TestSmallMatrixArith_7c()
{
    TestSmallMatrixArith_7<T,3,3>("3 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_7c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_7c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_7c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_7c<int>();
#endif
