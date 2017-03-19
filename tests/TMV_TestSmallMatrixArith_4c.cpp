
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4c()
{
    TestSmallMatrixArith_4<T,3,3>("3 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4c<int>();
#endif
