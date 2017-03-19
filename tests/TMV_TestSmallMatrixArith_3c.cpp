
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3c()
{
    TestSmallMatrixArith_3<T,3,3>("3 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3c<int>();
#endif
