
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1g()
{
    TestSmallMatrixArith_1<T,39,60>("39 60");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1g<int>();
#endif
