
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3g()
{
    TestSmallMatrixArith_3<T,39,60>("39 60");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3g<int>();
#endif
