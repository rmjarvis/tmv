
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4g()
{
    TestSmallMatrixArith_4<T,39,60>("39 60");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4g<int>();
#endif
