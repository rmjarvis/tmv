
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_7.h"

template <class T> void TestSmallMatrixArith_7g()
{
    TestSmallMatrixArith_7<T,39,60>("39 60");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_7g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_7g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_7g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_7g<int>();
#endif
