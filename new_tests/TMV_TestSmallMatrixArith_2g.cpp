
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2g()
{
    TestSmallMatrixArith_2<T,39,60>("39 60");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2g<int>();
#endif
