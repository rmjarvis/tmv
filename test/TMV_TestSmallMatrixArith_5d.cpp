
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5d()
{
    TestSmallMatrixArith_5<T,3,6,7>("3 6 7");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5d<int>();
#endif
