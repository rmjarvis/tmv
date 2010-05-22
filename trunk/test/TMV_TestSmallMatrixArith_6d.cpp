
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6d()
{
    TestSmallMatrixArith_6<T,3,6,7>("3 6 7");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6d<int>();
#endif
