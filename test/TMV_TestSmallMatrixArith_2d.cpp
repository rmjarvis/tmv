
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2d()
{
    TestSmallMatrixArith_2<T,3,6>("3 6");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2d<int>();
#endif
