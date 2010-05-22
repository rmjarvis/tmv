
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4d()
{
    TestSmallMatrixArith_4<T,3,6>("3 6");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4d<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4d<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4d<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4d<int>();
#endif
