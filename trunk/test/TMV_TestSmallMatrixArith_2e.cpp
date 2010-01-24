
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2e()
{
    TestSmallMatrixArith_2<T,1,10>("1 10");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2e<int>();
#endif
