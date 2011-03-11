
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3f()
{
    TestSmallMatrixArith_3<T,3,6>("3 6");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3f<int>();
#endif
