
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6f()
{
#if (XTEST & 2)
    TestSmallMatrixArith_6<T,2,1,12>("2 1 12");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6f<int>();
#endif
