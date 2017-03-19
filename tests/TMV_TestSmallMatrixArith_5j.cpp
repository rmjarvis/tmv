
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5j()
{
    TestSmallMatrixArith_5<T,39,60,49>("39 60 49");
#if XTEST & 2
    TestSmallMatrixArith_5<T,60,39,49>("60 39 49");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5j<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5j<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5j<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5j<int>();
#endif
