
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3e()
{
    TestSmallMatrixArith_3<T,1,10>("1 10");
#if XTEST & 2
    TestSmallMatrixArith_3<T,10,1>("10 1");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_3e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_3e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_3e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_3e<int>();
#endif
