
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6g()
{
    TestSmallMatrixArith_6<T,9,3,1>("9 3 1");
#if (XTEST & 2)
    TestSmallMatrixArith_6<T,3,9,1>("3 9 1");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6g<int>();
#endif
