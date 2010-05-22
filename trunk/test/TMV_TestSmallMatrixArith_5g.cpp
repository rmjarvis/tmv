
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5g()
{
#if (XTEST & 2)
    TestSmallMatrixArith_5<T,9,3,1>("9 3 1");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5g<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5g<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5g<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5g<int>();
#endif
