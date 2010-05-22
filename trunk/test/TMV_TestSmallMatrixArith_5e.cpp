
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5e()
{
#if (XTEST & 2)
    TestSmallMatrixArith_5<T,1,10,8>("1 10 8");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5e<int>();
#endif
