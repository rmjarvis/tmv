
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6e()
{
    TestSmallMatrixArith_6<T,1,10,8>("1 10 8");
#if (XTEST & 2)
    TestSmallMatrixArith_6<T,10,1,8>("10 1 8");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6e<int>();
#endif
