
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5i()
{
#if (XTEST & 2)
    TestSmallMatrixArith_5<T,12,4,2>("12 4 2");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5i<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5i<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5i<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5i<int>();
#endif
