
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6j()
{
#if (XTEST & 2)
    TestSmallMatrixArith_6<T,339,607,489>("339 607 489");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6j<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6j<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6j<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6j<int>();
#endif
