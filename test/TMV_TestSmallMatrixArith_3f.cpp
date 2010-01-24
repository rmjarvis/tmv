
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_3.h"

template <class T> void TestSmallMatrixArith_3f()
{
#ifdef XTEST
#if (XTEST & 2)
    TestSmallMatrixArith_3<T,339,607>("339 607");
#endif
#endif
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
