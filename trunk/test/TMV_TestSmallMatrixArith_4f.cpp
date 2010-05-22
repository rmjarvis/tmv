
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4f()
{
#if (XTEST & 2)
    TestSmallMatrixArith_4<T,339,607>("339 607");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4f<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4f<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4f<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4f<int>();
#endif
