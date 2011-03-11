
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6i()
{
    TestSmallMatrixArith_6<T,82,4,2>("82 4 2");
#if XTEST & 2
    TestSmallMatrixArith_6<T,4,82,2>("4 82 2");
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6i<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6i<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6i<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6i<int>();
#endif
