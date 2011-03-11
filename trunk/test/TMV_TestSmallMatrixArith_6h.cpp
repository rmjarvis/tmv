
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6h()
{
    TestSmallMatrixArith_6<T,2,2,79>("2 2 79");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6h<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6h<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6h<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6h<int>();
#endif
