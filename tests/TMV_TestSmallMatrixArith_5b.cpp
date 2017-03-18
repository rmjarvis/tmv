
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5b()
{
    TestSmallMatrixArith_5<T,2,2,2>("2 2 2");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5b<int>();
#endif
