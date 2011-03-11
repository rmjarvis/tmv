
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4b()
{
    TestSmallMatrixArith_4<T,2,2>("2 2");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4b<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4b<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4b<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4b<int>();
#endif
