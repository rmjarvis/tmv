
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5a()
{
    TestSmallMatrixArith_5<T,2,2,2>("2 2 2");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5a<int>();
#endif
