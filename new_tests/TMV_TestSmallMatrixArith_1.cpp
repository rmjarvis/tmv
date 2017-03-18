
#include "TMV_TestSmallMatrixArith_1.h"

template <class T> void TestSmallMatrixArith_1()
{
    TestSmallMatrixArith_1a<T>();
    TestSmallMatrixArith_1b<T>();
    TestSmallMatrixArith_1c<T>();
    TestSmallMatrixArith_1d<T>();
    TestSmallMatrixArith_1e<T>();
    TestSmallMatrixArith_1f<T>();
    TestSmallMatrixArith_1g<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_1<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_1<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_1<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_1<int>();
#endif
