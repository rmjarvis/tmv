
#include "TMV_TestSmallMatrixArith_7.h"

template <class T> void TestSmallMatrixArith_7a()
{
    TestSmallMatrixArith_7<T,1,1>("1 1");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_7a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_7a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_7a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_7a<int>();
#endif
