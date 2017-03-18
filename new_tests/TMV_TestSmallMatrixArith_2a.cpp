
#include "TMV_TestSmallMatrixArith_2.h"

template <class T> void TestSmallMatrixArith_2a()
{
    TestSmallMatrixArith_2<T,1,1>("1 1");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_2a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_2a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_2a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_2a<int>();
#endif
