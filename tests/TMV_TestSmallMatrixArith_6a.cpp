
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6a()
{
    TestSmallMatrixArith_6<T,1,1,1>("1 1 1");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6a<int>();
#endif
