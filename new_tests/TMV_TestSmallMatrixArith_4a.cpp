
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4a()
{
    TestSmallMatrixArith_4<T,1,1>("1 1");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4a<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4a<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4a<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4a<int>();
#endif
