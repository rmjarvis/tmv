
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_4.h"

template <class T> void TestSmallMatrixArith_4e()
{
    TestSmallMatrixArith_4<T,1,10>("1 10");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_4e<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_4e<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_4e<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_4e<int>();
#endif
