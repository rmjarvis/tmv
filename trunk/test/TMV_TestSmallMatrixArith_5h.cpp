
#define NONSQUARE
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5h()
{
    TestSmallMatrixArith_5<T,4,2,3>("4 2 3");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5h<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5h<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5h<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5h<int>();
#endif
