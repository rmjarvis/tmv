
#include "TMV_TestSmallMatrixArith_5.h"

template <class T> void TestSmallMatrixArith_5c()
{
    TestSmallMatrixArith_5<T,4,4,4>("4 4 4");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_5c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_5c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_5c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_5c<int>();
#endif
