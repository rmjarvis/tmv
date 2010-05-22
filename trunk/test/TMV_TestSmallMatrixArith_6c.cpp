
#include "TMV_TestSmallMatrixArith_6.h"

template <class T> void TestSmallMatrixArith_6c()
{
    TestSmallMatrixArith_6<T,4,4,4>("4 4 4");
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_6c<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_6c<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_6c<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_6c<int>();
#endif
