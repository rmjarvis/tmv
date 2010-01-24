
#include "TMV_TestSmallMatrixArith_7.h"

template <class T> void TestSmallMatrixArith_7()
{
    TestSmallMatrixArith_7a<T>();
    TestSmallMatrixArith_7b<T>();
    TestSmallMatrixArith_7c<T>();
    TestSmallMatrixArith_7d<T>();
    TestSmallMatrixArith_7e<T>();
    TestSmallMatrixArith_7f<T>();
}

#ifdef TEST_DOUBLE
template void TestSmallMatrixArith_7<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallMatrixArith_7<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallMatrixArith_7<long double>();
#endif
#ifdef TEST_INT
template void TestSmallMatrixArith_7<int>();
#endif
